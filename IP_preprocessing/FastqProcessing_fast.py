#####################################################################################################
# Author    HJ Han
# Date
# Editor    YHLee
# Revised   2021-06-22 ~
# Note      (1) Excluding the step of profiling quality scores of merged (assembled) files to increase the speed
#               -> Common.PESHandle_no_quality_check
#           (2) Add prefix to sample name to see sample names correctly in Excel
#               -> if first component of file name is an integer, 'f_' prefix is added
#
#####################################################################################################

import getopt
import logging
import multiprocessing
import os
import sys
from datetime import datetime
from threading import Thread

import glob


from utils.MultiprocessLogger import log_listener
from utils.PESHandle_no_quality_check import PESParser
from utils.ProcessUtil import read_config
from utils.SeqExtracter_newMP import SeqExtracter
from utils.FileUtil import read_file_fast


# Change label of the raw files for the high readability
def readable_labels(
    meta_file: str,
    gz_dir: str,
    col_before: str = "Sample_name",
    col_after: str = "Label",
):
    file_format = os.path.basename(meta_file).split(".")[-1]
    if file_format == "csv":
        delim = ","
    elif file_format == "tsv":
        delim = "\t"
    else:
        print("Not supported file format. Must be csv or tsv.")
        sys.exit()
    change_table = {}
    new_labels = {}
    header, data = read_file_fast(meta_file, delim)
    for d in data:
        # if d[col_before] in change_table:
        #     print(
        #         "Duplicated original label exists. Please check the original labels in metadata file."
        #     )
        #     sys.exit()
        if d[col_after] == '':
            continue
        if d[col_after] in new_labels:
            print(
                "Duplicate new label exists. Please check the new labels in metadata file."
            )
            sys.exit()
        change_table[d[col_before]] = d[col_after]
        new_labels[d[col_after]] = None

    ## (230418) by SB
    # to prevent incomplete label change 
    # ex) Correct: 'origin1->label1_BCR', 'origin10->label10_BCR' / Incomplete change: 'origin10'->'label1_BCR0'
    change_table = {k: v for k, v in sorted(change_table.items(), key=lambda x: len(x[0]), reverse=True)}

    change_cnt = 0
    whole_cnt = len(glob.glob(os.path.join(gz_dir, "*.fastq*")))
    for origin, new in change_table.items():
        infiles = glob.glob(os.path.join(gz_dir, "%s*.fastq*" % origin))
        for infile in infiles:
            file_suffix = os.path.basename(infile).replace(origin, "")
            outfile = os.path.join(gz_dir, "%s%s" % (new, file_suffix))
            os.system("mv %s %s" % (infile, outfile))
            change_cnt += 1
    logging.info("File name changed: %s/%s" % (change_cnt, whole_cnt))


# Return dictionary for the extractin of q-filtered csv into project dir -> label dir
def project_by_label(meta_file: str, gz_dir: str, col_label: str = "Label"):
    file_format = os.path.basename(meta_file).split(".")[-1]
    if file_format == "csv":
        delim = ","
    elif file_format == "tsv":
        delim = "\t"
    else:
        print("Not supported file format. Must be csv or tsv.")
        sys.exit()
    project_dict = {}
    header, data = read_file_fast(meta_file, delim)
    for d in data:
        project_name = d["Project"]
        label = d[col_label]
        
        if label == '':
            continue
        if len(glob.glob(os.path.join(gz_dir, "%s*.fastq*" % label))) < 2:
            continue
        if label in project_dict:
            print("Duplicate label was found. Please check the label names.")
            sys.exit()
        project_dict[label] = project_name
    return project_dict


def main():
    # Maximum process numbers for the system
    max_process = multiprocessing.cpu_count() * 0.8

    # default values for argument values of each option
    default_process = 5
    _dir, step, num_max_threads, platform, is_gz, meta = (
        None,
        None,
        default_process,
        None,
        "True",
        None,
    )
    py_name = os.path.basename(__file__)

    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "hd:s:p:j:g:m:",
            ["help", "dir=", "step=", "platform=", "threads=", "is_gz=", "meta_file="],
        )
    except getopt.GetoptError:
        print(
            "USAGE: %s [-h] [-d directory] [-s step] [-p platform] [-j threads] [-g is_gz] [-m meta_file]"
            % py_name
        )
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(
                """
USAGE: %s [-h] [-d directory] [-s step] [-p platform] [-j threads] [-g is_gz] [-m meta_file]

Optional arguments:
    -h, --help          show this help message and exit
    -d, --dir           working directory path (containing gz folder and config.conf file)
    -s, --step          step integer (1, 2). if None, run whole steps
                            1: Run PEAR (pear-end merger)
                            2: Apply Q-filtering, extract data (full_NT, readcount) and statistics
                        <default = None>
    -p, --platform      Illumina platform (nextseq, miseq, novaseq)
                            Apply length filtration by platform (nextseq: 150bp, miseq/novaseq: 250bp)
                        <default = None (no length filtration)>
    -j, --threads       number of threads
                        <default = 5>
    -g, --is_gz         whether fastq files are gz format or not (True / False)
                        <default = True>
    -m, --meta_file     path of metadata file
"""
                % py_name
            )
            sys.exit()
        elif opt in ("-d", "--dir"):
            _dir = arg
        elif opt in ("-s", "--step"):
            step = arg
        elif opt in ("-j", "--threads"):
            num_max_threads = int(arg)
            if num_max_threads > max_process:
                num_max_threads = max_process
        elif opt in ("-p", "--platform"):
            platform = arg
        elif opt in ("-g", "--gz_file"):
            is_gz = arg
            if is_gz not in ("True", "False"):
                is_gz = "True"
        elif opt in ("-m", "--meta_file"):
            meta = arg

    if not _dir:
        print(
            "USAGE: %s [-h] [-d directory] [-s step] [-p platform] [-j threads] [-g is_gz] [-m meta_file]"
            % py_name
        )
        sys.exit(2)
    else:
        _dir = os.path.normpath(_dir)

    # read config
    config_file = os.path.join(_dir, "config.conf")
    config = read_config(config_file)
    # dbname = os.path.basename(_dir).lower().replace('-', '_')

    # start logging thread
    log_queue = multiprocessing.Queue(-1)
    log_thread = Thread(target=log_listener, args=(log_queue, config_file), daemon=True)
    log_thread.start()

    # Process start
    start_time = datetime.now()

    # Step 1. PES Parsing
    if not step or "1" in step:
        # File name change for the readability (using metadata file)
        readable_labels(meta_file=meta, gz_dir=os.path.join(_dir, "gz"))
        parser = PESParser(_dir, config["PESParser"], is_gz)
        parser.run(num_max_threads)

    # Step 2. raw data csv & stats extraction
    if not step or "2" in step:
        try:
            config["SeqExtracter"]["platform"] = platform
        except:
            config["SeqExtracter"]["platform"] = ""
        project_by_label_dict = project_by_label(
            meta_file=meta, gz_dir=os.path.join(_dir, "gz")
        )
        extracter = SeqExtracter(
            _dir, config["SeqExtracter"], project_by_label_dict, log_queue
        )
        extracter.run(num_max_threads)

        # final output check
        merged_files = glob.glob(os.path.join(_dir, "merge", "*.assembled.fastq"))
        
        ## (230418) by SB
        # to correctly capture processed csv file in each file directory 
        csv_files = []
        for project_dir in os.listdir(os.path.join(_dir, "out")):
            for sample_dir in os.listdir(os.path.join(_dir, "out", project_dir)):
                sample_dir = os.path.join(_dir, "out", project_dir, sample_dir)
                csv_files.extend(glob.glob(os.path.join(sample_dir, "*.csv")))
            
        logging.info(
            "Csv output stat: %d/%d files are well extracted..."
            % (len(csv_files), len(merged_files))
        )

    logging.info("--- %s seconds elapsed ---" % (datetime.now() - start_time))
    log_queue.put_nowait(None)
    log_thread.join()


if __name__ == "__main__":
    main()
