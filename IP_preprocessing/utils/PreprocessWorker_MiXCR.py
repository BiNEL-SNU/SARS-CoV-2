#####################################################################################################
# Author    YHLee
# Date      2021-04-26 ~
# Editor
# Revised   2021-06-30 (YHLee),
# Note      user module for multiprocessing in pre-processing script.
#           version 1.1
#           applying MiXCR, not UMI processing
#           Merge whole case: rev_primer_trimming option
#####################################################################################################

import csv
import logging
import multiprocessing
import os
import re
from copy import deepcopy
import tempfile
import socket


from multiprocessing import Process
from datetime import datetime
from .MultiprocessLogger import QueueHandler
from .SeqUtil import Helper, is_nucleotide_or_protein
from .FileUtil import read_file_fast, write_file
from .Enum import ChainType
from .ProcessUtil import run_igblast_new

# for memory usage test
import psutil

# New import for blast_run
from Bio.Blast.Applications import NcbiblastnCommandline

# get host name
hostname = socket.gethostname()


class PreprocessWorker(object):
    """Do multiple works for preprocessing.
    step 1. Error correction
            1a) primer recognition -> return dict
            1b) apply MiXCR -> return error_correction file
    step 2. Isotype annotation (by BLAST)
    step 3. Ig-related Information annotation (by IgBLAST)
    step 4. functional reads filtration
    """

    def __init__(
        self,
        work_dir,
        sample_name_list,
        infile_suffix="",
        outfile_suffix="",
        fwd_primers=None,
        rev_primers=None,
        c_ref_file=None,
        evalue=0.001,
        step=None,
        multiprocess_log_queue=None,
        chain_type=None,
        file_type="tsv",
        debug=False,
        rev_primer_trimming=True,
    ):
        # Setup working directory
        if not os.path.isdir(work_dir):
            logging.info("Not valid working directory path...")
            raise RuntimeError
        else:
            self.dir = work_dir

        self.sample_list = sample_name_list
        self.log_queue = multiprocess_log_queue

        self.infile_suffix = infile_suffix
        self.outfile_suffix = outfile_suffix

        self.fwd_primers = fwd_primers
        self.rev_primers = rev_primers
        self.c_ref_file = c_ref_file
        self.evalue = evalue
        self.step = step

        self.chain_type = chain_type

        if file_type in ["tsv", "csv"]:
            self.file_type = file_type
        else:
            logging.info("Not supported file type... Must be csv or tsv")
            raise ValueError
        self.delim = "," if file_type == "csv" else "\t"

        self.initial_RAM_used = psutil.virtual_memory().used / 1024**3

        self.debug = debug
        self.rev_primer_trimming = rev_primer_trimming

    def run(self, num_process):
        set_num_process = max(2, num_process)
        if hostname == "JunhoLab" and set_num_process > 30:
            set_num_process = (
                30  # max number of process is limited to 30 (Junholab server)
            )
        manager = multiprocessing.Manager()

        # for the RAM usage check
        max_RAM_used = psutil.virtual_memory().used / 1024**3

        chunk_size_max = 20000000

        # copy the sample_list. Considering that whole steps will be done at once...
        samples = self.sample_list.copy()

        while True:
            # Do works according to step
            if len(samples) == 0:
                break
            self.sample_name_to_put = samples.pop(0)
            sample_dir = os.path.join(self.dir, self.sample_name_to_put)

            # case separation: Step 1, 2, 3, 4
            if self.step == 1:
                start_time_main = datetime.now()
                logging.info("Started Error Correction for %s", self.sample_name_to_put)

                in_dir = sample_dir
                out_dir = os.path.join(sample_dir, "1_error_correction")
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)
                target_file = os.path.join(in_dir, "%s.csv" % self.sample_name_to_put)
                out_file = os.path.join(
                    out_dir,
                    "%s_%s.%s"
                    % (self.sample_name_to_put, self.outfile_suffix, self.file_type),
                )

                ##################################################################################
                ######## sub-step: 1a. primer recognition ########################################
                ##################################################################################
                if True:
                    # 1a. primer recognition
                    sub_step = "a"
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    if len(self.fwd_primers) < 1 or len(self.rev_primers) < 1:
                        primer_recognized_dict = {}
                        with open(target_file, "r") as handle:
                            header = handle.readline().strip().split(",")
                            col_dict = {h: header.index(h) for h in header}
                            while True:
                                line = handle.readline()
                                if not line:
                                    break
                                d = line.strip().split(",")
                                seq = d[col_dict["full_NT"]]
                                rc = int(d[col_dict["readcount"]])
                                if seq not in primer_recognized_dict:
                                    primer_recognized_dict[seq] = {
                                        "full_NT": seq,
                                        "readcount": rc,
                                    }
                                else:
                                    primer_recognized_dict[seq]["readcount"] += rc

                    else:
                        # get chunk size according to num_process and put it to the queue
                        file_line_len = (
                            self.get_file_lines(target_file) - 1
                        )  # exclude header line (csv file)
                        chunk_size = int(file_line_len / set_num_process) + 1
                        if chunk_size > file_line_len:
                            chunk_size = file_line_len
                        elif chunk_size > chunk_size_max:
                            chunk_size = int(
                                chunk_size / (int(chunk_size / chunk_size_max) + 1)
                            )
                        chunk_start_pos_list = self.extract_chunk_start_pos(
                            file_line_len, chunk_size
                        )

                        target_queue = manager.Queue()
                        for chunk_range in chunk_start_pos_list:
                            target_queue.put(chunk_range)

                        # make queue to get result data from processes
                        result_queue = manager.Queue()

                        ############################## started MP #########
                        # start ErrorCorrecter processes (MP)
                        logging.info(
                            "Spawning %d ErrorCorrecter processes for step 1a.",
                            set_num_process,
                        )
                        extract_workers = []
                        for i in range(set_num_process):
                            p = ErrorCorrecter(
                                sub_step,
                                result_queue,
                                self.log_queue,
                                RAM_used_queue,
                                in_queue=target_queue,
                                chunk_size=chunk_size,
                                input_file=target_file,
                                _dir=in_dir,
                                fwd_primer_list=self.fwd_primers,
                                rev_primer_list=self.rev_primers,
                                rev_primer_trimming=self.rev_primer_trimming,
                            )

                            p.daemon = True
                            p.start()

                            extract_workers.append(p)
                            target_queue.put_nowait(None)

                        # wait until all extracting works done
                        for p in extract_workers:
                            p.join()
                        ############################## finished MP #########

                        logging.info("1a. Primer recognition was done.")
                        logging.info(
                            "--- %s seconds elapsed ---" % (datetime.now() - start_time)
                        )

                        # merge result into one dict
                        out_result_queue = p.get_result_queue()
                        primer_recognized_dict = self.merge_chunk_dict_type(
                            out_result_queue
                        )

                    # temporary file write for the comparison
                    if self.debug:
                        outfile_dir = sample_dir
                        output_file = os.path.join(
                            outfile_dir,
                            "%s_%s_primer_recognized.tsv"
                            % (self.sample_name_to_put, self.outfile_suffix),
                        )

                        w_header = ["full_NT", "readcount"]
                        with open(output_file, "w") as handle:
                            write_file(
                                handle,
                                w_header,
                                list(primer_recognized_dict.values()),
                                delim="\t",
                            )

                ##################################################################################
                ######## sub-step: 1b. apply MiXCR ###############################################
                ##################################################################################
                if True:
                    # 1b. MiXCR-based error correction
                    sub_step = "b"
                    start_time = datetime.now()

                    # start ErrorCorrecter processes (no MP)
                    logging.info(
                        "Applying MiXCR-based error correction using ErrorCorrecter processes started."
                    )

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    extract_workers = []
                    p = ErrorCorrecter(
                        sub_step,
                        result_queue,
                        self.log_queue,
                        RAM_used_queue,
                        input_dict=primer_recognized_dict,
                        evalue=self.evalue,
                    )
                    p.daemon = True
                    p.start()

                    extract_workers.append(p)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()

                    logging.info("1b. Applying MiXCR was done.")
                    logging.info(
                        "--- %s seconds elapsed ---" % (datetime.now() - start_time)
                    )

                    # merge result into one dict
                    out_result_queue = p.get_result_queue()
                    error_corrected_list = out_result_queue.get()

                    # get RAM_used_queue
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    # write result files
                    w_header = ["full_NT", "readcount"]
                    with open(out_file, "w") as handle:
                        write_file(
                            handle, w_header, error_corrected_list, delim=self.delim
                        )

                logging.info(
                    "Finished Error Correction for %s", self.sample_name_to_put
                )
                logging.info(
                    "--- %s seconds elapsed ---" % (datetime.now() - start_time_main)
                )

            elif self.step == 2:
                start_time_main = datetime.now()
                logging.info("Started Isotyping for %s", self.sample_name_to_put)

                in_dir = os.path.join(sample_dir, "1_error_correction")
                out_dir = os.path.join(sample_dir, "2_annotation")
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)

                if self.infile_suffix == "":
                    target_file = os.path.join(
                        in_dir, "%s.%s" % (self.sample_name_to_put, self.file_type)
                    )
                else:
                    target_file = os.path.join(
                        in_dir,
                        "%s_%s.%s"
                        % (self.sample_name_to_put, self.infile_suffix, self.file_type),
                    )
                if self.outfile_suffix == "":
                    out_file = os.path.join(
                        out_dir, "%s.%s" % (self.sample_name_to_put, self.file_type)
                    )
                else:
                    out_file = os.path.join(
                        out_dir,
                        "%s_%s.%s"
                        % (
                            self.sample_name_to_put,
                            self.outfile_suffix,
                            self.file_type,
                        ),
                    )

                # check if error-corrected file exists or not
                if not os.path.exists(target_file):
                    logging.info(
                        "Error corrected file does not exist... pass following processes"
                    )
                    continue

                ##################################################################################
                ######## sub-step: 2a. run blast on C gene reference #############################
                ##################################################################################
                if True:
                    sub_step = "a"
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    ############################## started MP #########
                    ## Applying the multiprocessing that the BLAST tool supports
                    logging.info(
                        "Do blast run (step 2a) of IsotypeExtractor with %d processes",
                        set_num_process,
                    )
                    extract_workers = []
                    p = IsotypeExtractor(
                        target_file,
                        sub_step,
                        self.log_queue,
                        RAM_used_queue,
                        out_dir=out_dir,
                        c_ref_file=self.c_ref_file,
                        num_process=set_num_process,
                        delim=self.delim,
                        rev_primer_trimming=self.rev_primer_trimming,
                    )
                    p.daemon = True
                    p.start()
                    extract_workers.append(p)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # RAM usage update -> not for sure... because it cannot monitor the inside the tool BLAST...
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info("2a. Run blast on C gene reference was done.")
                    logging.info(
                        "--- %s seconds elapsed ---" % (datetime.now() - start_time)
                    )

                ##################################################################################
                ######## sub-step: 2b. Parsing blast result file #################################
                ##################################################################################
                if True:
                    sub_step = "b"
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # blast_result out directory
                    blast_result_dir = os.path.join(out_dir + "temp_cblast", "result")

                    # get chunk size according to num_process and put it to the queue
                    blast_parsing_iterator_size = (
                        self.get_file_lines(target_file) - 1
                    )  # exclude header
                    chunk_size = int(blast_parsing_iterator_size / set_num_process) + 1
                    if chunk_size > blast_parsing_iterator_size:
                        chunk_size = blast_parsing_iterator_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(
                            chunk_size / (int(chunk_size / chunk_size_max) + 1)
                        )
                    chunk_start_pos_list = self.extract_chunk_start_pos(
                        blast_parsing_iterator_size, chunk_size
                    )

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    logging.info(
                        "Do parse isotype info (step 2b) of IsotypeExtracter with %d processes",
                        set_num_process,
                    )
                    extract_workers = []

                    for i in range(set_num_process):
                        p = IsotypeExtractor(
                            target_file,
                            sub_step,
                            self.log_queue,
                            RAM_used_queue,
                            out_dir=blast_result_dir,
                            in_queue=target_queue,
                            out_queue=result_queue,
                            chunk_size=chunk_size,
                            delim=self.delim,
                        )
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # combine all chunk result, sort following the line order of original input file (error-corrected file), and merge
                    out_result_queue = p.get_result_queue()
                    isotype_added_dict = self.merge_chunk_dict_type(out_result_queue)

                    # write file
                    w_header = ["full_NT", "isotype", "readcount"]
                    with open(out_file, "w") as handle:
                        write_file(
                            handle,
                            w_header,
                            list(isotype_added_dict.values()),
                            delim=self.delim,
                        )

                    # RAM usage update
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info("2b. Parse isotype information was done.")
                    logging.info(
                        "--- %s seconds elapsed ---" % (datetime.now() - start_time)
                    )

                logging.info("Finished Isotyping for %s", self.sample_name_to_put)
                logging.info(
                    "--- %s seconds elapsed ---" % (datetime.now() - start_time_main)
                )

            elif self.step == 3:
                start_time_main = datetime.now()
                logging.info("Started Annotation for %s", self.sample_name_to_put)

                in_dir = os.path.join(sample_dir, "2_annotation")
                out_dir = in_dir
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)

                if self.infile_suffix == "":
                    target_file = os.path.join(
                        in_dir, "%s.%s" % (self.sample_name_to_put, self.file_type)
                    )
                else:
                    target_file = os.path.join(
                        in_dir,
                        "%s_%s.%s"
                        % (self.sample_name_to_put, self.infile_suffix, self.file_type),
                    )
                if self.outfile_suffix == "":
                    out_file = os.path.join(
                        out_dir, "%s.%s" % (self.sample_name_to_put, self.file_type)
                    )
                else:
                    out_file = os.path.join(
                        out_dir,
                        "%s_%s.%s"
                        % (
                            self.sample_name_to_put,
                            self.outfile_suffix,
                            self.file_type,
                        ),
                    )

                # check if C gene annotated file exists or not
                if not os.path.exists(target_file):
                    logging.info(
                        "C gene annotated file does not exist... pass following processes"
                    )
                    continue

                # count line number of file and continue
                input_file_line = self.get_file_lines(target_file)
                if input_file_line == 1:
                    logging.info("Empty input file. pass the annotation step...")
                    continue

                ##################################################################################
                ######## sub-step: 3a. run IgBLAST  ##############################################
                ##################################################################################
                if True:
                    sub_step = "a"
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    ############################## started MP #########
                    ## Applying the multiprocessing that the BLAST tool supports
                    logging.info(
                        "Do IgBLAST run (step 3a) of InfoAnnotator with %d processes",
                        set_num_process,
                    )
                    extract_workers = []
                    p = InfoAnnotator(
                        target_file,
                        sub_step,
                        self.log_queue,
                        RAM_used_queue,
                        out_dir,
                        self.chain_type,
                        num_process=set_num_process,
                        delim=self.delim,
                    )
                    p.daemon = True
                    p.start()
                    extract_workers.append(p)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # RAM usage update -> not for sure... because it cannot monitor the inside the tool BLAST...
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info("3a. Run IgBLAST was done.")
                    logging.info(
                        "--- %s seconds elapsed ---" % (datetime.now() - start_time)
                    )

                ##################################################################################
                ######## sub-step: 3b. Parsing igblast result file ###############################
                ##################################################################################
                if True:
                    sub_step = "b"
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # igblast_result out directory
                    igblast_result_dir = os.path.join(out_dir, "temp_igblast")

                    # get chunk size according to num_process and put it to the queue
                    igblast_parsing_size = (
                        self.get_file_lines(target_file) - 1
                    )  # exclude header
                    chunk_size = int(igblast_parsing_size / set_num_process) + 1
                    if chunk_size > igblast_parsing_size:
                        chunk_size = igblast_parsing_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(
                            chunk_size / (int(chunk_size / chunk_size_max) + 1)
                        )
                    chunk_start_pos_list = self.extract_chunk_start_pos(
                        igblast_parsing_size, chunk_size
                    )

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    logging.info(
                        "Do info annotation (step 3b) of InfoAnnotator with %d processes",
                        set_num_process,
                    )
                    extract_workers = []

                    for i in range(set_num_process):
                        p = InfoAnnotator(
                            target_file,
                            sub_step,
                            self.log_queue,
                            RAM_used_queue,
                            igblast_result_dir,
                            self.chain_type,
                            in_queue=target_queue,
                            out_queue=result_queue,
                            chunk_size=chunk_size,
                            delim=self.delim,
                        )
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # combine all chunk result, sort following the line order of original input file (error-corrected file), and merge
                    out_result_queue = p.get_result_queue()
                    annotated_dict = self.merge_chunk_dict_type(out_result_queue)
                    # annotated_dict_origin_order = {k:v for k, v in sorted(annotated_dict.items(), key=lambda x:x[0], reverse=False)}
                    sorted_annotated_dict = {
                        k: v
                        for k, v in sorted(
                            annotated_dict.items(),
                            key=lambda x: x[1]["duplicate_count"],
                            reverse=True,
                        )
                    }

                    # write file
                    w_header = [
                        "duplicate_count",
                        "sequence",
                        "sequence_aa",
                        "sequence_alignment",
                        "germline_alignment",
                        "v_call",
                        "d_call",
                        "j_call",
                        "c_call",
                        "junction",
                        "junction_aa",
                        "cdr1",
                        "cdr2",
                        "cdr3",
                        "cdr1_aa",
                        "cdr2_aa",
                        "cdr3_aa",
                        "rev_comp",
                        "productive",
                        "v_cigar",
                        "d_cigar",
                        "j_cigar",
                        "v_alignment_length",
                        "v_alignment_mutation",
                    ]
                    if self.chain_type not in [
                        ChainType.HUMAN_HEAVY,
                        ChainType.MOUSE_HEAVY,
                        ChainType.RABBIT_HEAVY,
                        ChainType.HUMAN_BETA,
                        ChainType.HUMAN_DELTA,
                        ChainType.MOUSE_BETA,
                        ChainType.MOUSE_DELTA,
                    ]:
                        w_header.remove("d_call")
                        w_header.remove("d_cigar")
                    with open(out_file, "w") as handle:
                        write_file(
                            handle,
                            w_header,
                            list(sorted_annotated_dict.values()),
                            delim=self.delim,
                        )

                    # RAM usage update
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info("3b. Info annotation was done.")
                    logging.info(
                        "--- %s seconds elapsed ---" % (datetime.now() - start_time)
                    )

                logging.info("Finished Annotation for %s", self.sample_name_to_put)
                logging.info(
                    "--- %s seconds elapsed ---" % (datetime.now() - start_time_main)
                )

            elif self.step == 4:
                start_time_main = datetime.now()
                logging.info(
                    "Started Functionality check for %s", self.sample_name_to_put
                )

                in_dir = os.path.join(sample_dir, "2_annotation")
                out_dir = os.path.join(sample_dir, "3_functionality")
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)

                if self.infile_suffix == "":
                    target_file = os.path.join(
                        in_dir, "%s.%s" % (self.sample_name_to_put, self.file_type)
                    )
                else:
                    target_file = os.path.join(
                        in_dir,
                        "%s_%s.%s"
                        % (self.sample_name_to_put, self.infile_suffix, self.file_type),
                    )
                if self.outfile_suffix == "":
                    out_file = os.path.join(
                        out_dir, "%s.%s" % (self.sample_name_to_put, self.file_type)
                    )
                else:
                    out_file = os.path.join(
                        out_dir,
                        "%s_%s.%s"
                        % (
                            self.sample_name_to_put,
                            self.outfile_suffix,
                            self.file_type,
                        ),
                    )

                # check if annotated file exists or not
                if not os.path.exists(target_file):
                    logging.info(
                        "Annotated file does not exist... pass following processes"
                    )
                    continue

                # count line number of file and continue
                input_file_line = self.get_file_lines(target_file)
                if input_file_line == 1:
                    logging.info(
                        "Empty input file. pass the functional reads filtration step..."
                    )
                    continue

                if True:
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # get chunk size according to num_process and put it to the queue
                    target_file_size = (
                        self.get_file_lines(target_file) - 1
                    )  # exclude header
                    chunk_size = int(target_file_size / set_num_process) + 1
                    if chunk_size > target_file_size:
                        chunk_size = target_file_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(
                            chunk_size / (int(chunk_size / chunk_size_max) + 1)
                        )
                    chunk_start_pos_list = self.extract_chunk_start_pos(
                        target_file_size, chunk_size
                    )

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    logging.info(
                        "Do functional reads extraction (step 4) of functionalExtracter with %d processes",
                        set_num_process,
                    )
                    extract_workers = []

                    for i in range(set_num_process):
                        p = functionalExtracter(
                            target_file,
                            self.log_queue,
                            RAM_used_queue,
                            in_queue=target_queue,
                            out_queue=result_queue,
                            chunk_size=chunk_size,
                            delim=self.delim,
                        )
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # combine all chunk result, sort following the line order of readcount (duplicate_count), and merge
                    out_result_queue = p.get_result_queue()
                    functional_dict = self.merge_chunk_dict_type(out_result_queue)
                    sorted_functional_dict = {
                        k: v
                        for k, v in sorted(
                            functional_dict.items(),
                            key=lambda x: int(x[1]["duplicate_count"]),
                            reverse=True,
                        )
                    }

                    # set sequence_id & frequency
                    row_num = 1
                    rc_sum = sum(
                        [
                            int(d["duplicate_count"])
                            for d in sorted_functional_dict.values()
                        ]
                    )
                    for k in sorted_functional_dict:
                        sorted_functional_dict[k]["sequence_id"] = "%s-%d" % (
                            self.sample_name_to_put,
                            row_num,
                        )
                        sorted_functional_dict[k]["frequency"] = (
                            int(sorted_functional_dict[k]["duplicate_count"]) / rc_sum
                        )
                        row_num += 1

                    # write file
                    w_header = [
                        "sequence_id",
                        "duplicate_count",
                        "frequency",
                        "sequence",
                        "sequence_aa",
                        "sequence_alignment",
                        "germline_alignment",
                        "v_call",
                        "d_call",
                        "j_call",
                        "c_call",
                        "junction",
                        "junction_aa",
                        "cdr1",
                        "cdr2",
                        "cdr3",
                        "cdr1_aa",
                        "cdr2_aa",
                        "cdr3_aa",
                        "rev_comp",
                        "productive",
                        "v_cigar",
                        "d_cigar",
                        "j_cigar",
                        "v_alignment_length",
                        "v_alignment_mutation",
                    ]
                    if self.chain_type not in [
                        ChainType.HUMAN_HEAVY,
                        ChainType.MOUSE_HEAVY,
                        ChainType.RABBIT_HEAVY,
                        ChainType.HUMAN_BETA,
                        ChainType.HUMAN_DELTA,
                        ChainType.MOUSE_BETA,
                        ChainType.MOUSE_DELTA,
                    ]:
                        w_header.remove("d_call")
                        w_header.remove("d_cigar")

                    if self.chain_type in [
                        ChainType.CHICKEN_HEAVY,
                        ChainType.CHICKEN_LIGHT,
                    ]:
                        w_header = [
                            "sequence_id",
                            "duplicate_count",
                            "frequency",
                            "sequence",
                            "sequence_aa",
                            "v_call",
                            "j_call",
                            "cdr1",
                            "cdr2",
                            "cdr3",
                            "cdr1_aa",
                            "cdr2_aa",
                            "cdr3_aa",
                            "rev_comp",
                            "productive",
                            "v_alignment_length",
                            "v_alignment_mutation",
                        ]

                    with open(out_file, "w") as handle:
                        write_file(
                            handle,
                            w_header,
                            list(sorted_functional_dict.values()),
                            delim=self.delim,
                        )

                    # RAM usage update
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info("4. Functional reads extraction was done.")
                    logging.info(
                        "--- %s seconds elapsed ---" % (datetime.now() - start_time)
                    )

                logging.info(
                    "Finished Functionality check for %s", self.sample_name_to_put
                )
                logging.info(
                    "--- %s seconds elapsed ---" % (datetime.now() - start_time_main)
                )

            elif self.step == None:
                print('Need to set "step" parameter... (int type: 1, 2, 3, 4)')
                return ValueError

        logging.info("Step %d for all files are finished" % self.step)
        logging.info("Max RAM used : %.2f GB" % (max_RAM_used - self.initial_RAM_used))

    ### work functions for the ProcessWorker class
    def get_file_lines(self, target):
        count = 0
        with open(target, "r") as handle:
            for line in handle:
                count += 1
        return count

    def get_key_len(self, target_dict):
        count = len(target_dict)
        return count

    def extract_chunk_start_pos(self, target_len, chunk_size):
        chunk_start_list = []

        i = 0
        while True:
            if chunk_size * i > target_len:
                break

            line_start = chunk_size * i
            i += 1

            chunk_start_list.append(line_start)

        logging.info(
            "Size of divided chunks for target %s: %s"
            % (self.sample_name_to_put, chunk_size)
        )

        return chunk_start_list

    # merge chunk data according to the key
    # key ==
    #       1) trimmed sequence after primer recognition (primer recognized)
    #       2) lines in the original file (isotyping, annotation, functionality).
    def merge_chunk_dict_type(self, in_result_queue, rc_col_name="readcount"):
        # dict form
        # {key1: {}, key2: {}, ...}

        result_dict = {}
        while in_result_queue.qsize() > 0:
            result_chunk = in_result_queue.get()
            for k in result_chunk:
                if k not in result_dict:
                    result_dict[k] = result_chunk[k]
                else:
                    result_dict[k][rc_col_name] += result_chunk[k][rc_col_name]

        return result_dict


class ErrorCorrecter(Process):
    """Do error correction (UMI processing)

    1a. primer recognition : file -> dict
    1b. applying MiXCR: dict -> file (no MP)
    """

    def __init__(
        self,
        sub_step,
        out_queue,
        log_queue,
        RAM_used_queue,
        in_queue="",
        chunk_size="",
        input_file="",
        input_dict={},
        _dir=None,
        fwd_primer_list=None,
        rev_primer_list=None,
        evalue=0.001,
        rev_primer_trimming=True,
        chain_type=None,
    ):
        Process.__init__(self)

        self.input_file = input_file
        self.sub_step = sub_step

        self.out_chunk = {}

        self.out_queue = out_queue
        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.in_queue = in_queue
        self.chunk_size = chunk_size

        self.input_dict = input_dict
        self.dir = _dir
        self.fwd_primer_list = fwd_primer_list
        self.rev_primer_list = rev_primer_list

        self.evalue = evalue
        self.rev_primer_trimming = rev_primer_trimming

        self.chain_type = chain_type

    def run(self):
        if self.log_queue:
            logger = logging.getLogger()
            if len(logger.handlers) == 0:
                h = QueueHandler(self.log_queue)
                logger.addHandler(h)
                logger.setLevel(logging.DEBUG)

        if self.sub_step == "a":
            while True:
                # Get the work from the queue and expand the tuple
                chunk_start_pos = self.in_queue.get()
                if chunk_start_pos is None:
                    self.in_queue.task_done()
                    break

                self.primer_recognition(chunk_start_pos, self.input_file)

                tmp_RAM_used = psutil.virtual_memory().used / 1024**3
                self.RAM_used_queue.put(tmp_RAM_used)

                self.in_queue.task_done()

                # save chunk data into result queue
                self.out_queue.put(self.out_chunk)

        elif self.sub_step == "b":
            self.do_mixcr()

            tmp_RAM_used = psutil.virtual_memory().used / 1024**3
            self.RAM_used_queue.put(tmp_RAM_used)

            # save chunk data into result queue
            self.out_queue.put(self.out_chunk)

    # worker functions
    def primer_recognition(self, line_start, input_file_name):
        sHelper = Helper()

        # set the direction of primers into sense strand (plus/plus strand) (fwd, rev)
        forwardPrimerList = self.fwd_primer_list
        reversePrimerList = [
            sHelper.to_reverse_complement(reversePrimer)
            for reversePrimer in self.rev_primer_list
        ]

        # add reverse primer sequences containing 1 mismatch (regular expression pattern)
        reversePrimerMatchList = []
        for reversePrimer in reversePrimerList:
            reversePrimerMatchList.append(reversePrimer)
        for reversePrimer in reversePrimerList:
            for i in range(len(reversePrimer)):
                reversePrimerMatchList.append(
                    reversePrimer[:i] + reversePrimer[i + 1 :]
                )
                reversePrimerMatchList.append(
                    reversePrimer[:i] + "[ATGC]" + reversePrimer[i + 1 :]
                )

        # add forward primer sequences containing 1 mismatch
        forwardPrimerMatchList = []
        for forwardPrimer in forwardPrimerList:
            forwardPrimerMatchList.append(forwardPrimer)
        for forwardPrimer in forwardPrimerList:
            for i in range(len(forwardPrimer)):
                forwardPrimerMatchList.append(
                    forwardPrimer[:i] + forwardPrimer[i + 1 :]
                )
                forwardPrimerMatchList.append(
                    forwardPrimer[:i] + "[ATGC]" + forwardPrimer[i + 1 :]
                )

        # set the forward / reverse region to find pattern
        fwd_primer_max_len = max([len(fwd) for fwd in forwardPrimerList])
        rev_primer_max_len = max([len(rev) for rev in reversePrimerList])
        forwardTrimNum = fwd_primer_max_len + 5
        if self.chain_type == ChainType.CHICKEN_HEAVY:
            forwardTrimNum = fwd_primer_max_len + 50

        ### special case: dealing with human delta chain cases -> 11mer sequence was added
        if self.chain_type == ChainType.CHICKEN_HEAVY:
            reverseTrimNum = rev_primer_max_len + 50
        else:
            reverseTrimNum = rev_primer_max_len + 14 + 5

        recogSet = {}

        with open(os.path.join(input_file_name), "r") as f:
            reader = csv.reader(f)
            header = next(reader)

            if len(header) != 2:
                print("Input file column error")
                raise RuntimeError

            # go to start position
            for i in range(0, line_start):
                next(reader)

            line_count = 0
            accessed_line_num = line_start
            while True:
                try:
                    row = next(reader)
                except StopIteration:
                    break

                line_count += 1
                if line_count > self.chunk_size:
                    break

                sequence = row[0]
                readcount = int(row[1])

                success = False

                for reversePrimerMatch in reversePrimerMatchList:
                    reverseMatchResult = re.search(
                        reversePrimerMatch, sequence[-reverseTrimNum:]
                    )

                    if reverseMatchResult != None:
                        for forwardPrimerMatch in forwardPrimerMatchList:
                            forwardMatchResult = re.search(
                                forwardPrimerMatch, sequence[:forwardTrimNum]
                            )

                            if forwardMatchResult != None:
                                if reverseTrimNum - reverseMatchResult.span()[1] == 0:
                                    if self.rev_primer_trimming:
                                        newSequence = sequence[
                                            forwardMatchResult.span()[1] : -(
                                                reverseTrimNum
                                                - reverseMatchResult.span()[0]
                                            )
                                        ]
                                    else:
                                        newSequence = sequence[
                                            forwardMatchResult.span()[1] :
                                        ]
                                else:
                                    if self.rev_primer_trimming:
                                        newSequence = sequence[
                                            forwardMatchResult.span()[1] : -(
                                                reverseTrimNum
                                                - reverseMatchResult.span()[0]
                                            )
                                        ]  # erase primer regions
                                    else:
                                        newSequence = sequence[
                                            forwardMatchResult.span()[1] : -(
                                                reverseTrimNum
                                                - reverseMatchResult.span()[1]
                                            )
                                        ]  # erase fwd primer region
                                success = True

                                if newSequence not in recogSet:
                                    recogSet[newSequence] = {
                                        "full_NT": newSequence,
                                        "readcount": readcount,
                                    }
                                else:
                                    recogSet[newSequence]["readcount"] += readcount
                                accessed_line_num += 1
                                break
                            else:
                                continue
                        break

                # try primer recognition process for reverse complement case.
                if success == False:
                    sequence = sHelper.to_reverse_complement(sequence)

                    for reversePrimerMatch in reversePrimerMatchList:
                        reverseMatchResult = re.search(
                            reversePrimerMatch, sequence[-reverseTrimNum:]
                        )

                        if reverseMatchResult != None:
                            for forwardPrimerMatch in forwardPrimerMatchList:
                                forwardMatchResult = re.search(
                                    forwardPrimerMatch, sequence[:forwardTrimNum]
                                )

                                if forwardMatchResult != None:
                                    if (
                                        reverseTrimNum - reverseMatchResult.span()[1]
                                        == 0
                                    ):
                                        if self.rev_primer_trimming:
                                            newSequence = sequence[
                                                forwardMatchResult.span()[1] : -(
                                                    reverseTrimNum
                                                    - reverseMatchResult.span()[0]
                                                )
                                            ]
                                        else:
                                            newSequence = sequence[
                                                forwardMatchResult.span()[1] :
                                            ]
                                    else:
                                        if self.rev_primer_trimming:
                                            newSequence = sequence[
                                                forwardMatchResult.span()[1] : -(
                                                    reverseTrimNum
                                                    - reverseMatchResult.span()[0]
                                                )
                                            ]  # erase primer regions
                                        else:
                                            newSequence = sequence[
                                                forwardMatchResult.span()[1] : -(
                                                    reverseTrimNum
                                                    - reverseMatchResult.span()[1]
                                                )
                                            ]  # erase fwd primer region
                                    success = True
                                    if newSequence not in recogSet:
                                        recogSet[newSequence] = {
                                            "full_NT": newSequence,
                                            "readcount": readcount,
                                        }
                                    else:
                                        recogSet[newSequence]["readcount"] += readcount
                                    accessed_line_num += 1
                                    break
                                else:
                                    continue
                            break

        self.out_chunk = recogSet

    def do_mixcr(self, col_seq="full_NT"):
        data = list(self.input_dict.values())
        sample = "readcount"

        # Manipulating singleton-removed data
        logging.info("<<MiXCR>> Singleton removal started")
        sr_data = []
        for d in data:
            readcount = int(d[sample])
            if readcount > 1:
                sr_data.append(d)
        logging.info("---- %d -> %d clones" % (len(data), len(sr_data)))
        logging.info("<<MiXCR>> Singleton removal finished")

        logging.info("<<MiXCR>> Correcting errors started")
        sub_map = {"A": "GCTN", "G": "ACTN", "C": "AGTN", "T": "AGCN"}

        # generate dict
        map_seq_data = {}
        for d in sr_data:
            map_seq_data[d[col_seq]] = d
        remove_clones = set()
        logging.info("---- original %d clones" % len(sr_data))

        col_reads = sample
        num_remove = 0
        for d in sr_data:
            seq = d[col_seq]
            reads1 = d[col_reads]
            if reads1 * self.evalue < 2:
                pass
            else:
                # substitution
                for i in range(len(seq)):
                    for mut_base in sub_map[seq[i]]:
                        mut_seq = seq[:i] + mut_base + seq[i + 1 :]
                        if mut_seq in map_seq_data:
                            reads2 = map_seq_data[mut_seq][col_reads]
                            if reads1 * self.evalue > reads2:
                                # remove_clones.add(mut_seq)
                                map_seq_data[mut_seq][col_reads] = 0
                                num_remove += 1
        logging.info("---- read count correct %d clones" % num_remove)

        # remove corrected clones
        clone_filter_out = 0
        # for j in range(1):
        #     clone_filter_out.append([1])
        for seq in map_seq_data:
            # reads_list = []
            # reads_list += [map_seq_data[seq][sample]]
            # if reads_list in clone_filter_out:
            #     remove_clones.add(seq)
            reads = map_seq_data[seq][sample]
            if reads == clone_filter_out:
                remove_clones.add(seq)
        logging.info("---- remove zeros %d clones" % len(remove_clones))
        data_out = [d for d in sr_data if d[col_seq] not in remove_clones]

        logging.info("<<MiXCR>> Correct error finished")

        self.out_chunk = data_out  # data_out = [{}, {}, {}, ...]

    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue


class IsotypeExtractor(Process):
    """Do c gene alignment and isotype annotation

    2a. run blast : file -> file
    2b. parse blast run result and get isotype info : file -> dict
    """

    def __init__(
        self,
        input_file,
        sub_step,
        log_queue,
        RAM_used_queue,
        out_dir,
        c_ref_file=None,
        in_queue=None,
        out_queue=None,
        num_process=1,
        chunk_size=None,
        delim="\t",
        rev_primer_trimming=True,
    ):
        Process.__init__(self)

        self.input_file = input_file
        self.sub_step = sub_step

        self.out_chunk = {}

        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.out_dir = out_dir
        self.c_ref_file = c_ref_file
        self.in_queue = in_queue
        self.out_queue = out_queue

        self.num_process = num_process
        self.chunk_size = chunk_size

        self.delim = delim

        self.rev_primer_trimming = rev_primer_trimming

    def run(self):
        if self.log_queue:
            logger = logging.getLogger()
            if len(logger.handlers) == 0:
                h = QueueHandler(self.log_queue)
                logger.addHandler(h)
                logger.setLevel(logging.DEBUG)

        if self.sub_step == "a":
            self.run_blast(self.input_file)

            tmp_RAM_used = psutil.virtual_memory().used / 1024**3
            self.RAM_used_queue.put(tmp_RAM_used)

        elif self.sub_step == "b":
            ## load data into RAM: error-corrected files for loading sequence & readcount info
            targetHeader, targetData = read_file_fast(self.input_file, delim=self.delim)

            while True:
                # Get the work from the queue and expand the tuple
                chunk_start_pos = self.in_queue.get()
                if chunk_start_pos is None:
                    self.in_queue.task_done()
                    break

                # extract isotype hash table
                self.parse_isotype(chunk_start_pos, targetData)
                self.in_queue.task_done()

                # RAM usage check
                tmp_RAM_used = psutil.virtual_memory().used / 1024**3
                self.RAM_used_queue.put(tmp_RAM_used)

                # save chunk data into result queue
                self.out_queue.put(self.out_chunk)

    # worker functions
    def run_blast(self, input_file, **kwargs):
        chRefHeder, chRefList = read_file_fast(in_file=self.c_ref_file, delim=",")
        targetSetHeader, targetSeqList = read_file_fast(input_file, delim=self.delim)

        # directory setting
        if not self.out_dir:
            output_dir = "/".join(
                os.path.dirname(input_file).split("/")[:-1] + ["2_annotation"]
            )
        else:
            output_dir = self.out_dir
        blast_dir = output_dir + "temp_cblast"
        query_dir = os.path.join(blast_dir, "query")
        db_dir = os.path.join(blast_dir, "db")
        result_dir = os.path.join(blast_dir, "result")

        # make directories
        t_dirs = [blast_dir, query_dir, db_dir, result_dir]
        for t_dir in t_dirs:
            if not os.path.exists(t_dir):
                os.mkdir(t_dir)

        # set parameters
        if True:
            if "evalue" in kwargs:
                evalue = kwargs["evalue"]
            else:
                if self.rev_primer_trimming:
                    evalue = 0.1**5
                else:
                    evalue = 0.1**1

            if "gapopen" in kwargs:
                gapopen = kwargs["gapopen"]
            else:
                gapopen = 1

            if "gapextend" in kwargs:
                gapextend = kwargs["gapextend"]
            else:
                gapextend = 2

            if "word_size" in kwargs:
                word_size = kwargs["word_size"]
            else:
                word_size = 10

            if "num_alignments" in kwargs:
                num_alignments = kwargs["num_alignments"]
            else:
                num_alignments = 1

            query_list = [x["full_NT"] for x in targetSeqList]
            db_dict = {x[chRefHeder[0]]: x[chRefHeder[1]] for x in chRefList}
            result_file_name = "blasted_MiXCR.txt"

        result_file = os.path.join(result_dir, result_file_name)

        with open(os.path.join(query_dir, "query_MiXCR.fasta"), "w") as fasta_writer:
            for i, e_seq in enumerate(query_list):
                fasta_writer.write(">" + str(i) + "\n")
                fasta_writer.write(e_seq + "\n")

        with open(os.path.join(db_dir, "db_MiXCR.fasta"), "w") as fasta_writer:
            for name, e_seq in db_dict.items():
                fasta_writer.write(">" + name + "\n")
                fasta_writer.write(e_seq + "\n")

        # set blast_exe_path according to the server (host name)
        if hostname == "binel229":
            exe_path = "/home/team/IP-team/Tools/ncbi-blast-2.11.0+/bin"
        elif hostname == "JunhoLab":
            exe_path = "/Tools/ncbi-blast-2.7.1+/bin"

        # construct db
        format_cmd = "%s -in %s -dbtype nucl -input_type fasta -out %s" % (
            os.path.join(exe_path, "makeblastdb"),
            os.path.join(db_dir, "db_MiXCR.fasta"),
            os.path.join(db_dir, "db_MiXCR"),
        )
        os.system(format_cmd)

        # run blastn
        cline = NcbiblastnCommandline(
            num_threads=self.num_process,
            query=os.path.join(query_dir, "query_MiXCR.fasta"),
            db=os.path.join(db_dir, "db_MiXCR"),
            evalue=evalue,
            out=result_file,
            gapopen=gapopen,
            gapextend=gapextend,
            word_size=word_size,
            num_descriptions=num_alignments,
            num_alignments=num_alignments,
        )
        format_exe = os.path.join(exe_path, "blastn")
        os.system(format_exe + str(cline)[len("blastn") :])

    # make user defined blast parser for the time efficiency -> do not use NCBIStandalone module
    def parse_isotype(self, key_start, input_data_list):
        result_dict = {}
        blast_result_file = os.path.join(self.out_dir, "blasted_MiXCR.txt")

        with open(blast_result_file, "r") as handle:
            key_start_count = 0
            accessed_key_count = 0

            while True:
                info_accessed = 0
                line = handle.readline()
                if line == "":
                    break
                if line[: len("Query=")] == "Query=":
                    # go to start position
                    if key_start_count < key_start:
                        key_start_count += 1
                        continue

                    # stop working if the number of accessed key exceed the chunk_size
                    accessed_key_count += 1
                    if accessed_key_count > self.chunk_size:
                        break

                    access_num = int(line.split()[-1].strip())
                    seq = input_data_list[access_num]["full_NT"]
                    rc = input_data_list[access_num]["readcount"]

                    # get query length
                    handle.readline()
                    line = handle.readline()
                    if line[: len("Length=")] == "Length=":
                        query_len = int(line.split("=")[-1].strip())

                    query_start = 0
                    query_end = 0
                    sbjct_start = 0

                    # user defined blast parser (specific purpose: extract required info for isotyping)
                    while True:
                        line = handle.readline()
                        if "No hits found" in line:
                            break
                        elif line[: len(">")] == ">":
                            # sbjct_name = line.split()[-1].strip()[1:]
                            sbjct_name = line[1:].strip()
                            info_accessed += 1
                        elif line.strip()[: len("Identities")] == "Identities":
                            identities = (
                                line.split(",")[0].split("=")[-1].strip().split()[0]
                            )
                            matched, whole = identities.split("/")
                            info_accessed += 2
                        elif line.strip()[: len("Strand")] == "Strand":
                            direction = line.split("=")[-1].strip()
                            if direction != "Plus/Plus":
                                break
                            info_accessed += 1
                        elif line[: len("Query")] == "Query":
                            query_end = int(line.split()[-1].strip())
                            if query_start > 0:
                                continue
                            else:
                                query_start = int(line.split()[1].strip())
                                info_accessed += 1
                        elif line[: len("Sbjct")] == "Sbjct":
                            if sbjct_start > 0:
                                continue
                            sbjct_start = int(line.split()[1].strip())
                            info_accessed += 1

                        # end of alignment
                        elif (
                            line[: len("Effective search space used:")]
                            == "Effective search space used:"
                        ):
                            if query_end > 0:
                                info_accessed += 1
                            break

                    if info_accessed < 7:
                        continue

                    # pass not well aligned cases
                    numMismatches = int(whole) - int(matched)
                    # if int(whole) < 45 or int(whole) > 71 or int(numMismatches) > 5 or int(sbjct_start) > 10 or int(query_start) > 400:
                    #     continue

                    # if int(numMismatches) > 5 or int(sbjct_start) > 10 or (int(query_len) - int(query_end)) > 10:
                    #     continue

                    if int(numMismatches) > 5 or int(sbjct_start) > 10:
                        continue

                    # get isotype information and save
                    vdj_seq = seq[: query_start - sbjct_start]
                    # isotype = sbjct_name.split('*')[0][3:]
                    # isotype = sbjct_name
                    isotype = sbjct_name
                    if "TR" in isotype:
                        isotype = sbjct_name.split("_")[0]
                    new_key = "|".join([vdj_seq, isotype])
                    try:
                        result_dict[new_key]["readcount"] += rc
                    except KeyError:
                        result_dict[new_key] = {
                            "full_NT": vdj_seq,
                            "isotype": isotype,
                            "readcount": rc,
                        }

            self.out_chunk = result_dict

    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue


class InfoAnnotator(Process):
    """Do IgBLAST run and region annotation

    3a. run igblast : file -> file
    3b. parse igblast run result and get region info : file -> dict
    """

    def __init__(
        self,
        input_file,
        sub_step,
        log_queue,
        RAM_used_queue,
        out_dir,
        chain_type,
        in_queue=None,
        out_queue=None,
        num_process=1,
        chunk_size=None,
        delim="\t",
        add_keys={
            "sequence": "full_NT",
            "duplicate_count": "readcount",
            "c_call": "isotype",
        },
    ):
        Process.__init__(self)

        self.input_file = input_file
        self.sub_step = sub_step

        self.out_chunk = {}

        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.out_dir = out_dir
        self.chain_type = chain_type
        self.in_queue = in_queue
        self.out_queue = out_queue

        self.num_process = num_process
        self.chunk_size = chunk_size

        self.delim = delim

        self.add_keys = add_keys

    def run(self):
        if self.log_queue:
            logger = logging.getLogger()
            if len(logger.handlers) == 0:
                h = QueueHandler(self.log_queue)
                logger.addHandler(h)
                logger.setLevel(logging.DEBUG)

        if self.sub_step == "a":
            self.run_igblast()

            tmp_RAM_used = psutil.virtual_memory().used / 1024**3
            self.RAM_used_queue.put(tmp_RAM_used)

        elif self.sub_step == "b":
            ## load data into RAM: error-corrected files for loading sequence & readcount info
            targetHeader, targetData = read_file_fast(self.input_file, delim=self.delim)

            # # some parameters...
            # self.domain_system = 'imgt'

            while True:
                # Get the work from the queue and expand the tuple
                chunk_start_pos = self.in_queue.get()
                if chunk_start_pos is None:
                    self.in_queue.task_done()
                    break

                # extract isotype hash table
                self.annotate_info(chunk_start_pos, targetData)
                self.in_queue.task_done()

                # RAM usage check
                tmp_RAM_used = psutil.virtual_memory().used / 1024**3
                self.RAM_used_queue.put(tmp_RAM_used)

                # save chunk data into result queue
                self.out_queue.put(self.out_chunk)

    # worker functions → not use Helper class in SeqUtil.py. separate run_igblast and igblast_parser
    def run_igblast(self):
        ### parameter setting
        targetHeader, targetList = read_file_fast(self.input_file, delim=self.delim)
        sequence_list = [d["full_NT"] for d in targetList]

        blast_dir = os.path.join(os.path.join(self.out_dir, "temp_igblast"))
        if not os.path.exists(blast_dir):
            os.mkdir(blast_dir)

        igblast_name = os.path.basename(self.input_file)[: -len(".tsv")]
        file_blast = os.path.join(blast_dir, igblast_name + "_igblasted.txt")

        if self.chain_type in [
            ChainType.HUMAN_HEAVY,
            ChainType.HUMAN_LIGHT,
            ChainType.CHICKEN_HEAVY,
            ChainType.CHICKEN_LIGHT,
            ChainType.RAT_HEAVY,
            ChainType.RAT_LIGHT,
            ChainType.RABBIT_HEAVY,
            ChainType.RABBIT_KAPPA,
            ChainType.MOUSE_BALBC_HEAVY,
            ChainType.MOUSE_BALBC_LIGHT,
            ChainType.MOUSE_C57BL6_HEAVY,
            ChainType.MOUSE_C57BL6_LIGHT,
            ChainType.MOUSE_HEAVY,
            ChainType.MOUSE_LIGHT,
            ChainType.ALPACA_HEAVY,
        ]:
            seq_type = "Ig"
            # domain_system = 'kabat'
            domain_system = "imgt"
        elif self.chain_type in [
            ChainType.HUMAN_BETA,
            ChainType.HUMAN_ALPHA,
            ChainType.HUMAN_DELTA,
            ChainType.HUMAN_GAMMA,
            ChainType.MOUSE_BETA,
            ChainType.MOUSE_ALPHA,
            ChainType.MOUSE_GAMMA,
            ChainType.MOUSE_DELTA,
        ]:
            seq_type = "TCR"
            domain_system = "imgt"
        else:
            raise RuntimeError("chain type error")

        # check sequence type
        is_nucl = is_nucleotide_or_protein(sequence_list[0])
        if is_nucl == None or is_nucl == False:
            raise RuntimeError(
                "parameter sequence_list is not a nucleotide sequence list"
            )

        # make temporal query fasta file
        tp_query = open(os.path.join(blast_dir, "tmp.fasta"), "w")
        for i, seq in enumerate(sequence_list):
            tp_query.write(">%d\n" % i)
            if i == len(sequence_list) - 1:
                tp_query.write(seq)
            else:
                tp_query.write(seq + "\n")
        tp_query.close()

        # run igblast
        run_igblast_new(
            chain_type=self.chain_type,
            query=tp_query.name,
            out=file_blast,
            seq_type=seq_type,
            domain_system=domain_system,
            num_threads=self.num_process,
        )

        # erase tmp fasta file
        os.system("rm %s" % tp_query.name)

    def annotate_info(self, key_start, input_data_list):
        sHelper = Helper()
        igblast_name = os.path.basename(self.input_file)[: -len(".tsv")]
        igblast_result_file = os.path.join(
            self.out_dir, "%s_igblasted.txt" % igblast_name
        )

        ### parse igblast result txt file.
        igblast_parsed_dict = {}
        # igblast_parsed_list = []
        with open(igblast_result_file, "r") as igb_handle:
            key_start_count = 0
            accessed_key_count = 0

            # include BCR & TCR
            if self.chain_type in [
                ChainType.HUMAN_HEAVY,
                ChainType.HUMAN_LIGHT,
                ChainType.RAT_HEAVY,
                ChainType.RAT_LIGHT,
                ChainType.RABBIT_HEAVY,
                ChainType.RABBIT_KAPPA,
                ChainType.MOUSE_BALBC_HEAVY,
                ChainType.MOUSE_BALBC_LIGHT,
                ChainType.MOUSE_C57BL6_HEAVY,
                ChainType.MOUSE_C57BL6_LIGHT,
                ChainType.ALPACA_HEAVY,
                ChainType.MOUSE_HEAVY,
                ChainType.MOUSE_LIGHT,
                ChainType.HUMAN_BETA,
                ChainType.HUMAN_ALPHA,
                ChainType.HUMAN_DELTA,
                ChainType.HUMAN_GAMMA,
                ChainType.MOUSE_BETA,
                ChainType.MOUSE_ALPHA,
                ChainType.MOUSE_DELTA,
                ChainType.MOUSE_GAMMA,
            ]:
                # for the case that domain system is "imgt"
                while True:
                    line = igb_handle.readline()
                    if not line:
                        break

                    if line[: len("Query= ")] == "Query= ":
                        # go to start position
                        if key_start_count < key_start:
                            key_start_count += 1
                            continue

                        # stop working if the number of accessed key exceed the chunk_size
                        accessed_key_count += 1
                        if accessed_key_count > self.chunk_size:
                            break

                        row_num = int(line[len("Query= ") :].strip())

                        junction = ""
                        junction_aa = ""
                        junction_start = 0
                        junction_end = 0
                        CDR3_NT = ""
                        sequence_alignment = ""
                        seq_aa = ""

                        result = {
                            "query": line[len("Query= ") : -1],
                            "hit": False,
                            "v_call": "N/A",
                            "j_call": "N/A",
                            "cdr3": "N/A",
                            "cdr3_aa": "X",
                            "junction": "N/A",
                            "junction_aa": "X",
                            "rev_comp": "F",
                            "productive": "F",
                            "v_cigar": "",
                            "j_cigar": "",
                            "v_alignment_length": "N/A",
                            "v_alignment_mutation": "N/A",
                            "sequence_alignment": "N/A",
                            "germline_alignment": "N/A",
                        }
                        if self.chain_type in [
                            ChainType.HUMAN_HEAVY,
                            ChainType.RABBIT_HEAVY,
                            ChainType.MOUSE_C57BL6_HEAVY,
                            ChainType.MOUSE_BALBC_HEAVY,
                            ChainType.MOUSE_HEAVY,
                            ChainType.HUMAN_BETA,
                            ChainType.HUMAN_DELTA,
                            ChainType.MOUSE_BETA,
                            ChainType.MOUSE_DELTA,
                        ]:
                            result["d_call"] = "N/A"
                            result["d_cigar"] = ""

                        while True:
                            line = igb_handle.readline()

                            if "No hits found" in line:
                                result["hit"] = False
                                igblast_parsed_dict[row_num] = result
                                break

                            elif (
                                "Note that your query represents the minus strand"
                                in line
                            ):
                                result["rev_comp"] = "T"

                            elif (
                                line[: len("Sub-region sequence details")]
                                == "Sub-region sequence details"
                            ):
                                next_line_sp = igb_handle.readline().split("\t")
                                if next_line_sp[0] == "CDR3":
                                    result["hit"] = True
                                    CDR3_NT = next_line_sp[1]
                                    result["cdr3"] = CDR3_NT
                                    result["cdr3_aa"] = sHelper.translate(
                                        result["cdr3"]
                                    )
                                    result["cdr3_start"] = int(next_line_sp[3])
                                    junction_start = int(next_line_sp[3]) - 4
                                    if junction_start != 0:
                                        result["junction_start"] = junction_start
                                    result["cdr3_end"] = int(next_line_sp[4])
                                    junction_end = int(next_line_sp[4]) + 3
                                    if junction_end != 0:
                                        result["junction_end"] = junction_end
                                    result["fr4_from"] = int(next_line_sp[4])

                            elif line[: len("Alignments")] == "Alignments":
                                next_line_sp = igb_handle.readline()
                                next_line_sp = igb_handle.readline()
                                next_line_sp = igb_handle.readline().split(" ")
                                next_line_trimmed = []
                                for l in range(len(next_line_sp)):
                                    if next_line_sp[l] != "":
                                        next_line_trimmed.append(next_line_sp[l])
                                query_line = []
                                query_start_where = []
                                query_finish_where = []
                                V_germline_mismatch = []
                                D_germline_mismatch = []
                                J_germline_mismatch = []
                                V_mis_ratio = []
                                D_mis_ratio = []
                                J_mis_ratio = []
                                v_germ_start_where = []
                                v_germ_finish_where = []
                                d_germ_start_where = []
                                d_germ_finish_where = []
                                j_germ_start_where = []
                                j_germ_finish_where = []
                                sequence_a = []
                                dog = True
                                while dog:
                                    next_line_trimmed = []
                                    before_line_sp = next_line_sp
                                    next_line_sp = igb_handle.readline().split(" ")
                                    for l in range(len(next_line_sp)):
                                        if next_line_sp[l] != "":
                                            next_line_trimmed.append(next_line_sp[l])

                                    before_line_trimmed = []
                                    for l in range(len(before_line_sp)):
                                        if before_line_sp[l] != "":
                                            before_line_trimmed.append(
                                                before_line_sp[l]
                                            )

                                    findQu = re.compile("Query_")
                                    findQue = findQu.search(next_line_trimmed[0])

                                    findper = re.compile("%")
                                    if len(next_line_trimmed) >= 2:
                                        findperc = findper.search(next_line_trimmed[1])

                                    if findQue:
                                        query_start_where.append(next_line_trimmed[1])
                                        query_line.append(next_line_trimmed[2])
                                        query_finish_where.append(
                                            next_line_trimmed[3].replace("\n", "")
                                        )
                                        seq_a = ""
                                        for l in range(len(before_line_trimmed)):
                                            seq_a = seq_a + before_line_trimmed[l]
                                        sequence_a.append(seq_a.replace("\n", ""))

                                    elif next_line_trimmed[0] == "V" and findperc:
                                        if len(next_line_trimmed) >= 7:
                                            if len(V_mis_ratio) == 0:
                                                V_mis_ratio.append(
                                                    next_line_trimmed[1:3]
                                                )
                                            v_germ_start_where.append(
                                                next_line_trimmed[4]
                                            )
                                            V_germline_mismatch.append(
                                                next_line_trimmed[5]
                                            )
                                            v_germ_finish_where.append(
                                                next_line_trimmed[6].replace("\n", "")
                                            )

                                    elif next_line_trimmed[0] == "D" and findperc:
                                        if len(next_line_trimmed) >= 7:
                                            if len(D_mis_ratio) == 0:
                                                D_mis_ratio.append(
                                                    next_line_trimmed[1:3]
                                                )
                                            d_germ_start_where.append(
                                                next_line_trimmed[4]
                                            )
                                            D_germline_mismatch.append(
                                                next_line_trimmed[5]
                                            )
                                            d_germ_finish_where.append(
                                                next_line_trimmed[6].replace("\n", "")
                                            )

                                    elif next_line_trimmed[0] == "J" and findperc:
                                        if len(next_line_trimmed) >= 7:
                                            if len(J_mis_ratio) == 0:
                                                J_mis_ratio.append(
                                                    next_line_trimmed[1:3]
                                                )
                                            j_germ_start_where.append(
                                                next_line_trimmed[4]
                                            )
                                            J_germline_mismatch.append(
                                                next_line_trimmed[5]
                                            )
                                            j_germ_finish_where.append(
                                                next_line_trimmed[6].replace("\n", "")
                                            )

                                    elif next_line_trimmed[0] == "Lambda":
                                        dog = False
                                seq_aa = ""
                                for k in range(len(sequence_a)):
                                    seq_aa = seq_aa + sequence_a[k]
                                result["sequence_aa"] = seq_aa
                                result["query_line"] = query_line
                                for a in range(len(query_line)):
                                    sequence_alignment = (
                                        sequence_alignment + query_line[a]
                                    )
                                result["sequence_alignment"] = sequence_alignment
                                seq_without = ""
                                if len(sequence_alignment) != 0:
                                    for z in range(len(sequence_alignment)):
                                        if sequence_alignment[z] != "-":
                                            seq_without = (
                                                seq_without + sequence_alignment[z]
                                            )
                                if len(seq_without) != 0 and junction_start != 0:
                                    junction = seq_without[
                                        junction_start
                                        - int(query_start_where[0])
                                        + 1 : junction_end
                                        - int(query_start_where[0])
                                        + 1
                                    ]
                                if CDR3_NT != "":
                                    result["junction"] = junction
                                    result["junction_aa"] = sHelper.translate(junction)
                                else:
                                    result["junction"] = ""
                                    result["junction_aa"] = ""
                                v_finish_from_zero = 0
                                if len(v_germ_start_where) != 0:
                                    v_finish_from_zero = int(
                                        v_germ_finish_where[
                                            len(v_germ_finish_where) - 1
                                        ]
                                    ) - int(v_germ_start_where[0])
                                v_cigar_mismatch = []
                                v_cigar_mismatches = []
                                j_cigar_mismatches = []
                                d_cigar_mismatches = []
                                if (
                                    len(sequence_alignment) >= v_finish_from_zero + 1
                                ) and v_finish_from_zero != 0:
                                    for a in range(v_finish_from_zero + 1):
                                        if sequence_alignment[a] == "-":
                                            v_cigar_mismatch.append([a, "0"])
                                if len(v_germ_start_where) != 0:
                                    for a in range(len(v_germ_start_where)):
                                        for i in range(
                                            int(v_germ_finish_where[a])
                                            - int(v_germ_start_where[a])
                                            + 1
                                        ):
                                            if V_germline_mismatch[a][i] != ".":
                                                v_cigar_mismatch.append(
                                                    [
                                                        i
                                                        + int(v_germ_start_where[a])
                                                        - int(v_germ_start_where[0]),
                                                        V_germline_mismatch[a][i],
                                                    ]
                                                )
                                v_cigar_mismatch.sort(key=lambda x: x[0])
                                for kk in range(len(v_cigar_mismatch)):
                                    isit = 0
                                    for ll in range(len(v_cigar_mismatches)):
                                        if (
                                            v_cigar_mismatch[kk][0]
                                            == v_cigar_mismatches[ll][0]
                                        ):
                                            isit += 1
                                    if isit == 0:
                                        v_cigar_mismatches.append(v_cigar_mismatch[kk])

                                d_finish_from_zero = 0
                                if len(d_germ_start_where) != 0:
                                    d_finish_from_zero = int(
                                        d_germ_finish_where[
                                            len(d_germ_finish_where) - 1
                                        ]
                                    ) - int(d_germ_start_where[0])

                                d_cigar_mismatch = []
                                d_line_pos = []
                                if (len(d_germ_start_where)) != 0:
                                    for k in range(90):
                                        if D_germline_mismatch[0][k] != "-":
                                            d_line_po = k
                                            d_line_pos.append(k)
                                            break

                                seq_before_d_num = 0
                                if (
                                    d_finish_from_zero != 0
                                    and len(v_germ_finish_where) != 0
                                    and len(sequence_alignment) != 0
                                ):
                                    for q in range(len(v_germ_start_where) - 1):
                                        # seq_before_d_num : 90, 180, ...
                                        seq_before_d_num = seq_before_d_num + 90
                                    for a in range(len(d_cigar_mismatch)):
                                        if (
                                            sequence_alignment[
                                                seq_before_d_num
                                                + d_line_pos[0]
                                                + 90 * a
                                            ]
                                            == "-"
                                        ):
                                            d_cigar_mismatch.append(
                                                [
                                                    seq_before_d_num
                                                    + d_line_pos[0]
                                                    + 90 * a,
                                                    "0",
                                                ]
                                            )
                                if len(d_germ_start_where) == 1:
                                    for i in range(
                                        int(d_germ_finish_where[0])
                                        - int(d_germ_start_where[0])
                                    ):
                                        if (
                                            D_germline_mismatch[0][d_line_pos[0] + i]
                                            != "."
                                        ):
                                            d_cigar_mismatch.append(
                                                [
                                                    i
                                                    + seq_before_d_num
                                                    + d_line_pos[0],
                                                    D_germline_mismatch[0][
                                                        d_line_pos[0] + i
                                                    ],
                                                ]
                                            )
                                elif len(d_germ_start_where) >= 2:
                                    for a in range(90 - d_line_pos[0]):
                                        if (
                                            D_germline_mismatch[0][d_line_pos[0] + a]
                                            != "."
                                        ):
                                            d_cigar_mismatch.append(
                                                [
                                                    seq_before_d_num
                                                    + d_line_pos[0]
                                                    + a,
                                                    D_germline_mismatch[0][
                                                        d_line_pos[0] + a
                                                    ],
                                                ]
                                            )
                                    for a in range(
                                        int(d_germ_finish_where[1])
                                        - int(d_germ_start_where[0])
                                        + 1
                                        - 90
                                        + d_line_pos[0]
                                    ):
                                        if D_germline_mismatch[1][a] != ".":
                                            d_cigar_mismatch.append(
                                                [
                                                    seq_before_d_num + 90 + a,
                                                    D_germline_mismatch[1][a],
                                                ]
                                            )
                                d_cigar_mismatch.sort(key=lambda x: x[0])
                                for kk in range(len(d_cigar_mismatch)):
                                    isit = 0
                                    for ll in range(len(d_cigar_mismatches)):
                                        if (
                                            d_cigar_mismatch[kk][0]
                                            == d_cigar_mismatches[ll][0]
                                        ):
                                            isit += 1
                                    if isit == 0:
                                        d_cigar_mismatches.append(d_cigar_mismatch[kk])

                                j_finish_from_zero = 0
                                if len(j_germ_start_where) != 0:
                                    j_finish_from_zero = int(
                                        j_germ_finish_where[
                                            len(j_germ_finish_where) - 1
                                        ]
                                    ) - int(j_germ_start_where[0])

                                j_cigar_mismatch = []
                                j_line_pos = []
                                if (len(j_germ_start_where)) != 0:
                                    for k in range(90):
                                        if J_germline_mismatch[0][k] != "-":
                                            j_line_po = k
                                            j_line_pos.append(j_line_po)
                                            break

                                seq_before_j_num = 0
                                if (
                                    j_finish_from_zero != 0
                                    and len(v_germ_finish_where) != 0
                                ):
                                    if (len(d_germ_finish_where)) <= 1:
                                        for q in range(len(v_germ_start_where) - 1):
                                            # seq_before_j_num : 90, 180, ...
                                            seq_before_j_num = seq_before_j_num + 90
                                    elif (len(d_germ_finish_where)) == 2:
                                        for q in range(len(v_germ_start_where)):
                                            # seq_before_j_num : 90, 180, ...
                                            seq_before_j_num = seq_before_j_num + 90
                                    for a in range(len(j_cigar_mismatch)):
                                        if (
                                            sequence_alignment[
                                                seq_before_j_num
                                                + j_line_pos[0]
                                                + 90 * a
                                            ]
                                            == "-"
                                        ):
                                            j_cigar_mismatch.append(
                                                [
                                                    seq_before_j_num
                                                    + j_line_pos[0]
                                                    + 90 * a,
                                                    "0",
                                                ]
                                            )
                                if len(j_germ_start_where) == 1:
                                    for i in range(
                                        int(j_germ_finish_where[0])
                                        - int(j_germ_start_where[0])
                                    ):
                                        if (
                                            J_germline_mismatch[0][j_line_pos[0] + i]
                                            != "."
                                        ):
                                            j_cigar_mismatch.append(
                                                [
                                                    i
                                                    + seq_before_j_num
                                                    + j_line_pos[0],
                                                    J_germline_mismatch[0][
                                                        j_line_pos[0] + i
                                                    ],
                                                ]
                                            )
                                elif len(j_germ_start_where) >= 2:
                                    for a in range(90 - j_line_pos[0]):
                                        if (
                                            J_germline_mismatch[0][j_line_pos[0] + a]
                                            != "."
                                        ):
                                            j_cigar_mismatch.append(
                                                [
                                                    seq_before_j_num
                                                    + j_line_pos[0]
                                                    + a,
                                                    J_germline_mismatch[0][
                                                        j_line_pos[0] + a
                                                    ],
                                                ]
                                            )
                                    for pp in range(len(j_germ_start_where) - 1):
                                        for a in range(
                                            int(j_germ_finish_where[pp + 1])
                                            - int(j_germ_start_where[pp + 1])
                                            + 1
                                        ):
                                            if J_germline_mismatch[pp + 1][a] != ".":
                                                j_cigar_mismatch.append(
                                                    [
                                                        seq_before_j_num
                                                        + 90 * (pp + 1)
                                                        + a,
                                                        J_germline_mismatch[1][a],
                                                    ]
                                                )
                                j_cigar_mismatch.sort(key=lambda x: x[0])

                                for kk in range(len(j_cigar_mismatch)):
                                    isit = 0
                                    for ll in range(len(j_cigar_mismatches)):
                                        if (
                                            j_cigar_mismatch[kk][0]
                                            == j_cigar_mismatches[ll][0]
                                        ):
                                            isit += 1
                                    if isit == 0:
                                        j_cigar_mismatches.append(j_cigar_mismatch[kk])

                                germline_alignment = sequence_alignment
                                germline_alignment_one = ""
                                germline_alignment_two = ""
                                germline_alignment_three = ""
                                if (
                                    len(j_cigar_mismatches) != 0
                                    and len(germline_alignment) != 0
                                ):
                                    for p in range(len(j_cigar_mismatches)):
                                        germline_alignment = (
                                            germline_alignment[
                                                : j_cigar_mismatches[p][0]
                                            ]
                                            + j_cigar_mismatches[p][1]
                                            + germline_alignment[
                                                j_cigar_mismatches[p][0] + 1 :
                                            ]
                                        )
                                if (
                                    len(d_cigar_mismatches) != 0
                                    and len(germline_alignment) != 0
                                ):
                                    for p in range(len(d_cigar_mismatches)):
                                        germline_alignment = (
                                            germline_alignment[
                                                : d_cigar_mismatches[p][0]
                                            ]
                                            + d_cigar_mismatches[p][1]
                                            + germline_alignment[
                                                d_cigar_mismatches[p][0] + 1 :
                                            ]
                                        )
                                if (
                                    len(j_cigar_mismatch) != 0
                                    and len(germline_alignment) != 0
                                ):
                                    for p in range(len(j_cigar_mismatch)):
                                        germline_alignment = (
                                            germline_alignment[: j_cigar_mismatch[p][0]]
                                            + j_cigar_mismatch[p][1]
                                            + germline_alignment[
                                                j_cigar_mismatch[p][0] + 1 :
                                            ]
                                        )
                                the_string = ""
                                second_string = ""
                                if (
                                    (len(d_germ_start_where)) == 0
                                    and (len(j_germ_finish_where)) != 0
                                    and len(v_germ_finish_where) != 0
                                ):
                                    seq_before_j_num = 0
                                    if (
                                        j_finish_from_zero != 0
                                        and len(v_germ_finish_where) != 0
                                    ):
                                        if (len(d_germ_finish_where)) <= 1:
                                            for q in range(len(v_germ_start_where) - 1):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = seq_before_j_num + 90
                                        elif (len(d_germ_finish_where)) == 2:
                                            for q in range(len(v_germ_start_where)):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = seq_before_j_num + 90
                                    the_string = ""
                                    for u in range(
                                        seq_before_j_num
                                        + j_line_pos[0]
                                        - int(
                                            v_germ_finish_where[
                                                (len(v_germ_finish_where) - 1)
                                            ]
                                        )
                                        + int(v_germ_start_where[0])
                                    ):
                                        the_string = the_string + "-"
                                    if len(germline_alignment) != 0:
                                        germline_alignment = (
                                            germline_alignment[
                                                : int(
                                                    v_germ_finish_where[
                                                        len(v_germ_finish_where) - 1
                                                    ]
                                                )
                                                - int(v_germ_start_where[0])
                                            ]
                                            + the_string
                                            + germline_alignment[
                                                seq_before_j_num + j_line_pos[0] :
                                            ]
                                        )
                                elif (
                                    (len(d_germ_start_where)) != 0
                                    and len(j_germ_finish_where) != 0
                                    and len(v_germ_finish_where) != 0
                                ):
                                    seq_before_d_num = 0
                                    for q in range(len(v_germ_start_where) - 1):
                                        # seq_before_d_num : 90, 180, ...
                                        seq_before_d_num = (
                                            seq_before_d_num
                                            + int(v_germ_finish_where[0])
                                            - int(v_germ_start_where[0])
                                            + 1
                                        )

                                    seq_before_j_num = 0
                                    if (
                                        j_finish_from_zero != 0
                                        and len(v_germ_finish_where) != 0
                                    ):
                                        if (len(d_germ_finish_where)) <= 1:
                                            for q in range(len(v_germ_start_where) - 1):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = (
                                                    seq_before_j_num
                                                    + int(v_germ_finish_where[0])
                                                    - int(v_germ_start_where[0])
                                                    + 1
                                                )
                                        elif (len(d_germ_finish_where)) == 2:
                                            for q in range(len(v_germ_start_where)):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = (
                                                    seq_before_j_num
                                                    + int(v_germ_finish_where[0])
                                                    - int(v_germ_start_where[0])
                                                    + 1
                                                )

                                    the_string = ""
                                    second_string = ""
                                    for u in range(
                                        seq_before_d_num
                                        + d_line_pos[0]
                                        - int(
                                            v_germ_finish_where[
                                                (len(v_germ_finish_where) - 1)
                                            ]
                                        )
                                        + int(v_germ_start_where[0])
                                        - 1
                                    ):
                                        the_string = the_string + "-"
                                    for u in range(
                                        seq_before_j_num
                                        + j_line_pos[0]
                                        - seq_before_d_num
                                        - d_line_pos[0]
                                        - int(
                                            d_germ_finish_where[
                                                len(d_germ_finish_where) - 1
                                            ]
                                        )
                                        + int(d_germ_start_where[0])
                                        - 1
                                    ):
                                        second_string = second_string + "-"
                                    germline_alignment_one = germline_alignment[
                                        : int(
                                            v_germ_finish_where[
                                                len(v_germ_finish_where) - 1
                                            ]
                                        )
                                        - int(v_germ_start_where[0])
                                        + 1
                                    ]
                                    germline_alignment_two = germline_alignment[
                                        seq_before_d_num
                                        + d_line_pos[0] : seq_before_d_num
                                        + d_line_pos[0]
                                        + int(
                                            d_germ_finish_where[
                                                len(d_germ_finish_where) - 1
                                            ]
                                        )
                                        - int(d_germ_start_where[0])
                                        + 1
                                    ]
                                    germline_alignment_three = germline_alignment[
                                        seq_before_j_num + j_line_pos[0] :
                                    ]
                                    germline_alignment = (
                                        germline_alignment_one
                                        + the_string
                                        + germline_alignment_two
                                        + second_string
                                        + germline_alignment_three
                                    )

                                result["the_string"] = the_string
                                result["second_string"] = second_string
                                result["germline_alignment"] = germline_alignment
                                result["germline_alignment_1"] = germline_alignment_one
                                result["germline_alignment_2"] = germline_alignment_two
                                result[
                                    "germline_alignment_3"
                                ] = germline_alignment_three
                                result["j_line_pos"] = j_line_pos

                                v_cig = ""
                                cur_pos_is = -1
                                type_mis = 0
                                v_cig_num = 0
                                v_cigar = ""
                                if len(v_cigar_mismatches) >= 1:
                                    cur_pos_is = -1
                                    for e in range(len(v_cigar_mismatches)):
                                        if v_cigar_mismatches[e][1] == "0":
                                            v_cig_num = (
                                                int(v_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if v_cig_num == -1:
                                                v_cig = v_cig
                                            elif v_cig_num != 0:
                                                v_cig = (
                                                    v_cig
                                                    + str(v_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "D"
                                                )
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 1:
                                                v_cig = (
                                                    v_cig[: len(v_cig) - 2]
                                                    + str(
                                                        int(v_cig[len(v_cig) - 2]) + 1
                                                    )
                                                    + "D"
                                                )
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            else:
                                                v_cig = v_cig + str(1) + "D"
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            type_mis = 1
                                        elif v_cigar_mismatches[e][1] == "-":
                                            v_cig_num = (
                                                int(v_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if v_cig_num == -1:
                                                v_cig = v_cig
                                            elif v_cig_num != 0:
                                                v_cig = (
                                                    v_cig
                                                    + str(v_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "I"
                                                )
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 2:
                                                v_cig = (
                                                    v_cig[: len(v_cig) - 2]
                                                    + str(
                                                        int(v_cig[len(v_cig) - 2]) + 1
                                                    )
                                                    + "I"
                                                )
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            else:
                                                v_cig = v_cig + str(1) + "I"
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            type_mis = 2
                                        else:
                                            v_cig_num = (
                                                int(v_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if v_cig_num == -1:
                                                v_cig = v_cig
                                            elif v_cig_num != 0:
                                                v_cig = (
                                                    v_cig
                                                    + str(v_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "X"
                                                )
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 3:
                                                v_cig = (
                                                    v_cig[: len(v_cig) - 2]
                                                    + str(
                                                        int(v_cig[len(v_cig) - 2]) + 1
                                                    )
                                                    + "X"
                                                )
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            else:
                                                v_cig = v_cig + str(1) + "X"
                                                cur_pos_is = int(
                                                    v_cigar_mismatches[e][0]
                                                )
                                            type_mis = 3
                                if len(v_germ_finish_where) != 0:
                                    if (
                                        cur_pos_is
                                        != int(
                                            v_germ_finish_where[
                                                (len(v_germ_finish_where) - 1)
                                            ]
                                        )
                                        - int(v_germ_start_where[0])
                                        + 1
                                    ):
                                        v_cig_num = (
                                            int(
                                                v_germ_finish_where[
                                                    (len(v_germ_finish_where) - 1)
                                                ]
                                            )
                                            - int(v_germ_start_where[0])
                                            - cur_pos_is
                                        )
                                        v_cig = v_cig + str(v_cig_num) + "="
                                v_cigar = v_cig

                                d_cig = ""
                                type_mis = 0
                                d_cig_num = 0
                                d_cigar = ""
                                if (
                                    len(d_germ_finish_where) != 0
                                    and len(v_germ_finish_where) != 0
                                ):
                                    cur_pos_is = seq_before_d_num + d_line_pos[0] - 1
                                    for e in range(len(d_cigar_mismatches)):
                                        if d_cigar_mismatches[e][1] == "0":
                                            d_cig_num = (
                                                int(d_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if d_cig_num == -1:
                                                d_cig = d_cig
                                            elif d_cig_num != 0:
                                                d_cig = (
                                                    d_cig
                                                    + str(d_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "D"
                                                )
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 1:
                                                d_cig = (
                                                    d_cig[: len(d_cig) - 2]
                                                    + str(
                                                        int(d_cig[len(d_cig) - 2]) + 1
                                                    )
                                                    + "D"
                                                )
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            else:
                                                d_cig = d_cig + str(1) + "D"
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            type_mis = 1
                                        elif d_cigar_mismatches[e][1] == "-":
                                            d_cig_num = (
                                                int(d_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if d_cig_num == -1:
                                                d_cig = d_cig
                                            elif d_cig_num != 0:
                                                d_cig = (
                                                    d_cig
                                                    + str(d_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "I"
                                                )
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 2:
                                                d_cig = (
                                                    d_cig[: len(d_cig) - 2]
                                                    + str(
                                                        int(d_cig[len(d_cig) - 2]) + 1
                                                    )
                                                    + "I"
                                                )
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            else:
                                                d_cig = d_cig + str(1) + "I"
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            type_mis = 2
                                        else:
                                            d_cig_num = (
                                                int(d_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if d_cig_num == -1:
                                                d_cig = d_cig
                                            elif d_cig_num != 0:
                                                d_cig = (
                                                    d_cig
                                                    + str(d_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "X"
                                                )
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 3:
                                                d_cig = (
                                                    d_cig[: len(d_cig) - 2]
                                                    + str(
                                                        int(d_cig[len(d_cig) - 2]) + 1
                                                    )
                                                    + "X"
                                                )
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            else:
                                                d_cig = d_cig + str(1) + "X"
                                                cur_pos_is = int(
                                                    d_cigar_mismatches[e][0]
                                                )
                                            type_mis = 3
                                    if (
                                        cur_pos_is
                                        != int(
                                            d_germ_finish_where[
                                                (len(d_germ_finish_where) - 1)
                                            ]
                                        )
                                        - int(d_germ_start_where[0])
                                        + seq_before_d_num
                                        + d_line_pos[0]
                                        + 1
                                    ):
                                        d_cig_num = (
                                            int(
                                                d_germ_finish_where[
                                                    (len(d_germ_finish_where) - 1)
                                                ]
                                            )
                                            - int(d_germ_start_where[0])
                                            - cur_pos_is
                                            + seq_before_d_num
                                            + d_line_pos[0]
                                        )
                                        d_cig = d_cig + str(d_cig_num) + "="
                                    d_cigar = d_cig
                                elif len(d_germ_finish_where) == 0:
                                    d_cigar = ""

                                j_cig = ""
                                type_mis = 0
                                j_cig_num = 0
                                j_cigar = ""
                                if (len(j_line_pos)) != 0:
                                    cur_pos_is = seq_before_j_num + j_line_pos[0] - 1
                                    for e in range(len(j_cigar_mismatches)):
                                        if j_cigar_mismatches[e][1] == "0":
                                            j_cig_num = (
                                                int(j_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if j_cig_num == -1:
                                                j_cig = j_cig
                                            elif j_cig_num != 0:
                                                j_cig = (
                                                    j_cig
                                                    + str(j_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "D"
                                                )
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 1:
                                                j_cig = (
                                                    j_cig[: len(j_cig) - 2]
                                                    + str(
                                                        int(j_cig[len(j_cig) - 2]) + 1
                                                    )
                                                    + "D"
                                                )
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            else:
                                                j_cig = j_cig + str(1) + "D"
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            type_mis = 1
                                        elif j_cigar_mismatches[e][1] == "-":
                                            j_cig_num = (
                                                int(j_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if j_cig_num == -1:
                                                j_cig = j_cig
                                            elif j_cig_num != 0:
                                                j_cig = (
                                                    j_cig
                                                    + str(j_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "I"
                                                )
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 2:
                                                j_cig = (
                                                    j_cig[: len(j_cig) - 2]
                                                    + str(
                                                        int(j_cig[len(j_cig) - 2]) + 1
                                                    )
                                                    + "I"
                                                )
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            else:
                                                j_cig = j_cig + str(1) + "I"
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            type_mis = 2
                                        else:
                                            j_cig_num = (
                                                int(j_cigar_mismatches[e][0])
                                                - cur_pos_is
                                                - 1
                                            )
                                            if j_cig_num == -1:
                                                j_cig = j_cig
                                            elif j_cig_num != 0:
                                                j_cig = (
                                                    j_cig
                                                    + str(j_cig_num)
                                                    + "="
                                                    + str(1)
                                                    + "X"
                                                )
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            elif type_mis == 3:
                                                j_cig = (
                                                    j_cig[: len(j_cig) - 2]
                                                    + str(
                                                        int(j_cig[len(j_cig) - 2]) + 1
                                                    )
                                                    + "X"
                                                )
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            else:
                                                j_cig = j_cig + str(1) + "X"
                                                cur_pos_is = int(
                                                    j_cigar_mismatches[e][0]
                                                )
                                            type_mis = 3
                                    if (
                                        cur_pos_is
                                        != int(
                                            j_germ_finish_where[
                                                (len(j_germ_finish_where) - 1)
                                            ]
                                        )
                                        - int(j_germ_start_where[0])
                                        + seq_before_j_num
                                        + j_line_pos[0]
                                        + 1
                                    ):
                                        j_cig_num = (
                                            int(
                                                j_germ_finish_where[
                                                    (len(j_germ_finish_where) - 1)
                                                ]
                                            )
                                            - int(j_germ_start_where[0])
                                            - cur_pos_is
                                            + seq_before_j_num
                                            + j_line_pos[0]
                                        )
                                        j_cig = j_cig + str(j_cig_num) + "="
                                    j_cigar = j_cig

                                result["v_cigar"] = v_cigar
                                result["d_cigar"] = d_cigar
                                result["j_cigar"] = j_cigar
                                result["V_mis_ratio"] = V_mis_ratio
                                result["D_mis_ratio"] = D_mis_ratio
                                result["J_mis_ratio"] = J_mis_ratio
                                result["V_germline_mismatch"] = V_germline_mismatch
                                result["D_germline_mismatch"] = D_germline_mismatch
                                result["J_germline_mismatch"] = J_germline_mismatch
                                result["v_germ_start_where"] = v_germ_start_where
                                result["v_germ_finish_where"] = v_germ_finish_where
                                result["d_germ_start_where"] = d_germ_start_where
                                result["d_germ_finish_where"] = d_germ_finish_where
                                result["j_germ_start_where"] = j_germ_start_where
                                result["j_germ_finish_where"] = j_germ_finish_where
                                result["v_cigar_mismatches"] = v_cigar_mismatches
                                result["d_cigar_mismatches"] = d_cigar_mismatches
                                result["j_cigar_mismatches"] = j_cigar_mismatches

                            elif (
                                line[: len("Alignment summary")] == "Alignment summary"
                            ):
                                totalDistanceFromVGene = 0
                                totalAlignedLength = 0

                                while True:
                                    next_line_sp = igb_handle.readline().split("\t")
                                    region = (
                                        next_line_sp[0].replace("-IMGT", "").lower()
                                    )
                                    if region == "total":
                                        break
                                    if region == "FR1":
                                        fr1_start = int(next_line_sp[1]) - 1
                                        fr1_end = int(next_line_sp[2])
                                        if sHelper.RepresentsInt(next_line_sp[6]):
                                            fr1_gap = int(next_line_sp[6])
                                        else:
                                            fr1_gap = 0
                                        while (fr1_end - fr1_start + fr1_gap) % 3 != 0:
                                            fr1_start += 1
                                        result["fr1_from"] = fr1_start
                                        result["fr1_to"] = fr1_end
                                        result["fr1_length"] = int(next_line_sp[3])
                                        result["fr1_matches"] = int(next_line_sp[4])
                                        result["fr1_mismatches"] = int(next_line_sp[5])
                                        result["fr1_gaps"] = int(next_line_sp[6])
                                        totalAlignedLength += result["fr1_length"]
                                        totalDistanceFromVGene += (
                                            result["fr1_mismatches"]
                                            + result["fr1_gaps"]
                                        )
                                    else:
                                        try:
                                            result[region + "_from"] = (
                                                int(next_line_sp[1]) - 1
                                            )
                                            result[region + "_to"] = int(
                                                next_line_sp[2]
                                            )
                                            result[region + "_length"] = int(
                                                next_line_sp[3]
                                            )
                                            result[region + "_matches"] = int(
                                                next_line_sp[4]
                                            )
                                            result[region + "_mismatches"] = int(
                                                next_line_sp[5]
                                            )
                                            result[region + "_gaps"] = int(
                                                next_line_sp[6]
                                            )
                                            totalAlignedLength += result[
                                                region + "_length"
                                            ]
                                            totalDistanceFromVGene += (
                                                result[region + "_mismatches"]
                                                + result[region + "_gaps"]
                                            )
                                        except ValueError:
                                            continue

                                result["v_alignment_length"] = totalAlignedLength
                                result["v_alignment_mutation"] = totalDistanceFromVGene

                            elif (
                                line[: len("V-(D)-J rearrangement summary")]
                                == "V-(D)-J rearrangement summary"
                            ):
                                next_line_sp = igb_handle.readline().split("\t")
                                if self.chain_type in [
                                    ChainType.HUMAN_HEAVY,
                                    ChainType.RABBIT_HEAVY,
                                    ChainType.MOUSE_C57BL6_HEAVY,
                                    ChainType.MOUSE_BALBC_HEAVY,
                                    ChainType.MOUSE_HEAVY,
                                    ChainType.HUMAN_BETA,
                                    ChainType.HUMAN_DELTA,
                                    ChainType.MOUSE_BETA,
                                    ChainType.MOUSE_DELTA,
                                ]:
                                    result["v_call"] = next_line_sp[0]
                                    result["d_call"] = next_line_sp[1]
                                    result["j_call"] = next_line_sp[2]
                                    result["vj_frame"] = next_line_sp[5].split("-")[0]
                                    result["productive"] = next_line_sp[6]

                                elif self.chain_type in [
                                    ChainType.HUMAN_LIGHT,
                                    ChainType.RABBIT_KAPPA,
                                    ChainType.MOUSE_C57BL6_LIGHT,
                                    ChainType.MOUSE_BALBC_LIGHT,
                                    ChainType.MOUSE_LIGHT,
                                    ChainType.HUMAN_ALPHA,
                                    ChainType.HUMAN_GAMMA,
                                    ChainType.MOUSE_ALPHA,
                                    ChainType.MOUSE_GAMMA,
                                ]:
                                    result["v_call"] = next_line_sp[0]
                                    result["j_call"] = next_line_sp[1]
                                    result["vj_frame"] = next_line_sp[4].split("-")[0]
                                    result["productive"] = next_line_sp[5]

                            # End of one alignment result
                            elif (
                                line[: len("Effective search space used:")]
                                == "Effective search space used:"
                            ):
                                igblast_parsed_dict[row_num] = result
                                # igblast_parsed_list.append(result)
                                break

        ### combine annotation information with input_data_list (isotyped input file)
        # get range_from and range_to
        range_from = key_start
        range_to = max(igblast_parsed_dict.keys())

        # extract target sequence list
        target_list = input_data_list[range_from : (range_to + 1)]

        # get output dictionary
        result_dict = {}
        w_row_num = key_start
        for target, igblast_parsed_data in zip(
            target_list, igblast_parsed_dict.values()
        ):
            tmp_dict = {k: target[v] for k, v in self.add_keys.items()}
            # tmp_dict = {'sequence': target['full_NT'], 'duplicate_count': target['readcount'], 'c_call': target['isotype']}
            igblast_parsed_data.update(tmp_dict)

            ## check not annotated cases and fill it with default info (N/A or blank)
            if "sequence_aa" not in igblast_parsed_data:
                igblast_parsed_data["sequence_aa"] = ""
            if igblast_parsed_data["productive"] == "Yes":
                igblast_parsed_data["productive"] = "T"
            else:
                igblast_parsed_data["productive"] = "F"

            # add cdr1, 2
            sequence = target["full_NT"]
            if igblast_parsed_data["rev_comp"] == "T":
                sequence = sHelper.to_reverse_complement(sequence)
            if "fr4_from" in igblast_parsed_data:
                igblast_parsed_data["fr4"] = sequence[
                    igblast_parsed_data["fr4_from"] : igblast_parsed_data["fr4_from"]
                    + 33
                ]
            for r in ["fr1", "fr2", "fr3", "cdr1", "cdr2"]:
                if ("%s_from" % r) in igblast_parsed_data and (
                    "%s_to" % r
                ) in igblast_parsed_data:
                    igblast_parsed_data[r] = sequence[
                        igblast_parsed_data["%s_from" % r] : igblast_parsed_data[
                            "%s_to" % r
                        ]
                    ]
                    igblast_parsed_data[r + "_aa"] = sHelper.translate(
                        igblast_parsed_data[r]
                    )
                else:
                    igblast_parsed_data[r] = "N/A"
                    igblast_parsed_data[r + "_aa"] = "X"

            result_dict[w_row_num] = igblast_parsed_data
            w_row_num += 1

        self.out_chunk = result_dict

    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue


class functionalExtracter(Process):
    """Do extraction of functional reads"""

    def __init__(
        self,
        input_file,
        log_queue,
        RAM_used_queue,
        in_queue=None,
        out_queue=None,
        chunk_size=None,
        delim="\t",
    ):
        Process.__init__(self)

        self.input_file = input_file

        self.out_chunk = {}

        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.in_queue = in_queue
        self.out_queue = out_queue

        self.chunk_size = chunk_size

        self.delim = delim

    def run(self):
        if self.log_queue:
            logger = logging.getLogger()
            if len(logger.handlers) == 0:
                h = QueueHandler(self.log_queue)
                logger.addHandler(h)
                logger.setLevel(logging.DEBUG)

        # logging.info("Started run")
        while True:
            # Get the work from the queue and expand the tuple
            chunk_start_pos = self.in_queue.get()
            if chunk_start_pos is None:
                self.in_queue.task_done()
                break

            # extract functional reads
            self.extract_functional(chunk_start_pos)
            self.in_queue.task_done()

            # RAM usage check
            tmp_RAM_used = psutil.virtual_memory().used / 1024**3
            self.RAM_used_queue.put(tmp_RAM_used)

            # save chunk data into result queue
            self.out_queue.put(self.out_chunk)

        # logging.info("Finished run")

    # work function
    def extract_functional(self, key_start):
        # get output dictionary
        result_dict = {}
        w_row_num = key_start

        with open(self.input_file, "r") as igb_handle:
            key_start_count = 0
            accessed_key_count = 0

            header = (
                igb_handle.readline().strip().split(self.delim)
            )  # extract header line
            col_seq_aa = header.index("sequence_aa")
            col_v_call = header.index("v_call")
            col_j_call = header.index("j_call")

            col_cdr1_aa = header.index("cdr1_aa")
            col_cdr2_aa = header.index("cdr2_aa")
            col_cdr3_aa = header.index("cdr3_aa")

            col_rev_comp = header.index("rev_comp")
            col_productive = header.index("productive")

            while True:
                line = igb_handle.readline()
                if not line:
                    break

                # go to start position
                if key_start_count < key_start:
                    key_start_count += 1
                    continue

                # stop working if the number of accessed key exceed the chunk_size
                accessed_key_count += 1
                if accessed_key_count > self.chunk_size:
                    break

                data_list = line.strip().split(self.delim)
                seq_aa = data_list[col_seq_aa]
                v_call = data_list[col_v_call]
                j_call = data_list[col_j_call]

                cdr1_aa = data_list[col_cdr1_aa]
                cdr2_aa = data_list[col_cdr2_aa]
                cdr3_aa = data_list[col_cdr3_aa]

                rev_comp = data_list[col_rev_comp]
                productive = data_list[col_productive]

                # full aa stop codon filtering
                try:
                    seq_aa.index("*")
                    continue
                except ValueError:
                    pass

                # full aa frame filtering
                try:
                    seq_aa.index("X")
                    continue
                except ValueError:
                    pass

                # gene annotation filtering
                if v_call == "N/A" or j_call == "N/A":
                    continue

                # cdr1,2,3 annotation flitering
                if cdr1_aa == "N/A" or cdr2_aa == "N/A" or cdr3_aa == "N/A":
                    continue

                # cdr1,2,3 stop codon filtering
                try:
                    cdr1_aa.index("*")
                    continue
                except ValueError:
                    pass
                try:
                    cdr2_aa.index("*")
                    continue
                except ValueError:
                    pass
                try:
                    cdr3_aa.index("*")
                    continue
                except ValueError:
                    pass

                # cdr1,2,3 frame filtering
                try:
                    cdr1_aa.index("X")
                    continue
                except ValueError:
                    pass
                try:
                    cdr2_aa.index("X")
                    continue
                except ValueError:
                    pass
                try:
                    cdr3_aa.index("X")
                    continue
                except ValueError:
                    pass

                # rev_comp filtering: must be 'F'
                if rev_comp == "T":
                    continue

                # productive filtering: must be 'T'
                if productive == "F":
                    continue

                result_dict[w_row_num] = {k: v for k, v in zip(header, data_list)}
                w_row_num += 1

        self.out_chunk = result_dict

    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue
