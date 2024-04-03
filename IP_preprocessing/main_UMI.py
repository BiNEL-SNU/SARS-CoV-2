#####################################################################################################
# Author
# Date
# Editor    YHLee
# Revised   2022-12-05
# Note      Preprocessing code from raw csv to functionality csv
#           Including UMI processing step
#           Target> samples prepared for NovaSeq
#           human BCR, TCR
#           mouse BCR, TCR ...
#
#####################################################################################################

import os
import logging
from datetime import datetime

from utils.Enum import ChainType
from utils.PreprocessWorker_for_NovaSeq import PreprocessWorker
from utils.FileUtil import read_file_fast
from utils.CommandParseUtil import get_opts, get_argv, check_work_dir, check_threads

BASE_NAME = os.path.basename(__file__)


# Set ChainType by organism and chain_type
def define_ChainType(organism, chaintype):
    if organism == "human":
        if chaintype == "BCRH":
            return ChainType.HUMAN_HEAVY
        elif chaintype in ["BCRK", "BCRL"]:
            return ChainType.HUMAN_LIGHT
        elif chaintype == "TRB":
            return ChainType.HUMAN_BETA
        elif chaintype == "TRA":
            return ChainType.HUMAN_ALPHA
        elif chaintype == "TRD":
            return ChainType.HUMAN_DELTA
        elif chaintype == "TRG":
            return ChainType.HUMAN_GAMMA
    elif organism == "mouse":
        if chaintype == "BCRH":
            return ChainType.MOUSE_HEAVY
        elif chaintype in ["BCRK", "BCRL"]:
            return ChainType.MOUSE_LIGHT
        elif chaintype == "TRB":
            return ChainType.MOUSE_BETA
        elif chaintype == "TRA":
            return ChainType.MOUSE_ALPHA
        elif chaintype == "TRD":
            return ChainType.MOUSE_DELTA
        elif chaintype == "TRG":
            return ChainType.MOUSE_GAMMA


# main
def main():
    # Get argument values from the command, by each option
    opts = get_opts(py_name=BASE_NAME)
    argv_by_opt = get_argv(parsed_opts=opts, py_name=BASE_NAME)
    check_work_dir(working_dir=argv_by_opt["working_directory"], py_name=BASE_NAME)
    argv_by_opt["num_threads"] = check_threads(argv_by_opt["num_threads"])

    # Define the rev_trimming
    if argv_by_opt["organism"] == "human" and argv_by_opt["chain_type"] == "BCRH":
        rev_trimming = True
    else:
        rev_trimming = False

    # Define chainType (Enum)
    chain_type = define_ChainType(
        organism=argv_by_opt["organism"],
        chaintype=argv_by_opt["chain_type"],
    )

    # Define file debugging mode
    if argv_by_opt["debug"] == "True":
        file_debug = True
    elif argv_by_opt["debug"] == "False":
        file_debug = False
    else:
        file_debug = True

    # temporary for logging
    logging.basicConfig(
        format="[%(asctime)s] %(processName)s %(funcName)s\t[%(levelname)s]\t%(message)s",
        level=logging.INFO,
    )
    logging.info("Script started")

    # Set list of forward & reverse primer sequences
    fwd_file = argv_by_opt["fwd_primer_file"]
    rev_file = argv_by_opt["rev_primer_file"]
    header, data = read_file_fast(fwd_file, delim=",")
    forwardPrimerList = [d["sequence"] for d in data]
    header, data = read_file_fast(rev_file, delim=",")
    reversePrimerList = [d["sequence"] for d in data]

    # Set working directory and sample list
    working_dir = argv_by_opt["working_directory"]
    sampleList = os.listdir(working_dir)
    logging.info("The number of target samples: %d" % len(sampleList))

    # Check the script start time
    start_time = datetime.now()

    ### do work in multiprocessing manner.
    """
    All functions needs following parameters:
        - working_dir : base directory where all sample folders exists
        - sampleList : list of sample name
        - step : 1, 2, 3, or 4 according to what you want to do
            1 = error correction using UMI processing
            2 = isotyping
            3 = annotation
            4 = functionality check
        - file_type : csv or tsv. default=tsv
        - infile_suffix / outfile_suffix : according to the name of files, which is used or constructed at the process
            e[num] : error correction
            a[num] : isotyping & annotation
            f[num] : functionality check
    """
    # 1. Run error correction = primer recognition -> UMI merge -> sub-clustering -> length filtering -> consensus sequence extraction
    # mandatory input: forward primer list, reverse primer list, step=1
    if not argv_by_opt["step"] or argv_by_opt["step"] == "1":
        # for sub-clustering
        UMI_subcluster_dist_thresh = 5
        # for length filtering
        len_filter_hamming_thresh = 1
        # for consensus sequence extraction
        consensus_thresh = 0.6

        error_correcter = PreprocessWorker(
            work_dir=working_dir,
            sample_name_list=sampleList,
            fwd_primers=forwardPrimerList,
            rev_primers=reversePrimerList,
            dist_thresh=UMI_subcluster_dist_thresh,
            hamming_thresh=len_filter_hamming_thresh,
            vote_thresh=consensus_thresh,
            step=1,
            file_type="tsv",
            outfile_suffix="e1",
            debug=file_debug,
            rev_primer_trimming=rev_trimming,
        )
        error_correcter.run(argv_by_opt["num_threads"])

    # 2. c gene annotation (Isotyping)
    # mandatory input: c gene reference file
    if not argv_by_opt["step"] or argv_by_opt["step"] == "2":
        isotype_finder = PreprocessWorker(
            work_dir=working_dir,
            sample_name_list=sampleList,
            c_ref_file=argv_by_opt["c_gene_file"],
            step=2,
            file_type="tsv",
            infile_suffix="e1",
            outfile_suffix="e1_a0",
            rev_primer_trimming=rev_trimming,
        )
        isotype_finder.run(argv_by_opt["num_threads"])

    # 3. Annotation using IgBLAST -> need to written in user function, for simplicity.
    # mandatory input: chain type
    if not argv_by_opt["step"] or argv_by_opt["step"] == "3":
        Annotator = PreprocessWorker(
            work_dir=working_dir,
            sample_name_list=sampleList,
            chain_type=chain_type,
            step=3,
            file_type="tsv",
            infile_suffix="e1_a0",
            outfile_suffix="e1_a1",
        )
        Annotator.run(argv_by_opt["num_threads"])

    # 4. Functional reads filtration
    # mandatory input: chain_type
    if not argv_by_opt["step"] or argv_by_opt["step"] == "4":
        functional_extracter = PreprocessWorker(
            work_dir=working_dir,
            sample_name_list=sampleList,
            chain_type=chain_type,
            step=4,
            file_type="tsv",
            infile_suffix="e1_a1",
            outfile_suffix="e1_a1_f1",
        )
        functional_extracter.run(argv_by_opt["num_threads"])

    logging.info("--- %s seconds elapsed ---" % (datetime.now() - start_time))
    logging.info("Script finished")


# run main function
if __name__ == "__main__":
    main()
