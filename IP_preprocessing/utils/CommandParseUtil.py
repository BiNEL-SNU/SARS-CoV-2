#####################################################################################################
# Author    Yonghee Lee
# Date      2022-12-05
# Editor
# Revised
# Note      Functions related to command-line parsing
#
#####################################################################################################

import os
import logging
import multiprocessing
from datetime import datetime
import sys, getopt

MIN_THREAD = 5
MAX_THREAD = int(multiprocessing.cpu_count() * 0.8)
OPT_DICT = {
    "d": "working_directory",
    "s": "step",
    "j": "num_threads",
    "c": "chain_type",
    "organism": "organism",
    "fwd_primer": "fwd_primer_file",
    "rev_primer": "rev_primer_file",
    "c_gene_file": "c_gene_file",
    "debug": "debug",
}


# Print usage
def print_usage(py_name):
    print(
        """
Usage: %s [-h, --help] [-d path] [-s step] [-j threads] [-c chaintype] [--organism organism] [--fwd_primer path] [--rev_primer path] [--c_gene_file path] [--debug boolean]
"""
        % py_name
    )


# Print detailed descriptions about options
def print_descriptions():
    print(
        """
Optional arguments:
    -h, --help  show this help message and exit
    -d          working directory path
    -s          step integer (1 ~ 4). if None, run whole steps
                <default = None>
    -j          number of threads
                <default = 5>
    -c          chain type (BCRH, BCRK, BCRL, TRB, TRA, TRD or TRG)
                <default = BCRH>

    --organism      organism (human, mouse, rabbit, ...)
                    <default = human>
    --fwd_primer    forward primer file
    --rev_primer    reverse primer file
    --c_gene_file   c gene file
    --debug         Save intermediete outputs for debugging (True / False)
                    <default = True>
"""
    )


# Define options
def define_getopt_params():
    short_opts = [opt for opt in OPT_DICT if len(opt) < 2]
    long_opts = [opt for opt in OPT_DICT if len(opt) > 1]
    short_opt_param = "h" + ":".join(short_opts) + ":"
    long_opt_param = ["help"] + [opt + "=" for opt in long_opts]
    return short_opt_param, long_opt_param


# Try parsing options from the command
def get_opts(py_name):
    short_opt_param, long_opt_param = define_getopt_params()
    try:
        parsed_opts, parsed_args = getopt.getopt(
            sys.argv[1:], short_opt_param, long_opt_param
        )
    except getopt.GetoptError:
        print_usage(py_name)
        sys.exit(2)
    return parsed_opts


# Define default values for each option -> dict
def default_argv():
    argv_by_opt = {opt: None for opt in OPT_DICT.values()}
    argv_by_opt["num_threads"] = MIN_THREAD
    argv_by_opt["organism"] = "human"
    argv_by_opt["chain_type"] = "BCRH"
    argv_by_opt["debug"] = "True"
    return argv_by_opt


# Get argument values by options into dictionary
def get_argv(parsed_opts, py_name):
    argv_by_opt = default_argv()
    for opt, arg in parsed_opts:
        if opt in ("-h", "--help"):
            print_usage(py_name)
            print_descriptions()
            sys.exit(0)
        else:
            opt = opt.replace("-", "")
            argv_by_opt[OPT_DICT[opt]] = arg
    return argv_by_opt


# Check working directory
def check_work_dir(working_dir, py_name):
    if not working_dir:
        print("\nMust specify the working directory path.")
        print_usage(py_name)
        sys.exit(2)


# Check if number of threads exceeds maximum threads number, and return
def check_threads(num_threads):
    if int(num_threads) > MAX_THREAD:
        num_threads = MAX_THREAD
    return int(num_threads)


if __name__ == "__main__":
    opts = get_opts(py_name=os.path.basename(__file__))
    argv_by_opt = get_argv(parsed_opts=opts, py_name=os.path.basename(__file__))
    check_work_dir(
        working_dir=argv_by_opt["working_directory"], py_name=os.path.basename(__file__)
    )
    argv_by_opt["num_threads"] = check_threads(argv_by_opt["num_threads"])
