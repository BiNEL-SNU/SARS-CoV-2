import configparser
import getopt
import glob
import gzip
import logging
import os
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .ProcessUtil import run_pear


class PESParser(object):
    """Parsing NGS 'P'aired-'E'nd 'S'equencing Data.

    Handling of given data follows the steps below.
    List up g_target gz file pairs and pass them to new process to conduct following steps.
    """

    def __init__(self, _dir, config, is_gz: "str" = "T"):
        """"""
        # Setup working directory
        if not os.path.isdir(_dir):
            raise RuntimeError
        else:
            self.dir = _dir

        config["directory"] = self.dir
        try:
            self.primer_forward_length = int(config["primer_forward_length"])
        except KeyError:
            self.primer_forward_length = 0
        try:
            self.primer_reverse_length = int(config["primer_reverse_length"])
        except KeyError:
            self.primer_reverse_length = 0
        self.config = config

        self.is_gz = is_gz

    def run(self, num_threads):
        """"""
        logging.info("Started run")

        working_dir = os.path.join(self.dir, "merge"), os.path.join(self.dir, "stat")
        for directory in working_dir:
            if not os.path.exists(directory):
                os.makedirs(directory)

        # list up pairs
        target_matches = self._find_all_pair()

        logging.info("Started work on %d targets" % (len(target_matches)))

        name_overlap_check = {}
        for m in target_matches:
            sampleName = m.group(1)
            # print(sampleName)
            sample_list = sampleName.split("_")
            if len(sample_list) >= 3 and sample_list[-1][0] == "S":
                sampleID = "_".join(sample_list[:-1])
            elif (
                len(sample_list) >= 3
                and sample_list[-2][0] == "S"
                and sample_list[-1][0] == "L"
            ):
                sampleID = "_".join(sample_list[:-2])
            elif len(sample_list) == 2 and sample_list[-1][0] == "S":
                sampleID = sample_list[0]
            else:
                sampleID = "_".join(sample_list)

            if sampleID not in name_overlap_check:
                name_overlap_check[sampleID] = 1
            else:
                sampleID_original = sampleID
                sampleID = sampleID + "_add%d" % name_overlap_check[sampleID]
                name_overlap_check[sampleID_original] += 1

            if "Undetermined" in sampleName:
                logging.info("pass the PEMerge of the sample Undetermined")
                continue

            target = (
                sampleID,
                m.string[: m.start(2)] + "R1" + m.string[m.end(2) :],
                m.string[: m.start(2)] + "R2" + m.string[m.end(2) :],
            )
            assem_file = self._merge_pair(target, num_threads)
            # self._quailty_check('stat', target, assem_file)

        logging.info("Finshed run")

    def _find_all_pair(self):
        """Find out all gz files

        :return: found file name list with %s instead of 'R1' and 'R2', without .gz
        """
        target_matches = []

        if self.is_gz == "True":
            gz_files = glob.glob(os.path.join(self.dir, "gz", "*.fastq.gz"))
            gz_files = [[os.path.basename(file), False] for file in gz_files]

            # p = re.compile(r'(.*)_S[0-9]+_L[0-9]+_(R1)[a-zA-Z0-9_-]*.fastq.gz')
            p = re.compile(r"(.*)_(R1)[a-zA-Z0-9_-]*.fastq.gz")
            # p = re.compile(r'(.*)_S[0-9]+_(R1)[a-zA-Z0-9_-]*.fastq.gz')
            for i, (file, flag) in enumerate(gz_files):
                if flag:
                    continue
                # file name is match to R1 format
                m = p.fullmatch(file)
                if m:
                    r2 = file[: m.start(2)] + "R2" + file[m.end(2) :]
                    # find if corresponding R2 file exists
                    if gz_files.count([r2, False]) == 1:
                        j = gz_files.index([r2, False])
                        gz_files[i][1] = True
                        gz_files[j][1] = True
                        target_matches.append(m)

        elif self.is_gz == "False":
            gz_files = glob.glob(os.path.join(self.dir, "gz", "*.fastq"))
            gz_files = [[os.path.basename(file), False] for file in gz_files]

            # p = re.compile(r'(.*)_S[0-9]+_L[0-9]+_(R1)[a-zA-Z0-9_-]*.fastq.gz')
            p = re.compile(r"(.*)_(R1)[a-zA-Z0-9_-]*.fastq")
            # p = re.compile(r'(.*)_S[0-9]+_(R1)[a-zA-Z0-9_-]*.fastq.gz')
            for i, (file, flag) in enumerate(gz_files):
                if flag:
                    continue
                # file name is match to R1 format
                m = p.fullmatch(file)
                if m:
                    r2 = file[: m.start(2)] + "R2" + file[m.end(2) :]
                    # find if corresponding R2 file exists
                    if gz_files.count([r2, False]) == 1:
                        j = gz_files.index([r2, False])
                        gz_files[i][1] = True
                        gz_files[j][1] = True
                        target_matches.append(m)

        target_with_size = [
            (
                m,
                os.stat(os.path.join(self.dir, "gz", m.string)).st_size
                + os.stat(
                    os.path.join(
                        self.dir,
                        "gz",
                        m.string[: m.start(2)] + "R2" + m.string[m.end(2) :],
                    )
                ).st_size,
            )
            for m in target_matches
        ]

        target_with_size.sort(key=lambda x: x[1], reverse=True)

        return [t[0] for t in target_with_size]

    def _merge_pair(self, target, num_threads) -> str:
        """Merge the paired-end sequence reads using PEAR - Paired-End reAd mergeR

        PEAR url: http://sco.h-its.org/exelixis/web/software/pear/
        """

        if not os.path.exists(os.path.join(self.dir, "merge")):
            os.makedirs(os.path.join(self.dir, "merge"))

        out_file = os.path.join(self.dir, "merge", target[0])
        if out_file[-1] == "-" or out_file[-1] == "_":
            out_file = out_file[:-1]
        out_file_assem = Path(out_file + ".assembled.fastq")

        if out_file_assem.is_file() and os.stat(str(out_file_assem)).st_size > 0:
            logging.info(
                "%s file already exists. Pass merge process."
                % os.path.basename(out_file)
            )
            return str(out_file_assem)
        else:
            # Merge forward and reverse consensus data
            logging.info("Start pear.exe on %s." % os.path.basename(out_file))

            fastq_r1 = os.path.join(self.dir, "gz", target[1])
            fastq_r2 = os.path.join(self.dir, "gz", target[2])

            if run_pear(fastq_r1, fastq_r2, out_file, num_threads) == None:
                return None

        return str(out_file_assem)

    def _quailty_check(self, directory: str, target, assem_file: str):
        stat_file = os.path.join(
            self.dir,
            directory,
            os.path.basename(assem_file.replace(".assembled.fastq", ".csv")),
        )
        if Path(stat_file).is_file():
            logging.info(
                "%s file already exists. Pass _quailty_check process."
                % os.path.basename(stat_file)
            )
            return
        else:
            handle = open(assem_file, "r")
            total = 0
            count_min = [0] * 41  # 40 ~ 0
            count_median = [0] * 41  # 40 ~ 0
            cumul_min = [0] * 41  # 40 ~ 0
            cumul_median = [0] * 41  # 40 ~ 0
            min_l, max_l, sum_l = 99999, 0, 0
            while True:
                seq_id = handle.readline().strip()[1:]
                sequence = handle.readline().strip()
                handle.readline()
                score = [ord(c) - 33 for c in handle.readline().strip()]
                score = score[
                    self.primer_forward_length : len(score) - self.primer_reverse_length
                ]
                if seq_id == "":
                    break

                total += 1

                if len(score) > 0:
                    count_min[40 - min(score)] += 1
                    count_median[40 - int(np.median(score))] += 1

                if len(sequence) > max_l:
                    max_l = len(sequence)
                if len(sequence) < min_l:
                    min_l = len(sequence)
                sum_l += len(sequence)
            handle.close()

            if total == 0:
                logging.info(
                    "%s file is empty.. Pass _quality_check process."
                    % os.path.basename(stat_file)
                )
                return

            fastq_r1 = os.path.join(self.dir, "gz", target[1])
            fastq_r2 = os.path.join(self.dir, "gz", target[2])
            with gzip.open(fastq_r1, "r") as handle:
                r1_base_count = 0
                r1_score_sum = 0
                while True:
                    seq_id = handle.readline().strip()[1:].decode("ascii")
                    sequence = handle.readline().strip().decode("ascii")
                    handle.readline()
                    score = [
                        ord(c) - 33 for c in handle.readline().strip().decode("ascii")
                    ]
                    score = score[
                        self.primer_forward_length : len(score)
                        - self.primer_reverse_length
                    ]
                    if seq_id == "":
                        break

                    r1_base_count += len(score)
                    r1_score_sum += sum(score)

            with gzip.open(fastq_r2, "r") as handle:
                r2_base_count = 0
                r2_score_sum = 0
                while True:
                    seq_id = handle.readline().strip()[1:].decode("ascii")
                    sequence = handle.readline().strip().decode("ascii")
                    handle.readline()
                    score = [
                        ord(c) - 33 for c in handle.readline().strip().decode("ascii")
                    ]
                    score = score[
                        self.primer_forward_length : len(score)
                        - self.primer_reverse_length
                    ]
                    if seq_id == "":
                        break

                    r2_base_count += len(score)
                    r2_score_sum += sum(score)

            with open(stat_file, "w") as f:
                f.write("Basic statistics\n")
                f.write("File Name,%s\n" % os.path.basename(assem_file))
                f.write("Total Sequences (merged),%d\n" % total)
                f.write("Sequence Min Length (merged),%d\n" % min_l)
                f.write("Sequence Max Length (merged),%d\n" % max_l)
                f.write("Sequence Average Length (merged),%f\n" % (sum_l / total))
                f.write("R1 Average Q-score,%f\n" % (r1_score_sum / r1_base_count))
                f.write("R2 Average Q-score,%f\n" % (r2_score_sum / r2_base_count))
                f.write(
                    "Total Average Q-score,%f\n\n"
                    % ((r1_score_sum + r2_score_sum) / (r1_base_count + r2_base_count))
                )

                f.write("Q-score density\n")
                f.write(
                    "Q-score,min,min(%),min(cumulative),min(cumulative,%),median,median(%),median(cumulative),"
                    + "median(cumulative,%)\n"
                )
                sum_min, sum_median = 0, 0
                for i in range(len(count_min)):
                    sum_min += count_min[i]
                    sum_median += count_median[i]
                    cumul_min[i] = sum_min
                    cumul_median[i] = sum_median
                    f.write(
                        "%d,%d,%.2f,%d,%.2f,%d,%.2f,%d,%.2f\n"
                        % (
                            40 - i,
                            count_min[i],
                            count_min[i] / total * 100,
                            sum_min,
                            sum_min / total * 100,
                            count_median[i],
                            count_median[i] / total * 100,
                            sum_median,
                            sum_median / total * 100,
                        )
                    )

            stat_figure = stat_file.replace(".csv", ".png")
            fig, ax = plt.subplots()
            ind = np.arange(41)
            width = 0.4
            r1 = ax.bar(ind, [c / total for c in cumul_min], width, color="c")
            r2 = ax.bar(
                ind + width, [c / total for c in cumul_median], width, color="y"
            )

            # add some text for labels, title and axes ticks
            ax.set_ylabel("Cumulative Frequency")
            title = "Cumulative Frequency Distribution of Min / Median Q-Value\n"
            title += "Sample: %s, Total Sequences: %d" % (
                os.path.basename(assem_file).replace(".assembled.fastq", ""),
                total,
            )
            ax.set_title(title)
            ax.set_xticks(np.arange(0, 41, 5) + width)
            ax.set_xticklabels(["%.1f" % float(i) for i in range(40, -1, -5)])

            ax.legend((r1[0], r2[0]), ("Min", "Median"), loc="lower right")
            plt.savefig(stat_figure)
