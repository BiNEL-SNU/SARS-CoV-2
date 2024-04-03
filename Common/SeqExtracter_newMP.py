import glob
import logging
import multiprocessing
import os
import re
# from ast import literal_eval as make_tuple
from collections import defaultdict
from multiprocessing import Process
from Common.MultiprocessLogger import QueueHandler
# from Bio.Data import CodonTable
from Common.SeqUtil import Helper
from Common.FileUtil import read_file, write_file

# for memory usage test
import psutil

class SeqExtracter(object):
    """ Extracting  from given set of files, multiprocessing each file into file chunk
          - list up g_target assembled fastq files and pass them to ExtractWorker process through the queue.
          - one file at a time, split into file chunk and spawn to multiple processes
    """

    def __init__(self, _dir, config, multiprocess_log_queue=None):
        """"""
        # Setup working directory
        if not os.path.isdir(_dir):
            raise RuntimeError
        else:
            self.dir = _dir

        self.config = config
        self.log_queue = multiprocess_log_queue

        self.platform = config['platform']

        self.initial_RAM_used = psutil.virtual_memory().used / 1024 ** 3


    def run(self, num_threads):
        # for the RAM usage check
        max_RAM_used = psutil.virtual_memory().used / 1024 ** 3

        # list up all merged files (targets)
        targets = self._find_all_targets()

        # # (211213. temporary!) for Behcet alpha chain data (B4_t1_A)
        # targets = [os.path.join(self.dir, 'merge', 'B4_t1_A.assembled.fastq')]

        # run remain steps in parallel
        manager = multiprocessing.Manager()

        # RAM_used queue
        RAM_used_queue = manager.Queue()


        # set filtering length according to NGS platform (miseq / novaseq : 250, nextseq : 150)
        if self.platform == 'miseq' or self.platform == 'novaseq':
            self.filter_len = 250
            # self.filter_len = 200
            logging.info('Platform: %s... Setting filtering length as 250' % self.platform)
        elif self.platform == 'nextseq':
            self.filter_len = 150
            logging.info('Platform: %s... Setting filtering length as 150' % self.platform)
        else:
            logging.info('Not supporting NGS platform... Setting filtering length as 0')
            self.filter_len = 0


        # g_target queueing, within chunk size
        while True:
            # set queue for getting data from ExtractWorker processes
            raw_queue = manager.Queue()
            stat_queue = manager.Queue()

            if len(targets) == 0:
                break

            chunk_size_max = 20000000

            target_to_put = targets[0]

            # initialization of csv dict. for each target file.
            self.rawCSV = {}
            self.statCSV = {}

            # self.target_name = os.path.basename(target_to_put).split('.')[0]

            # (210520 -> 210622). for dealing with the problem of sample name in csv files
            tmp_name = os.path.basename(target_to_put).split('.')[0]
            tmp_name_prefix = tmp_name.replace('_', '-').split('-')[0]
            try:
                int(tmp_name_prefix)
                self.target_name = 'f_' + tmp_name
            except ValueError:
                self.target_name = tmp_name

            logging.info("Started extract_target %s", self.target_name)

            # Get readcount from discarded and unassembled files
            self.count_unmerged = self._get_read_count_residue(target_to_put)

            targets.remove(target_to_put)

            self.num_process = max(1, num_threads)

            # get line number of file and get chunk start & end position list
            file_line_len = self.get_file_lines(target_to_put)
            chunk_size = (int(file_line_len / (self.num_process * 4)) + 1) * 4
            if chunk_size < 4:
                chunk_size = 4
            elif chunk_size > chunk_size_max:
                chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
            chunk_start_pos_list = self.extract_chunk_start_pos(file_line_len, chunk_size)


            # save chunk start & end position data into input queue
            target_queue = manager.Queue()
            for chunk_range in chunk_start_pos_list:
                target_queue.put(chunk_range)

            logging.info("Spawning %d ExtractWorker processes.", self.num_process)

            extract_workers = []
            for i in range(self.num_process):
                p = ExtractWorker(self.dir, target_to_put, target_queue, chunk_size,
                                  self.config, self.log_queue, self.filter_len,
                                  raw_queue, stat_queue, RAM_used_queue)
                p.daemon = True
                p.start()

                extract_workers.append(p)
                target_queue.put_nowait(None)

                # self.update_chunk_data()

            # wait until all extracting works done
            for p in extract_workers:
                p.join()

            # merge into one dict
            out_raw_queue = p.get_raw_queue()
            out_stat_queue = p.get_stat_queue()

            self.merge_chunk_data(out_raw_queue, out_stat_queue)

            logging.info("Finished extract_target %s", self.target_name)

            # get RAM_used_queue
            out_RAM_used_queue = p.get_RAM_used_queue()
            while out_RAM_used_queue.qsize() > 0:
                RAM_used = out_RAM_used_queue.get()
                if RAM_used > max_RAM_used :
                    max_RAM_used = RAM_used

            # write result files
            self.write_csv()

        logging.info("Extraction for all files are finished")
        logging.info('Max RAM used : %.2f GB' % (max_RAM_used - self.initial_RAM_used))


    def _find_all_targets(self):
        """
        Find out all assembled files ordered in file size (decrease).
        """
        target_files = glob.glob(os.path.join(self.dir, 'merge', '*.assembled.fastq'))
        target_with_size = [(f, os.stat(f).st_size) for f in target_files]
        target_with_size.sort(key=lambda x: x[1], reverse=True)

        return [t[0] for t in target_with_size]


    def get_file_lines(self, target):
        count = 0
        with open(target, 'r') as handle:
            for line in handle:
                count += 1
        return count


    def extract_chunk_start_pos(self, target_line_len, chunk_size):
        chunk_start_list = []

        i = 0
        while True:
            if chunk_size * i > target_line_len:
                break

            line_start = chunk_size * i
            i += 1

            chunk_start_list.append(line_start)

        logging.info('Chunk size for target %s, distributed into %d processes: %d' % (self.target_name, self.num_process, chunk_size))

        return chunk_start_list


    def merge_chunk_data(self, in_raw_queue, in_stat_queue):
        # for raw csv
        while in_raw_queue.qsize() > 0:
            raw_chunk = in_raw_queue.get()
            for k in raw_chunk:
                if k not in self.rawCSV:
                    self.rawCSV[k] = raw_chunk[k]
                else:
                    self.rawCSV[k]['readcount'] += raw_chunk[k]['readcount']

        # for stat csv
        while in_stat_queue.qsize() > 0:
            stat_chunk = in_stat_queue.get()
            if len(self.statCSV) == 0:
                self.statCSV = stat_chunk
            else:
                self.statCSV['merged_reads'] += stat_chunk['merged_reads']
                self.statCSV['q-filtered_reads'] += stat_chunk['q-filtered_reads']


    def write_csv(self):
        raw_header = ['full_NT', 'readcount']
        stat_header = ['sample_name', 'qfilter',
                       'total_reads', 'merged_reads',
                       'q-filtered_reads', 'unique_q-filtered_reads',
                       'merge_eff', 'q-filter_eff', 'whole_eff']

        # extra info setting for statCSV.
        count_merged = self.statCSV['merged_reads']
        count_nfiltered = self.statCSV['q-filtered_reads']
        count_total = count_merged + self.count_unmerged

        if count_total == 0:
            self.statCSV['merge_eff'] = 'N/A'
            self.statCSV['q-filter_eff'] = 'N/A'
            self.statCSV['whole_eff'] = 'N/A'
        elif count_merged == 0:
            self.statCSV['merge_eff'] = 0
            self.statCSV['q-filter_eff'] = 'N/A'
            self.statCSV['whole_eff'] = 0
        else:
            self.statCSV['merge_eff'] = count_merged / count_total * 100
            self.statCSV['q-filter_eff'] = count_nfiltered / count_merged * 100
            self.statCSV['whole_eff'] = count_nfiltered / count_total * 100

        self.statCSV['sample_name'] = self.target_name
        self.statCSV['total_reads'] = count_total
        self.statCSV['unique_q-filtered_reads'] = len(self.rawCSV)

        # write raw data file
        out_dir = os.path.join(self.dir, 'out')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        with open(os.path.join(out_dir, '%s.csv' % self.target_name), 'w') as handle:
            write_file(handle, raw_header, list(self.rawCSV.values()))


        # write whole statistics file
        stat_dir = os.path.join(self.dir, 'stat')
        if not os.path.exists(stat_dir):
            os.mkdir(stat_dir)

        stat_outfile = os.path.join(stat_dir, 'Statistics.csv')
        if not os.path.isfile(stat_outfile):
            with open(stat_outfile, 'w') as handle:
                write_file(handle, stat_header, [self.statCSV])
        else:
            with open(stat_outfile, 'a+') as handle:
                handle.writelines([','.join([str(self.statCSV[h]) for h in stat_header]) + '\n'])


    @staticmethod
    def _get_read_count_residue(target):
        """Get counts of discarded and unassembled reads
        """

        discarded = target.replace('.assembled.', '.discarded.')
        unassembled = target.replace('.assembled.', '.unassembled.forward.')

        count = 0
        d_files = glob.glob(discarded)
        u_files = glob.glob(unassembled)
        if len(d_files) == 1:
            records = open(d_files[0], 'r')
            try:
                for _ in records:
                    count += 1
                count = int(count / 4)
            except ValueError as e:
                logging.info(e)
        if len(u_files) == 1:
            records = open(u_files[0], 'r')
            try:
                for _ in records:
                    count += 1
                count = int(count / 4)
            except ValueError as e:
                logging.info(e)

        target_file = os.path.basename(target)
        split = re.split('[-_]', target_file)
        library_name = split[0]
        if len(split) > 1:
            extra_info = '%s' % split[1]
        else:
            extra_info = ''
        logging.info("%s-%s Unmerged reads: %d", library_name, extra_info, count)
        return count


# noinspection PyTypeChecker
class ExtractWorker(Process):
    """Extract the amino acid sequence of HCDR3 region from given assembled sequence file

    Handling of given data follows the bellowing steps.
    First, filter by quality score and length. (following the option: qfilter & platform)
    Finally, write down full sequence and statistics file.
    """

    def __init__(self, _dir, target_file, in_queue, chunk_size, config, log_queue, filter_len, raw_queue, stat_queue, RAM_used_queue):
        Process.__init__(self)

        self.dir = _dir
        self.helper = Helper()
        self.target = target_file
        self.in_queue = in_queue
        self.chunk_size = chunk_size
        self.log_queue = log_queue

        self.raw_queue = raw_queue
        self.stat_queue = stat_queue

        self.filter_text = config['filter_condition']
        filter_m = re.match('q(\d+)p(\d+)', self.filter_text)
        self.filter_q, self.filter_p = int(filter_m.group(1)), float(filter_m.group(2)) / 100
        self.run_type = config['run_type']

        self.rawCSV_chunk = {}
        self.statCSV_chunk = {}

        self.filter_len = filter_len

        self.RAM_used_queue = RAM_used_queue


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
            self.extract_sequence(chunk_start_pos)

            tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
            self.RAM_used_queue.put(tmp_RAM_used)

            self.in_queue.task_done()

            # save chunk data into result queue
            self.raw_queue.put(self.rawCSV_chunk)
            self.stat_queue.put(self.statCSV_chunk)

        # logging.info("Finished run")


    def extract_sequence(self, chunk_start_pos):
        """ Extract g_target sequences(nucleotide) for each NGS merged read and put them into the result queue. """

        count_merged = 0
        count_nfiltered = 0
        sequence_hash = defaultdict(int)


        # read each line from the chunk start to end
        handle = open(self.target, 'r')

        # go to chunk_start position
        for i in range(0, chunk_start_pos):
            handle.readline()

        try:
            count = 0
            while True:
                seq_id = handle.readline().strip()[1:]
                sequence = handle.readline().strip()
                handle.readline()
                score = [ord(c) - 33 for c in handle.readline().strip()]

                count += 4

                if seq_id == '':
                    break
                if count > self.chunk_size:
                    break

                count_merged += 1

                if len(sequence) < self.filter_len:
                    continue

                p_q = sum(1 for q in score if q >= self.filter_q) / len(score)
                # filter out low-quality sequences
                if p_q < self.filter_p:
                    pass
                else:
                    # filter out sequences with N nucleotide
                    if 'N' in sequence.upper():
                        pass
                    else:
                        count_nfiltered += 1
                        sequence_hash[sequence] += 1

        except ValueError as e:
            logging.error(e)

        handle.close()

        self.rawCSV_chunk = {}
        for k, v in sequence_hash.items():
            self.rawCSV_chunk[k] = {'full_NT': k, 'readcount': v}
        self.statCSV_chunk = {'qfilter': self.filter_text, 'merged_reads': count_merged, 'q-filtered_reads': count_nfiltered}

        # logging.info("Finished extract on chunk %d", chunk_num)


    def get_raw_queue(self):
        return self.raw_queue

    def get_stat_queue(self):
        return self.stat_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue
