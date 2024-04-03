import glob
import logging
import multiprocessing
import os
import re
from ast import literal_eval as make_tuple
from collections import defaultdict
from multiprocessing import Process
from Common.MultiprocessLogger import QueueHandler
from Bio.Data import CodonTable
from Common.SeqUtil import Helper
from Common.FileUtil import read_file, write_file

# for memory usage test
import psutil

class SeqExtracter(object):
    """Extracting amino acid sequence of HCDR3 region from given set of files

    First list up g_target assembled fastq file and pass them to ExtractWorker process through the queue.
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

        self.initial_RAM_used = psutil.virtual_memory().used / 1024 ** 3


    def run(self, num_threads):
        """"""
        logging.info("Started run")

        batch_size = 15
        batch_num = 1

        if num_threads > batch_size:
            num_threads = batch_size

        # list up pairs
        targets = self._find_all_targets()

        # run remain steps in parallel
        manager = multiprocessing.Manager()

        # for the test...
        RAM_used_queue = manager.Queue()
        max_RAM_used = psutil.virtual_memory().used / 1024 ** 3


        # g_target queueing, within batch size
        while True:
            if len(targets) == 0:
                break

            batch_targets = targets[:batch_size]
            targets = targets[batch_size:]
            logging.info('Batch number %d started...' % batch_num)

            target_queue = manager.Queue()
            for job in batch_targets:
                target_queue.put(job)
            extract_workers = []
            num_process = max(1, num_threads - 1)
            logging.info("Spawning %d ExtractWorker processes.", num_process)
            for i in range(num_process):
                p = ExtractWorker(self.dir, target_queue, self.config, self.log_queue, RAM_used_queue)
                p.daemon = True
                p.start()
                extract_workers.append(p)
                target_queue.put_nowait(None)

            # wait until all extracting works done
            for p in extract_workers:
                p.join()

            # get RAM_used_queue
            out_RAM_used_queue = p.get_RAM_used_queue()

            while out_RAM_used_queue.qsize() > 0:
                RAM_used = out_RAM_used_queue.get()
                if RAM_used > max_RAM_used:
                    max_RAM_used = RAM_used

            logging.info('Batch number %d finished...' % batch_num)
            batch_num += 1

        logging.info("Finished run")

        logging.info('Max RAM used : %.2f GB' % (max_RAM_used - self.initial_RAM_used))


    def _find_all_targets(self):
        """
        Find out all assembled files ordered in file size (decrease).
        """
        target_files = glob.glob(os.path.join(self.dir, 'merge', '*.assembled.fastq'))
        target_with_size = [(f, os.stat(f).st_size) for f in target_files]
        target_with_size.sort(key=lambda x: x[1], reverse=True)

        return [t[0] for t in target_with_size]


# noinspection PyTypeChecker
class ExtractWorker(Process):
    """Extract the amino acid sequence of HCDR3 region from given assembled sequence file

    Handling of given data follows the bellowing steps.
    First, filter by quality score and length. (min score=20, min length = 150)
    Second, find primer region and extract scFv sequence.
    Third, extract HCDR3 region.
    Finally, write down full sequence and statistics file.
    """

    def __init__(self, _dir, in_queue, config, log_queue, RAM_used_queue):
        Process.__init__(self)

        self.dir = _dir
        self.helper = Helper()
        self.in_queue = in_queue
        self.log_queue = log_queue

        self.platform = config['platform']
        self.filter_text = config['filter_condition']
        filter_m = re.match('q(\d+)p(\d+)', self.filter_text)
        self.filter_q, self.filter_p = int(filter_m.group(1)), float(filter_m.group(2)) / 100
        self.run_type = config['run_type']

        self.rawCSV_each = {}
        self.statCSV_each = {}

        self.RAM_used_queue = RAM_used_queue

        try:
            re_trans = re.compile("\(.*?,.*?\)")
            trans_groups = re_trans.findall(config['translation'])
            standard_table = CodonTable.unambiguous_dna_by_id[1]
            self.codon_table = standard_table.forward_table
            for triplet in standard_table.stop_codons:
                self.codon_table[triplet] = '*'
            for tu in trans_groups:
                triplet, aminoacid = make_tuple(tu)
                self.codon_table[triplet] = aminoacid
        except KeyError:
            self.codon_table = None

    def run(self):
        if self.log_queue:
            logger = logging.getLogger()
            if len(logger.handlers) == 0:
                h = QueueHandler(self.log_queue)
                logger.addHandler(h)
                logger.setLevel(logging.DEBUG)

        logging.info("Started run")
        while True:
            # Get the work from the queue and expand the tuple
            target = self.in_queue.get()
            if target is None:
                self.in_queue.task_done()
                break
            self.extract_sequence(target)

            tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
            self.RAM_used_queue.put(tmp_RAM_used)

            self.in_queue.task_done()

        self.extract_rawCSV_and_statCSV()

        logging.info("Finished run")


    def extract_sequence(self, target, chunk_size=200000):
        """Extract g_target sequences(nucleotide) for each NGS merged read and put them into the result queue.
        ;
        """
        target_file = os.path.basename(target)

        logging.info("Started extract_target %s", target_file)

        sample_name = target_file.split('.')[0]
        split = re.split('[-_]', sample_name)
        library_name = split[0]
        if len(split) > 1:
            extra_info = '_'.join(split[1:])
        else:
            extra_info = ''

        # Get read count from discarded and unassembled files
        count_unmerged = self._get_read_count_residue(target)

        count_merged = 0
        count_qfiltered = 0
        count_nfiltered = 0
        sequence_hash = defaultdict(int)

        # set filtering length according to NGS platform (miseq / novaseq : 250, nextseq : 150)
        if self.platform == 'miseq' or self.platform == 'novaseq':
            filter_len = 250
            logging.info('Platform: %s... Setting filtering length as 250' % self.platform)
        elif self.platform == 'nextseq':
            filter_len = 150
            logging.info('Platform: %s... Setting filtering length as 150' % self.platform)
        else:
            logging.info('Not supporting NGS platform... Setting filtering length as 0')
            filter_len = 0

        handle = open(target, 'r')
        try:
            while True:
                # read sequence
                seq_id = handle.readline().strip()[1:]
                sequence = handle.readline().strip()
                handle.readline()
                score = [ord(c) - 33 for c in handle.readline().strip()]

                if seq_id == '':
                    break

                count_merged += 1

                if len(sequence) < filter_len:
                    continue

                if count_merged % (chunk_size * 5) == 0:
                    count_total = count_merged + count_unmerged
                    logging.info('Processing %s-%s, %d hit in %d reads (quality & length filtering)', library_name, extra_info, count_qfiltered, count_total)
                    logging.info('Processing %s-%s, %d hit in %d reads (n-filtering)', library_name, extra_info, count_nfiltered, count_qfiltered)

                p_q = sum(1 for q in score if q >= self.filter_q) / len(score)
                # filter out low-quality sequences
                if p_q < self.filter_p:
                    pass
                else:
                    count_qfiltered += 1

                    # filter out sequences with N nucleotide
                    if 'N' in sequence.upper():
                        pass
                    else:
                        count_nfiltered += 1
                        sequence_hash[sequence] += 1

        except ValueError as e:
            logging.error(e)

        handle.close()

        self.rawCSV_each[sample_name] = []
        for k, v in sequence_hash.items():
            self.rawCSV_each[sample_name].append({'full_NT': k, 'readcount': v})

        count_clones = len(sequence_hash)
        count_total = count_merged + count_unmerged
        logging.info('Done %s-%s, %d hit in %d reads (quality & length filtering)', library_name, extra_info, count_qfiltered, count_total)
        logging.info('Done %s-%s, %d hit in %d reads (n-filtering)', library_name, extra_info, count_nfiltered, count_qfiltered)

        # tmp_dict = {'sample_name': sample_name, 'qfilter': self.filter_text,
        #             'count_total': count_total, 'count_merged': count_merged, 'count_unmerged': count_unmerged,
        #             'count_qfiltered': count_qfiltered, 'percent_qfiltered': count_qfiltered/count_total*100,
        #             'count_nfiltered': count_nfiltered, 'percent_nfiltered': count_nfiltered/count_qfiltered*100,
        #             'count_unique_seq_nt': count_clones}

        tmp_dict = {'sample_name': sample_name, 'qfilter': self.filter_text,
                    'total_reads': count_total, 'merged_reads': count_merged,
                    'functional_reads': count_nfiltered, 'unique_functional_reads': count_clones}
        if count_total == 0:
            tmp_dict['quality_filter_eff'] = 'N/A'
            tmp_dict['whole_eff'] = 'N/A'
        elif count_total > 0 and count_merged == 0:
            tmp_dict['quality_filter_eff'] = 'N/A'
            tmp_dict['whole_eff'] = count_nfiltered / count_total * 100
        else:
            tmp_dict['quality_filter_eff'] = count_nfiltered / count_merged * 100
            tmp_dict['whole_eff'] = count_nfiltered / count_total * 100

        self.statCSV_each[sample_name] = tmp_dict

        logging.info("Finished extract_target %s", target_file)


    def extract_rawCSV_and_statCSV(self):
        raw_header = ['full_NT', 'readcount']
        stat_header = ['sample_name', 'qfilter',
                       'total_reads', 'merged_reads',
                       'functional_reads', 'unique_functional_reads',
                       'quality_filter_eff', 'whole_eff']

        out_dir = os.path.join(self.dir, 'out')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # write each raw data files according to sample_id
        for sample_name in self.rawCSV_each:
            with open(os.path.join(out_dir, '%s.csv' % sample_name), 'w') as handle:
                write_file(handle, raw_header, self.rawCSV_each[sample_name])

        # write whole statistics file
        stat_dir = os.path.join(self.dir, 'stat')
        if not os.path.exists(stat_dir):
            os.mkdir(stat_dir)

        stat_outfile = os.path.join(stat_dir, 'Statistics.csv')
        if not os.path.isfile(stat_outfile):
            with open(stat_outfile, 'w') as handle:
                write_file(handle, stat_header, list(self.statCSV_each.values()))
        else:
            with open(stat_outfile, 'a+') as handle:
                handle.writelines([','.join([str(x[h]) for h in stat_header]) + '\n' for x in list(self.statCSV_each.values())])


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

    def get_RAM_used_queue(self):
        return self.RAM_used_queue



def sorting(primer_info, sequence_list, sample_name, id_prefix):
    header_p, data_p = read_file(primer_info)
    header_s, data_s = read_file(sequence_list)

    helper = Helper()

    seq_map = {}
    primer_names = []
    for pr in data_p:
        if pr['sample'] == sample_name:
            primer_names.append(pr['name'])

    total_reads = 0
    seq_map['total'] = {'sequence': 'sum', 'sequence AA': '-', 'total': 0}
    for pn in primer_names:
        seq_map['total'][pn] = 0

    for sd in data_s:
        seq = sd['target_sequence']
        if seq not in seq_map:
            seq_map[seq] = {'sequence': seq, 'sequence AA': helper.translate(seq[1:]), 'total': 0}
            for pn in primer_names:
                seq_map[seq][pn] = 0
        seq_map[seq][sd['primer_name']] = int(sd['read_count'])
        seq_map['total'][sd['primer_name']] += int(sd['read_count'])
        seq_map[seq]['total'] += int(sd['read_count'])
        seq_map['total']['total'] += int(sd['read_count'])
        total_reads += int(sd['read_count'])

    seq_data = list(seq_map.values())
    seq_data.sort(key=lambda x: x['total'], reverse=True)

    digit = len(str(len(seq_data)))
    for i, sd in enumerate(seq_data):
        # for pn in primer_names:
        #    sd[pn + ' (PPM)'] = sd[pn] / total_reads * 1000000
        sd['clone_id'] = id_prefix + ('%0' + str(digit) + 'd') % i

    seq_header = ['clone_id', 'sequence', 'sequence AA', 'total'] + [pr for pr in primer_names]
    # list(itertools.chain(*[[pr, pr + ' (PPM)'] for pr in primer_names]))
    with open(sequence_list.replace('.csv', '_sorted.csv'), 'w') as f:
        write_file(f, seq_header, seq_data)


def fullseq(primer_info, sorted_list, sample_name):
    header_p, data_p = read_file(primer_info)
    header_s, data_s = read_file(sorted_list)

    seq_data = []
    primer_set = {}
    for pr in data_p:
        if pr['sample'] == sample_name:
            primer_set[pr['name']] = primer_set


if __name__ == '__main__':
    sorting(
        r'D:\Celemics\SyntheticBiologyTeam\Clustering\Data\180205_NMO_Mi\primer_info.csv',
        r'D:\Celemics\SyntheticBiologyTeam\Clustering\Data\180205_NMO_Mi\analysis\vK_all.csv',
        'K', '180205_NMO_vK-'
    )
