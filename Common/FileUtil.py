import ast
import csv
import glob
import io
import os
import re

from typing import List, Dict, Tuple
from xlsxwriter.workbook import Workbook

TableData = List[Dict[str, str]]
TableHeader = List[str]


# read fasta files into dictionary: key=header, value=sequence
def read_fasta(inputFile: str):
    # edit reference files (multi-line sequences into one line)
    original_cnt = 0
    fasta_dict = {}
    with open(inputFile, 'r') as f:
        lines = f.readlines()
    head = ''
    for i in range(0, len(lines)):
        if i == len(lines) - 1:
            seq += lines[i].strip('\n')
            fasta_dict[head] = seq
        else:
            if lines[i][0] == '>':
                original_cnt += 1
                if head != '':
                    fasta_dict[head] = seq
                head = lines[i].strip('\n')[1:]
                seq = ''
            else:
                seq += lines[i].strip('\n')

    return fasta_dict


def read_fasta_v2(in_file: str) -> Dict[str, str]:
    resultSeqList = []
    with open(in_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                if line[0] != '>':
                    print('File format error!')
                    raise RuntimeError
                tempSeq = {'name':line[1:-1], 'sequence':''}
            elif line[0] == '>':
                resultSeqList.append(tempSeq)
                tempSeq = {'name':line[1:-1], 'sequence':''}
            else:
                tempSeq['sequence'] += line[:-1].upper()
        resultSeqList.append(tempSeq)
    return resultSeqList


def read_file(in_file: str, delim: str = ',', header: List[str] = None) -> Tuple[TableHeader, TableData]:
    # noinspection PyBroadException
    data = []
    with open(in_file) as f:
        lines = f.readlines()
        f.close()

        if header is None:
            header = lines[0].strip().split(delim)
            for i in range(len(header)):
                header[i] = header[i].replace('"', '')
                header[i] = header[i].replace("'", '')
            data_start = 1
        else:
            data_start = 0

        for line in lines[data_start:]:
            args = line.strip().split(delim)
            if len(header) == len(args):
                line_data = dict()
                for i in range(len(header)):
                    if args[i] == 'NULL':
                        line_data[header[i]] = 'None'
                    else:
                        try:
                            line_data[header[i]] = ast.literal_eval(args[i].replace("'", '').replace('"', ""))
                        except ValueError:
                            line_data[header[i]] = args[i].replace("'", '').replace('"', "")
                        except SyntaxError:
                            line_data[header[i]] = args[i].replace("'", '').replace('"', "")
                data.append(line_data)
    return header, data


def read_file_new(in_file: str, delim: str = ',', header: List[str] = None) -> Tuple[TableHeader, TableData]:
    # noinspection PyBroadException
    data = []
    with open(in_file) as f:
        lines = f.readlines()
        f.close()
        # return [], []

        if header is None:
            header = lines[0].strip().split(delim)
            for i in range(len(header)):
                header[i] = header[i].replace('"', '')
                header[i] = header[i].replace("'", '')
            data_start = 1
        else:
            data_start = 0

        for line in lines[data_start:]:
            args = line.strip().split(delim)
            if len(header) == len(args):
                line_data = dict()
                for i in range(len(header)):
                    if args[i] == 'NULL':
                        line_data[header[i]] = 'None'
                    else:
                        if len(re.findall('[-_]', args[i].replace("'", '').replace('"', ""))) > 0:
                            line_data[header[i]] = args[i].replace("'", '').replace('"', "")
                        else:
                            try:
                                line_data[header[i]] = ast.literal_eval(args[i].replace("'", '').replace('"', ""))
                            except ValueError:
                                line_data[header[i]] = args[i].replace("'", '').replace('"', "")
                            except SyntaxError:
                                line_data[header[i]] = args[i].replace("'", '').replace('"', "")
                data.append(line_data)
    return header, data


def read_file_fast(in_file: str, delim: str = ',', header: List[str] = None) -> Tuple[TableHeader, TableData]:
    # noinspection PyBroadException
    data = []
    with open(in_file, encoding='cp949', mode='r') as f:
        if header is None:
            header = f.readline().strip().replace('"', '').replace("'", "").split(delim)

        numeric_keys = {}
        while True:
            line = f.readline()
            if not line:
                break
            # args = line.strip().split(delim)
            if len(line.split('"')) > 1 and len(line.split('"')) % 2 == 1:
                ind = 1
                splitted_line = line.split('"')
                while True:
                    splitted_line[ind] = splitted_line[ind].replace(delim, ' ')
                    ind += 2
                    if ind > len(splitted_line) - 1:
                        break
                line = ''.join(splitted_line)
            elif len(line.split("'")) > 1 and len(line.split("'")) % 2 == 1:
                ind = 1
                splitted_line = line.split("'")
                while True:
                    splitted_line[ind] = splitted_line[ind].replace(delim, ' ')
                    ind += 2
                    if ind > len(splitted_line) - 1:
                        break
                line = ''.join(splitted_line)

            args = line.split(delim)
            args[-1] = args[-1].strip()
            if len(header) == len(args):
                line_data = {k.replace('"', '').replace("'", ""):v.replace('"', '').replace("'", "") for (k, v) in zip(header, args)}
                if len(numeric_keys) == 0:
                    for k in line_data:
                        if '-' in line_data[k] or '_' in line_data[k]:
                            if 'E-' in line_data[k] or 'e-' in line_data[k]:
                                pass
                            else:
                                continue
                        try:
                            line_data[k] = int(line_data[k])
                            numeric_keys[k] = None
                        except ValueError:
                            try:
                                line_data[k] = float(line_data[k])
                                numeric_keys[k] = None
                            except ValueError:
                                if line_data[k] in ['', '-', 'None', 'N/A']:
                                    numeric_keys[k] = None
                                continue

                else:
                    for k in numeric_keys:
                        try:
                            line_data[k] = int(line_data[k])
                        except:
                            try:
                                line_data[k] = float(line_data[k])
                            except:
                                continue

                data.append(line_data)

    return header, data


def read_header(in_file: str, delim: str = ',') -> TableHeader:
    # noinspection PyBroadException
    data = []
    with open(in_file) as f:
        line = f.readline()
        f.close()

        header = line.strip().split(delim)
        for i in range(len(header)):
            header[i] = header[i].replace('"', '')
            header[i] = header[i].replace("'", '')

    return header


def write_file(handle: io.TextIOWrapper, header: TableHeader, data: TableData, delim: str=',') -> None:
    handle.write(delim.join(header) + '\n')
    handle.writelines([delim.join([str(x[h]) for h in header]) + '\n' for x in data])


def write_file_no_header(handle: io.TextIOWrapper, header: TableHeader, data: TableData, delim: str=',') -> None:
    handle.writelines([delim.join([str(x[h]) for h in header]) + '\n' for x in data])


def write_file_uncond(out_file: str, header: TableHeader, data: TableData, delim: str=',') -> None:
    with open(out_file, 'w') as g:
        g.write(delim.join(header) + '\n')
        g.writelines([delim.join([str(x[h]) for h in header]) + '\n' for x in data])


def read_file_new_tsv(in_file: str, delim: str = '\t', header: List[str] = None) -> Tuple[TableHeader, TableData]:
    # noinspection PyBroadException
    data = []
    with open(in_file) as f:
        lines = f.readlines()
        f.close()

        if header is None:
            header = lines[0].strip().split(delim)
            for i in range(len(header)):
                header[i] = header[i].replace('"', '')
                header[i] = header[i].replace("'", '')
            data_start = 1
        else:
            data_start = 0

        for line in lines[data_start:]:
            args = line.strip().split(delim)
            if len(header) == len(args):
                line_data = dict()
                for i in range(len(header)):
                    if args[i] == 'NULL':
                        line_data[header[i]] = 'None'
                    else:
                        if len(re.findall('[-_]', args[i].replace("'", '').replace('"', ""))) > 0:
                            line_data[header[i]] = args[i].replace("'", '').replace('"', "")
                        else:
                            try:
                                line_data[header[i]] = ast.literal_eval(args[i].replace("'", '').replace('"', ""))
                            except ValueError:
                                line_data[header[i]] = args[i].replace("'", '').replace('"', "")
                            except SyntaxError:
                                line_data[header[i]] = args[i].replace("'", '').replace('"', "")
                data.append(line_data)
    return header, data


def write_file_tsv(out_file: str, header: TableHeader, data: TableData) -> None:
    """

    :rtype: object
    """
    if os.path.exists(out_file):
        print('File already exist!')
        return
    with open(out_file, 'w') as g:
        g.write('\t'.join(header) + '\n')
        g.writelines(['\t'.join([str(x[h]) for h in header]) + '\n' for x in data])


def split_file(in_file: str, lines_for_file: int, delim: str = ','):
    header, data = read_file(in_file, delim)

    for i in range(int((len(data) - 1) / lines_for_file) + 1):
        split_data = data[i * lines_for_file:(i + 1) * lines_for_file]
        with open(in_file.replace('.csv', '_%d.csv' % (i + 1)), 'w') as f:
            write_file(f, header, split_data)


def union_files_in_directories(target_dirs: List[str], output_dir: str):
    map_files = {}
    for td in target_dirs:
        files = glob.glob(os.path.join(td, '*.*'))
        for f in files:
            fname = os.path.basename(f)
            if fname not in map_files:
                map_files[fname] = [f, ]
            else:
                map_files[fname].append(f)

    for fname, files_to_union in map_files.items():
        with open(os.path.join(output_dir, fname), 'w') as outfile:
            for fname in files_to_union:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)


def fastq_batch_generator(handle: io.TextIOWrapper, batch_size: int):
    """Returns sequence lists of length batch_size.

    This can be used on any handle, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a handle function, and it returns lists of the
    entries from the supplied handle.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = [handle.readline() for _ in range(4)]
            except StopIteration:
                entry = None
                break
            if entry[0] == '':
                # End of file
                entry = None
                break
            batch.append(entry)
        if batch:
            yield batch


def csv_to_excel(in_files, out_file, sheet_names=None, format={}, **kwargs):
    """
    Combine multiple (or single) .csv files into a .xlsx file.
    
    :param in_files: csv files to be combined into one xlsx file
    :param out_file: output xlsx file full path
    :param sheet_names: sheet name for each csv file. If None (default), sheet names will be set as corresponding file name.  
    :return: 
    """

    if type(in_files) is str:
        in_files = [in_files, ]
    if type(sheet_names) is str:
        sheet_names = [sheet_names, ]
    if sheet_names == None:
        sheet_names = [os.path.basename(f) for f in in_files]
    column_settings = [None for _ in in_files]
    if 'column_settings' in kwargs:
        column_settings = kwargs['column_settings']
        if len(column_settings) != len(sheet_names):
            column_settings = [None for _ in in_files]
    highlight_rows = [None for _ in in_files]
    if 'highlight_rows' in kwargs:
        highlight_rows = kwargs['highlight_rows']
        if len(highlight_rows) != len(sheet_names):
            highlight_rows = [None for _ in in_files]

    workbook = Workbook(out_file, {'constant_memory': True})
    workbook.use_zip64()
    fm = workbook.add_format(format)
    fm_highlight = workbook.add_format({**format, "bg_color": "#92D050"})
    for f, sn, cs_list, hr_list in zip(in_files, sheet_names, column_settings, highlight_rows):
        worksheet = workbook.add_worksheet(sn)
        if cs_list != None:
            for cs in cs_list:
                worksheet.set_column(*cs)
        with open(f, 'rt') as f:
            reader = csv.reader(f)
            for r, row in enumerate(reader):
                for c, col in enumerate(row):
                    try:
                        val = ast.literal_eval(col.replace("'", '').replace('"', ""))
                    except ValueError:
                        val = col.replace("'", '').replace('"', "")
                    except SyntaxError:
                        val = col.replace("'", '').replace('"', "")

                    cell_format = get_format(col)
                    if hr_list != None and r in hr_list:
                        worksheet.write(r, c, detach_format(col), fm_highlight)
                    elif len(cell_format) > 0:
                        cell_fm = workbook.add_format({**format, **cell_format})
                        worksheet.write(r, c, detach_format(col), cell_fm)
                    else:
                        worksheet.write(r, c, val, fm)

    workbook.close()


def detach_format(s: str) -> str:
    d = s
    for m in re.finditer(r'{[^{}]*:[^{}]*}', s):
        d = d.replace(m.group(0), '')
    return d


def get_format(s: str) -> str:
    f = {}
    for m in re.finditer(r'{([^{}]*):([^{}]*)}', s):
        try:
            val = ast.literal_eval(m.group(2).replace("'", '').replace('"', ""))
        except ValueError:
            val = m.group(2).replace("'", '').replace('"', "")
        except SyntaxError:
            val = m.group(2).replace("'", '').replace('"', "")
        f[m.group(1)] = val
    return f
