import configparser
import glob
import logging
import os
import platform
import subprocess
import socket

# import pymysql

from .Enum import ChainType
from Common.FileUtil import read_file_fast, write_file


def setup_database(sqlpath: str, dbname: str):
    """Create and setup a database schema to store g_target data."""

    logging.info('Started setup_database')

    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user='biteam',
                                 password='Antibody54321',
                                 port=3307,
                                 charset='utf8mb4',
                                 max_allowed_packet=1 * 1024 * 1024 * 1024,
                                 cursorclass=pymysql.cursors.DictCursor)

    # Create DB
    try:
        with connection.cursor() as cursor:
            cursor.execute("CREATE SCHEMA %s" % dbname)
            connection.commit()
    except Exception as e:
        logging.info('Pass database setup process.: %s', str(e))
        return
    finally:
        pass

    # Reconnect to the g_target database
    connection = pymysql.connect(host='localhost',
                                 db=dbname,
                                 user='biteam',
                                 password='Antibody54321',
                                 port=3307,
                                 charset='utf8mb4',
                                 max_allowed_packet=1 * 1024 * 1024 * 1024,
                                 cursorclass=pymysql.cursors.DictCursor)

    # Run setup sql scripts
    script_files = os.path.join(sqlpath, 'setup_*.sql')
    setup_sqls = glob.glob(script_files)
    for sql in setup_sqls:
        with open(sql, 'r') as f:
            sql_commands = f.read().split(';')
            for command in sql_commands:
                try:
                    with connection.cursor() as cursor:
                        cursor.execute(command)
                        connection.commit()
                except pymysql.err.InternalError:
                    pass
                finally:
                    pass

    logging.info('Finished setup_database')


def read_config(config_file: str) -> configparser.ConfigParser:
    directory = os.path.dirname(config_file)
    with open(config_file) as cf:
        config_str = cf.read()
        config_str = config_str.replace('$logfile', "r'%s'" % os.path.join(os.path.abspath(directory), 'log.log'))
    with open(config_file, 'w') as cf:
        cf.write(config_str)

    # read config
    config = configparser.ConfigParser()
    config.read(config_file)

    return config


def run_pear(forward: str, reverse: str, output: str, threads: int = 1) -> str:
    """Merge the paired-end sequence reads using PEAR - Paired-End reAd mergeR

    PEAR url: http://sco.h-its.org/exelixis/web/software/pear/
    """
    if platform.system() == 'Windows':
        pear_cmd = r'pear -f %s -r %s -o %s -j %d' % (forward, reverse, output, threads)
        process = subprocess.Popen(pear_cmd, stdout=subprocess.PIPE)
    elif platform.system() == 'Linux':
        script_directory = os.path.dirname(os.path.abspath(__file__))
        pear_cmd = ['pear', '-f', forward, '-r', reverse, '-o', output, '-j', str(threads)]
        process = subprocess.Popen(pear_cmd, stdout=subprocess.PIPE)
        process.communicate('Antibody54321\n')
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()

    return 1


def run_clustal_omega(input: str, output: str) -> str:
    if platform.system() == 'Windows':
        clustalo_cmd = r'clustalo -i %s -o %s --force' % (input, output)
        process = subprocess.Popen(clustalo_cmd, stdout=subprocess.DEVNULL)
    elif platform.system() == 'Linux':
        # script_directory = os.path.dirname(os.path.abspath(__file__))
        # clustalo_path = os.path.join(script_directory, 'clustalo-1.2.4-Ubuntu-x86_64')
        clustalo_path = r'/Tools/clustalo'
        # clustalo_cmd = ['sudo', '-S', clustalo_path, '-i', input, '-o', output, '--force']
        clustalo_cmd = [clustalo_path, '-i', input, '-o', output, '--force']
        process = subprocess.Popen(clustalo_cmd, stdout=subprocess.DEVNULL)
        # process.communicate('Antibody54321\n')
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()


def run_igblast(chain_type: ChainType, query: str, out: str, seq_type: str = 'Ig', domain_system: str = 'imgt', **kwargs):
    igblast_db_basepath = '/Tools/ncbi-igblast-1.8.0/bin/database'
    igblast_human_db_path = os.path.join(igblast_db_basepath, 'human')
    igblast_chicken_db_path = os.path.join(igblast_db_basepath, 'chicken')
    igblast_rat_db_path = os.path.join(igblast_db_basepath, 'rat')
    igblast_rabbit_db_path = os.path.join(igblast_db_basepath, 'rabbit')
    igblast_mouse_db_path = os.path.join(igblast_db_basepath, 'mouse')

    if chain_type in [ChainType.HUMAN_HEAVY, ChainType.HUMAN_LIGHT,
                      ChainType.HUMAN_BETA, ChainType.HUMAN_ALPHA, ChainType.HUMAN_DELTA, ChainType.HUMAN_GAMMA]:
        if seq_type == 'TCR':
            germline_db_V = os.path.join(igblast_human_db_path, r'imgt_human_trv_functional')
            germline_db_J = os.path.join(igblast_human_db_path, r'imgt_human_trj_functional')
            germline_db_D = os.path.join(igblast_human_db_path, r'imgt_human_trd_functional')
        else:
            germline_db_V = os.path.join(igblast_human_db_path, r'imgt_human_igv_functional')
            germline_db_J = os.path.join(igblast_human_db_path, r'imgt_human_igj_functional')
            germline_db_D = os.path.join(igblast_human_db_path, r'imgt_human_igd_functional')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s -ig_seqtype %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/human_gl.aux',
                                       query, out, domain_system, seq_type) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism human"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
    elif chain_type in [ChainType.CHICKEN_HEAVY, ChainType.CHICKEN_LIGHT]:
        germline_db_V = os.path.join(igblast_chicken_db_path, r'imgt_chicken_v')
        germline_db_J = os.path.join(igblast_chicken_db_path, r'imgt_chicken_j')
        germline_db_D = os.path.join(igblast_chicken_db_path, r'imgt_chicken_d')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/chicken_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism chicken"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
    elif chain_type in [ChainType.RAT_HEAVY, ChainType.RAT_LIGHT]:
        germline_db_V = os.path.join(igblast_rat_db_path, r'imgt_rattus_novegicus_igv')
        germline_db_J = os.path.join(igblast_rat_db_path, r'imgt_rattus_novegicus_igj')
        germline_db_D = os.path.join(igblast_rat_db_path, r'imgt_rattus_novegicus_igd')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/rat_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism rat"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
    elif chain_type in [ChainType.RABBIT_HEAVY, ChainType.RABBIT_KAPPA]:
        germline_db_V = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igv_functional')
        germline_db_J = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igj_functional')
        germline_db_D = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igd_functional')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/rabbit_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism rabbit"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
    elif chain_type in [ChainType.MOUSE_BALBC_HEAVY, ChainType.MOUSE_BALBC_LIGHT, ChainType.MOUSE_C57BL6_HEAVY, ChainType.MOUSE_C57BL6_LIGHT,
                        ChainType.MOUSE_HEAVY, ChainType.MOUSE_LIGHT,
                        ChainType.MOUSE_ALPHA, ChainType.MOUSE_BETA, ChainType.MOUSE_GAMMA, ChainType.MOUSE_DELTA]:
        if seq_type == 'TCR':
            germline_db_V = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trv_functional')
            germline_db_J = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trj_functional')
            germline_db_D = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trd_functional')
        elif seq_type == 'Ig':
            germline_db_V = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_igv_functional')
            germline_db_J = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_igj_functional')
            germline_db_D = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_igd_functional')


        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s -ig_seqtype %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/mouse_gl.aux',
                                       query, out, domain_system, seq_type) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism mouse"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
        # for the test
        cmd += ' -evalue 1'

    if platform.system() == 'Windows':
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    elif platform.system() == 'Linux':
        igblast_path = os.path.join('/Tools', 'ncbi-igblast-1.8.0', 'bin', 'igblastn')
        linux_cmd = [igblast_path] + cmd.split(' ')[1:]
        process = subprocess.Popen(linux_cmd, stdout=subprocess.DEVNULL)
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()


def run_igblast_new(chain_type: ChainType, query: str, out: str, seq_type: str = 'Ig', domain_system: str = 'kabat', **kwargs):
    hostname = socket.gethostname()
    if hostname == 'JunhoLab':
        igblast_db_basepath = '/Tools/ncbi-igblast-1.17.1/database'
    elif hostname == 'binel229':
        igblast_db_basepath = '/home/team/IP-team/Tools/ncbi-igblast-1.17.1/database'
    igblast_human_db_path = os.path.join(igblast_db_basepath, 'human')
    igblast_chicken_db_path = os.path.join(igblast_db_basepath, 'chicken')
    igblast_rat_db_path = os.path.join(igblast_db_basepath, 'rat')
    igblast_rabbit_db_path = os.path.join(igblast_db_basepath, 'rabbit')
    igblast_mouse_db_path = os.path.join(igblast_db_basepath, 'mouse')
    # igblast_alpaca_db_path = os.path.join(igblast_db_basepath, 'alpaca')

    if chain_type in [ChainType.HUMAN_HEAVY, ChainType.HUMAN_LIGHT,
                      ChainType.HUMAN_BETA, ChainType.HUMAN_ALPHA, ChainType.HUMAN_DELTA, ChainType.HUMAN_GAMMA]:
        if seq_type == 'TCR':
            germline_db_V = os.path.join(igblast_human_db_path, r'imgt_human_trv_functional')
            germline_db_J = os.path.join(igblast_human_db_path, r'imgt_human_trj_functional')
            germline_db_D = os.path.join(igblast_human_db_path, r'imgt_human_trd_functional')
        else:
            germline_db_V = os.path.join(igblast_human_db_path, r'imgt_human_igv_functional')
            germline_db_J = os.path.join(igblast_human_db_path, r'imgt_human_igj_functional')
            germline_db_D = os.path.join(igblast_human_db_path, r'imgt_human_igd_functional')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s -ig_seqtype %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/human_gl.aux',
                                       query, out, domain_system, seq_type) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism human"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
    elif chain_type in [ChainType.CHICKEN_HEAVY, ChainType.CHICKEN_LIGHT]:
        germline_db_V = os.path.join(igblast_chicken_db_path, r'imgt_chicken_v')
        germline_db_J = os.path.join(igblast_chicken_db_path, r'imgt_chicken_j')
        germline_db_D = os.path.join(igblast_chicken_db_path, r'imgt_chicken_d')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/chicken_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism chicken"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
    elif chain_type in [ChainType.RAT_HEAVY, ChainType.RAT_LIGHT]:
        germline_db_V = os.path.join(igblast_rat_db_path, r'imgt_rattus_novegicus_igv')
        germline_db_J = os.path.join(igblast_rat_db_path, r'imgt_rattus_novegicus_igj')
        germline_db_D = os.path.join(igblast_rat_db_path, r'imgt_rattus_novegicus_igd')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/rat_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism rat"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
    elif chain_type in [ChainType.RABBIT_HEAVY, ChainType.RABBIT_KAPPA]:
        germline_db_V = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igv_functional')
        germline_db_J = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igj_functional')
        germline_db_D = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igd_functional')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/rabbit_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism rabbit"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))
    elif chain_type in [ChainType.MOUSE_BALBC_HEAVY, ChainType.MOUSE_BALBC_LIGHT, ChainType.MOUSE_C57BL6_HEAVY, ChainType.MOUSE_C57BL6_LIGHT,
                        ChainType.MOUSE_HEAVY, ChainType.MOUSE_LIGHT,
                        ChainType.MOUSE_ALPHA, ChainType.MOUSE_BETA, ChainType.MOUSE_GAMMA, ChainType.MOUSE_DELTA]:
        if seq_type == 'TCR':
            germline_db_V = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trv_functional')
            germline_db_J = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trj_functional')
            germline_db_D = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trd_functional')
        elif seq_type == 'Ig':
            germline_db_V = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_igv_functional')
            germline_db_J = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_igj_functional')
            germline_db_D = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_igd_functional')

        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s -ig_seqtype %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/mouse_gl.aux',
                                       query, out, domain_system, seq_type) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism mouse"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    if platform.system() == 'Windows':
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    elif platform.system() == 'Linux':
        if hostname == 'JunhoLab':
            igblast_path = os.path.join('/Tools', 'ncbi-igblast-1.17.1', 'bin', 'igblastn')
        elif hostname == 'binel229':
            igblast_path = os.path.join('/home', 'team', 'IP-team', 'Tools', 'ncbi-igblast-1.17.1', 'bin', 'igblastn')
        linux_cmd = [igblast_path] + cmd.split(' ')[1:]
        process = subprocess.Popen(linux_cmd, stdout=subprocess.DEVNULL)
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()


def run_IMPre(query: str, seq_type: str = 'Ig', chain: str = 'heavy', **kwargs):
    IMPre_outPath = '/Tools/IMPre-master/Out_temp'
    if 'vn' in kwargs:
        vn = kwargs['vn']
    else:
        if seq_type == 'Ig':
            vn = 300
        else:
            vn = 200
    if 'v_seed' in kwargs:
        v_seed = kwargs['v_seed']
    else:
        v_seed = 200
    if 'iter_num' in kwargs:
        iter_num = kwargs['iter_num']
    else:
        iter_num = ''
    if len(str(iter_num)) > 0:
        IMPre_outPath_final = os.path.join(IMPre_outPath, '%s_vn=%d_vseed=%d_iter=%d' %(chain, vn, v_seed, iter_num))
    else:
        IMPre_outPath_final = os.path.join(IMPre_outPath, '%s_vn=%d_vseed=%d' % (chain, vn, v_seed))
    if not os.path.exists(IMPre_outPath_final):
        os.mkdir(IMPre_outPath_final)
    if seq_type == 'Ig':
        if chain == 'heavy':
            cmd_s1 = r'IMPre.pl -i %s -n %s -o %s -vm 50 -jm 65 -v_seed_m 20 -vn %d -v_seed %d -jf_ave 2'\
                     % (query, seq_type, IMPre_outPath_final, vn, v_seed)
        elif chain=='light':
            cmd_s1 = r'IMPre.pl -i %s -n %s -o %s -vm 50 -jm 65 -v_seed_m 10 -vn %d -v_seed %d -jf_ave 2'\
                     % (query, seq_type, IMPre_outPath_final, vn, v_seed)
        else:
            logging.info('Wrong chain : For Ig, chain=[heavy, light]')
    elif seq_type == 'TCR':
        if chain == 'beta':
            cmd_s1 = r'IMPre.pl -i %s -n %s -o %s -vm 40 -jm 60 -v_seed_m 20 -vn 200 -jf_ave 0.5' % (query, seq_type, IMPre_outPath_final)
        elif chain == 'alpha':
            cmd_s1 = r'IMPre.pl -i %s -n %s -o %s -vm 40 -jm 60 -v_seed_m 10 -vn 200 -jf_ave 0.5' % (query, seq_type, IMPre_outPath_final)
        else:
            logging.info('Wrong chain : For TCR, chain=[beta, alpha]')
    else:
        logging.info('Wrong seq_type : Ig(default) or TCR')

    # for key, value in kwargs.items():
    #     cmd_s1 += ' -%s %s' % (key, str(value))

    if platform.system() == 'Windows':
        win_cmd1 = 'perl ' + cmd_s1
        win_cmd2 = 'sh Execute_all.sh'
        process = subprocess.Popen(win_cmd1, stdout=subprocess.PIPE)
        process.wait()
        process = subprocess.Popen(win_cmd2, stdout=subprocess.PIPE)
    elif platform.system() == 'Linux':
        IMPre_path = os.path.join('/Tools', 'IMPre-master', 'IMPre.pl')
        linux_cmd1 = ['perl', IMPre_path] + cmd_s1.split(' ')[1:]
        IMPre_path = os.path.join(IMPre_outPath_final, 'Execute_all.sh')
        linux_cmd2 = ['sh', IMPre_path]
        process = subprocess.Popen(linux_cmd1, stdout=subprocess.DEVNULL)
        process.wait()
        process = subprocess.Popen(linux_cmd2, stdout=subprocess.DEVNULL)
    else:
        logging.info('Not supporting OS')
        return None

    logging.info('%s' % cmd_s1)

    process.wait()


def functionality_check(sample_name, targetFile, writeFile, col_rc='duplicate_count'):
    filetype = targetFile[-len('csv'):]
    if filetype == 'csv':
        delim = ','
    elif filetype == 'tsv':
        delim = '\t'
    else:
        print('Unsupported file type... must be csv or tsv')
        return

    targetHeader, targetList = read_file_fast(in_file=targetFile, delim=delim)
    writeList = []

    for target in targetList:
        # full aa stop codon filtering
        try:
            target['sequence_aa'].index('*')
            continue
        except ValueError:
            pass

        # full aa frame filtering
        try:
            target['sequence_aa'].index('X')
            continue
        except ValueError:
            pass

        # gene annotation filtering
        if target['v_call'] == 'N/A' or target['j_call'] == 'N/A':
            continue

        # cdr1,2,3 annotation flitering
        if target['cdr1_aa'] == 'N/A' or target['cdr2_aa'] == 'N/A' or target['cdr3_aa'] == 'N/A':
            continue

        # cdr1,2,3 stop codon filtering
        try:
            target['cdr1_aa'].index('*')
            continue
        except ValueError:
            pass
        try:
            target['cdr2_aa'].index('*')
            continue
        except ValueError:
            pass
        try:
            target['cdr3_aa'].index('*')
            continue
        except ValueError:
            pass

        # cdr1,2,3 frame filtering
        try:
            target['cdr1_aa'].index('X')
            continue
        except ValueError:
            pass
        try:
            target['cdr2_aa'].index('X')
            continue
        except ValueError:
            pass
        try:
            target['cdr3_aa'].index('X')
            continue
        except ValueError:
            pass

        # rev_comp filtering: must be 'F'
        if target['rev_comp'] == 'T':
            continue

        # productive filtering: must be 'T'
        if target['productive'] == 'F':
            continue

        writeList.append(target)

    # set sequence_id
    sorted_write_list = sorted(writeList, key=lambda x:x[col_rc], reverse=True)
    row_num = 1
    for d in sorted_write_list:
        d['sequence_id'] = '%s-%d' % (sample_name, row_num)
        row_num += 1

    # write file
    w_header = ['sequence_id'] + targetHeader

    with open(writeFile, 'w') as handle:
        write_file(handle=handle, header=w_header, data=sorted_write_list, delim=delim)

run_clustal_omega(input="/home/Backup3/namphilkim/lr2_s3_v.fasta", output="/home/Backup3/namphilkim/lr2_s3_v_out.fasta")