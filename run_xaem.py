#!/usr/bin/env python

import pandas as pd
import os
import sys
import argparse
import logging
import datetime
from time import asctime
import configparser
import subprocess
import threading
import queue


def args_parse():
    parser = argparse.ArgumentParser( #prog=__file__,
                                     usage='Auto pipline for XAEM',
                                     description='',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='''
    Note:
         1. index trascript ref  
         2. zcat multi source fastq data to one fastq and rename to samle name
         3. generate eqclass for each sample
         4. Create count matrix and AEM update X_beta for multiple samples in sample dir

         ''')
    parser.add_argument('-i','--infile',
                        dest='infile', metavar='',
                        required=True,
                        help=('File for sample information, including '
                              "(sample name\\tdata source name\\tfq1\\tfq2)"))
    parser.add_argument('-t','--transcript',metavar='',dest='transcript')                 
    parser.add_argument('-c','--config',metavar='',
                        dest='config')
    parser.add_argument('-f','--force',action='store_true',
                         dest='force')     
    parser.add_argument('-o','--outdir',
                        dest='outdir', metavar='',
                        default='./workdir')
    parser.add_argument('--xaem-dir',
                        dest='xaem_dir')
    parser.add_argument('--xaem-index',
                        dest='xaem_index')
    parser.add_argument('--x-matrix',
                        dest='x_matrix')

    args = parser.parse_args()
    return args


def checkfq(fqfile):
    suffixs = ('fq.gz','fq','fastq','fastq.gz','fa.gz','fa','fasta','fasta.gz')
    if not fqfile.endswith(suffixs):
        logging.error(f'This is not an fq/fa file,suffix should be {";".join(suffixs)}')


def read_sampleinfo(args):
    infi = args.infile
    sample_info = []
    with open(infi) as f:
        for line in f:
            if line.strip() == '' or  line.startswith('#'):
                continue
            sample, lib, fq1, fq2  = line.strip().split()[0:4]
            sample_info.append([sample, lib, os.path.abspath(fq1), os.path.abspath(fq2)])
    df = pd.DataFrame(sample_info,columns = ['sample','lib','fq1','fq2'])
    return df


def index_ref(step=1):
    trascript = config.get('xaem','transcript_fa') if config.get('xaem','transcript_fa') else f'{bindir}/ref/refseq_38/transcript.fa.gz'
    logging.warning(f'Parameter: transcript_fa  is {trascript}')
    cmd = ''
    outfa = f'{outdir}/ref/{os.path.basename(trascript).replace(".gz","")}'
    index_dir =  f'{outdir}/ref/TxIndexer_idx'

    if trascript.endswith('gz'):
        cmd += f'gunzip -c {trascript} >{outfa}\n'
    else:
        cmd += f'ln -fs  {trascript} {outfa}\n'
    
    cmd += f'{xaem_dir}/bin/TxIndexer -t {outfa} --out {index_dir} \n'
    outf = open(f'{shelldir}/Step{step}.index_fa.sh','w')
    outf.write(cmd)
    return f'{index_dir}', f'{shelldir}/Step{step}.index_fa.sh'


def get_eqclass(df,step=2,TxIndexer_idx=None):
    shell_lst = []
    for sample,val in df.groupby('sample'):
        cmd = ''
        fq1_lst = list(val['fq1'])
        fq2_lst = list(val['fq2'])
        sample_fq1 = f'{seqdir}/{sample}_1.fq.gz'
        sample_fq2 = f'{seqdir}/{sample}_2.fq.gz'

        if len(fq1_lst) == 1:
            cmd += f'ln -fs {fq1_lst[0]} {sample_fq1}\n'
            cmd += f'ln -fs {fq1_lst[0]} {sample_fq2}\n'
        elif len(fq1_lst) > 1:
            cmd += f'zcat {" ".join(fq1_lst)} | gzip -cf > {sample_fq1} &\n'
            cmd += f'zcat {" ".join(fq2_lst)} | gzip -cf > {sample_fq2} &\n'
            cmd += 'wait \n'
        cmd += f"""{xaem_dir}/bin/XAEM \\
        -i {TxIndexer_idx} \\
        -l IU \\
        -1 <(gunzip -c {sample_fq1}) \\
        -2 <(gunzip -c {sample_fq2}) \\
        -p 2 \\
        -o {outdir}/results/{sample}\n"""


        outf = open(f'{shelldir}/Step{step}.gen_eqclass_{sample}.sh','w')
        outf.write(cmd)

        shell_lst.append(f'{shelldir}/Step{step}.gen_eqclass_{sample}.sh')
    return shell_lst


def count_matrix(step=3):
    x_matrix = config.get('xaem','x_matrix') if config.get('xaem','x_matrix') else f'{bindir}/ref/refseq_38/X_matrix.RData'
    logging.warning(f'Parameter: x_matrix is {x_matrix}')
    cmd = f"Rscript {xaem_dir}/R/Create_count_matrix.R workdir={outdir}/results core=8 design.matrix={x_matrix}\n"
    cmd += f"""Rscript {xaem_dir}/R/AEM_update_X_beta.R \\
        workdir={resdir} \\
        core={config.getint('xaem','update_cpu')} \\
        design.matrix={x_matrix} \\
        merge.paralogs={config.getboolean('xaem','merge.paralogs')} \\
        isoform.method={config.get('xaem','isoform.method')} \\
        remove.ycount={config.getboolean('xaem','remove.ycount')}"""
    outf = open(f'{shelldir}/Step{step}.matrix_samples.sh','w')
    
    outf.write(cmd)
    return f'{shelldir}/Step{step}.matrix_samples.sh'


def get_all_shells():
    """
    write all shells and initialize the jobs status
    """
    shell_info = []
    step_n = 1

    ## index
    if args.xaem_index:
        TxIndexer_idx =  os.path.abspath(args.xaem_index)
    else:
        TxIndexer_idx,cmd = index_ref(step = step_n)
        shell_info.append([cmd,step_n,'index'])
        step_n = step_n + 1

    ## eqclass
    shell_info.extend([[i,step_n,'eqclass'] for i in get_eqclass(df_sample,step=step_n,TxIndexer_idx=TxIndexer_idx)])
    step_n = step_n + 1

    ## count matrix and update beta
    cmd  = count_matrix(step=step_n)
    shell_info.append([cmd,step_n,'matrix'])
     
    ## initialize job status
    df_shell_info =  pd.DataFrame(shell_info)
    df_shell_info['status'] = 'Ready'
    df_shell_info.columns =['shell','step','name','status']
    df_shell_info.to_csv(status_fi,index=False,sep='|')


class myThread(threading.Thread):
    def __init__(self, q, df_status):
        threading.Thread.__init__(self)
        self.q = q
        self.df_status = df_status

    def run(self):
        while True:
            try:
                cmd = self.q.get(timeout=2)
                job = Job_status(cmd,self.df_status)
                self.df_status = job.df_status
            except:
                break

class Job_status(object):
    def __init__(self, sh_nm, df_status):
        self.sh_nm = sh_nm
        self.df_status = df_status
        self.run()

    def run(self):
        logger.warning(f"Running: {self.sh_nm} ")
        res = subprocess.run(f'bash  {self.sh_nm} 1>{self.sh_nm}.stdout 2>{self.sh_nm}.stderr',
                            shell=True)
        if res.returncode == 0:
            logger.warning(f"Success: {self.sh_nm} ")
            self.change_status('Success')
        else:
            logger.warning(f"Error: Place check {self.sh_nm}.stderr")
            self.change_status('Error')
    
    def change_status(self, staus_item):
        threadLock.acquire()
        self.df_status.loc[self.df_status['shell'] == self.sh_nm,'status'] = staus_item
        threadLock.release()

def is_success(df, step_name):
    df_step = df[df['name'] == step_name]
    if df_step.shape[0] == (df_step['status'] == 'Success').sum():
        return True
    else:
        return False


threadLock = threading.Lock()
FORMAT = '%(asctime)0s %(message)s'
log_file = f'{datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")}.sQTL_XAEM.info.log'
logging.basicConfig(format=FORMAT,filename=log_file)
print(f'Project starting\nFor logging information,please check {log_file}')
logger = logging.getLogger()


args = args_parse()
df_sample = read_sampleinfo(args)
bindir = os.path.split(os.path.realpath(sys.argv[0]))[0]
outdir = os.path.abspath(args.outdir)
seqdir = f'{outdir}/seqData'
resdir = f'{outdir}/results'
shelldir = f'{outdir}/shell'
refdir = f'{outdir}/ref'


for item in [outdir,seqdir,resdir,shelldir,refdir]:
    if not os.path.exists(item): os.mkdir(item)

config = configparser.ConfigParser()
if args.config:
    config.read(args.config,encoding="utf-8")
else:
    config.read(f'{bindir}/config.ini', encoding="utf-8")

for each in config['xaem']:
    logging.warning(f'Parameter: {each} is {config.get("xaem",each)}')

xaem_dir = config.get('xaem','xaem_dir') if config.get('xaem','xaem_dir')  else f'{bindir}/XAEM/XAEM-binary-0.1.1-cq' 
logging.warning(f'Parameter: xaem_dir is {xaem_dir}')

status_fi = f'{shelldir}/JOB.Status'

logging.warning('Project starting')
if not os.path.exists(status_fi):
    get_all_shells()
    logger.warning(f'Finished generate all jobs in {shelldir}')
else:
    if args.force:
        logger.warning('Force to restart all jobs')
        get_all_shells()
        logger.warning(f'Finished generate all jobs in {shelldir}')
    else:
        logger.warning(f'{status_fi} exists. It is will continue to run unfinished jobs. if you want to restart the project, pleace add --force parameter')

df_status = pd.read_csv(status_fi,sep='|')
status_dict = df_status.groupby('status').size().to_dict()
logger.warning(f'There are {df_status.shape[0]} jobs, {status_dict}')

def write_status(df):
    df.to_csv(status_fi,sep='|',index=False)
    df[df['status'] == 'Error'].to_csv(status_fi + '.Error',sep='|',index=False)

def single_job_run(sh_nm, df):
    job = Job_status(sh_nm,df)
    df_status = job.df_status
    write_status(df_status)
    return df_status

### step1 index reference, not required
if df_status[df_status['name'] == 'index'].shape[0] == 1:
    if not is_success(df_status,'index'):
        shell_lst = df_status[df_status['name'] == 'index']['shell'].to_list()
        df_status = single_job_run(shell_lst[0],df_status)

if is_success(df_status,'index'):
    logger.warning('Index reference finished, Starting eqclass')
else:
    logger.warning(f'Index Error, please check {shell_lst[0]}.stderr ')


### step 2 eqclass, multiple batch to run jobs, each batch contains some jobs，
thread_n = int(config.getint('xaem','eqclass_cpu')/2)
count = 1
while not is_success(df_status,'eqclass') and count<=2:
    df_status_eqclass = df_status.query("name == 'eqclass'")
    shell_lst = df_status_eqclass.query("status != 'Success'")['shell'].to_list()
    
    ## 
    workQueue = queue.Queue(len(shell_lst))
    for fi in shell_lst:
        workQueue.put(fi)
    threads = []
    if len(shell_lst)<thread_n:
        thread_n = len(shell_lst)
    for i in range(thread_n):
        thread = myThread(workQueue,df_status)
        thread.start()
        df_status = thread.df_status
        threads.append(thread)
    for thread in threads:
        thread.join()
    count = count + 1
    write_status(df_status)

### step3, need all eqclass finished
if is_success(df_status,'eqclass'):
    df_status_eqclass = df_status.query("name == 'eqclass'")
    success_lst = df_status_eqclass.query("status == 'Success'")['shell'].to_list()
    logger.warning(f'{len(success_lst)} eqclass finished, Starting matrix')
    shell_lst = df_status[df_status['name'] == 'matrix']['shell'].to_list()
    df_status = single_job_run(shell_lst[0],df_status)
else:
    error_lst = df_status_eqclass.query("status != 'Error'")['shell'].to_list()
    logger.warning(f'There are {len(error_lst)} error, please check {status_fi}.Error carefully !!!')

if is_success(df_status,'matrix'):
    logger.warning(f'All job finished ')
else:
    logger.warning(f'Matrix not finished. Please check {shell_lst[0]} carefully !!!')

