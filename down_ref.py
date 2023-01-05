
#!/usr/bin/env python
import os
import argparse
import subprocess

import gzip

def check_integrity(local_file):
    if local_file.endswith('.gz'):
        with gzip.open(local_file, 'rb') as f:
            try:
                f.seek(-1, os.SEEK_END)
            except:
                return -1
    return 0


BASE_PATH = os.path.abspath(os.path.split(os.path.abspath(__file__))[0])
DATA_PATH = os.path.join(BASE_PATH, "ref")


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

def get_transcript(db): 
    return os.path.join(DATA_PATH, f"{db}/transcript.fa.gz")

def get_x_matrix(db): 
    return os.path.join(DATA_PATH, f"{db}/X_matrix.RData")

def gunzip_flag():
    if args.force:
        return '-f'
    return ''

def download_transcript_db(data_path,url):

    cmd = (
        f'cd {data_path} && '
        f'wget -c {url} -O tmp_transcript.fa.gz'
    )
    res = subprocess.run(cmd, shell=True)
    if res.returncode == 0:
        try:
            os.system(f'cd {data_path} && mv tmp_transcript.fa.gz transcript.fa.gz')
            print(f"Success: success download {url}, and save at {data_path}/transcript.fa.gz")
        except:
            print(f'mv {data_path}/tmp_transcript.fa.gz to {data_path}/transcript.fa.gz error' )
    else:
        print(f"Error: fail download {url}, and delete at {data_path}/transcript.fa.gz")

def download_xmatrix_db(data_path,url):

    cmd = (
        f'cd {data_path} && '
        f'wget  -c {url} -O tmp_X_matrix.RData'
    )
    res = subprocess.run(cmd, shell=True)
    if res.returncode == 0:
        try:
            os.system(f'cd {data_path} && mv tmp_X_matrix.RData X_matrix.RData')
            print(f"Success: success download {url}, and save at {data_path}/X_matrix.RData")
        except:
            print(f'mv {data_path}/tmp_X_matrix.RData to {data_path}/tmp_X_matrix.RData error' )
    else:
        print(f"Error: fail download {url}, and delete at {data_path}/X_matrix.RData")


def args_parse():
    parser = argparse.ArgumentParser( usage='Download referece data for XAEM',
                                     description='',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='''
    Note:
        Download reference data from public
         '''
    )
    parser.add_argument('-f', action="store_true", dest='force',
                        help='forces download even if the files exist')

    parser.add_argument('-db', '--db',dest='db',choices=('refseq_38', 'refseq_37','gencode_38','gencode_37'),
                        help='reference transcript,default=refseq_38',default='gencode_38'
    
    )


    args = parser.parse_args()
    return args


transcript_url={
    'refseq_38':'https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna.fna.gz',
    'gencode_38':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz'
}

xmatrix_url={
    'gencode_38': 'https://github.com/ZhengCQ/iGTEx/releases/download/iGTEx_XAEM_v0.1.1/gencode_38.X_matrix.RData'
}

args = args_parse()

db = args.db

is_force = args.force

if not os.path.exists(f'{DATA_PATH}/{db}'):
    os.mkdir(f'{DATA_PATH}/{db}')
    

print(f'Starting download {db} transcript fasta...')
if args.force or not os.path.exists(get_transcript(db)) or check_integrity(get_transcript(db)) == -1:
    print(f'Downloading reference {db} data at {DATA_PATH}/{db}...')
    download_transcript_db(f'{DATA_PATH}/{db}', transcript_url[db])
else:
    print(f'reference {db} have exists at {get_transcript(db)}')

print(f'Starting download {db} x_matrix data')
if args.force or not os.path.exists(get_x_matrix(db)):
    print(f'Downloading X matrix {db} data at {DATA_PATH}/{db}...')
    download_xmatrix_db(f'{DATA_PATH}/{db}', xmatrix_url[db])
else:
    print(f'reference {db} have exists at {get_x_matrix(db)}')
