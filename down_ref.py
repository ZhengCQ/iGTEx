
#!/usr/bin/env python
import os
import argparse
import subprocess


BASE_PATH = os.path.abspath(os.path.split(os.path.abspath(__file__))[0])
DATA_PATH = os.path.join(BASE_PATH, "ref")


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

def get_transcript(db): 
    return os.path.join(DATA_PATH, f"{db}/transcript.fa.gz")

def gunzip_flag():
    if args.force:
        return '-f'
    return ''

def download_transcript_db(data_path):
    ncbi_ref_url = 'https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/110/GCF_000001405.40_GRCh38.p14/'
    fi_name  = 'GCF_000001405.40_GRCh38.p14_rna.fna.gz'

    url = ncbi_ref_url + fi_name
    cmd = (
        f'cd {data_path} && '
        f'wget  {url} -O transcript.fa.gz'
    )
    res = subprocess.run(cmd, shell=True)
    if res.returncode == 0:
        print(f"Success: success download {url}, and save at {data_path}/transcript.fa.gz")
    else:
        print(f"Error: fail download {url}, and delete at {data_path}/transcript.fa.gz")


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

    parser.add_argument('--db',dest='db',choices=('refseq_38', 'refseq_37','genocode_38','genocode_37'),
                        help='reference transcript,default=refseq_38',default='refseq_38'
    
    )


    args = parser.parse_args()
    return args




args = args_parse()

db = args.db
is_force = args.force

print(f'Starting download {db} transcript fasta...')
if args.force or not os.path.exists(get_transcript(db)):
    print(f'Downloading reference {db} data at {DATA_PATH}/{db}...')
    download_transcript_db(f'{DATA_PATH}/{db}')
else:
    print(f'reference {db} have exists at {get_transcript(db)}')

