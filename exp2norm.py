import argparse
import os
import pandas as pd

import numpy as np
import scipy.stats as stats
import statsmodels.api as sm

import sys


def args_parse():
    parser = argparse.ArgumentParser( #prog=__file__,
                                     usage='isoform expression transform to norm for QTL',
                                     description='',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='''

         ''')
    parser.add_argument('-i','--isoform',
                        dest='isoform', metavar='',
                        required=True,
                        help=('isoform expression file'))
    parser.add_argument('--ref',dest='ref', choices=('refseq_38','gencode_38'),
                        help='reference transcript,default=gencode_38', default='gencode_38')
    parser.add_argument('--covariates', dest='covariates',
                         help='covariates file') 
    parser.add_argument('--isoform-ratio',dest='isoform_ratio', action='store_true',help='calcualte the isoform expression splice ratio') 
    parser.add_argument('--prefix',dest='prefix',
                        help='prefix name, recommnd tissue name',default='ExpNorm')   
    parser.add_argument('-o','--outdir',
                        dest='outdir', metavar='')
    parser.add_argument('--tissue',dest='tissue',
                        help='GTEx tissue name')
    args = parser.parse_args()
    return args



def filtered_isoform(df_exp, df_anno):
    # # 过滤掉在所有表达量值都为0的isoform
    # noexp = ((self.df_exp == 0).sum(axis=1) == self.df_exp.shape[1])

    # tpm < 0.1 在20%样本中存在
    noexp = ((df_exp >= 0.1).sum(axis=1) <= df_exp.shape[1]*0.2)         
    print(f'There are {noexp.sum()} isoforms are exclued due to the tpm >= 0.1 less than 20% samples')


    df_exp  = df_exp[~noexp]

    # 过滤掉一个基因一个isoform的情况
    df_exp = df_anno[['gene_id']].merge(
                                df_exp,
                                left_index=True,
                                right_index=True)

    ## gene
    df_gene_exp = df_exp.groupby('gene_id').sum()
    print(f'There are {df_gene_exp.shape[0]} genes for {df_gene_exp.shape[1]} samples')

    
    ## isoform
    muti_isf_gene = [idx for idx,val in df_exp.groupby('gene_id').size().items() if val>1]
    df_iso_exp = df_exp[df_exp['gene_id'].isin(muti_isf_gene)]
    print(f'There are {df_iso_exp.shape[0]} isoforms for {df_gene_exp.shape[1]} samples')

    return df_gene_exp, df_iso_exp

class CallNorm(object):
    def __init__(self, exp, cov,isratio):
        self.df_exp = exp
        self.df_cov = cov
        self.isratio = isratio
        self.main()

    def zscore(self, x):
        x = pd.Series(x)
        return stats.norm.ppf((x.rank()-0.5)/(~pd.isnull(x)).sum())
    
    def norm(self, df):
        df_in = df.copy().iloc[:,1:]
        df_in.fillna(0,inplace=True)
        df_in.loc[:,:] = [i for i in map(self.zscore, df_in.values)]
        return df_in.T

    def splice_ratio(self):
        df_splice_ratio = (self.df_exp.iloc[:,1:]/self.df_exp.iloc[:,1:].groupby(self.df_exp['gene_id']).transform('sum'))
        return self.norm(df_splice_ratio)
    
    def lm_covariates(self,exp_mtx, covariates, phenos=None,covs=None):
        #exp_mtx矩阵sample,isoform/gene
        #covariates矩阵sample,covariate
        #返回 表型矩阵sample,phenoratio
        if not phenos:
            phenos = list(exp_mtx.columns)
        if not covs:
            covs = list(covariates.columns)
        
        #合并，目的将两个表格数据按照样本对应，并去除不存在covs的样本
        print(f'There are {exp_mtx.shape[0]} samples in expresstion matrix')
        print(f'There are {covariates.shape[0]} samples in covariates matrix')
        exp_mtx = exp_mtx.merge(covariates,
                                      left_index=True,
                                      right_index=True)
            
        print(f'There are {exp_mtx.shape[0]} samples include after combine covariates information')

        #lm,y为isoform x Nsamples，x为covs x Nsamples
        resid_infos = [sm.OLS(np.array(exp_mtx[i]), 
                        sm.add_constant(np.array(exp_mtx[covs]))).fit().resid 
                 for i in phenos]
        #resid_infos = [i for i in parallel_applymap(lambda y:sm.OLS(y, sm.add_constant(np.array(splice_rate[covs]))).fit().resid,
        #                                              splice_rate[isoforms].T.values)]
            
        df_phe = pd.DataFrame(resid_infos,index=phenos,columns=exp_mtx.index).T
        return df_phe
    
    def main(self):
        if self.isratio:
            print('Starting splice ratio and norm...')
            self.df_norm = self.splice_ratio()
        else:
            print('Starting norm ...')
            self.df_norm = self.norm(self.df_exp)

        print('Sarting lm ...')
        self.df_pheo = self.lm_covariates(self.df_norm,self.df_cov)
        print('Job finished')



args = args_parse()

bindir = os.path.split(os.path.realpath(sys.argv[0]))[0]

refdb = args.ref
isoform_fi = args.isoform
tissue = args.tissue
isratio =  args.isoform_ratio
if isratio:
    print ("The --isoform-ratio para set to True, this will calcualte each isoform ratio as phenotype")
else:
    print ("The --isoform-ratio para set to False, this will use each isoform abundance as phenotype")

# out dir 
if  args.outdir:
    outdir = os.path.abspath(args.outdir)
else:
    outdir = os.path.abspath(os.path.dirname(isoform_fi))

# make path
for item in [outdir]:
    if not os.path.exists(item): os.mkdir(item)


if args.prefix:
    prefix = args.prefix
else:
    if args.tissue:
        prefix = tissue
    else:
        prefix = ''

if args.covariates:
    cov_fi = args.covariates
else:
    if args.tissue:
        cov_fi = f'{bindir}/data/covariates/GTEx_Analysis_v8_sQTL_covariates_forProject/{tissue}.v8.sqtl_covariates.txt.gz'
    else:
        sys.exit('covariates file missed,please set --covariates or set --tissue in GTEx_v8')

    outfi = f'{outdir}/{prefix}.isoform_splice_ratio.tsv'
    print(f'{prefix} output isoform splice ratio file is: {outfi}')
#表达量数据 XAEM
print(f'input isoform data: {isoform_fi}')

if isoform_fi.endswith('RData'):
    os.system(f'Rscript {bindir}/isoform_rdata2exp.R inRdata={isoform_fi}')
    isoform_fi = isoform_fi.replace(".RData", "_tpm.tsv")
    print(f'Trans isoforms data to: {isoform_fi}')

df_exp = pd.read_csv(isoform_fi,sep='\t',index_col=0)
#注释数据 来自gencode gtf
df_anno = pd.read_csv(f'{bindir}/ref/{refdb}/transcript_gene_info.tsv.gz',sep='\t').set_index('transcript_id')
# covariates
df_cov = pd.read_csv(cov_fi,sep='\t',index_col=0)



def write_file(df,outf):
    ## 避免plink在处理样本的时候报错
    print (f'Starting to write {outf}')
    df.index = ['-'.join(i.split('-')[0:2]) for i in df.index]
    df.index.name = 'IID'
    df.to_csv(outf,sep='\t')
    print (f'Finished write {outf}')


## 过滤后得到gene与isoform 的表达量
df_gene, df_iso = filtered_isoform(df_exp,df_anno)

## isoform norm and lm with covariates
### isratio true: calculated isoform splice ratio
### isratio False: get isoform abundance
res_isf = CallNorm(df_iso,df_cov.T,isratio=isratio)

if isratio:
    outfi_iso = f'{outdir}/{prefix}.isoform_splice_ratio.tsv'
    print(f'{prefix} output isoform splice ratio file is: {outfi_iso}')
else:
    outfi_iso = f'{outdir}/{prefix}.isoform_abundance.tsv'
    print(f'{prefix} output isoform isoform_abundance file is: {outfi_iso}')

write_file(res_isf.df_pheo, outfi_iso)


## gene norm and  lm with covariates
outfi_gene = f'{outdir}/{prefix}.gene_abundance.tsv'
print(f'{prefix} output gene isoform_abundance file is: {outfi_gene}')
if os.path.exists(outfi_gene):
    print(f'Warings: The {outfi_gene} eixits, ignore writing')
else:
    res_gene = CallNorm(df_gene,df_cov.T,isratio=False)
    write_file(res_gene.df_pheo, outfi_gene)
