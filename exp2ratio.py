import argparse
import os
import pandas as pd

import numpy as np
import scipy.stats as stats
import statsmodels.api as sm

import sys


def args_parse():
    parser = argparse.ArgumentParser( #prog=__file__,
                                     usage='isoform expression transform to splice ratio for QTL',
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
    parser.add_argument('--prefix',dest='prefix',
                        help='prefix name, recommnd tissue name')   
    parser.add_argument('-o','--outdir',
                        dest='outdir', metavar='')
    parser.add_argument('--tissue',dest='tissue',
                        help='GTEx tissue name')
    args = parser.parse_args()
    return args



class CalIR(object):
    def __init__(self, exp, anno, cov):
        self.df_exp = exp
        self.df_anno = anno
        self.df_cov = cov
        self.main()
        

    def filtered(self):
        # 过滤掉在所有表达量值都为0的isoform
        noexp = ((self.df_exp == 0).sum(axis=1) == self.df_exp.shape[1])
        print(f'There are {noexp.sum()} isoforms are exclued due to expression equal to 0 in all samples')

        self.df_exp  = self.df_exp[~noexp]

        # 过滤掉一个基因一个isoform的情况
        self.df_exp = self.df_anno[['gene_id']].merge(
                                 self.df_exp,
                                 left_index=True,
                                 right_index=True)

        muti_isf_gene = [idx for idx,val in self.df_exp.groupby('gene_id').size().items() if val>1]
        self.df_exp = self.df_exp[self.df_exp['gene_id'].isin(muti_isf_gene)]

    def zscore(self, x):
        x = pd.Series(x)
        return stats.norm.ppf((x.rank()-0.5)/(~pd.isnull(x)).sum())
    
    def splce_ratio(self):
        df_splice_ratio = (self.df_exp.iloc[:,1:]/self.df_exp.iloc[:,1:].groupby(self.df_exp['gene_id']).transform('sum'))
        df_splice_ratio.fillna(0,inplace=True)
        df_splice_ratio.loc[:,:] = [i for i in map(self.zscore, df_splice_ratio.values)]
        return df_splice_ratio.T
    
    def lm_covariates(self,splice_rate, covariates, isoforms=None,covs=None):
        #splice_rate矩阵sample,isoform
        #covariates矩阵sample,covariate
        #返回 表型矩阵sample,pheno
        if not isoforms:
            isoforms = list(splice_rate.columns)
        if not covs:
            covs = list(covariates.columns)
        
        #合并，目的将两个表格数据按照样本对应，并去除不存在covs的样本
        print(f'There are {splice_rate.shape[0]} samples in expresstion matrix')
        splice_rate = splice_rate.merge(covariates,
                                      left_index=True,
                                      right_index=True)
        print(f'There are {splice_rate.shape[0]} samples include after combine covariates information')

        #lm,y为isoform x Nsamples，x为covs x Nsamples
        resid_infos = [sm.OLS(np.array(splice_rate[i]), 
                        sm.add_constant(np.array(splice_rate[covs]))).fit().resid 
                 for i in isoforms]
        
        #resid_infos = [i for i in parallel_applymap(lambda y:sm.OLS(y, sm.add_constant(np.array(splice_rate[covs]))).fit().resid,
        #                                              splice_rate[isoforms].T.values)]
            
        df_phe = pd.DataFrame(resid_infos,index=isoforms,columns=splice_rate.index).T
        return df_phe
        
    
    def main(self):
        print('Starting filter ...')
        self.filtered()
        print('Starting splice ratio ...')
        self.df_sr = self.splce_ratio()
        print('Sarting lm ...')
        self.df_pheo = self.lm_covariates(self.df_sr,self.df_cov)
        print('Job finished')

args = args_parse()

bindir = os.path.split(os.path.realpath(sys.argv[0]))[0]

refdb = args.ref
isoform_fi = args.isoform
tissue = args.tissue

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


res = CalIR(df_exp,df_anno,df_cov.T)

res.df_pheo.index = ['-'.join(i.split('-')[0:2]) for i in res.df_pheo.index]

res.df_pheo.index.name = 'IID'
outfi = f'{outdir}/{prefix}.isoform_splice_ratio.tsv'
print(f'{prefix} output isoform splice ratio file is: {outfi}')
res.df_pheo.to_csv(outfi,sep='\t')







    
    
    


        




