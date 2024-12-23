# iGTEx
This project describes an automated pipeline to quantify isoform expressions for the GTEx V8 RNA sequencing data. The protocol uses [XAEM](https://github.com/WenjiangDeng/XAEM), a powerful method for isoform expression estimation across multiple samples. You also can find the detailed information about XAEM [website](https://www.meb.ki.se/sites/biostatwiki/xaem) and in the published paper in [Bioinformatics](https://academic.oup.com/bioinformatics/article/36/3/805/5545974).

## Prerequisites
```
R (recommended version >= 3.5.1)
Python (recommended version >= 3.7)
```


## Installation

#### Step 1: Prepare using conda
```
conda create -n iGTEx python=3.10 r=4.2 -c conda-forge
conda activate iGTEx
pip install pandas
R -e "install.packages(c('foreach','doParallel'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
```

#### Step 2: Setup tools
Install the iGTEx_XAEM tool for this protocol via the bash command:
```
wget https://github.com/ZhengCQ/iGTEx/archive/refs/tags/iGTEx_XAEM_v0.1.2.zip
unzip iGTEx_XAEM_v0.1.2.zip
ln -fs iGTEx-iGTEx_XAEM_v0.1.2 iGTEx_XAEM
```


#### Step 3: Download the annotation reference
Run the following commands to download the reference annotating the transcripts:
```
cd /path/to/iGTEx_XAEM
python down_ref.py
```

(Optional) To down a particular reference with gencode hg38, use:
```
python down_ref.py -db gencode_38
```

## Example
An example is prepared in the project **Example** folder, executable as
```
cd /path/to/iGTEx_XAEM/Example
sh run_example.sh 
```

## Isoform Estimation using GTEx data

XAEM performs better when multiple samples of similar data type are considered. In the GTEx V8 data, for **each tissue**, we create a project directory:
```
mkdir -p /path/to/Tissue1
cd /path/to/Tissue1
```
#### Input files
In `/path/to/Tissue1`, create a file `/path/to/Tissue1/infastq_lst.tsv` listing the FASTQ input files. The file is a tab-delimited text file with 4 columns: `Sample name`, `Source name`, `FASTQ file name for paired-end read 1`, and `FASTQ file name for paired-end read 2`. `Source name` indicates the batch or sequencing library of the sample, so that the same sample may correspond to more than one sources. A standard example, where each sample has only a single batch or multiple batches, is given as `/path/to/iGTEx_XAEM/Example/infastq_lst.tsv`:

##### single batch
```
sample4 S0007   S0007_1.fg.gz   S0007_2.fg.gz
sample5 S0008   S0008_1.fg.gz   S0008_2.fg.gz
```
##### multiple batches
```
sample1 S0001   S0001_1.fg.gz   S0001_2.fg.gz
sample1 S0002   S0002_1.fg.gz   S0002_2.fg.gz
sample2 S0003   S0003_1.fg.gz   S0003_2.fg.gz
sample2 S0004   S0004_1.fg.gz   S0004_2.fg.gz
sample3 S0005   S0005_1.fg.gz   S0005_2.fg.gz
sample3 S0006   S0006_1.fg.gz   S0006_2.fg.gz
```



#### Run XAEM 
XAEM can be easily run with:
```
python /path/to/iGTEx_XAEM/run_xaem.py -i /path/to/Tissue1/infastq_lst.tsv
```
(Optional) To specify a particular reference with gencode hg38, use:
```
--ref gencode_38
```
(Optional) To specify a particular output directory, use:
```
-o /path/to/Tissue1_output_directory
```
(Optional) Further customized configuration of XAEM can be setup by:
```
-c /path/to/Tissue1_config.ini
```
An example of the `config.ini` file can be found in `/path/to/iGTEx_XAEM/`.



#### Calculate isoform ratio  
Isoform ratio can be easily run with
```
python /path/to/iGTEx_XAEM/exp2ratio.py -i /path/to/XAEM_isoform_expression.RData
```
(Optional) To specify a particular reference with gencode hg38, use:
```
--ref gencode_38
```
(Optional) The default output directory as same as input isoform file. To specify a particular output directory, use:
```
-o /path/to/Tissue1_output_directory
```
(Optional) To specify a covariates file for linear regression:
```
--covariates /path/to/covariates_file
```
(For GTEx) For GTEx v8, we have prepared the covariates information for calculation, you can specify a tissue name:
```
--tissue Brain_Amygdala 
```
(For GTEx) For GTEx v8 demo:
```
python /path/to/iGTEx_XAEM/exp2ratio.py -i /path/to/Brain_Amygdala_XAEM_isoform_expression.RData --ref gencode_38 --tissue Brain_Amygdala
```

#### sQTL THISTLE
##### BOD file
```
/path/to/osca --efile /path/to/TissueName/isoform_splice_ratio.tsv --gene-expression --make-bod --no-fid --out TissueName 
/path/to/osca --befile TissueName --update-opi /path/to/iGTEx_XAEM/ref/gencode_38/anno_gene_info.opi
```
##### sQTL
```
/path/to/osca --sqtl --bfile /path/to/Genotype/BED_All/TissueName_Genotype --befile TissueName --maf 0.05 --call 0.85 --cis-wind 1000 --thread-num 10 --task-num 1 --task-num 1 --task-id 1 --to-smr --bed /path/to/iGTEx_XAEM/ref/gencode_38/anno_gene_info.bed --out sQTL_results/TissueName
```



