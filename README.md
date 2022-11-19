# aXAEM

This is the autopipeline for [XAEM](https://github.com/WenjiangDeng/XAEM). You also can find the detailed instructions from XAEM [website](https://www.meb.ki.se/sites/biostatwiki/xaem).

### Prerequisites
```
R (recommended version >= 3.5.1)
Python (recommended version >= 3.7)

```


### Installing
```
git https://github.com/ZhengCQ/aXAEM.git
```
or 
```
wget https://github.com/ZhengCQ/aXAEM/archive/refs/tags/v1.0.0_beta.zip
unzip v1.0.0_beta.zip
ln -fs aXAEM-1.0.0_beta aXAEM
```
#### R Dependencies
```
install.packages("foreach")
install.packages("parallel")
install.packages("doParallel")
```

#### python Dependencies
```
pip instrall pandas
```

### Download for the annotation reference
```
cd /path/to/aXAEM
python down_ref.py
```

### run example
```
cd /path/to/aXAEM/Example
sh run_example.sh 
```

### run project
#### creat your project directory
```
mkdir -p /path/to/project
cd /path/to/project
```
#### Preparation for input files
Demo: /path/to/aXAEM/Example/infastq_lst.tsv
```
sample1 S0001   S0001_1.fg.gz   S0001_2.fg.gz
sample1 S0002   S0002_1.fg.gz   S0002_2.fg.gz
sample2 S0003   S0003_1.fg.gz   S0003_2.fg.gz
sample2 S0004   S0004_1.fg.gz   S0004_2.fg.gz
sample3 S0005   S0005_1.fg.gz   S0005_2.fg.gz
sample3 S0006   S0006_1.fg.gz   S0006_2.fg.gz
sample4 S0007   S0007_1.fg.gz   S0007_2.fg.gz
sample5 S0008   S0008_1.fg.gz   S0008_2.fg.gz
```
#SampleName\tsource_name\tfastq read1\tfastq read2

#### Perform XAEM 
##### Default parameter
```
python /path/to/aXAEM/run_xaem.py -i infastq_lst.tsv
```

##### Custom parameter
```
cp /path/to/aXAEM/config.ini config_custom.ini
python run_xaem.py -i infastq_lst.tsv -c config_custom.ini
```

##### target output directory

```
python run_xaem.py -i infastq_lst.tsv -o our_anlysis_dir
```

