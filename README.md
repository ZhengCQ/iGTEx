# sQTL_XAEM

### Download for the annotation reference
```
python down_ref.py
```
### Preparation for input files

#### one sample with one source fqstaq file:
|#SampleName|LibraryName|Fastq1|Fastq2|
| --- |--- |--- | --- |
| sample2 | sample2 | \*1.fataq.gz | \*2.fataq.gz|

#### one sample with multiple source fqstaq files:
|#SampleName|LibraryName|Fastq1|Fastq2|
| --- |--- |--- | --- |
| sample1 | source_name1 | \*1.fataq.gz | \*2.fataq.gz|
| sample1 | source_name2 | \*1.fataq.gz | \*2.fataq.gz|


#### fasta file

|#SampleName|LibraryName|Fasta|
| --- |--- |--- |
| sample3 | sample3 | \*.fasta.gz |

Note: header columns is not necessary，for display only. Demo file can be find here:./Example/infq.lst


### Perform XAEM 
#### Default parameter
```
python run_xaem.py -i infq.lst
```

#### Custom parameter
```
cp config.ini config_custom.ini
python run_xaem.py -i infq.lst -c config_custom.ini
```

#### target output directory

```
python run_xaem.py -i infq.lst -o our_anlysis_dir
```