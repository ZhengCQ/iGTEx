
osca=/path/to/osca
bfile=/path/to/plink_Genotype
befile=/path/to/BODfile

$osca --sqtl \
--bfile $bfile 
--befile  $befile
--maf 0.05 
--call 0.85 
--cis-wind 1000 
—thread-num 10 
--task-num 1 --task-num 1 
--task-id 1 --to-smr 
--bed ../ref/gencode_38/anno_gene_info.bed 
—out Muscle_Skeletal