osca=/path/to/osca
isoform_ratio_file=/path/to/isoform_ratio_file
mkdir BOD_files
$osca --efile $isoform_ratio_file --gene-expression --make-bod --no-fid --out BOD_files/Muscle_Skeletal
$osca --befile BOD_files/Muscle_Skeletal --update-opi ../ref/gencode_38/anno_gene_info.opi
