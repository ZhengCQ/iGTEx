scriptdir="/path/to/R"
workdir='NULL'


args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")


for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	  if (res[1]=="inRdata") inRdata=as.character(res[2])
}


options(stringsAsFactors=FALSE)
load(inRdata)

write.table(file=gsub(".RData", "_tpm.tsv",inRdata) ,isoform_tpm,sep='\t',quote=FALSE)
