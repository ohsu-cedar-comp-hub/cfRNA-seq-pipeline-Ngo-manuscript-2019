files = snakemake@input

#make a empty dataframe, ncol= number of samples +1, nrow= number of rows in gene count file
all = data.frame(matrix(ncol = length(files) + 2, nrow=nrow(read.table(files[[1]], header = F, sep = "\t"))))

#get row names
samp1=read.table(file = files[[1]], header = F, sep = "\t")
genelist=samp1[,1:2]
genelist=as.data.frame(genelist)
colnames(genelist)=c("ensembl_gene_id","gene_id")
head(genelist)

#combine counts from files
for(i in 1:length(files)){
  counts = read.table(file = files[[i]], header = F, stringsAsFactors  = F, sep = "\t")
  print(files[[i]])
  samp=gsub("_htseq_exon_count.txt","",files[[i]])
  
  colnames(counts)=c("ensembl_gene_id","gene_id",samp)
  
  if(i==1){
    all=merge(genelist,counts)  
  } 
  else
  {all=merge(all,counts) }
  
}

write.table(all, file = snakemake@output[[1]], quote = F, sep = "\t", row.names = F)
