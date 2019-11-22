library(Glimma)
library(limma)
library(DESeq2)

condition = snakemake@params[['condition']]
cat(sprintf(c('Condition: ',condition,'\n')))
Type <- snakemake@config[['linear_model']]

title = snakemake@params[["contrast"]]

contrast = c(condition, snakemake@params[["contrast"]])
rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

ens2geneID <- snakemake@config[['ens2geneID']]

out_path = file.path(getwd(),'results','diffexp')
dir.create(out_path)
print(out_path)
rds = readRDS(rds)
groups.df = as.data.frame(colData(rds))


#### by contrasts
contrasts_to_plot = resultsNames(rds)
res <- results(rds, contrast=contrast)
res$padj[is.na(res$padj)] = 1

rnaseq = as.data.frame(counts(rds, normalized=T))

## Replace ensemble id's with gene id's
gene_id = read.delim(ens2geneID)

## Remove unique identifier .xx from heatmap data
rownames(res) <- sub("\\.[0-9]*", "", rownames(res))
iv <- match(rownames(res), gene_id$ensembl_gene_id)
head(gene_id[iv,])

res$GeneID  <- paste(gene_id[iv, "external_gene_name"])

res <- res[order(res$padj),]
# Remove duplicated geneIDs (this will, be default, remove the ones with the higher p-values, as we have ordered by p-value above
res <- res[!duplicated(res$GeneID),]

# Also subset your rnaseq counts table by these ensembl IDs so all of the data matches
rownames(rnaseq) <- sub("\\.[0-9]*", "", rownames(rnaseq))
rnaseq <- rnaseq[rownames(rnaseq) %in% rownames(res),]
# Now match gene symbols to these ensembl IDS
iv <- match(rownames(rnaseq), gene_id$ensembl_gene_id)
head(gene_id[iv,])
# Add gene symbols as rownames
rownames(rnaseq)  <- paste(gene_id[iv, "external_gene_name"])

# Now remove old ensembl IDs from results and paste gene symbol there
rownames(res) <- res$GeneID
res$GeneID <- NULL

genes = as.data.frame(row.names(res))
colnames(genes) = 'GeneID'

status_frame = res[,c('log2FoldChange','padj')]
status_frame['status'] = 0
status_frame$padj[is.na(status_frame$padj)] = 1
status_frame[status_frame$padj<0.05 & status_frame$log2FoldChange < 0 ,'status'] = -1
status_frame[status_frame$padj<0.05 & status_frame$log2FoldChange > 0 ,'status'] = 1

title = paste(title[1],'vs',title[2],sep='-')

glMDPlot(res, anno=genes, status=status_frame$status, samples=colnames(rnaseq), 
         counts=log2(rnaseq + 0.0001), groups=groups.df[[Type]], main=strsplit(res@elementMetadata$description[2],': ')[[1]][2], 
         transform=F, side.ylab='Log2-expression',launch=FALSE,side.main='GeneID', html = paste(title,'ma_plot',sep='.'),path=out_path)

## Volcano plot
glXYPlot(x=res$log2FoldChange, y=-log10(res$pvalue), xlab="logFC", ylab="logodds",path=out_path,
         status=status_frame$status, launch=FALSE,counts=log2(rnaseq + 0.0001), groups=groups.df[[Type]], anno=genes,main=strsplit(res@elementMetadata$description[2],': ')[[1]][2],html = paste(title,'volcano_plot',sep='.'))

