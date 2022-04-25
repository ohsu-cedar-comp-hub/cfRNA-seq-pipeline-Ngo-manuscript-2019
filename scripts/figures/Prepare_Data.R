# Prepare data for further analysis
library(dplyr)
setwd("/Users/breesheyroskams-hieter/Desktop/cfRNA/manuscript/revised_paper_with_validation/tables/")

ensCounts <- read.table("counts_table_updated.txt", header = TRUE, stringsAsFactors = FALSE)
metadata <- read.csv("PP_metadata_keep_FINAL_updated.csv", stringsAsFactors = FALSE)

# Filter ensebl counts
ensCounts <- ensCounts[,colnames(ensCounts) %in% metadata$PP_ID]
ensCounts <- cbind(rownames(ensCounts), ensCounts)
colnames(ensCounts)[1] <- "gene"
rownames(ensCounts) <- NULL

## Keep ENS IDs
old <- ensCounts
## Replace ensemble id's with gene id's
gene_id = read.delim("/Users/breesheyroskams-hieter/Desktop/cfRNA/biomart_ensembl_geneid.txt")

## Remove unique identifier .xx from gencode IDs
ensCounts$gene <- sub("\\.[0-9]*", "", ensCounts$gene)
iv <- match(ensCounts$gene, gene_id$ensembl_gene_id)
head(gene_id[iv,])

## Replace ensemble id with gene name
ensCounts$gene <- paste(gene_id[iv, "external_gene_name"])

## Write table output with new gene names
write.table(ensCounts, file = "counts_tables_geneID_updated.txt", 
            sep="\t", row.names = FALSE, quote=F)

# Generate RPM
tmp <- sweep(ensCounts[,2:67], 2, colSums(ensCounts[,2:67]), '/') * 1e6
RPM <- data.frame(ensCounts$gene, tmp)
colnames(RPM)[1] <- "gene"

# Export to a table
write.table(RPM, file = "RPM_updated.txt", 
            sep="\t", row.names = FALSE, quote=F)

# Write RPMs for ens IDs
tmp <- sweep(old[,2:67], 2, colSums(ensCounts[,2:67]), '/') * 1e6
RPM <- data.frame(old$gene, tmp)
colnames(RPM)[1] <- "gene"
RPM$gene <- sub("\\.[0-9]*", "", RPM$gene)

# Export to a table
write.table(RPM, file = "RPM_updated_ensIDs.txt", 
            sep="\t", row.names = FALSE, quote=F)
