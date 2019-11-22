library("DESeq2")
library("ggplot2")
library("pheatmap")
library("dplyr")
library("vsn")
library("RColorBrewer")
library("genefilter")

cat(sprintf(c('Working directory',getwd())))

cat(sprintf('Setting parameters'))

pca_plot <- snakemake@output[['pca']]
cat(sprintf(c('PCA plot: ',pca_plot)))

labels <- snakemake@params[['pca_labels']]
cat(sprintf(c('PCA Labels: ',labels)))

sd_mean_plot <- snakemake@output[['sd_mean_plot']]
cat(sprintf(c('SD Mean plot: ',sd_mean_plot,'\n')))

distance_plot <- snakemake@output[['distance_plot']]
cat(sprintf(c('Distance plot: ',distance_plot,'\n')))

heatmap_plot <- snakemake@output[['heatmap_plot']]
cat(sprintf(c('Heatmap Plot: ', heatmap_plot, '\n')))

rds_out <- snakemake@output[['rds']]
cat(sprintf(c('RDS Output: ', rds_out, '\n')))

rld_out <- snakemake@output[['rld_out']]
cat(sprintf(c('RLD Output: ', rld_out, '\n')))

counts <- snakemake@input[['counts']]
cat(sprintf(c('Counts table: ', counts, '\n')))

metadata <- snakemake@params[['samples']]
cat(sprintf(c('Metadata: ', metadata, '\n')))

sampleID <- snakemake@params[['sample_id']]
cat(sprintf(c('Sample ID: ', sampleID, '\n')))

Type <- snakemake@params[['linear_model']]
cat(sprintf(c('Linear Model: ', Type, '\n')))

group <- snakemake@params[['LRT']]
cat(sprintf(c('Subsetted group: ', group, '\n')))

plot_cols <- snakemake@config[['meta_columns_to_plot']]
subset_cols = names(plot_cols)

# color palette
colors <- snakemake@params[['colors']]
discrete <- snakemake@params[['discrete']]

# function to grab the ggplot2 colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Dir <- "results/diffexp/group/"

md <- read.delim(file=metadata, sep = "\t", stringsAsFactors = FALSE)
md <- md[order(md[sampleID]),]

# Read in counts table
cts <- read.table(counts, header=TRUE, row.names=1, sep="\t", check.names=F)
cts <- cts[,order(colnames(cts))]

# Put sample IDs as rownames of metadata
rownames(md) <- md[[sampleID]]
md[[sampleID]] <- NULL

# Ensure that we subset md to have exactly the same samples as in the counts table
md <- md[colnames(cts),]
dim(md)

# Check
stopifnot(rownames(md)==colnames(cts))

# Define colours based on number of Conditions
if(colors[[1]] !='NA' & discrete[[1]] =='NA'){
    if (brewer.pal.info[colors[[1]],]$maxcolors >= length(unique(md[[Type]]))) {
        pal <- brewer.pal(length(unique(md[[Type]])),name=colors[[1]])
    } 
} else if(discrete[[1]] != 'NA' & length(discrete)==length(unique(md[[Type]]))){
        pal <- unlist(discrete)
} else {
        pal <- gg_color_hue(length(unique(md[[Type]])))
}

# Create dds object from counts data and correct columns
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=md,
                              design= as.formula(paste('~',Type)))

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) >= 1, ]

# Likelihood Ratio test to look at differential expression across ALL types, and not just pairs of types (contrast)
dds.lrt <- DESeq(dds, test="LRT", reduced=~1)
res.lrt <- results(dds.lrt, cooksCutoff = Inf, independentFiltering=FALSE)
head(res.lrt)

# Obtain normalized counts
rld <- rlog(dds.lrt, blind=FALSE)

# Pairwise PCA Plot
pcaData <- plotPCA(rld, intgroup=labels[[1]], returnData=TRUE)

pdf(pca_plot)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes_string("PC1", "PC2", color=labels[[1]])) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  scale_colour_manual(values=pal)
dev.off()

# Pairwise PCA Plot with more than one PCA parameter
if (length(labels)>1) {
  pca_plot2 <- sub("$","twoDimensional_pca_plot.pdf", Dir)
  pcaData <- plotPCA(rld, intgroup=c(labels[[1]], labels[[2]]), returnData=TRUE)
  pdf(pca_plot2, 5, 5)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes_string("PC1", "PC2", color=labels[[1]], shape=labels[[2]])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    scale_colour_manual(values=pal)
  dev.off()
}

# SD mean plot
pdf(sd_mean_plot)
meanSdPlot(assay(rld))
dev.off()

# Heatmap of distances
pdf(distance_plot)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- colnames(rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, fontsize=5, scale="row",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# Heatmap across all samples
# List top 50 genes for group comparisons
topGenes <- head(order(res.lrt$padj), 50)

# Extract topGenes from rld object
plot <- assay(rld)[topGenes,] #for 2+ types

# Generate data frame with samples as the rownames and single colData as the first row
# Default when we subset creates an incompatible dataframe so this is a check
df <- as.data.frame(colData(rld))
if (length(subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
  names(annot) <- c("SampleID", subset_cols[1])
  rownames(annot) <- annot[[sampleID]]
  annot[[sampleID]] <- NULL
} else {
  annot <- df[,subset_cols]
}

filt <- plot[apply(plot, MARGIN = 1, FUN = function(x) sd(x) != 0),]

pdf(heatmap_plot)
pheatmap(filt, cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes across all samples"))
dev.off()

saveRDS(dds, file=rds_out)
saveRDS(rld, file=rld_out)

group <- as.vector(group)

# If LRT group has been specified, run the analysis for that group
if (length(group)>0) {
  md <- read.delim(file=metadata, sep = "\t", stringsAsFactors = FALSE)
  md <- md[order(md[sampleID]),]
  cts <- read.table(counts, header=TRUE, row.names=1, sep="\t")
  cts <- cts[,order(colnames(cts))]
  md <- md[md[[Type]] %in% group,]
  rownames(md) <- md[[sampleID]]
  md[[sampleID]] <- NULL
  keep <- colnames(cts)[colnames(cts) %in% rownames(md)]
  cts <- cts[, keep]
  dim(cts)
  md <- md[colnames(cts),]
  dim(md)

  dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=md,
                              design= as.formula(paste('~',Type)))
  dds <- dds[ rowSums(counts(dds)) >= 1, ]
  dds.lrt <- DESeq(dds, test="LRT", reduced=~1)
  res.lrt <- results(dds.lrt, cooksCutoff = Inf, independentFiltering=FALSE)
  rld <- rlog(dds.lrt, blind=FALSE)
    
  # Pairwise PCA Plot
  pdf(sub("$", "subsetted_pca_plot.pdf", Dir), 5, 5)
  plotPCA(rld, intgroup=labels[[1]])
  dev.off()
  # Pairwise PCA Plot with more than one PCA parameter
  if (length(labels)>1) {
    pcaData <- plotPCA(rld, intgroup=c(labels[[1]], labels[[2]]), returnData=TRUE)
    pdf(sub("$", "subsetted_twoDimensional_pca_plot.pdf", Dir), 5, 5)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes_string("PC1", "PC2", color=labels[[1]], shape=labels[[2]])) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed()
    dev.off()
  }
  
  # Heatmap
  topGenes <- head(order(res.lrt$padj), 50)
  # Extract topGenes from rld object
  plot <- assay(rld)[topGenes,] #for 2+ types

  df <- as.data.frame(colData(rld))
  if (length(subset_cols)==1) {
    annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
    names(annot) <- c("SampleID", subset_cols[1])
    rownames(annot) <- annot[[sampleID]]
    annot[[sampleID]] <- NULL
  } else {
    annot <- df[,subset_cols]
  }

  pdf(sub("$", "subsetted_heatmap.pdf", Dir), 5, 5)
  pheatmap(assay(rld)[topGenes,], cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes across selected samples"))
  dev.off()
}
