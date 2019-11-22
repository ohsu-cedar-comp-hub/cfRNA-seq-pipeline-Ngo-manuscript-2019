library("DESeq2")
library("reshape2")
library("cowplot")
library("limma")
library("vsn")
library("genefilter")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("pheatmap")
library("hexbin")

# output files
MDS_out <- snakemake@output[['mds_plot']]
MDS_table <- snakemake@output[['mds_table']]
heatmap_out <- snakemake@output[['heatmap_plot']]
sd_out <- snakemake@output[['sd_plot']]
normCounts_out <- snakemake@output[['rlogCounts_plot']]
normCounts_fac <- snakemake@output[['rlogCounts_fac_plot']]
rawCounts_out <- snakemake@output[['counts_plot']]
rawCounts_fac <- snakemake@output[['counts_fac_plot']]

# parameters
sampleID <- snakemake@params[['sample_id']]
Type = snakemake@params[['linear_model']]
plot_cols <- snakemake@config[['meta_columns_to_plot']]
subset_cols = names(plot_cols)

# color palette
colors <- snakemake@params[['colors']]
discrete <- snakemake@params[['discrete']]

# DESeq2 objects
rld <- snakemake@input[['rld']]
dds <- snakemake@input[['rds']]

rld <- readRDS(rld)
dds <- readRDS(dds)

# function to grab the ggplot2 colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

rawCounts <- counts(dds, normalized=FALSE)
md <- as.data.frame(colData(rld))
md$SampleID <- rownames(md)

if(colors[[1]] !='NA' & discrete[[1]] =='NA'){
    if (brewer.pal.info[colors[[1]],]$maxcolors >= length(unique(md[[Type]]))) {
        pal <- brewer.pal(length(unique(md[[Type]])),name=colors[[1]])
    } 
} else if(discrete[[1]] != 'NA' & length(discrete)==length(unique(md[[Type]]))){
        pal <- unlist(discrete)
} else {
        pal <- gg_color_hue(length(unique(md[[Type]])))
}
    

df1 <- melt(rawCounts) %>%
  dplyr::rename(Gene=Var1) %>%
  dplyr::rename(SampleID=Var2) %>%
  dplyr::rename(counts=value)

iv <- match(df1$SampleID, md$SampleID)
df1$Condition <- paste(md[iv,][[Type]])
df1$SampleID <- factor(df1$SampleID, levels=unique(md$SampleID))

# aesthetic for plots
dodge <- position_dodge(width = 0.6)
theme_update(plot.title = element_text(hjust = 0.5))

p1 <- ggplot(data=df1, mapping=aes(x=SampleID, y=counts, fill=Condition)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.2, outlier.colour=NA, position = dodge, color="gray28") +
  scale_y_log10() +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(hjust=1, angle=45, size=6))

# width of pdf to ensure all sampleIDs are visible when exported to pdf
# This was generated with a use case of 16 samples and a width of 7 fitting well, the +8 is to account for the margins
if (nrow(md)<8) {
  width <- 6
} else {
  width <- 7/24*(nrow(md)+8)
}

# raw counts boxplot
pdf(rawCounts_out, width, 5)
print({
  p1
})
dev.off()

# faceted by condition
p2 <- ggplot(data=df1, mapping=aes(x=SampleID, y=counts, fill=Condition)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.2, outlier.colour=NA, position = dodge, color="gray28") +
  scale_y_log10() +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(hjust=1, angle=45, size=4)) +
  facet_wrap(~Condition)

pdf(rawCounts_fac, 2*width, 5)
print({
  plot_grid(p1, p2)
})
dev.off()

# Run same analysis for log2-transformed normalized counts
df2 <- melt(assay(rld)) %>%
  dplyr::rename(Gene=Var1) %>%
  dplyr::rename(SampleID=Var2) %>%
  dplyr::rename(normCounts=value)

# Add Condition information to this dataframe
iv <- match(df2$SampleID, md$SampleID)
df2$Condition <- paste(md[iv,][[Type]])
df2$SampleID <- factor(df2$SampleID, levels=unique(md$SampleID))

p1 <- ggplot(data=df2, mapping=aes(x=SampleID, y=normCounts, fill=Condition)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.2, outlier.colour=NA, position = dodge, color="gray28") +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(hjust=1, angle=45, size=6)) +
  ylab("regularized log expression")

# raw counts boxplot
pdf(normCounts_out, width, 5)
print({
  p1
})
dev.off()

# faceted by condition
p2 <- ggplot(data=df2, mapping=aes(x=SampleID, y=normCounts, fill=Condition)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.2, outlier.colour=NA, position = dodge, color="gray28") +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(hjust=1, angle=45, size=4)) +
  facet_wrap(~Condition) +
  ylab("regularized log expression")

pdf(normCounts_fac, 2*width, 5)
print({
  plot_grid(p1, p2)
})
dev.off()

# Standard deviation vs. mean
ntd <- normTransform(dds)

pdf(sd_out)
meanSdPlot(assay(ntd))
dev.off()

# Generate annotation column for heatmap
if (length(subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(md), paste(md[[subset_cols[1]]])))
  names(annot) <- c("SampleID", subset_cols[1])
  rownames(annot) <- annot$SampleID
  annot$SampleID <- NULL
} else {
  annot <- md[,subset_cols]
}

# Remove rows with no standard deviation so that pheatmap executes properly
filt <- assay(rld)[apply(assay(rld), MARGIN = 1, FUN = function(x) sd(x) != 0),]

hm <- pheatmap(filt, show_rownames=F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", clustering_method = "average", 
         annotation_col = annot, scale = "row", 
         main="Unsupervised heatmap of all gene counts across samples",
         fontsize_row=4, fontsize_col=6, fontsize=8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

save_pheatmap_pdf <- function(x, filename) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(hm, heatmap_out)

# use plotMA function from limma, then extract data from this variable to plot with ggplot2
p <- plotMDS(assay(rld), top = 1000)
df <- data.frame(x=p$x, y=p$y, name=names(p$x))
iv <- match(df$name, md$SampleID)
df$Condition <- paste(md[iv,][[Type]])

pdf(MDS_out)
ggplot(data=df, mapping=aes(x=x,y=y)) +
  geom_point(size=3, colour = "black", show.legend = TRUE) +
  geom_point(aes(color=Condition), size=2.2) +
  scale_colour_manual(values=pal) +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  ggtitle("MDS Plot")
dev.off()

# Export the table for MDS
write.table(df, file=MDS_table, sep="\t", quote=F, row.names=FALSE)

