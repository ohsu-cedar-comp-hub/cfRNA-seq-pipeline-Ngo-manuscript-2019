# This is a script to generate figures for Thuy Ngo's manuscript on cfRNA
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)

setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/figures/supplementary/")

# Read in important results tables
anno <- get(load("/Users/roskamsh/Desktop/CEDAR/anno/hg38.Ens_90.biomaRt.geneAnno.Rdata"))
counts <- read.csv("../../tables/RPM_without_dates_updated.csv", sep=",", stringsAsFactors = F, row.names = 1)
coverage <- read.delim("../../tables/read_coverage_updated.txt", stringsAsFactors = F)
md <- read.csv("../../tables/PP_metadata_keep_FINAL_updated.csv", stringsAsFactors = F)
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

rownames(coverage) <- coverage$Sample
rownames(md) <- md$PP_ID

## Format data ##
counts$gene <- as.character(counts$gene)

# Consolidate duplicate gene names
datam <- melt(counts,id="gene")
datac <- dcast(datam,gene~variable,fun.aggregate = mean)
rownames(datac) <- datac[,1]
counts <- datac[,-1]

# Check
stopifnot(colnames(counts)==rownames(md) & rownames(md)==rownames(coverage))

# Define colours
colourHD <- colors[colors$Status=="HD",]$Colour
colourLuCa <- colors[colors$Status=="LuCa",]$Colour
colourMM <- colors[colors$Status=="MM",]$Colour
colourMGUS <- colors[colors$Status=="MGUS",]$Colour
colourCirr <- colors[colors$Status=="Cirr",]$Colour
colourHCC <- colors[colors$Status=="HCC",]$Colour

# New graph theme
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    #panel.border = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(color='black', size=18),
    #axis.text.x=element_text(color='black',size=8),
    panel.grid.major.y = element_line(colour = "black"),
    panel.grid.minor.y = element_line(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_text(angle=45, hjust=1, colour="black"),
    axis.text.y=element_text(colour="black"),
    plot.title=element_text(size=18, face="bold", hjust = 0.5)
  )

dodge <- position_dodge(width = 0.6)
theme_update(plot.title = element_text(hjust = 0.5))

#################################################
# UMI Counts across all samples
#################################################

# Calculate UMI count for each sample
md$libSize <- colSums(counts)

md$Status <- factor(md$Status, levels = c("HD","LuCa","MGUS","MM","Cirr","HCC"))

##################################################
# Intron and Exon fraction
##################################################

# Import intron and exon fraction results into md table
iv <- match(md$PP_ID, coverage$Sample)
md$exon_fraction <- round(coverage[iv,]$Exon, 3)
md$intron_fraction <- round(coverage[iv,]$Intron, 3)
md$intergenic_fraction <- round(coverage[iv,]$Intergenic, 3)

# Now visualize in a different way, retaining the intergenic fraction as well, representing as a stacked barplot
df1 <- dplyr::select(md, intron_fraction)
df1$attribute <- "intronFraction"
names(df1) <- c("Value","attribute")
df1$Sample <- rownames(df1)
df2 <- dplyr::select(md, exon_fraction)
df2$attribute <- "exonFraction"
names(df2) <- c("Value","attribute")
df2$Sample <- rownames(df2)
df3 <- dplyr::select(md, intergenic_fraction)
df3$attribute <- "intergenicFraction"
names(df3) <- c("Value","attribute")
df3$Sample <- rownames(df3)

# Bind these three dataframes together, and then add in other information about each sample for plotting
forPlot <- rbind(df1,df2,df3)
iv <- match(forPlot$Sample, md$PP_ID)
head(md[iv,])
forPlot$Status <- paste(md[iv,]$Status)
forPlot$CEDAR_ID <- paste(md[iv,]$CEDAR_ID)

# Order dataframe by type
forPlot <- forPlot[order(match(forPlot$Status, c("HD","LuCa","MGUS","MM","Cirr","HCC"))),]
forPlot$attribute <- factor(forPlot$attribute, levels=c("intergenicFraction","intronFraction","exonFraction"))

fractionBar <- ggplot(forPlot, aes(x=Sample, y=Value, fill=attribute)) +
  geom_bar(stat="identity") +
  ylab("Fraction of each gene attribute") +
  xlab("Sample") +
  scale_fill_brewer( palette = "BuPu", labels = c("Intergenic fraction","Intron fraction","Exon fraction") ) +
  theme(text = element_text(family = "Arial"),
        plot.margin = unit(c(0.3,0,0.3,0),"cm"),
        axis.text.x=element_text(angle = 90, size = 5),
        axis.text.y=element_text(size=5),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
fractionBar

legend <- cowplot::get_legend(fractionBar)

fractionBar <- fractionBar + theme(legend.position = "none")

pdf("QC/Legend.Stacked.barPlot.byFraction.pdf", 5, 2, useDingbats = FALSE)
grid.newpage()
grid.draw(legend)
dev.off()

# Barplots of different QC metrics
md <- md[order(match(md$Status, c("HD","LuCa","MGUS","MM","Cirr","HCC"))),]
md$PP_ID <- factor(md$PP_ID, levels = md$PP_ID)
inputReadsBar <- ggplot(data=md, mapping=aes(x=PP_ID,y=total_input_reads,fill=Status)) +
  geom_bar(stat="identity") +
  ylab("Number of input reads") +
  xlab("Sample") +
  scale_fill_manual(values = c(colourHD, colourLuCa, colourMGUS, colourMM, colourCirr, colourHCC)) +
  theme(text = element_text(family = "Arial"),
        plot.margin = unit(c(0.3,0,0.3,0),"cm"),
        axis.text.x=element_text(angle = 90, size = 5),
        axis.text.y=element_text(size=5),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
inputReadsBar

legend <- cowplot::get_legend(inputReadsBar)

inputReadsBar <- inputReadsBar + theme(legend.position = "none")

pdf("QC/Legend.num.input.reads.pdf", 5, 2, useDingbats = FALSE)
grid.newpage()
grid.draw(legend)
dev.off()


################################################
# Annotation of biotypes for gene counts
################################################

# Remove unique gene identifiers from counts table
counts <- counts[-1,]
fD <- as.data.frame(rownames(counts))
names(fD) <- "ensID"
rownames(fD) <- fD$ensID
iv <- match(rownames(fD), anno$external_gene_name)
fD$biotype <- anno[iv,]$gene_biotype

# generate new dataframe expr which includes gene information
# here presence or absence of gene expression is denoted by a 1/0
# Stacked barplot with all other biotypes marked as "other"
expr <- counts
expr <- expr[rowSums(expr)>0,]
expr <- expr!=0
expr[expr=="FALSE"] <- 0
expr <- as.data.frame(expr)
expr$ensID <- rownames(expr)

# Add geneID and biotype information to this matrix
m <- match(rownames(expr), rownames(fD))
expr <- left_join(expr, fD[m,])
rownames(expr) <- expr$ensID

# Melt this matrix so it can be plotted
geneTable <- melt(expr, id=c("ensID","biotype"))
geneTable <- geneTable[geneTable$value>0,]

# Let's filter out some biotypes so that we don't have such a messy graph
ordered <- as.data.table(geneTable)[, .N, by = biotype][order(-N)]
# Grab the first 6 options
filtBio <- as.vector(ordered[6:1,]$biotype)
other <- as.vector(ordered[7:nrow(ordered)]$biotype)

# Filter for top 6 biotypes of interest
# Assign all biotypes that are not in the top 6 to be defined as "other"
geneTable$biotype[geneTable$biotype %in% other] <- "other"
geneTable$biotype <- factor(geneTable$biotype, levels=c("other", filtBio)) # factor so that our graph is in a comprehensive order

# Sum counts across all biotype/sampleID combinations
sampleTable <- as.data.table(geneTable)[, .N, by= .(variable)]
geneTable <- as.data.table(geneTable)[, .N, by = .(biotype, variable)]

# Calculate a fraction for each biotype per sample
s <- match(geneTable$variable, sampleTable$variable)
geneTable$Frac <- geneTable$N/sampleTable[s,]$N

# Add in Type information
c <- match(geneTable$variable, rownames(md))
geneTable$Type <- paste(md[c,]$Status)

biotypeBar2 <- ggplot(geneTable, aes(x=variable, y=Frac, fill=biotype)) +
  geom_bar(stat="identity") +
  ylab("Fraction of each biotype") +
  xlab("Sample") +
  scale_fill_brewer( palette = "YlGnBu" ) +
  theme(text = element_text(family = "Arial"),
        plot.margin = unit(c(0.3,0,0.3,0),"cm"),
        axis.text.x=element_text(angle = 90, size = 5),
        axis.text.y=element_text(size=5),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
biotypeBar2

legend <- cowplot::get_legend(biotypeBar2)

biotypeBar2 <- biotypeBar2 + theme(legend.position = "none")

pdf("QC/Legend.biotype.barplot.pdf", 5, 2, useDingbats = FALSE)
grid.newpage()
grid.draw(legend)
dev.off()

pdf("QC/QC_barplots_grid.pdf", 6.5, 7.5, useDingbats = FALSE)
cowplot::plot_grid(inputReadsBar, fractionBar, biotypeBar2, nrow = 3)
dev.off()
