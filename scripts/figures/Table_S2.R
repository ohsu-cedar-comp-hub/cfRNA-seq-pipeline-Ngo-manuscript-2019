# This is a script to read in and obtain the quantitative information for general stats
library(dplyr)
library(reshape2)
library(data.table)

setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/tables")

#######################################
# General stats
#######################################

counts <- read.csv("RPM_without_dates_updated.csv",sep=",", row.names = 1)
md <- read.csv("PP_metadata_keep_FINAL_updated.csv",sep = ",", stringsAsFactors = FALSE)
anno <- get(load("/Users/roskamsh/Desktop/CEDAR/anno/hg38.Ens_90.biomaRt.geneAnno.Rdata"))

# Filter for only samples of interest
rownames(md) <- md$PP_ID
md$total_unique_reads <- as.numeric(md$total_unique_reads)

stopifnot(rownames(md)==colnames(counts)[2:ncol(counts)])

# General stats across ALL samples
allSampStats <- as.data.frame(mean(md$total_unique_reads))
names(allSampStats) <- "Mean_Unique_Reads"
rownames(allSampStats) <- "All_Samples"
allSampStats$Mean_Exon_Fraction <- mean(md$exon_fraction)

# Consolidate duplicate gene names
datam <- melt(counts,id="gene")
datac <- dcast(datam,gene~variable,fun.aggregate = mean)
rownames(datac) <- datac[,1]
counts <- datac[,-1]

####################
# Biotype breakdown
####################

# Remove unique gene identifiers from counts table
fD <- as.data.frame(rownames(counts))
names(fD) <- "ensID"
rownames(fD) <- fD$ensID
iv <- match(rownames(fD), anno$external_gene_name)
fD$biotype <- anno[iv,]$gene_biotype

# generate new dataframe expr which includes gene information
# here presence or absence of gene expression is denoted by a 1/0
expr <- counts
expr <- expr[rowSums(expr)>0,]
expr <- expr!=0
expr[expr=="FALSE"] <- 0
expr <- as.data.frame(expr)
expr$ensID <- rownames(expr)

# Add the total number of annotated features to allSampleStats
allSampStats$Number_Annotated_Features <- nrow(expr)

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
c <- match(geneTable$variable, md$PP_ID)
geneTable$Type <- paste(md[c,]$Status)

# Save fraction protein coding for each sample 
temp <- filter(geneTable, biotype=="protein_coding")
md$protein_coding_frac <- temp$Frac

# Calculate mean protein coding fraction acorss all samples
allSampStats$Mean_Protein_Coding_Fraction <- mean(temp$Frac)

## Fill in this information for each Type
temp <- as.data.frame(group_by(md, Status) %>%
  summarise(Mean_Unique_Reads=mean(total_unique_reads),
            Mean_Exon_Fraction=mean(exon_fraction),
            Mean_Intron_Fraction=mean(intron_fraction),
            Mean_Intergenic_Fraction=mean(intergenic_fraction),
            Mean_Protein_Coding_Fraction=mean(protein_coding_frac)))

iv <- match(typeData$Type, temp$Status)
typeData <- cbind(typeData, temp[iv,]$Mean_Unique_Reads, temp[iv,]$Mean_Exon_Fraction, temp[iv,]$Mean_Intron_Fraction, 
                  temp[iv,]$Mean_Intergenic_Fraction, temp[iv,]$Mean_Protein_Coding_Fraction)
names(typeData) <- c("Type","Diagnosis","Mean Unique Reads","Mean Exon Fraction","Mean Intron Fraction",
                     "Mean Intergenic Fraction","Mean Protein Coding Fraction")

cols <- c("Mean Exon Fraction","Mean Intron Fraction","Mean Intergenic Fraction","Mean Protein Coding Fraction")

# Round to two decimal places
for (i in 1:length(cols)) {
  col <- cols[i]
  typeData[[col]] <- sprintf("%.2f", typeData[[col]])
}

## Export as table
# Create a metadata table for information on all samples
perSampleStats <- as.data.frame(cbind(as.character(md$CEDAR.Number), as.character(md$Status), as.character(md$Gender), as.character(md$Date_of_RNA_Extraction), 
                                      as.character(md$Date_of_Lib_Prep), md$number_of_reads, md$total_unique_reads, md$exon_fraction, md$intron_fraction, 
                                      md$intergenic_fraction, md$protein_coding_frac))
names(perSampleStats) <- c("Sample ID", "Diagnosis", "Gender", "Date of RNA Extraction", "Date of Library Prep", "Total Number of Reads", 
                           "Number of Unique Reads", "Exon Fraction", "Intron Fraction", "Intergenic Fraction", "Fraction Protein Coding")
# Order the same way as the typeData
m <- match(perSampleStats$Diagnosis, typeData$Type)
head(perSampleStats[m,])
perSampleStats <- perSampleStats[order(m),]

# Round the values of the allSampStats to 2 decimal places for better viewing
cols <- c("Mean_Intron_Fraction", "Mean_Protein_Coding_Fraction", "Mean_Intergenic_Fraction", "Mean_Exon_Fraction")
for (i in 1:length(cols)) {
  col <- cols[i]
  allSampStats[[col]] <- sprintf("%.2f", allSampStats[[col]])
}

# Create a dataframe for all sample stats that include an attribute row name
allSampPlot <- as.data.frame(cbind(gsub("_"," ",names(allSampStats)), t(allSampStats)))
names(allSampPlot) <- c("Attribute","All Samples")

# Remove all but two decimal points from columns in PerSampleStats
# These are factored, so we need to refer to the values here
cols <- c("Exon Fraction", "Intron Fraction", "Intergenic Fraction", "Fraction Protein Coding")

# Reduce all fractions to only two decimal places
for (i in 1:length(cols)) {
  col <- cols[i]
  perSampleStats[[col]] <- sprintf("%.2f", as.numeric(levels(perSampleStats[[col]])[perSampleStats[[col]]]))
}

write.table(perSampleStats, file="perSampleStats.csv", sep = ",", quote = F, row.names = FALSE)

library(gridExtra)
library(cowplot)

p1 <- tableGrob(allSampPlot, rows = NULL)
p2 <- tableGrob(perSampleStats, rows = NULL)
p3 <- tableGrob(typeData, rows = NULL)

cairo_pdf("All_Sample_Statistics.pdf",height=2, width=10)
plot_grid(p1)
dev.off()

cairo_pdf("By_Type_Statistics.pdf",height=2, width=18)
plot_grid(p3)
dev.off()

cairo_pdf("Per_Sample_Statistics.pdf",height=21, width=20)
plot_grid(p2)
dev.off()

write.table(perSampleStats, file="AllSamples.Statistics.txt", sep="\t", quote=F, row.names=FALSE)
