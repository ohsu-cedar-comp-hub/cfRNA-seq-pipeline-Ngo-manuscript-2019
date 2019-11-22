# This a script to make the figures for the permutation test (Figure S3)

library(DESeq2)
library(dplyr)
library(ggplot2)

# Set Working Directory
setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/figures/supplementary/permutationTest")

# Read in previously run data

GeneratePermutationTest <- function(Results, Counts, Metadata, Baseline, Target, cutoff = 0.01) {
  numdiff <- read.csv(Results, row.names=1, header=TRUE)
  names(numdiff) <- "NumDiffGenes"
  
  # Generate subdata 
  md <- read.csv(file=Metadata, stringsAsFactors = FALSE)
  cutoff <- 0.01
  
  # Extract only the Types that we want in further analysis & only the PP_ID and Status informative columns
  md <- filter(md, Status == Baseline | Status == Target)
  
  # Read in counts table
  subdata <- read.table(Counts, header=TRUE, row.names=1, sep="\t")
  
  # Keep only the PP_IDs of the types we have chosen in the metadata table above
  rownames(md) <- md$PP_ID
  md$PP_ID <- NULL
  keep <- colnames(subdata)[colnames(subdata) %in% rownames(md)]
  subdata <- subdata[, keep]
  dim(subdata)
  
  # Obtain the number of genes that meet padj<0.01 for reference line in histogram
  dim(md)
  
  # Order md and subdata so that they are in the same order, so that they can be properly read into DESeq
  md <- md[order(rownames(md)),]
  subdata <- subdata[,order(colnames(subdata))]
  
  dds_design <- "Status"
  
  dds <- DESeqDataSetFromMatrix(countData=subdata,
                                colData=md,
                                design= as.formula(paste('~',dds_design)))
  
  # Remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) >= 1, ]
  
  # Normalization and pre-processing
  dds <- DESeq(dds)
  
  # Extract results and the number of significant genes with padj<0.01
  results = results(dds, contrast = c("Status", Target, Baseline), independentFiltering = FALSE,cooksCutoff = Inf)
  numSig <- sum(results$padj < cutoff, na.rm=TRUE)
  
  # we want xlim to be close to numSig but this messes up the zero binning, so remove these values before plotting
  if (Target=="HCC") {
    numdiff[numdiff$NumDiffGenes>20,] <- 0
  } else if (Target=="Cirr") {
    numdiff[numdiff$NumDiffGenes>20,] <- 0
  } else if (Target=="MGUS") {
    numdiff[numdiff$NumDiffGenes>10,] <- 0
  } else if (Target=="MM") {
    numdiff[numdiff$NumDiffGenes>100,] <- 0
  }
  
  numdiff$Actual <- numSig
  
  p <- ggplot(numdiff, aes(x=NumDiffGenes)) +
    geom_histogram(bins=100) +
    theme_bw() +
    geom_vline(data=numdiff, mapping=aes(xintercept = numSig, color = "Correct Labels"), 
               linetype="longdash", size=0.3, show.legend = T) +
    scale_color_manual(values = "gray75", name = "Number of \nDE genes") +
    ggtitle(paste(Baseline, "vs", Target)) +
    xlab("Number of significant \ngenes") +
    ylab("Frequency") +
    theme(text = element_text(family = "Arial"),
          aspect.ratio=1,
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.title = element_text(size=8, hjust = 0.5),
          panel.border = element_blank(), axis.line = element_line(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  legend <- cowplot::get_legend(p)
  
  p <- p + theme(legend.position = "none")

  pdf(paste0("Legend", ".pdf"), 5, 2, useDingbats = FALSE)
  grid.newpage()
  grid.draw(legend)
  dev.off()
  
  return(p)
}

LuCa_perm <- GeneratePermutationTest(Results = "../../../tables/permutation/LuCa-vs-HD.number.diff.genes.csv", 
                                     Counts = "../../../tables/counts_table_updated.txt",
                                     Metadata = "../../../tables/PP_metadata_keep_FINAL_updated.csv",
                                     Baseline = "HD", Target = "LuCa")
MM_perm <- GeneratePermutationTest(Results = "../../../tables/permutation/MM-vs-HD.number.diff.genes.csv", 
                                   Counts = "../../../tables/counts_table_updated.txt",
                                   Metadata = "../../../tables/PP_metadata_keep_FINAL_updated.csv",
                                   Baseline = "HD", Target = "MM")
HCC_perm <- GeneratePermutationTest(Results = "../../../tables/permutation/HCC-vs-HD.number.diff.genes.csv", 
                                    Counts = "../../../tables/counts_table_updated.txt",
                                    Metadata = "../../../tables/PP_metadata_keep_FINAL_updated.csv",
                                    Baseline = "HD", Target = "HCC")

pdf("permutationTest_grid.pdf", 8.3, 2.6, useDingbats = FALSE)
cowplot::plot_grid(LuCa_perm, MM_perm, HCC_perm, nrow = 1)
dev.off()

# Check specificity of our DESeq2 results
LuCa_numdiff <- read.csv(file = "../../../tables/permutation/LuCa-vs-HD.number.diff.genes.csv", row.names=1, header=TRUE)
names(LuCa_numdiff) <- c("NumDiffGenes","Actual")
LuCa_specificity <- 1 - nrow(LuCa_numdiff[LuCa_numdiff$NumDiffGenes>0 & LuCa_numdiff$NumDiffGenes != 1820,]) / 500

MM_numdiff <- read.csv(file = "../../../tables/permutation/MM-vs-HD.number.diff.genes.csv", row.names=1, header=TRUE)
names(MM_numdiff) <-c("NumDiffGenes","Actual")
MM_specificity <- 1 - nrow(MM_numdiff[MM_numdiff$NumDiffGenes>0 & MM_numdiff$NumDiffGenes != 106,]) / 500

HCC_numdiff <- read.csv(file = "../../../tables/permutation/HCC-vs-HD.number.diff.genes.csv", row.names=1, header=TRUE)
names(HCC_numdiff) <-c("NumDiffGenes","Actual")
HCC_specificity <- 1 - nrow(HCC_numdiff[HCC_numdiff$NumDiffGenes>0 & HCC_numdiff$NumDiffGenes != 12,]) / 500
