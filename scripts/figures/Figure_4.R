## qPCR boxplots
library(ggplot2)
library(data.table)
library(tidyr)
library(reshape2)
library(useful)
library(ggpubr)
library(dplyr)

setwd("/Users/breesheyroskams-hieter/Desktop/cfRNA/manuscript/revised_paper_with_validation/figures/qPCR")

metadata <- read.csv("../../tables/PP_metadata_keep_FINAL_updated.csv")
PP_qPCR <- read.table("../../tables/qPCR_data_pilot_cohort_filt.txt", header = T, stringsAsFactors = F, sep = "\t")
Targets <- read.table("../../tables/Target_guide.txt", stringsAsFactors = F, header = T)
counts_table_updated <- read.delim("../../tables/counts_table_updated.txt")
RPM_updated <- read.csv("../../tables/RPM_without_dates_updated.csv", stringsAsFactors = F, row.names = 1)
biomart_ensembl_geneid <- read.delim("../../../../biomart_ensembl_geneid.txt")
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

## read in LVQ results
LVQ_HCCvsHD <- readRDS(file="../../tables/LVQ/HD_HCC_importance.rds")
LVQ_MMvsHD <- readRDS(file="../../tables/LVQ/HD_MM_importance.rds")
HCC.vs.HD_top5lvq <- rownames(LVQ_HCCvsHD$importance[order(LVQ_HCCvsHD$importance$HD, decreasing = TRUE),])[1:5]
MM.vs.HD_top5lvq <- rownames(LVQ_MMvsHD$importance[order(LVQ_MMvsHD$importance$HD, decreasing = TRUE),])[1:5]

## Assign all missing values as 40
PP_qPCR[is.na(PP_qPCR$CT),3] <- 40

# Add group information
iv <- match(PP_qPCR$Sample, metadata$PP_ID)
PP_qPCR$Group <- metadata[iv,]$Status
iv <- match(PP_qPCR$Target_Name, Targets$Target_Name)
PP_qPCR$Target_type <- Targets[iv,]$Target_Type

# Rename all non-patient samples to a consistent naming structure
PP_qPCR[which(PP_qPCR$Sample %like% "NTC"),4] <- "NTC"
PP_qPCR[which(PP_qPCR$Sample =='(+)1'),4] <- "Positive Control"
PP_qPCR[which(PP_qPCR$Sample =='(+)2'),4] <- "Positive Control"
PP_qPCR[which(PP_qPCR$Sample =='(+)3'),4] <- "Positive Control"
PP_qPCR[which(PP_qPCR$Sample =='(+)4'),4] <- "Positive Control"
PP_qPCR[which(PP_qPCR$Sample =='(+)5'),4] <- "Positive Control"
PP_qPCR[is.na(PP_qPCR$Group),4] <- 'Control'

# Create separate dataframe for ACTB and B2M
PP_ACTB <- PP_qPCR[which(PP_qPCR$Target_Name == 'ACTB'),]
PP_B2M <- PP_qPCR[which(PP_qPCR$Target_Name == "B2M"),]

# Add in ACTB and B2M values to original dataframe
iv <- match(PP_qPCR$Sample, PP_ACTB$Sample)
PP_qPCR$ACTB <- PP_qPCR[iv,]$CT
iv <- match(PP_qPCR$Sample, PP_B2M$Sample)
PP_qPCR$B2M <- PP_qPCR[iv,]$CT

# Calculate deltaCT values
PP_qPCR$delta_ACTB <- PP_qPCR$CT - PP_qPCR$ACTB
PP_qPCR$delta_B2M <- PP_qPCR$CT - PP_qPCR$B2M

## Create separate dataframes for each type
PP_HCC <- PP_qPCR[which(PP_qPCR$Target_type == "HCC-LVQ"),]
PP_MM <- PP_qPCR[which(PP_qPCR$Target_type == "MM-LVQ"),]

## plot boxplots for pairwise comparisons of HD-HCC and HD-MM
GenerateBoxplot <- function(data, target, baseline, genelist, colors) {
  # Filter for types and targets
  filt <- data[data$Group %in% c(target, baseline) & data$Target_Name %in% genelist,]
  filt$Group <- factor(filt$Group, levels = c(baseline, target)) 
  # Colors
  colourBaseline <- colors[colors$Status==baseline,]$Colour
  colourTarget <- colors[colors$Status==target,]$Colour
  
  my_comparisons <- list(c(baseline,target))
  
  p <- ggplot(filt, mapping = aes(x = Group, y = CT, color = Group)) +
    geom_boxplot() + 
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 5),
          legend.position = "none") + 
    facet_wrap(~Target_Name, scales = 'free_y', ncol = 5) + 
    stat_compare_means(label = "p.signif", method = "t.test", 
                       comparisons = my_comparisons, 
                       show.legend = FALSE, size = 2, hide.ns = TRUE) +
    scale_y_reverse(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))),
                     expand = c(0.1,0)) + 
    scale_color_manual(breaks=c(baseline,target),
                       values=c(colourBaseline,colourTarget))
  
  return(p)
}

HCC_boxplot <- GenerateBoxplot(data = PP_qPCR, target = "HCC", baseline = "NC", genelist = HCC.vs.HD_top5lvq, colors = colors)
MM_boxplot <- GenerateBoxplot(data = PP_qPCR, target = "MM", baseline = "NC", genelist = MM.vs.HD_top5lvq, colors = colors)

pdf("qRT_PCR_top5LVQ_boxplots.pdf", 4, 4.5, useDingbats = FALSE)
cowplot::plot_grid(MM_boxplot, HCC_boxplot, nrow = 2)
dev.off()

# to run the correlation plot we first need to create a matrix
PP_qPCR1 <- PP_qPCR[which(!PP_qPCR$Group %like% "NTC"),]
PP_qPCR1 <- PP_qPCR1[which(!PP_qPCR1$Group %like% "Positive"),]
PP_target_set <- unique(PP_qPCR1$Target_Name)
RPM_subset <- RPM_updated[which(RPM_updated$gene %in% PP_target_set),]
rownames(RPM_subset) <- RPM_subset[,1]
RPM_subset <- data.matrix(RPM_subset[,-1])
logRPM <- log2(RPM_subset + 1)

# making matrices for corrplot
PP_qpcr_raw <- PP_qPCR1[,1:3]
PP_qpcr_raw <- pivot_wider(PP_qpcr_raw, names_from = Target_Name, values_from = CT)

# Transpose and reformat
t_qpcr <- transpose(PP_qpcr_raw)
rownames(t_qpcr) <- colnames(PP_qpcr_raw)
colnames(t_qpcr) <- t_qpcr[1,]
t_qpcr <- t_qpcr[-1,]
t_qpcr$gene <- rownames(t_qpcr)

# Melt for plotting purposes
mlt_qpcr <- reshape2::melt(t_qpcr, id.vars = "gene", variable.name = "Sample", value.name = "CT")
logRPM <- data.frame(logRPM)
logRPM$gene <- rownames(logRPM)
mlt_RPM <- reshape2::melt(logRPM, id.vars = "gene", variable.name = "Sample", value.name = "log2RPM")
logRPM_qpcr <- merge(mlt_qpcr, mlt_RPM, by = c("Sample", "gene"))
logRPM_qpcr$CT <- as.numeric(logRPM_qpcr$CT)
logRPM_qpcr$log2RPM <- as.numeric(logRPM_qpcr$log2RPM)

## Remove log2RPM values that are zero or CT values that are 40
keep <- setdiff(1:nrow(logRPM_qpcr), rownames(logRPM_qpcr[logRPM_qpcr$CT==40 | logRPM_qpcr$log2RPM==0,]))
filt <- logRPM_qpcr[keep,]

## Remove CT values > 28
filt <- filt[filt$CT < 28,]

logRPM_qpcr_plot <- ggplot(filt, aes(x = log2RPM, y = CT)) + 
  stat_cor(method = "pearson", show.legend = FALSE, label.y = 30) + 
  geom_point(show.legend = FALSE, alpha = 0.7) + 
  xlab("RNA-Seq log2(RPM + 1)") + 
  ylab("RT-qPCR Ct")

pdf("qRT_PCR_RPM_corrPlot.pdf", 4, 4)
print(logRPM_qpcr_plot)
dev.off()
