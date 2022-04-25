# Script to generate plots in Figure 1, excluding the LOOCV tables, which are made separately

# Load libraries
library(reshape2)
library(ggplot2)
library(MASS)
library(pROC)
library(randomForest)
library(caret)
library(lattice)
library(expss)
library(cowplot)
library(grid)
library(gridExtra)
library(ggpubr)

# Set working directory
setwd("/Users/breesheyroskams-hieter/Desktop/cfRNA/manuscript/revised_paper_with_validation/figures/figure1/")

############################################
# Import Data
############################################

data <- read.csv("../../tables/RPM_without_dates_updated.csv",sep=",", row.names = 1)
metadata <- read.table("../../tables/combined_metadata_with_stage.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

# DE from DESeq2
DE_HCCvsHD <- read.table(file="../../tables/DE/HCC-vs-HD.diffexp.geneID.tsv",header = TRUE,stringsAsFactors = F)
DE_MMvsHD <- read.table(file="../../tables/DE/MM-vs-HD.diffexp.geneID.tsv",header = TRUE,stringsAsFactors = F)

# LVQ genes from lvq script
LVQ_HCCvsHD <- readRDS(file="../../tables/LVQ/HD_HCC_importance.rds")
LVQ_MMvsHD <- readRDS(file="../../tables/LVQ/HD_MM_importance.rds")

colors <- read.delim("../../tables/colors_by_stage.txt", header = T)

############################################
# Format Data
############################################
data$gene <- as.character(data$gene)

# Consolidate duplicate gene names
datam <- reshape2::melt(data,id="gene")
datac <- reshape2::dcast(datam,gene~variable,fun.aggregate = mean)
# Remove empty first row
datac <- datac[-1,]
rownames(datac) <- datac[,1]
datac <- datac[,-1]

## Select only pilot samples
metadata <- metadata[metadata$Set=="Pilot",]
metadata <- metadata[order(metadata$SeqID),]
datac <- datac[,order(colnames(datac))]

#######################################################
## LDA boxplots split by Stage, Etiology and Sex
#######################################################

# initialize variables
Data=datac
Metadata=metadata
Baseline="NC"
Target="HCC"

GenerateROCforSexStageEtiology <- function(Data, Metadata, Baseline, Target, colors) {
  if (Target=="HCC") {
    DElist=DE_HCCvsHD
    LVQlist=LVQ_HCCvsHD
  } else if (Target=="MM") {
    DElist=DE_MMvsHD
    LVQlist=LVQ_MMvsHD
  }
  Colors = colors
  
  plots <- list()
  ######### LOOCV LDA and DE ###################
  sample_sel <- Metadata[Metadata$Status==Baseline|Metadata$Status==Target,]$SeqID
  datac_sel <- Data[,colnames(Data) %in% sample_sel]
  names(datac_sel) <- make.names(names(datac_sel))
  
  metadata_info <- Metadata[Metadata$SeqID %in% sample_sel,]
  rownames(metadata_info)=metadata_info$SeqID
  metadata_sel=metadata_info[,c(4,5)]
  
  # Select LVQ genes
  genes_sel <- rownames(LVQlist$importance[order(LVQlist$importance$HD, decreasing = TRUE),])[1:10]
  select_datac= which((row.names(datac)) %in% (genes_sel))
  
  data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_train <- data_train[,-1]
  data_train$Status <- factor(data_train$Status, levels = c(Baseline, Target))
  table(data_train$Status)
  
  names(data_train) <- make.names(names(data_train))
  lda.fit = lda(Status ~., data_train)
  lda.pred=predict(lda.fit,data_train)
  table(data_train$Status,lda.pred$class)
  
  LD1_lvq <- data.frame(LDA_lvq = lda.pred$x,Status = data_train$Status)
  LD1_lvq$SeqID <- rownames(LD1_lvq)
  LD1_lvq <- reshape2::melt(LD1_lvq, by=c("Status","SeqID"))
  
  ## Add in stage, etiology and type information
  iv <- match(LD1_lvq$SeqID, metadata_info$SeqID)
  LD1_lvq <- data.frame(LD1_lvq, metadata_info[iv,c("Sex","Stage","Etiology_group","Age")])
  colnames(LD1_lvq)[3:4] <- c("Model","LD1")
  LD1_lvq$Model <- "LDA_lvq"
  
  LD1_both <- LD1_lvq
  LD1_both$Model <- plyr::mapvalues(LD1_both$Model, c("LDA_lvq"), c("LDA/LVQ"))
  
  ############ Boxplot of linear discrimination ###############
  ## Define comparisons for statistical test
  basecomp <- "na"
  stagecomp <- unique(LD1_both$Stage)[!(unique(LD1_both$Stage) %in% basecomp)]
  stagecomp <- stagecomp[order(stagecomp)]
  # Create combined stage column
  # for HCC - B,C,D
  # for MM - II,III
  LD1_both$combined_stage <- ifelse(LD1_both$Stage %in% as.vector(stagecomp[2:length(stagecomp)]), 
                                    paste(stagecomp[2:length(stagecomp)], collapse = ","), 
                                    ifelse(LD1_both$Stage == stagecomp[1], stagecomp[1], "NC"))
  
  # Now generate my_comparisons list for ggplot command
  my_comparisons <- list()
  if (Target == "MM") {
    my_comparisons <- list(c("NC","I"),c("NC","II,III"))
  } else if (Target == "HCC") {
    my_comparisons <- list(c("NC","A"),c("NC","B,C,D"))
  }
  
  LD1_both$combined_stage <- factor(LD1_both$combined_stage, 
                                    levels = c("NC",as.vector(stagecomp[1]),
                                               paste(stagecomp[2:length(stagecomp)], collapse = ",")))
  
  ## Thuy only wants LDA_LVQ set
  LD1_both <- LD1_both[LD1_both$Model=="LDA/LVQ",]
  
  # Colours for combined stages and either HCC or MM
  if (Target == "HCC") {
    collabs <- c("na","A","B,C,D")
    colvals <- c(colors[colors$Stage=="na",]$Colour, 
                 colors[colors$Stage=="A",]$Colour,
                 colors[colors$Stage=="B_C_D",]$Colour)
  } else if (Target == "MM") {
    collabs <- c("na","I","II,II")
    colvals <- c(colors[colors$Stage=="na",]$Colour, 
                 colors[colors$Stage=="I",]$Colour,
                 colors[colors$Stage=="II_III",]$Colour)
  }
  ## CHECK
  ggplot(data = LD1_both, mapping=aes(x=Status,y=LD1, colour = Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4) +
    ylab("Linear \nDiscrimination") +
    xlab("Stage") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          legend.position = "none") 
  
  bstage <- ggplot(data = LD1_both, mapping=aes(x=combined_stage,y=LD1, colour = combined_stage)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4) +
    ylab("Linear \nDiscrimination") +
    theme_bw() +
    ylim(c(min(LD1_both$LD1) - 1, max(LD1_both$LD1 + 2))) +
    scale_colour_manual(values = colvals, labels = collabs) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          legend.position = "none") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 2)
  
  ## Create new ID for Sex that links Status to Sex
  LD1_both$Status_Sex <- ifelse(LD1_both$Status==Target, LD1_both$Sex,
                                Baseline)
  my_comparisons <- list(c(Baseline, "M"),
                         c(Baseline, "F"))
  
  # Colors for sex with either target
  collabs <- c(Baseline, paste(Target, "M", sep = "_"), paste(Target, "F", sep = "_"))
  colvals <- c(colors[colors$Stage=="na",]$Colour, "#807CFF", "#F28333") 
  
  LD1_both$Status_Sex <- factor(LD1_both$Status_Sex, levels = c(Baseline, "M", "F"))
  
  bsex <- ggplot(data = LD1_both, mapping=aes(x=Status_Sex,y=LD1, colour = Status_Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4) +
    ylab("Linear \nDiscrimination") +
    theme_bw() +
    scale_colour_manual(values = colvals, labels = collabs) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          legend.position = "none") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 2)
  
  if (Target=="HCC") {
    # Colours
    collabs <- c("NC","NASH","HCV+","Other")
    colvals <- c(colors[colors$Stage=="na",]$Colour,"blue","red","dark green") 
    
    # Comparisons
    my_comparisons <- list(c("NC","NASH"),c("NC","Other"))
    
    # Subset data for anova test
    sub_dat <- LD1_both[LD1_both$Etiology_group %in% c("NASH","HCV+","Other"),]
    # Calculate one-way anova on this subsetted dataframe
    one.way <- aov(LD1 ~ Etiology_group, data = sub_dat)
    pval <- round(summary(one.way)[[1]][["Pr(>F)"]][[1]], 2)
    
    LD1_both$Etiology_group <- plyr::mapvalues(LD1_both$Etiology_group, c("na"),c("NC"))
    LD1_both$Etiology_group <- factor(LD1_both$Etiology_group, levels = c("NC","NASH","HCV+","Other"))
    
    betiology <- ggplot(data = LD1_both, mapping=aes(x=Etiology_group,y=LD1, colour = Etiology_group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.4) +
      ylab("Linear \nDiscrimination") +
      scale_colour_manual(values = colvals, labels = collabs) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            legend.position = "none") +
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 2, test = "wilcox")
    
    plots[["LD1_boxplot_by_etiology"]] <- betiology
  }
  
  plots[["LD1_boxplot_by_stage"]] <- bstage
  plots[["LD1_boxplot_by_sex"]] <- bsex

  return(plots)
}

HCC_plots <- GenerateROCforSexStageEtiology(Data = datac, Metadata = metadata, 
                                            Baseline = "NC", Target = "HCC", colors = colors)
MM_plots <- GenerateROCforSexStageEtiology(Data = datac, Metadata = metadata, 
                                            Baseline = "NC", Target = "HCC", colors = colors)

## Export plots separately at same size as ROC curves
## HCC
pdf("../ROC/pilot_HCC_LD1_boxplots_sex.pdf", 2, 2)
print(HCC_plots$LD1_boxplot_by_sex)
dev.off()

pdf("../ROC/pilot_HCC_LD1_boxplots_age.pdf", 2, 2)
print(HCC_plots$LD1_boxplot_by_age)
dev.off()

pdf("../ROC/pilot_HCC_LD1_boxplots_stage.pdf", 2, 2)
print(HCC_plots$LD1_boxplot_by_stage)
dev.off()

pdf("../ROC/pilot_HCC_LD1_boxplots_etiology.pdf", 2.8, 2.8)
print(HCC_plots$LD1_boxplot_by_etiology)
dev.off()

## MM
pdf("../ROC/pilot_MM_LD1_boxplots_sex.pdf", 2, 2)
print(MM_plots$LD1_boxplot_by_sex)
dev.off()

pdf("../ROC/pilot_MM_LD1_boxplots_age.pdf", 2, 2)
print(MM_plots$LD1_boxplot_by_age)
dev.off()

pdf("../ROC/pilot_MM_LD1_boxplots_stage.pdf", 2, 2)
print(MM_plots$LD1_boxplot_by_stage)
dev.off()
