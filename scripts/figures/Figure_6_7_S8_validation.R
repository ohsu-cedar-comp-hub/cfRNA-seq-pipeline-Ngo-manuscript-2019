## LD1 boxplots by sex, stage, etiology and age with validation cohort tested on training data

library(reshape2)
library(tsne)
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

val_data <- read.table("../../tables/RPM_validation.txt", sep="\t", stringsAsFactors = F, header = T)
metadata <- read.table("../../tables/combined_metadata_with_stage.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
colors <- read.delim("../../tables/colors_validation.txt", stringsAsFactors = FALSE)

# DE from DESeq2
DE_HCCvsHD <- read.table(file="../../tables/DE/HCC-vs-HD.diffexp.geneID.tsv",header = TRUE)
DE_MMvsHD <- read.table(file="../../tables/DE/MM-vs-HD.diffexp.geneID.tsv",header = TRUE)

# LVQ genes from lvq script
LVQ_HCCvsHD <- readRDS(file="../../tables/LVQ/HD_HCC_importance.rds")
LVQ_MMvsHD <- readRDS(file="../../tables/LVQ/HD_MM_importance.rds")

# Colors
colors <- read.delim("../../tables/colors_by_stage.txt", header = T)

############################################
# Format Data
############################################

val_data <- val_data[!is.na(val_data$gene),]
val_data$gene <- as.character(val_data$gene)

# Consolidate duplicate gene names
val_datam <- reshape2::melt(val_data,id="gene")
val_datac <- reshape2::dcast(val_datam,gene~variable,fun.aggregate = mean)
rownames(val_datac) <- val_datac[,1]
val_datac <- val_datac[,-1]

## Select only pilot samples
val_metadata <- metadata[metadata$Set=="Validation",]
val_metadata <- val_metadata[order(val_metadata$SeqID),]
val_datac <- val_datac[,order(colnames(val_datac))]

### PILOT DATA
## Read in pilot cohort tables to generate model and test validation samples with it
pilot_data <- read.csv("../../tables/RPM_without_dates_updated.csv",sep=",", row.names = 1)
pilot_colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)
pilot_metadata <- metadata[metadata$Set=="Pilot",]

## Format data
pilot_data$gene <- as.character(pilot_data$gene)

# Consolidate duplicate gene names
pilot_datam <- reshape2::melt(pilot_data,id="gene")
pilot_datac <- reshape2::dcast(pilot_datam,gene~variable,fun.aggregate = mean)
pilot_datac <- pilot_datac[-1,]
rownames(pilot_datac) <- pilot_datac[,1]
pilot_datac <- pilot_datac[,-1]

#######################################################
## LDA boxplots split by Stage, Etiology and Sex
#######################################################

GenerateROCforSexStageEtiology <- function(Data_pilot, Data_val, Metadata_pilot, Metadata_val, Baseline, Target, colors) {
  if (Target=="HCC") {
    DElist=DE_HCCvsHD
    LVQlist=LVQ_HCCvsHD
  } else if (Target=="MM") {
    DElist=DE_MMvsHD
    LVQlist=LVQ_MMvsHD
  }
  
  plots <- list()
  ######### LOOCV LDA and DE ###################
  sample_sel <- Metadata_pilot[Metadata_pilot$Status==Baseline|Metadata_pilot$Status==Target,]$SeqID
  datac_sel <- Data_pilot[,colnames(Data_pilot) %in% sample_sel]
  #make genename become legal name in R
  names(datac_sel) <- make.names(names(datac_sel))
  
  metadata_info <- Metadata_pilot[Metadata_pilot$SeqID %in% sample_sel,]
  rownames(metadata_info)=metadata_info$SeqID
  
  ## Order datac_sel and metadata_info the same
  datac_sel <- datac_sel[,order(colnames(datac_sel))]
  metadata_info <- metadata_info[order(metadata_info$SeqID),]
  metadata_sel=metadata_info[,c(4,5)]
  
  #####################################
  ## Configure testing data the same
  #####################################
  sample_sel <- Metadata_val[Metadata_val$Status==Baseline|Metadata_val$Status==Target,]$SeqID
  val_datac_sel <- Data_val[,colnames(Data_val) %in% sample_sel]
  names(val_datac_sel) <- make.names(names(val_datac_sel))
  
  val_metadata_info <- Metadata_val[Metadata_val$SeqID %in% sample_sel,]
  rownames(val_metadata_info)=val_metadata_info$SeqID
  
  ## Order datac_sel and metadata_info the same
  val_datac_sel <- val_datac_sel[,order(colnames(val_datac_sel))]
  val_metadata_info <- val_metadata_info[order(val_metadata_info$SeqID),]
  val_metadata_sel <- val_metadata_info[,c(4,5)]
  
  # Check
  stopifnot(colnames(val_datac_sel)==rownames(val_metadata_sel))
  
  ########## LOOCV LDA and lvq ###################
  ## LVQ and pilot dataset to make training data
  genes_sel <- rownames(LVQlist$importance[order(LVQlist$importance$HD, decreasing = TRUE),])[1:10]
  
  select_datac= which((row.names(Data_pilot)) %in% (genes_sel))
  
  data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_train <- data_train[,-1]
  data_train$Status <- factor(data_train$Status, levels = c(Baseline, Target))
  table(data_train$Status)
  names(data_train) <- make.names(names(data_train))
  
  ## Now configure testing data
  select_datac= which((row.names(val_datac_sel)) %in% (genes_sel))
  
  data_test <- cbind(val_metadata_sel,t(log(val_datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_test <- data_test[,-1]
  data_test$Status <- factor(data_test$Status, levels = c(Baseline, Target))
  table(data_test$Status)
  
  # Fit LDA model
  lda.fit = lda(Status ~., data_train)
  lda.pred=predict(lda.fit,data_test) # test on validation set
  table(data_test$Status,lda.pred$class)
  
  LD1_lvq <- data.frame(LDA_lvq = lda.pred$x,Status = data_test$Status)
  LD1_lvq$SeqID <- rownames(LD1_lvq)
  LD1_lvq <- reshape2::melt(LD1_lvq, by=c("Status","SeqID"))
  ## Add in stage, etiology and type information
  iv <- match(LD1_lvq$SeqID, val_metadata_info$SeqID)
  LD1_lvq <- data.frame(LD1_lvq, val_metadata_info[iv,c("Sex","Stage","Etiology_group","Age")])
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
    my_comparisons <- list(c("NC","A"),c("NC","B,C"))
  }
  
  LD1_both$combined_stage <- factor(LD1_both$combined_stage, 
                                    levels = c("NC",as.vector(stagecomp[1]),
                                               paste(stagecomp[2:length(stagecomp)], collapse = ",")))
  
  # Colours for combined stages and either HCC or MM
  if (Target == "HCC") {
    collabs <- c("NC","A","B,C,D")
    colvals <- c(colors[colors$Stage=="na",]$Colour, 
                 colors[colors$Stage=="A",]$Colour,
                 colors[colors$Stage=="B_C_D",]$Colour)
  } else if (Target == "MM") {
    collabs <- c("NC","I","II,II")
    colvals <- c(colors[colors$Stage=="na",]$Colour, 
                 colors[colors$Stage=="I",]$Colour,
                 colors[colors$Stage=="II_III",]$Colour)
  }
  
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
  
  ## merged ID for Sex that links Status to Sex
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
          #aspect.ratio = 1,
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
    collabs <- c("na","HCV+","Other")
    colvals <- c(colors[colors$Stage=="na",]$Colour,"red","dark green") 
    
    # Comparisons
    my_comparisons <- list(c("NC","HCV+"),c("NC","Other"))
    
    # Subset data for anova test
    sub_dat <- LD1_both[LD1_both$Etiology_group %in% c("HCV+","Other"),]
    # Calculate one-way anove on this subsetted dataframe
    one.way <- aov(LD1 ~ Etiology_group, data = sub_dat)
    pval <- round(summary(one.way)[[1]][["Pr(>F)"]][[1]], 2)
    
    # Factor both
    LD1_both$Etiology_group <- plyr::mapvalues(LD1_both$Etiology_group, c("na"),c("NC"))
    LD1_both$Etiology_group <- factor(LD1_both$Etiology_group, levels = c("NC","HCV+","Other"))
    
    betiology <- ggplot(data = LD1_both, mapping=aes(x=Etiology_group,y=LD1, colour = Etiology_group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.4) +
      ylab("Linear \nDiscrimination") +
      scale_colour_manual(values = colvals, labels = collabs) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
            #aspect.ratio = 1,
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
  
  ## What about by age as a continuous variable?
  bage <- ggplot(data = LD1_both, mapping=aes(x=Age,y=LD1, colour = Status)) +
    geom_jitter(width = 0.2, size = 1.5) +
    ylab("Linear \nDiscrimination") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.title.y = element_text(size = 8),
          legend.position = "none") 
  
  plots[["LD1_boxplot_by_stage"]] <- bstage
  plots[["LD1_boxplot_by_sex"]] <- bsex
  plots[["LD1_boxplot_by_age"]] <- bage
  
  return(plots)
}

HCC_plots <- GenerateROCforSexStageEtiology(Data_pilot=pilot_datac, Data_val=val_datac, 
                                            Metadata_pilot=pilot_metadata, Metadata_val=val_metadata,
                                            Baseline = "NC", Target = "HCC", colors = colors)
MM_plots <- GenerateROCforSexStageEtiology(Data_pilot=pilot_datac, Data_val=val_datac, 
                                            Metadata_pilot=pilot_metadata, Metadata_val=val_metadata,
                                            Baseline = "NC", Target = "MM", colors = colors)

## Export plots separately at same size as ROC curves
## HCC
pdf("../ROC/validation_HCC_LD1_boxplots_sex.pdf", 2, 2)
print(HCC_plots$LD1_boxplot_by_sex)
dev.off()

pdf("../ROC/validation_HCC_LD1_boxplots_age.pdf", 2, 2)
print(HCC_plots$LD1_boxplot_by_age)
dev.off()

pdf("../ROC/validation_HCC_LD1_boxplots_stage.pdf", 2, 2)
print(HCC_plots$LD1_boxplot_by_stage)
dev.off()

pdf("../ROC/validation_HCC_LD1_boxplots_etiology.pdf", 2.8, 2.8)
print(HCC_plots$LD1_boxplot_by_etiology)
dev.off()

## MM
pdf("../ROC/validation_MM_LD1_boxplots_sex.pdf", 2, 2)
print(MM_plots$LD1_boxplot_by_sex)
dev.off()

pdf("../ROC/validation_MM_LD1_boxplots_age.pdf", 2, 2)
print(MM_plots$LD1_boxplot_by_age)
dev.off()

pdf("../ROC/validation_MM_LD1_boxplots_stage.pdf", 2, 2)
print(MM_plots$LD1_boxplot_by_stage)
dev.off()
