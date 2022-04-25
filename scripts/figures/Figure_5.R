# Script to generate plots in Figure 1, excluding the LOOCV tables, which are made separately

# Load libraries
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
library(glmnet)
library(pheatmap)
library(dplyr)

# Set working directory
setwd("/Users/breesheyroskams-hieter/Desktop/cfRNA/manuscript/revised_paper_with_validation/figures/validation/")

############################################
# Import Data
############################################
data <- read.table("../../tables/RPM_validation.txt", sep="\t", stringsAsFactors = F, header = T)
metadata <- read.table("../../tables/validation_metadata.txt",sep = "\t", stringsAsFactors = FALSE, header = T)
colors <- read.delim("../../tables/colors_validation.txt", stringsAsFactors = FALSE)

# DE from DESeq2
DE_HCCvsHD <- read.table(file="../../tables/DE/HCC-vs-HD.diffexp.geneID.tsv",header = TRUE)
DE_MMvsHD <- read.table(file="../../tables/DE/MM-vs-HD.diffexp.geneID.tsv",header = TRUE)

# LVQ genes from lvq script
LVQ_HCCvsHD <- readRDS(file="../../tables/LVQ/HD_HCC_importance.rds")
LVQ_MMvsHD <- readRDS(file="../../tables/LVQ/HD_MM_importance.rds")

############################################
# Format Data
############################################
data <- data[!is.na(data$gene),]
data$gene <- as.character(data$gene)

# Consolidate duplicate gene names
datam <- melt(data,id="gene")
datac <- dcast(datam,gene~variable,fun.aggregate = mean)
rownames(datac) <- datac[,1]
datac <- datac[,-1]

# Filter for good QC samples
datac <- datac[,metadata$Seq_ID]

## Select different gene sets
## select DE gene lists for MM and HCC
DE_MMvsHD_genes <- DE_MMvsHD[DE_MMvsHD$padj<0.01 & abs(DE_MMvsHD$log2FoldChange)>0,]$GeneID
DE_HCCvsHD_genes <- DE_HCCvsHD[DE_HCCvsHD$padj<0.01 & abs(DE_HCCvsHD$log2FoldChange)>0,]$GeneID
LVQ_MMvsHD_genes <- rownames(LVQ_MMvsHD$importance[order(LVQ_MMvsHD$importance$HD, decreasing = TRUE),])[1:10]
LVQ_HCCvsHD_genes <- rownames(LVQ_HCCvsHD$importance[order(LVQ_HCCvsHD$importance$HD, decreasing = TRUE),])[1:10]


############################################
#  PCA of genes selected by variance
############################################

# Import Arial type font
library(extrafont)
extrafont::font_import()
extrafont::fonts()
extrafont::loadfonts()

## Summary of good quality samples
sample_counts <- metadata %>%
  group_by(Group) %>%
  summarise(num_samples = n())

######################################################
# ROC Curves and boxplots for Linear Discrimination
######################################################

## Read in pilot cohort tables to generate model and test validation samples with it
pilot_data <- read.csv("../../tables/RPM_without_dates_updated.csv",sep=",", row.names = 1)
pilot_metadata <- read.csv("../../tables/PP_metadata_keep_FINAL_updated.csv",sep = ",", stringsAsFactors = FALSE)
pilot_colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

## Combined metadata
comb_md <- read.table("../../tables/validation_and_pilot_combined_metadata.txt", header = T, sep = "\t", stringsAsFactors = F)
## Format data
pilot_data$gene <- as.character(pilot_data$gene)

# Consolidate duplicate gene names
pilot_datam <- melt(pilot_data,id="gene")
pilot_datac <- dcast(pilot_datam,gene~variable,fun.aggregate = mean)
rownames(pilot_datac) <- pilot_datac[,1]
pilot_datac <- pilot_datac[,-1]

GenerateROCforAllFourModels <- function(Data, Data_val, Metadata, pilot_Baseline, pilot_Target, validation_Baseline,
                                        validation_Target, BaselineName, TargetName,
                                        DElist, LVQlist, pilot_Colors, validation_Colors) {
  plots <- list()
  ######### LOOCV LDA and DE ###################
  sample_sel <- Metadata[Metadata$cohort=="pilot" & Metadata$group==pilot_Baseline|Metadata$group==pilot_Target,]$id
  datac_sel <- Data[,colnames(Data) %in% sample_sel]
  #make genename become legal name in R
  names(datac_sel) <- make.names(names(datac_sel))
  
  metadata_sel <- Metadata[Metadata$id %in% sample_sel,]
  rownames(metadata_sel)=metadata_sel$id
  metadata_sel=metadata_sel[,c("group","combined_groups")]
  
  # Or Select genes based on DE results
  select_datac= which((row.names(Data)) %in% (DElist))
  
  ## order metadata_sel and datac_sel the same
  metadata_sel <- metadata_sel[order(rownames(metadata_sel)),]
  datac_sel <- datac_sel[,order(colnames(datac_sel))]
  
  stopifnot(colnames(datac_sel)==rownames(metadata_sel))
  
  data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_train <- data_train[,-2]
  data_train$group <- factor(data_train$group, levels = c(pilot_Baseline, pilot_Target))
  table(data_train$group)
  #make genename become legal name in R
  names(data_train) <- make.names(names(data_train))
  
  #####################################
  ## Configure testing data the same
  #####################################
  sample_sel <- Metadata[Metadata$cohort=="validation" & Metadata$group==validation_Baseline|Metadata$group==validation_Target,]$id
  val_datac_sel <- Data_val[,colnames(Data_val) %in% sample_sel]
  #make genename become legal name in R
  names(val_datac_sel) <- make.names(names(val_datac_sel))
  
  val_metadata_sel <- Metadata[Metadata$id %in% colnames(val_datac_sel),]
  rownames(val_metadata_sel)=val_metadata_sel$id
  val_metadata_sel=val_metadata_sel[,c("group","combined_groups")]
  
  # Or Select genes based on DE results
  val_select_datac= which((row.names(Data_val)) %in% (DElist))
  
  ## order val_metadata_sel and val_datac_sel the same
  val_metadata_sel <- val_metadata_sel[order(rownames(val_metadata_sel)),]
  val_datac_sel <- val_datac_sel[,order(colnames(val_datac_sel))]
  
  # Check
  stopifnot(colnames(val_datac_sel)==rownames(val_metadata_sel))
  
  data_test <- cbind(val_metadata_sel,t(log(val_datac_sel[val_select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_test <- data_test[,-2]
  data_test$group <- factor(data_test$group, levels = c(validation_Baseline, validation_Target))
  table(data_test$group)
  #make genename become legal name in R
  names(data_test) <- make.names(names(data_test))
  
  ## How about running it on the whole dataframe?
  model_lda <- lda(group~.,data_train)
  p2 <- predict(model_lda,data_test)$class
  
  pred_out = data.frame(sample = rownames(data_test), pred = as.character(p2), actual_class = data_test$group)
  
  dir.create("Figure1/LDA_DE/", recursive = T)
  filename = paste0("Figure1/LDA_DE/LDA_", validation_Baseline,"_vs_",validation_Target,"_DE0.01.csv")
  write.csv(pred_out, file=filename)
  tabLDA <- table(Predict = pred_out[,2], Actual = pred_out[,3])
  print("LDA DE LOOCV table")
  print(tabLDA)
  print(sum(diag(tabLDA))/sum(tabLDA))
  
  LDA_DE_LOOCV <- as.data.frame(tabLDA)
  
  lda.fit = lda(group ~., data_train)
  #plot(lda.fit)
  lda.pred=predict(lda.fit,data_test)
  table(data_test$group,lda.pred$class)
  LD1_DE <- data.frame(lda.pred$x,data_test$group)
  names(LD1_DE) <- c("LDA_DE", "Status")
  LD1_DE <- melt(LD1_DE, by="Status")
  
  # Save LD1_DE for plotting later
  model.roc.LDA.DE <-roc(data_test$group,lda.pred$posterior[,2])
  print(pROC::auc(model.roc.LDA.DE))
  LDA_DE_AUC <- round(as.numeric(strsplit(as.character(pROC::auc(model.roc.LDA.DE)), " ")[[1]][1]), 2)
  
  ######### LOOCV Random Forest and DE #############
  set.seed(5)
  model_rf <- randomForest(group~.,data=data_train)
  
  gini <- varImp(model_rf)
  gini <- cbind(row.names(gini),gini)
  gini_sort <- gini[rev(order(gini$Overall)),]
  
  pred <- predict(model_rf,data_test)
  pred_out_rf <- data.frame(as.data.frame(pred), actual_class = data_test$group)
  
  dir.create("Figure1/RF_DE/", recursive = T)
  filename = paste0("Figure1/RF_DE/RF_",validation_Baseline,"_vs_",validation_Target,"_DE0.01.csv")
  write.csv(pred_out_rf, file=filename)
  tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
  print("RF DE LOOCV table")
  print(tabrf)
  print(sum(diag(tabrf))/sum(tabrf))
  
  RF_DE_LOOCV <- as.data.frame(tabrf)
  
  # Run ROC for random forest/DE model
  set.seed(5)
  model_train <- randomForest(group~.,data=data_train)
  pred <- predict(model_train, data_test, type = "prob")

  model.roc.RF.DE<-roc(pred_out_rf$actual_class, pred[,2])
  
  print(pROC::auc(model.roc.RF.DE))
  RF_DE_AUC <- round(as.numeric(strsplit(as.character(pROC::auc(model.roc.RF.DE)), " ")[[1]][1]), 2)
  
  ########## LOOCV LDA and lvq ###################
  genes_sel <- rownames(LVQlist$importance[order(LVQlist$importance$HD, decreasing = TRUE),])[1:10]
  
  select_datac= which((row.names(Data)) %in% (genes_sel))
  
  ## order val_metadata_sel and val_datac_sel the same
  metadata_sel <- metadata_sel[order(rownames(metadata_sel)),]
  datac_sel <- datac_sel[,order(colnames(datac_sel))]
  
  # Check
  stopifnot(colnames(datac_sel)==rownames(metadata_sel))
  
  data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_train <- data_train[,-2]
  data_train$group <- factor(data_train$group, levels = c(pilot_Baseline, pilot_Target))
  table(data_train$group)
  # make genename become legal name in R
  names(data_train) <- make.names(names(data_train))
  
  ## Make new testing dataframe
  val_select_datac= which((row.names(Data_val)) %in% (genes_sel))
  
  ## order val_metadata_sel and val_datac_sel the same
  val_metadata_sel <- val_metadata_sel[order(rownames(val_metadata_sel)),]
  val_datac_sel <- val_datac_sel[,order(colnames(val_datac_sel))]
  
  # Check
  stopifnot(colnames(val_datac_sel)==rownames(val_metadata_sel))
  
  data_test <- cbind(val_metadata_sel,t(log(val_datac_sel[val_select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_test <- data_test[,-2]
  data_test$group <- factor(data_test$group, levels = c(validation_Baseline, validation_Target))
  table(data_test$group)
  #make genename become legal name in R
  names(data_test) <- make.names(names(data_test))
  
  ## How about running it on the whole dataframe?
  model_lda <- lda(group~.,data_train)
  p2 <- predict(model_lda,data_test)$class
  
  pred_out = data.frame(sample = rownames(data_test), pred = as.character(p2), actual_class = data_test$group)
  
  dir.create("Figure1/LDA_LVQ/", recursive = T)
  filename = paste0("Figure1/LDA_LVQ/LDA_",validation_Baseline,"_vs_",validation_Target,"_lvq.csv")
  write.csv(pred_out, file=filename)
  
  tabLDA <- table(Predict = pred_out[,2], Actual = pred_out[,3])
  print("LDA lvq LOOCV table")
  print(tabLDA)
  print(sum(diag(tabLDA))/sum(tabLDA))
  
  LDA_LVQ_LOOCV <- as.data.frame(tabLDA)
  
  lda.fit = lda(group ~., data_train)
  #plot(lda.fit)
  lda.pred=predict(lda.fit,data_test)
  table(data_test$group,lda.pred$class)
  LD1_lvq <- data.frame(lda.pred$x,data_test$group)
  names(LD1_lvq) <- c("LDA_LVQ", "Status")
  LD1_lvq <- melt(LD1_lvq, by="Status")
  
  LD1_both <- rbind(LD1_DE, LD1_lvq)
  names(LD1_both) <- c("Status","Model","Value")
  
  LD1_both$Model <- plyr::mapvalues(LD1_both$Model, c("LDA_DE", "LDA_LVQ"), c("LDA/DE","LDA/LVQ"))
  
  # Add new type assignments to match pilot cohort
  LD1_both$Group <- plyr::mapvalues(LD1_both$Status, c(validation_Target), c(pilot_Target))
  
  # Find colours for baseline and target
  colourBaselinePilot <- pilot_Colors[pilot_Colors$Status==pilot_Baseline,]$Colour
  colourTargetPilot <- pilot_Colors[pilot_Colors$Status==pilot_Target,]$Colour
  
  ############ Boxplot of linear discrimination ###############
  my_comparisons <- list(c(as.vector(unique(LD1_both$Group))))
  
  LD1_both$Group <- factor(LD1_both$Group, levels = c("NC", pilot_Target))
  
  # Find colours for baseline and target
  colourBaseline <- validation_Colors[validation_Colors$Status=="NC",]$Colour
  colourTarget <- validation_Colors[validation_Colors$Status==pilot_Target,]$Colour
  
  ## Remove LD/DE
  filt <- LD1_both[LD1_both$Model=="LDA/LVQ",]
  
  b <- ggplot(data = filt, mapping=aes(x=Group,y=Value, colour = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4) +
    ylab("Linear \nDiscrimination") +
    facet_grid(~Model) +
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25,0.15,0.2,0.15),"cm"),
          #aspect.ratio = 1,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.title.y = element_text(size = 8),
          legend.position = "none") +
    scale_color_manual(breaks=c("NC", pilot_Target),
                       values=c(colourBaseline, colourTarget)) +
    scale_y_continuous(breaks = scales::pretty_breaks(), 
                       limits = c(min(LD1_both$Value), max(LD1_both$Value) + 2)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 2)
  
  b2 <- ggplot(data = filt, mapping=aes(x=Group,y=Value, colour = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4) +
    ylab("Linear \nDiscrimination") +
    facet_grid(~Model) +
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5),
          #plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
          #aspect.ratio = 1,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.title.y = element_text(size = 8),
          legend.position = "none") +
    scale_color_manual(breaks=c("NC", pilot_Target),
                       values=c(colourBaseline, colourTarget)) +
    scale_y_continuous(breaks = scales::pretty_breaks(), 
                       limits = c(min(LD1_both$Value), max(LD1_both$Value) + 2)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.label", size = 2)

  
  plots[["boxplot"]] <- b
  plots[["boxplot_with_pvalue"]] <- b2
  
  model.roc.LDA.lvq <-roc(data_test$group,lda.pred$posterior[,2])
  print(pROC::auc(model.roc.LDA.lvq))
  LDA_lvq_AUC <- round(as.numeric(strsplit(as.character(pROC::auc(model.roc.LDA.lvq)), " ")[[1]][1]), 2)
  
  ######### LOOCV Random Forest and lvq  ###################
  set.seed(5)
  model_rf <- randomForest(group~.,data=data_train)
  
  gini <- varImp(model_rf)
  gini <- cbind(row.names(gini),gini)
  gini_sort <- gini[rev(order(gini$Overall)),]
  
  pred <- predict(model_rf,data_test)
  pred_out_rf <- data.frame(as.data.frame(pred), actual_class = data_test$group)
  
  dir.create("Figure1/RF_LVQ/", recursive = T)
  filename = paste0("Figure1/RF_LVQ/RF_",validation_Baseline,"_vs_",validation_Target,"_lvq.csv")
  write.csv(pred_out_rf, file=filename)
  tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
  print("RF DE LOOCV table")
  print(tabrf)
  print(sum(diag(tabrf))/sum(tabrf))
  
  RF_LVQ_LOOCV <- as.data.frame(tabrf)
  
  # Run ROC for random forest/DE model
  set.seed(5)
  model_train <- randomForest(group~.,data=data_train)
  pred <- predict(model_train, data_test, type = "prob")
  
  model.roc.RF.lvq<-roc(pred_out_rf$actual_class, pred[,2])
  
  print(pROC::auc(model.roc.RF.lvq))
  RF_lvq_AUC <- round(as.numeric(strsplit(as.character(pROC::auc(model.roc.RF.lvq)), " ")[[1]][1]), 2)
  ##############################################################
  
  # Plot all four ROC curves on top of each other
  g <- ggroc(data = list(LDA_lvq = model.roc.LDA.lvq, RF_lvq = model.roc.RF.lvq)
             , size = 0.7, alpha = 0.8, legacy.axes = TRUE) +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme_bw() +
    #geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5),
          #plot.margin = unit(c(0,0.15,0.65,0.15),"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1,
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.box.background = element_rect(),
          #legend.margin = unit(c(0,0,0,0),"cm"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 5)) +
    scale_colour_manual(values = c("#BFAE0A","blue"),
                        labels = c(paste0("LDA/LVQ (AUC: ", LDA_lvq_AUC, ")"),
                                   paste0("RF/LVQ (AUC: ", RF_lvq_AUC, ")")))
  legend <- get_legend(g)
  
  plots[["ROC"]] <- g + theme(legend.position = "none")
  plots[["legend"]] <- legend
  
  pdf(paste0("ROC_legend-", "NC", "-vs-", pilot_Target, ".pdf"), 5, 2)
  grid.newpage()
  grid.draw(legend)
  dev.off()
  
  return(plots)
}

ROC_MM <- GenerateROCforAllFourModels(Data=pilot_datac, Data_val=datac, Metadata=comb_md, pilot_Baseline="NC", pilot_Target="MM", 
                                       validation_Baseline="NC",validation_Target="MMV",
                                       DElist=DE_MMvsHD_genes, LVQlist=LVQ_MMvsHD, pilot_Colors = pilot_colors,validation_Colors = colors)
ROC_HCC <- GenerateROCforAllFourModels(Data=pilot_datac, Data_val=datac, Metadata=comb_md, pilot_Baseline="NC", pilot_Target="HCC", 
                                       validation_Baseline="NC",validation_Target="HCCT_A",
                                       DElist=DE_HCCvsHD_genes, LVQlist=LVQ_HCCvsHD, pilot_Colors = pilot_colors,validation_Colors = colors)

pdf("ROC_curve_and_LDA_boxplot.pdf", 3.5, 4, useDingbats = FALSE)
cowplot::plot_grid(ROC_MM$boxplot, ROC_MM$ROC,
                   ROC_HCC$boxplot, ROC_HCC$ROC,
                   rel_widths = c(1,1.3),ncol = 2)
dev.off()
