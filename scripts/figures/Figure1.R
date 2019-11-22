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

# Set working directory
setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/figures/figure1/")

############################################
# Import Data
############################################

data <- read.csv("../../tables/RPM_without_dates_updated.csv",sep=",", row.names = 1)
metadata <- read.csv("../../tables/PP_metadata_keep_FINAL_updated.csv",sep = ",", stringsAsFactors = FALSE)
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

# DE from DESeq2
DE_LuCavsHD <- read.table(file="../../tables/DE/LuCa-vs-HD.diffexp.geneID.tsv",header = TRUE,stringsAsFactors = F)
DE_HCCvsHD <- read.table(file="../../tables/DE/HCC-vs-HD.diffexp.geneID.tsv",header = TRUE,stringsAsFactors = F)
DE_MMvsHD <- read.table(file="../../tables/DE/MM-vs-HD.diffexp.geneID.tsv",header = TRUE,stringsAsFactors = F)

# LVQ genes from lvq script
LVQ_LuCavsHD <- readRDS(file="../../tables/LVQ/HD_LuCa_importance.rds")
LVQ_HCCvsHD <- readRDS(file="../../tables/LVQ/HD_HCC_importance.rds")
LVQ_MMvsHD <- readRDS(file="../../tables/LVQ/HD_MM_importance.rds")

############################################
# Format Data
############################################

data$gene <- as.character(data$gene)

# Consolidate duplicate gene names
datam <- melt(data,id="gene")
datac <- dcast(datam,gene~variable,fun.aggregate = mean)
rownames(datac) <- datac[,1]
datac <- datac[,-1]

############################################
#  PCA of genes selected by variance
############################################

# Import Arial type font
library(extrafont)
extrafont::font_import()
extrafont::fonts()
extrafont::loadfonts()

# PCA function
GeneratePCA_n_genes <- function(Data, Metadata, Baseline, Target, Colors, n=500) {
  plots <- list()
  sample_sel <- Metadata[Metadata$Status==Baseline|Metadata$Status==Target,]$PP_ID
  datac_sel <- Data[,colnames(Data) %in% sample_sel]
  # Make genename become legal name in R
  names(datac_sel) <- make.names(names(datac_sel))
  # Transform data to log scale
  genecount <- log(datac_sel+1,2)
  # Which values in genecount are greater than 1.5?
  thresh <- genecount  > 1.5
  # This produces a logical matrix with TRUEs and FALSEs
  head(thresh)
  # Summary of how many TRUEs there are in each row
  table(rowSums(thresh))
  # we would like to keep genes that have at least 2 TRUES in each row of thresh
  keep <- rowSums(thresh) >= 5
  # Subset the rows of countdata to keep the more highly expressed genes
  counts.keep <- genecount[keep,]
  summary(keep)
  dim(counts.keep)
  # We estimate the variance for each row in the logcounts matrix
  var_genes <- apply(genecount, 1, var)
  head(var_genes)
  # Get the gene names for the top 500 most variable genes
  select_var <- names(sort(var_genes, decreasing=TRUE))[1:n]
  head(select_var)
  # Subset logcounts matrix
  highly_variable_lcpm <- genecount[select_var,]
  dim(highly_variable_lcpm)
  head(highly_variable_lcpm)
  Genes_wo_hit <- apply(highly_variable_lcpm,1,function(x) all(x==0))
  genecount_nozero <- highly_variable_lcpm[!Genes_wo_hit,]
  genecount_PCA <- t(genecount_nozero)
  
  # Run PCA
  PCA_out <- prcomp(genecount_PCA, scale. = TRUE)
  metadata_sel <- Metadata[Metadata$PP_ID %in% sample_sel,]
  rownames(metadata_sel)=metadata_sel$PP_ID
  metadata_sel=metadata_sel[,c(4,5)]
  Groups <- metadata_sel$Status
  # Generate the summary stats for prcomp
  eigs <- PCA_out$sdev^2
  summary <- rbind(
    SD = sqrt(eigs),
    Proportion = eigs/sum(eigs),
    Cumulative = cumsum(eigs)/sum(eigs))
  # Find the proprtion of PC1 and PC2
  propPC1 <- round(as.numeric(summary[2,1]) * 100)
  propPC2 <- round(as.numeric(summary[2,2]) * 100)
  # Find colours for baseline and target
  colourBaseline <- Colors[Colors$Status==Baseline,]$Colour
  colourTarget <- Colors[Colors$Status==Target,]$Colour
  # Plot PCA
  my.pca <- data.frame(PCA_out$x[,1:3], Groups)
  my.pca$Groups <- factor(my.pca$Groups, levels = c(Baseline, Target))
  
  p <- ggplot(data = my.pca, mapping = aes(x = PC1,y = PC2, colour = Groups)) +
    geom_point(size = 2) +
    theme_bw() +
    xlab(paste0("PC-1 (", propPC1, "%)")) +
    ylab(paste0("PC-2 (", propPC2, "%)")) +
    scale_color_manual(breaks=c(Baseline, Target),
                       values=c(colourBaseline, colourTarget)) +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.15,0.15,0.15,0.15),"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1,
          legend.title = element_blank(),
          legend.position="bottom",
          legend.direction="horizontal",
          legend.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 5)) 
  
  legend <- cowplot::get_legend(p)
  
  pdf(paste0("../figure1/Legend-", Target, ".pdf"), 5, 2, useDingbats = FALSE)
  grid.newpage()
  grid.draw(legend)
  dev.off()
  
  plots[["plot"]] <- p + theme(legend.position = "none")
  return(plots)
}

LuCa_PCA <- GeneratePCA_n_genes(Data = datac, Metadata = metadata, Baseline = "HD", Target = "LuCa", 
                    Colors = colors, n = 500)
HCC_PCA <- GeneratePCA_n_genes(Data = datac, Metadata = metadata, Baseline = "HD", Target = "HCC", 
                    Colors = colors, n = 500)
MM_PCA <- GeneratePCA_n_genes(Data = datac, Metadata = metadata, Baseline = "HD", Target = "MM", 
                    Colors = colors, n = 500)

######################################################
# ROC Curves and boxplots for Linear Discrimination
######################################################

GenerateROCforAllFourModels <- function(Data, Metadata, Baseline, Target, BaselineName, TargetName,
                                        DElist, LVQlist, Colors) {
  plots <- list()
  ######### LOOCV LDA and DE ###################
  sample_sel <- Metadata[Metadata$Status==Baseline|Metadata$Status==Target,]$PP_ID
  datac_sel <- Data[,colnames(Data) %in% sample_sel]
  #make genename become legal name in R
  names(datac_sel) <- make.names(names(datac_sel))
  
  metadata_sel <- Metadata[Metadata$PP_ID %in% sample_sel,]
  rownames(metadata_sel)=metadata_sel$PP_ID
  metadata_sel=metadata_sel[,c(4,2)]
  
  # Or Select genes based on DE results
  sig_genes=DElist[DElist$padj<0.01 & abs(DElist$log2FoldChange)>0,]
  select_datac= which((row.names(Data)) %in% (sig_genes$GeneID))
  
  data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_train <- data_train[,-2]
  data_train$Status <- factor(data_train$Status, levels = c(Baseline, Target))
  table(data_train$Status)
  #make genename become legal name in R
  names(data_train) <- make.names(names(data_train))
  
  pred_out <- data.frame(sample = character(), pred = character(), actual_class = character())
  for(i in 1:nrow( data_train)){
    training <- data_train[-i,]
    testing <- data_train[i,]
    model_lda <- lda(Status~.,training)
    p2 <- predict(model_lda,testing)$class
    pred_out <- rbind(pred_out,cbind(rownames(testing),data.frame(p2),testing$Status))
  }
  
  dir.create("../../tables/LDA_DE/")
  filename = paste0("../../tables/LDA_DE/LDA_", Baseline,"_vs_",Target,"_DE0.01.csv")
  write.csv(pred_out, file=filename)
  tabLDA <- table(Predict = pred_out[,2], Actual = pred_out[,3])
  print("LDA DE LOOCV table")
  print(tabLDA)
  print(sum(diag(tabLDA))/sum(tabLDA))
  
  LDA_DE_LOOCV <- as.data.frame(tabLDA)
  
  lda.fit = lda(Status ~., data_train)
  #plot(lda.fit)
  lda.pred=predict(lda.fit,data_train)
  table(data_train$Status,lda.pred$class)
  LD1_DE <- data.frame(lda.pred$x,data_train$Status)
  names(LD1_DE) <- c("LDA_DE", "Status")
  LD1_DE <- melt(LD1_DE, by="Status")
  
  # Save LD1_DE for plotting later
  model.roc.LDA.DE <-roc(data_train$Status,lda.pred$posterior[,2])
  print(pROC::auc(model.roc.LDA.DE))
  LDA_DE_AUC <- round(as.numeric(strsplit(as.character(pROC::auc(model.roc.LDA.DE)), " ")[[1]][1]), 2)
  
  ######### LOOCV Random Forest and DE #############
  pred_out_rf <- data.frame(sample = character(), pred = character(), actual_class = character())
  
  for(i in 1:nrow(data_train)){
    
    training <- data_train[-i,]
    validation <- data_train[i,]
    
    set.seed(5)
    model_rf <- randomForest(Status~.,data=training)
    
    gini <- varImp(model_rf)
    gini <- cbind(row.names(gini),gini)
    gini_sort <- gini[rev(order(gini$Overall)),]
    
    pred <- predict(model_rf,validation)
    pred_out_rf <- rbind(pred_out_rf,cbind(data.frame(pred),validation$Status))
    
  }
  
  dir.create("../../tables/RF_DE/")
  filename = paste("../../tables/RF_DE/RF_HD_vs_LG_DE0.01.csv")
  write.csv(pred_out_rf, file=filename)
  tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
  print("RF DE LOOCV table")
  print(tabrf)
  print(sum(diag(tabrf))/sum(tabrf))
  
  RF_DE_LOOCV <- as.data.frame(tabrf)
  
  # Run ROC for random forest/DE model
  set.seed(5)
  model_train <- randomForest(Status~.,data=data_train)
  
  model.roc.RF.DE<-roc(data_train$Status,model_train$votes[,2])
  print(pROC::auc(model.roc.RF.DE))
  RF_DE_AUC <- round(as.numeric(strsplit(as.character(pROC::auc(model.roc.RF.DE)), " ")[[1]][1]), 2)
  
  ########## LOOCV LDA and lvq ###################
  genes_sel <- rownames(LVQlist$importance[order(LVQlist$importance[[Baseline]], decreasing = TRUE),])[1:10]
  
  select_datac= which((row.names(datac)) %in% (genes_sel))
  
  data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_train <- data_train[,-2]
  data_train$Status <- factor(data_train$Status, levels = c(Baseline, Target))
  table(data_train$Status)
  
  # make genename become legal name in R
  names(data_train) <- make.names(names(data_train))
  pred_out <- data.frame(sample = character(), pred = character(), actual_class = character())
  for(i in 1:nrow( data_train)){
    training <- data_train[-i,]
    testing <- data_train[i,]
    model_lda <- lda(Status~.,training)
    p2 <- predict(model_lda,testing)$class
    pred_out <- rbind(pred_out,cbind(rownames(testing),data.frame(p2),testing$Status))
  }
  
  dir.create("../../tables/LDA_LVQ/")
  filename = paste0("../../tables/LDA_LVQ/LDA_",Baseline,"_vs_",Target,"_lvq.csv")
  write.csv(pred_out, file=filename)
  
  tabLDA <- table(Predict = pred_out[,2], Actual = pred_out[,3])
  print("LDA lvq LOOCV table")
  print(tabLDA)
  print(sum(diag(tabLDA))/sum(tabLDA))
  
  LDA_LVQ_LOOCV <- as.data.frame(tabLDA)
  
  lda.fit = lda(Status ~., data_train)
  #plot(lda.fit)
  lda.pred=predict(lda.fit,data_train)
  table(data_train$Status,lda.pred$class)
  
  LD1_lvq <- data.frame(lda.pred$x,data_train$Status)
  names(LD1_lvq) <- c("LDA_LVQ", "Status")
  LD1_lvq <- melt(LD1_lvq, by="Status")
  
  LD1_both <- rbind(LD1_DE, LD1_lvq)
  names(LD1_both) <- c("Status","Model","Value")
  
  LD1_both$Model <- plyr::mapvalues(LD1_both$Model, c("LDA_DE", "LDA_LVQ"), c("LDA/DE","LDA/LVQ"))
  
  # Find colours for baseline and target
  colourBaseline <- Colors[Colors$Status==Baseline,]$Colour
  colourTarget <- Colors[Colors$Status==Target,]$Colour
  
  ############ Boxplot of linear discrimination ###############
  my_comparisons <- list(c(as.vector(unique(LD1_both$Status))))
  
  LD1_both$Status <- factor(LD1_both$Status, levels = c(Baseline, Target))
  
  b <- ggplot(data = LD1_both, mapping=aes(x=Status,y=Value, colour = Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4) +
    ylab("Linear \nDiscrimination") +
    facet_grid(~Model) +
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
          #aspect.ratio = 1,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          legend.position = "none") +
    scale_color_manual(breaks=c(Baseline, Target),
                       values=c(colourBaseline, colourTarget)) +
    scale_y_continuous(breaks = scales::pretty_breaks(), 
                       limits = c(min(LD1_both$Value), max(LD1_both$Value) + 2)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 2)
  
  b2 <- ggplot(data = LD1_both, mapping=aes(x=Status,y=Value, colour = Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4) +
    ylab("Linear \nDiscrimination") +
    facet_grid(~Model) +
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25,0.15,0.52,0.15),"cm"),
          #aspect.ratio = 1,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          legend.position = "none") +
    scale_color_manual(breaks=c(Baseline, Target),
                       values=c(colourBaseline, colourTarget)) +
    scale_y_continuous(breaks = scales::pretty_breaks(), 
                       limits = c(min(LD1_both$Value), max(LD1_both$Value) + 2)) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.label", size = 2)

  
  plots[["boxplot"]] <- b
  plots[["boxplot_with_pvalue"]] <- b2
  
  model.roc.LDA.lvq<-roc(data_train$Status,lda.pred$posterior[,2])
  print(pROC::auc(model.roc.LDA.lvq))
  LDA_lvq_AUC <- round(as.numeric(strsplit(as.character(pROC::auc(model.roc.LDA.lvq)), " ")[[1]][1]), 2)
  
  ######### LOOCV Random Forest and lvq  ###################
  pred_out_rf <- data.frame(sample = character(), pred = character(), actual_class = character())
  
  for(i in 1:nrow(data_train)){
    
    training <- data_train[-i,]
    validation <- data_train[i,]
    
    set.seed(5)
    model_rf <- randomForest(Status~.,data=training)
    
    gini <- varImp(model_rf)
    gini <- cbind(row.names(gini),gini)
    gini_sort <- gini[rev(order(gini$Overall)),]
    
    pred <- predict(model_rf,validation)
    pred_out_rf <- rbind(pred_out_rf,cbind(data.frame(pred),validation$Status))
    
  }
  
  dir.create("../../tables/RF_LVQ/")
  filename = paste0("../../tables/RF_LVQ/RF_",Baseline,"_vs_",Target,"_lvq.csv")
  write.csv(pred_out_rf, file=filename)
  
  tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
  print("RF lvq LOOCV table")
  print(tabrf)
  print(sum(diag(tabrf))/sum(tabrf))
  
  RF_LVQ_LOOCV <- as.data.frame(tabrf)
  
  # ROC on random Forest model
  set.seed(5)
  model_train <- randomForest(Status~.,data=data_train)
  model.roc.RF.lvq<-roc(data_train$Status,model_train$votes[,2])
  print(pROC::auc(model.roc.RF.lvq))
  RF_lvq_AUC <- round(as.numeric(strsplit(as.character(pROC::auc(model.roc.RF.lvq)), " ")[[1]][1]), 2)
  
  # Plot all four ROC curves on top of each other
  g <- ggroc(data = list(LDA_DE = model.roc.LDA.DE, LDA_lvq = model.roc.LDA.lvq, RF_DE = model.roc.RF.DE, RF_lvq = model.roc.RF.lvq)
             , size = 0.7, alpha = 0.8, legacy.axes = TRUE) +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.15,0.15,0.15,0.15),"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1,
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.box.background = element_rect(),
          #legend.margin = unit(c(0,0,0,0),"cm"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 5)) +
    scale_colour_manual(values = c("#15C2FF","#BFAE0A","#7A4924","blue"),
                        labels = c(paste0("LDA/DE (AUC: ", LDA_DE_AUC, ")"), paste0("LDA/LVQ (AUC: ", LDA_lvq_AUC, ")"),
                                   paste0("RF/DE (AUC: ", RF_DE_AUC, ")"), paste0("RF/LVQ (AUC: ", RF_lvq_AUC, ")")))
  
  legend <- get_legend(g)
  
  plots[["ROC"]] <- g + theme(legend.position = "none")
  plots[["legend"]] <- legend
  
  pdf(paste0("ROC_legend-", Target, ".pdf"), 5, 2)
  grid.newpage()
  grid.draw(legend)
  dev.off()
  
  return(plots)
}

ROC_LuCa <- GenerateROCforAllFourModels(Data=datac, Metadata=metadata, Baseline="HD", Target="LuCa", 
                            DElist=DE_LuCavsHD, LVQlist=LVQ_LuCavsHD, Colors = colors)
ROC_MM <- GenerateROCforAllFourModels(Data=datac, Metadata=metadata, Baseline="HD", Target="MM", 
                            DElist=DE_MMvsHD, LVQlist=LVQ_MMvsHD, Colors = colors)
ROC_HCC <- GenerateROCforAllFourModels(Data=datac, Metadata=metadata, Baseline="HD", Target="HCC", 
                            DElist=DE_HCCvsHD, LVQlist=LVQ_HCCvsHD, Colors = colors)

pdf("../figure1/Figure1_A_I.pdf", 5.5, 4.5, useDingbats = FALSE)
cowplot::plot_grid(LuCa_PCA$plot, ROC_LuCa$boxplot, ROC_LuCa$ROC,
          MM_PCA$plot, ROC_MM$boxplot, ROC_MM$ROC,
          HCC_PCA$plot, ROC_HCC$boxplot, ROC_HCC$ROC,
          ncol = 3, rel_widths = c(1, 1.75, 1))
dev.off()

pdf("../figure1/LD1_boxplots_w_pvalue.pdf", 8, 2)
cowplot::plot_grid(ROC_LuCa$boxplot_with_pvalue, ROC_MM$boxplot_with_pvalue, ROC_HCC$boxplot_with_pvalue, ncol = 3)
dev.off()
