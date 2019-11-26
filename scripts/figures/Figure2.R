# Figure 2 for cfRNA manuscript 2019 and Table S4
# Load libraries
library(reshape2)
library(tsne)
library(ggplot2)
library(MASS)
library(devtools)
library(ggord)
library(plotly)
library(plyr)
library(scatterplot3d)
library(randomForest)
library(caret)
library(irr)
library(ggpubr)
library(cowplot)
library(plotly)

# Set working directory
setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/figures/figure2/")

############################################
# Import Data
############################################

data <- read.csv("../../tables/RPM_without_dates_updated.csv",sep=",", row.names = 1)
metadata <- read.csv("../../tables/PP_metadata_keep_FINAL_updated.csv",sep = ",", stringsAsFactors = FALSE)
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

DE_LuCavsHD <- read.table(file="../../tables/DE/LuCa-vs-HD.diffexp.geneID.tsv",header = TRUE)
DE_HCCvsHD <- read.table(file="../../tables/DE/HCC-vs-HD.diffexp.geneID.tsv",header = TRUE)
DE_CirrvsHD <- read.table(file="../../tables/DE/Cirr-vs-HD.diffexp.geneID.tsv",header = TRUE)
DE_MMvsHD <- read.table(file="../../tables/DE/MM-vs-HD.diffexp.geneID.tsv",header = TRUE)
DE_MGUSvsHD <- read.table(file="../../tables/DE/MGUS-vs-HD.diffexp.geneID.tsv",header = TRUE)

DE_LuCavsHCC <- read.table(file="../../tables/DE/LuCa-vs-HCC.diffexp.geneID.tsv",header = TRUE)
DE_LuCavsMM <- read.table(file="../../tables/DE/LuCa-vs-MM.diffexp.geneID.tsv",header = TRUE)
DE_HCCvsMM <- read.table(file="../../tables/DE/MM-vs-HCC.diffexp.geneID.tsv",header = TRUE)

# LVQ
LVQ_LuCavsHD <- readRDS(file="../../tables/LVQ/HD_LuCa_importance.rds")
LVQ_HCCvsHD <- readRDS(file="../../tables/LVQ/HD_HCC_importance.rds")
LVQ_CirrvsHD <- readRDS(file="../../tables/LVQ/Cirr_HD_importance.rds")
LVQ_MMvsHD <- readRDS(file="../../tables/LVQ/HD_MM_importance.rds")
LVQ_MGUSvsHD <- readRDS(file="../../tables/LVQ/HD_MGUS_importance.rds")

LVQ_LuCavsMM <- readRDS(file="../../tables/LVQ/LuCa_MM_importance.rds")
LVQ_HCCvsMM <- readRDS(file="../../tables/LVQ/HCC_MM_importance.rds")
LVQ_MGUSvsMM <- readRDS(file="../../tables/LVQ/MGUS_MM_importance.rds")
LVQ_CirrvsHCC <- readRDS(file="../../tables/LVQ/Cirr_HCC_importance.rds")
LVQ_HCCvsLuCa <- readRDS(file="../../tables/LVQ/HCC_LuCa_importance.rds")

############################################
# Format Data
############################################

data$gene <- as.character(data$gene)

# Consolidate duplicate gene names
datam <- melt(data,id="gene")
datac <- dcast(datam,gene~variable,fun.aggregate = mean)
rownames(datac) <- datac[,1]
datac <- datac[,-1]

# Define colours
colourHD <- colors[colors$Status=="HD",]$Colour
colourHCC <- colors[colors$Status=="HCC",]$Colour
colourLuCa <- colors[colors$Status=="LuCa",]$Colour
colourMM <- colors[colors$Status=="MM",]$Colour

# Select samples for Healthy and 3 cancers
sample_sel <- metadata[metadata$Status=="HD"|metadata$Status=="LuCa"|metadata$Status=="HCC"|metadata$Status=="MM",]$PP_ID
datac_sel <- datac[,colnames(datac) %in% sample_sel]
# make genename legal in R
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[metadata$PP_ID %in% sample_sel,]
row.names(metadata_sel)=metadata_sel$PP_ID
metadata_sel=metadata_sel[,c(4,2)]

# Select custom gene set
LuCa.vs.HD_top5lvq <- rownames(LVQ_LuCavsHD$importance[order(LVQ_LuCavsHD$importance$HD, decreasing = TRUE),])[1:5]
HCC.vs.HD_top5lvq <- rownames(LVQ_HCCvsHD$importance[order(LVQ_HCCvsHD$importance$HD, decreasing = TRUE),])[1:5]
MM.vs.HD_top5lvq <- rownames(LVQ_MMvsHD$importance[order(LVQ_MMvsHD$importance$HD, decreasing = TRUE),])[1:5]
MM.vs.HCC_top5lvq <- rownames(LVQ_HCCvsMM$importance[order(LVQ_HCCvsMM$importance$MM, decreasing = TRUE),])[1:5]
MM.vs.LuCa_top5lvq <- rownames(LVQ_LuCavsMM$importance[order(LVQ_LuCavsMM$importance$MM, decreasing = TRUE),])[1:5]
LuCa.vs.HCC_top5lvq <- rownames(LVQ_HCCvsLuCa$importance[order(LVQ_HCCvsLuCa$importance$HCC, decreasing = TRUE),])[1:5]

# Import Arial type font
library(extrafont)
extrafont::font_import()
extrafont::fonts()
extrafont::loadfonts()

GenerateBoxPlotOfTop5LVQGenes <- function(Data, Metadata, LVQlist, Comparison, Colors) {
  # Generate boxplots for top 5 lvq genes for across all types
  select_datac= which((row.names(Data)) %in% (LVQlist))
  data_sub <- cbind(Metadata,t(log(Data[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_sub <- data_sub[,-2]
  data_plot <- melt(data_sub, id="Status")
  
  # Colors
  colourHD <- Colors[Colors$Status=="HD",]$Colour
  colourHCC <- Colors[Colors$Status=="HCC",]$Colour
  colourLuCa <- Colors[Colors$Status=="LuCa",]$Colour
  colourMM <- Colors[Colors$Status=="MM",]$Colour
  
  my_comparisons <- list(c("HD","LuCa"),c("HD","MM"),c("HD","HCC"))
  data_plot$Status <- factor(data_plot$Status, levels = c("HD","LuCa","MM","HCC"))
  
  g <- ggplot(data = data_plot, mapping = aes(x = Status, y = value, colour = Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size=0.4) +
    facet_wrap(~variable, nrow=1, scale="free") +
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 5),
          legend.position = "none") +
    ylab("Counts (log2(RPM+1))") +
    scale_color_manual(breaks=c("HD","LuCa","MM","HCC"),
                       values=c(colourHD, colourLuCa, colourMM, colourHCC)) +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))),
                       expand = c(0.1,0)) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 2)

  return(g)
}

LuCa_boxplot <- GenerateBoxPlotOfTop5LVQGenes(Data = datac_sel, Metadata = metadata_sel, LVQlist = LuCa.vs.HD_top5lvq, 
                                              Comparison = "HD-vs-LuCa", Colors = colors)
MM_boxplot <- GenerateBoxPlotOfTop5LVQGenes(Data = datac_sel, Metadata = metadata_sel, LVQlist = MM.vs.HD_top5lvq, 
                                            Comparison = "HD-vs-MM", Colors = colors)
HCC_boxplot <- GenerateBoxPlotOfTop5LVQGenes(Data = datac_sel, Metadata = metadata_sel, LVQlist = HCC.vs.HD_top5lvq, 
                                             Comparison = "HD-vs-HCC", Colors = colors)

pdf("Figure2_A_C.pdf", 4, 5, useDingbats = FALSE)
cowplot::plot_grid(LuCa_boxplot, MM_boxplot, HCC_boxplot, nrow = 3)
dev.off()

GenerateLegend <- function(Data, Metadata, LVQlist, Colors) {
  # Generate legend for top 5 lvq genes for across all types
  select_datac= which((row.names(Data)) %in% (LVQlist))
  data_sub <- cbind(Metadata,t(log(Data[select_datac,]+1,2))) #transform data to log scale and attach with metadata
  data_sub <- data_sub[,-2]
  data_plot <- melt(data_sub, id="Status")
  
  # Colors
  colourHD <- Colors[Colors$Status=="HD",]$Colour
  colourHCC <- Colors[Colors$Status=="HCC",]$Colour
  colourLuCa <- Colors[Colors$Status=="LuCa",]$Colour
  colourMM <- Colors[Colors$Status=="MM",]$Colour
  
  my_comparisons <- list(c("HD","HCC"),c("HD","LuCa"),c("HD","MM"))
  data_plot$Status <- factor(data_plot$Status, levels = c("HD","LuCa","MM","HCC"))

  g <- ggplot(data = data_plot, mapping = aes(x = Status, y = value, colour = Status)) +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    facet_wrap(~variable, nrow=1, scale="free") +
    scale_color_manual(values=c(colourHD, colourLuCa, colourMM, colourHCC)) +
    theme(text = element_text(family = "Arial"),
          legend.text = element_text(size=6),
          legend.title = element_text(size=8))
  
  legend <- cowplot::get_legend(g)
  
  pdf(paste0("Legend-Boxplot.pdf"), 5, 2, useDingbats = FALSE)
  grid.newpage()
  grid.draw(legend)
  dev.off()
}

GenerateLegend(Data = datac_sel, Metadata = metadata_sel, LVQlist = LuCa.vs.HD_top5lvq, Colors = colors)

####################################
# Multiclass discrimination
####################################
# Select top 5 most important genes from all pairwise comparisons through LVQ analysis
genes_sel = c(LuCa.vs.HD_top5lvq,HCC.vs.HD_top5lvq,MM.vs.HD_top5lvq,MM.vs.HCC_top5lvq,MM.vs.LuCa_top5lvq,LuCa.vs.HCC_top5lvq)

# Index the rows in our data which include these genes
select_datac= which((row.names(datac)) %in% (genes_sel))

# Transform data to log scale and attach with metadata
data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) 
data_train <- data_train[,-2]
data_train$Status <- factor(data_train$Status, levels = c("HD","LuCa","MM","HCC"))
table(data_train$Status)

names(data_train) <- make.names(names(data_train))

# Linear discriminant analysis
linear <- lda(Status~.,data_train)

p <- predict(linear,data_train)
p

p1 <- predict(linear,data_train)$class
tab1 <- table(Predicted = p1, Actual = data_train$Status)
tab1
sum(diag(tab1))/sum(tab1)

Groups <- data_train$Status
my.lda <- data.frame(p$x[,1:3], Groups)

# Trained model on our training data
# Visualization of first three components of LDA
pl <- plot_ly(my.lda, x = ~LD1, y = ~LD2, z = ~LD3, color = ~Groups,colors = c(colourHD, colourLuCa, colourMM, colourHCC)) 
pl

dir.create("../../tables/LDA_multiclass/")
filename = paste("../../tables/LDA_multiclass/LDA_multiclass.csv")
write.csv(my.lda, file=filename)

# Write out my.lda with colours attached to feed into our 3d plot python script "Figure2_3D_plot.py"
my.lda$Colour <- plyr::mapvalues(my.lda$Groups, c("HD","LuCa","MM","HCC"), c(colourHD, colourLuCa, colourMM, colourHCC))
rownames(my.lda) <- 0:(nrow(my.lda)-1)
my.lda <- my.lda[order(my.lda$Groups),]
write.table(my.lda, file="../../tables/LDA_multiclass/LDA_table_for_3D_plot.txt", quote = F, sep="\t")

#####################################################################
# LDA LOOCV HD LuCa HCC MM using top 5 lvq genes pairwise  
#####################################################################

# Prediction based on model
pred_out <- data.frame(sample = character(), pred = character(), actual_class = character())
for(i in 1:nrow(data_train)){
  training <- data_train[-i,]
  testing <- data_train[i,]
  model_lda <- lda(Status~.,training)
  p2 <- predict(model_lda,testing)$class
  pred_out <- rbind(pred_out,cbind(rownames(testing),data.frame(p2),testing$Status))
}

tab1 <- table(Predict = pred_out[,2], Actual = pred_out[,3])
tab1
sum(diag(tab1))/sum(tab1)

#####################################################################
# Random Forest LOOCV HD LuCa HCC MM using top 5 lvq genes pairwise  
#####################################################################
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

dir.create("../../tables/RF_multiclass/")
filename = paste("../../tables/RF_multiclass/RF_HD.LG.LCx.MM_lvqtop5.csv")
write.csv(pred_out_rf, file=filename)

tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
tabrf
sum(diag(tabrf))/sum(tabrf)

########################################################################
# LOOCV Cohen's Kappa and Confusion Matrix HD HCC LuCa MM
########################################################################

rownames(pred_out) <- rownames(pred_out_rf)
LOOCV <- data.frame(pred_out, pred_out_rf)
LOOCV$rownames.testing. <- NULL
LOOCV$testing.Status <- NULL
names(LOOCV) <- c("LDA", "RF", "Actual")

library(irr)
# Cohen Kappa RF vs Actual
RF_Act <- kappa2(LOOCV[,c(2,3)], "squared")$value #Kappa = 0.827

## Temporary code for generating barplots
df <- as.data.frame(cbind("LOOCV RF \nvs. Actual", RF_Act))

names(df) <- c("Comparison","Cohen's Kappa")
df$Comparison <- as.character(df$Comparison)
df$`Cohen's Kappa` <- as.numeric(paste(df$`Cohen's Kappa`))

p <- ggplot(data = df, mapping = aes(x = Comparison, y = `Cohen's Kappa`)) +
  geom_bar(stat = "identity", width = 0.8) +
  ylim(0,1) +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 8),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

pdf("Cohens_Kappa_Barplot.pdf", 1.8, 2.3, useDingbats = FALSE)
print(p)
dev.off()
