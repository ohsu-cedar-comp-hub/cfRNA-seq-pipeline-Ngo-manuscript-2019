#Load Libraries 
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

# Set the working directory
setwd("/Users/breesheyroskams-hieter/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/figures/figure3/")

############################################
# Import Data
############################################
data <- read.csv("../../tables/RPM_without_dates_updated.csv",sep=",", row.names = 1)
metadata <- read.csv("../../tables/PP_metadata_keep_FINAL_updated.csv",sep = ",", stringsAsFactors = FALSE)
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

# LVQ genes from lvq script
LVQ_CirrvsHCC <- readRDS(file="../../tables/LVQ/Cirr_HCC_importance.rds")
LVQ_CirrvsHD <- readRDS(file="../../tables/LVQ/Cirr_HD_importance.rds")
LVQ_HCCvsHD <- readRDS(file="../../tables/LVQ/HD_HCC_importance.rds")

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

###############################################
######## HD Cirr HCC and lvq top10 genes ######
# Select sample subset
sample_sel <- which(metadata$Status=='NC'|metadata$Status=='Cirr'|metadata$Status=='HCC')
datac_sel <- datac[,sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[sample_sel,]
row.names(metadata_sel)=metadata_sel$PP_ID
metadata_sel=metadata_sel[,c(4,2)]

# Select custom gene set
HCC.vs.HD_top10lvq <- rownames(LVQ_HCCvsHD$importance[order(LVQ_HCCvsHD$importance[['HD']], decreasing = TRUE),])[1:10]
Cirr.vs.HD_top10lvq <- rownames(LVQ_CirrvsHD$importance[order(LVQ_CirrvsHD$importance[['HD']], decreasing = TRUE),])[1:10]

genes_sel = c(HCC.vs.HD_top10lvq,Cirr.vs.HD_top10lvq)
select_datac= which( (row.names(datac)) %in% (genes_sel))

# Format training data
data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
data_train <- data_train[,-2]
str(data_train)
data_train$Status <- as.factor(data_train$Status)
table(data_train$Status)
names(data_train) <- make.names(names(data_train))

# Linear discriminant analysis
linear <- lda(Status~.,data_train)
linear
linear$prior
linear$counts
attributes(linear)
p <- predict(linear,data_train)
p

p1 <- p$class
tab1 <- table(Predicted = p1, Actual = data_train$Status)
tab1
sum(diag(tab1))/sum(tab1)

### LDA LOOCV ###
pred_out <- data.frame(sample = character(), pred = character(), actual_class = character())
for(i in 1:nrow( data_train)){
  training <- data_train[-i,]
  testing <- data_train[i,]
  model_lda <- lda(Status~.,training)
  p2 <- predict(model_lda,testing)$class
  pred_out <- rbind(pred_out,cbind(rownames(testing),data.frame(p2),testing$Status))
}

dir.create("LDA_LVQ/")
filename = paste0("LDA_LVQ/HD_Cirr_HCC_lvq.csv")
write.csv(pred_out, file=filename)

tabLDA <- table(Predict = pred_out[,2], Actual = pred_out[,3])
tabLDA
sum(diag(tabLDA))/sum(tabLDA)

# RF LOOCV model
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

filename = paste("RF_LVQ/RF_LOOCV_HD_HCC.Cirr.use_lvqtop10.csv")
write.csv(pred_out_rf, file=filename)

tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
tabrf
sum(diag(tabrf))/sum(tabrf)

########## HD, Cirr lvq genes ###################
# LDA LOOCV
sample_sel <- which(metadata$Status=='NC'|metadata$Status=='Cirr')
datac_sel <- datac[,sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[sample_sel,]
row.names(metadata_sel)=metadata_sel$PP_ID
metadata_sel=metadata_sel[,c(4,2)]

# Select custom gene set
Cirr.vs.HD_top10lvq <- rownames(LVQ_CirrvsHD$importance[order(LVQ_CirrvsHD$importance[['HD']], decreasing = TRUE),])[1:10]

genes_sel = c(Cirr.vs.HD_top10lvq)
select_datac= which( (row.names(datac)) %in% (genes_sel))

# Format training data
data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
data_train <- data_train[,-2]
str(data_train)
data_train$Status <- as.factor(data_train$Status)
table(data_train$Status)
names(data_train) <- make.names(names(data_train))

pred_out <- data.frame(sample = character(), pred = character(), actual_class = character())
for(i in 1:nrow( data_train)){
  training <- data_train[-i,]
  testing <- data_train[i,]
  model_lda <- lda(Status~.,training)
  p2 <- predict(model_lda,testing)$class
  pred_out <- rbind(pred_out,cbind(rownames(testing),data.frame(p2),testing$Status))
}

filename = paste("LDA_LVQ/LDA_LOOCV_HD.Cirr.use_lvqtop10.csv")
write.csv(pred_out, file=filename)

tabLDA <- table(Predict = pred_out[,2], Actual = pred_out[,3])
tabLDA
sum(diag(tabLDA))/sum(tabLDA)

# RF LOOCV
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

filename = paste("RF_LVQ/RF_LOOCV_HD.Cirr.use_lvqtop10.csv")
write.csv(pred_out_rf, file=filename)

tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
tabrf
sum(diag(tabrf))/sum(tabrf)

########## HCC vs Cirr top10 lvq genes ###################
# LDA LOOCV
sample_sel <- which(metadata$Status=='HCC'|metadata$Status=='Cirr')
datac_sel <- datac[,sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[sample_sel,]
row.names(metadata_sel)=metadata_sel$PP_ID
metadata_sel=metadata_sel[,c(4,2)]

# Select custom gene set
HCC.vs.Cirr_top10lvq <- rownames(LVQ_CirrvsHCC$importance[order(LVQ_CirrvsHCC$importance[['Cirr']], decreasing = TRUE),])[1:10]

genes_sel = c(HCC.vs.Cirr_top10lvq)
select_datac= which( (row.names(datac)) %in% (genes_sel))

data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
data_train <- data_train[,-2]
str(data_train)
data_train$Status <- as.factor(data_train$Status)
table(data_train$Status)
names(data_train) <- make.names(names(data_train))

pred_out <- data.frame(sample = character(), pred = character(), actual_class = character())
for(i in 1:nrow( data_train)){
  training <- data_train[-i,]
  testing <- data_train[i,]
  model_lda <- lda(Status~.,training)
  p2 <- predict(model_lda,testing)$class
  pred_out <- rbind(pred_out,cbind(rownames(testing),data.frame(p2),testing$Status))
}

filename = paste("LDA_LVQ/LDA_LOOCV_HCC.Cirr.use_lvqtop10.csv")
write.csv(pred_out, file=filename)

tabLDA <- table(Predict = pred_out[,2], Actual = pred_out[,3])
tabLDA
sum(diag(tabLDA))/sum(tabLDA)

# RF LOOCV
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

filename = paste("RF_LVQ/RF_LOOCV_HCC.Cirr.use_lvqtop10.csv")
write.csv(pred_out_rf, file=filename)

tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
tabrf
sum(diag(tabrf))/sum(tabrf)
