# Re-run LOOCV for healthy, pre-cancer and cancer types - MM, MGUS and HD
library(MASS)
library(reshape2)
library(ggplot2)
library(devtools)
library(randomForest)

# Set the working directory
setwd("/Users/breesheyroskams-hieter/Desktop/cfRNA/manuscript/revised_paper_with_validation/figures/figure3/")

############# Import Data ####################
data <- read.csv("../../tables/RPM_without_dates_updated.csv",header=TRUE, sep=",", row.names = 1)
metadata <- read.csv("../../tables/PP_metadata_keep_FINAL_updated.csv",sep = ",", stringsAsFactors = FALSE)
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

# LVQ genes from lvq script
LVQ_MGUSvsMM <- readRDS(file="../../tables/LVQ/MGUS_MM_importance.rds")
LVQ_MGUSvsHD <- readRDS(file="../../tables/LVQ/NC_MGUS_importance.rds")
LVQ_MMvsHD <- readRDS(file="../../tables/LVQ/HD_MM_importance.rds")

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
## Remove PP084 from analysis (misclassified patient)
metadata <- metadata[metadata$PP_ID != "PP084",]
datac <- datac[,colnames(datac) != "PP084",]

########  HD MGUS MM and lvq top10 genes ########
# Subset for types of interest
sample_sel <- which(metadata$Status=='NC'|metadata$Status=='MGUS'|metadata$Status=='MM')
datac_sel <- datac[,sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[sample_sel,]
row.names(metadata_sel)=metadata_sel$PP_ID
metadata_sel=metadata_sel[,c(4,2)]

# Select custom gene set (LVQ gene set)
MM.vs.HD_top10lvq <- rownames(LVQ_MMvsHD$importance[order(LVQ_MMvsHD$importance[['HD']], decreasing = TRUE),])[1:10]
MG.vs.HD_top10lvq <- rownames(LVQ_MGUSvsHD$importance[order(LVQ_MGUSvsHD$importance[['NC']], decreasing = TRUE),])[1:10]

genes_sel = c(MM.vs.HD_top10lvq,MG.vs.HD_top10lvq)
select_datac= which( (row.names(datac)) %in% (genes_sel))

# Format training data
data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) 
data_train <- data_train[,-2]
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

Groups <- data_train$Status
my.lda <- data.frame(p$x[,1:2], Groups)

pred_out <- data.frame(sample = character(), pred = character(), actual_class = character())
for(i in 1:nrow( data_train)){
  training <- data_train[-i,]
  testing <- data_train[i,]
  model_lda <- lda(Status~.,training)
  p2 <- predict(model_lda,testing)$class
  pred_out <- rbind(pred_out,cbind(rownames(testing),data.frame(p2),testing$Status))
}

dir.create("LDA_LVQ/")
filename = paste0("LDA_LVQ/HD_MGUS_MM_lvq.csv")
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

dir.create("RF_LVQ/")
filename = paste("RF_LVQ/RF_LOOCV_HD_MM.MG.use_lvqtop10.csv")
write.csv(pred_out_rf, file=filename)

tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
tabrf
sum(diag(tabrf))/sum(tabrf)

########## NC, MGUS lvq genes #############
# LDA LOOCV
sample_sel <- which(metadata$Status=='NC'|metadata$Status=='MGUS')
datac_sel <- datac[,sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[sample_sel,]
row.names(metadata_sel)=metadata_sel$PP_ID
metadata_sel=metadata_sel[,c(4,2)]

# Select custom gene set
MG.vs.HD_top10lvq <- rownames(LVQ_MGUSvsHD$importance[order(LVQ_MGUSvsHD$importance[['NC']], decreasing = TRUE),])[1:10]
genes_sel = c(MG.vs.HD_top10lvq)
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

filename = paste("LDA_LVQ/LDA_LOOCV_HD.MG.use_lvqtop10.csv")
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

filename = paste("RF_LVQ/RF_LOOCV_HD.MG.use_lvqtop10.csv")
write.csv(pred_out_rf, file=filename)

tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
tabrf
sum(diag(tabrf))/sum(tabrf)

##########  MM.vs.MGUS_top10lvq genes ###################
# LDA LOOCV
sample_sel <- which(metadata$Status=='MM'|metadata$Status=='MGUS')
datac_sel <- datac[,sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[sample_sel,]
row.names(metadata_sel)=metadata_sel$PP_ID
metadata_sel=metadata_sel[,c(4,2)]

# Select custom gene set
MM.vs.MG_top10lvq <- rownames(LVQ_MGUSvsMM$importance[order(LVQ_MGUSvsMM$importance[['MGUS']], decreasing = TRUE),])[1:10]
genes_sel = MM.vs.MG_top10lvq
select_datac= which( (row.names(datac)) %in% (genes_sel))

# Format training data
data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) #transform data to log scale and attach with metadata
data_train <- data_train[,-2]
str(data_train)
data_train$Status <- as.factor(data_train$Status)
table(data_train$Status)
names(data_train) <- make.names(names(data_train))

# LOOCV
pred_out <- data.frame(sample = character(), pred = character(), actual_class = character())
for(i in 1:nrow( data_train)){
  training <- data_train[-i,]
  testing <- data_train[i,]
  model_lda <- lda(Status~.,training)
  p2 <- predict(model_lda,testing)$class
  pred_out <- rbind(pred_out,cbind(rownames(testing),data.frame(p2),testing$Status))
}

filename = paste("LDA_LVQ/LDA_LOOCV_MM.MG.use_lvqtop10.csv")
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

filename = paste("RF_LVQ/RF_LOOCV_MM.MG.use_lvqtop10.csv")
write.csv(pred_out_rf, file=filename)

tabrf <- table(Predict = pred_out_rf[,1], Actual = pred_out_rf[,2])
tabrf
sum(diag(tabrf))/sum(tabrf)
