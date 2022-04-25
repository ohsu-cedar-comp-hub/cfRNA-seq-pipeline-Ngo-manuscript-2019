# Generate LVQ results for all pairwise comparisons
# Load libraries
library(reshape2)
library(caret)
library(glmnet)
library(mlbench)
library(psych)

# Set working directory
setwd("/Users/breesheyroskams-hieter/Desktop/cfRNA/manuscript/revised_paper_with_validation/")

############################################
# Import Data
############################################

data <- read.csv("tables/RPM_without_dates_updated.csv",sep=",", row.names = 1)
metadata <- read.csv("tables/PP_metadata_keep_FINAL_updated.csv",sep = ",", stringsAsFactors = FALSE)

############################################
# Format data 
############################################

# Format gene names as characters
data$gene <- as.character(data$gene)

# Consolidate duplicate gene names
datam <- melt(data,id="gene")
datac <- dcast(datam,gene~variable,fun.aggregate = mean)
# Set rownames as gene symbols
rownames(datac) <- datac[,1]
datac <- datac[,-1]

## Remove PP084 from analysis
metadata <- metadata[metadata$PP_ID != "PP084",]
datac <- datac[,colnames(datac) != "PP084"]

############################################
# Feature Selection by lvq using caret 
############################################

dir.create("tables/LVQ/", recursive = TRUE)

RunLVQ <- function(Baseline, Target, Data, Metadata) {
  sample_sel <- Metadata[Metadata$Status==Baseline|Metadata$Status==Target,]$PP_ID
  datac_sel <- Data[,colnames(Data) %in% sample_sel]
  #make genename become legal name in R
  names(datac_sel) <- make.names(names(datac_sel))
  
  metadata_sel <- Metadata[Metadata$PP_ID %in% sample_sel,]
  rownames(metadata_sel)=metadata_sel$PP_ID
  metadata_sel=metadata_sel[,c(4,5)]
  
  #transform data to log scale
  genecount <- log(datac_sel+1,2)
  # Which values in genecount are greater than 0.5?
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
  
  # Create training data with metadata attached
  data_train <- cbind(metadata_sel,t(counts.keep)) 
  data_train <- data_train[,-2]
  str(data_train)
  data_train$Status <- as.factor(data_train$Status)
  table(data_train$Status)
  # Generate legal names in R
  names(data_train) <- make.names(names(data_train))
  
  #############************** Rank Features By Importance using lvq in caret *************########################
  # ensure results are repeatable
  set.seed(7)
  # prepare training scheme
  control <- trainControl(method="repeatedcv", number=10, repeats=3)
  # train the model
  model <- train(Status~., data=data_train, method="lvq", preProcess="scale", trControl=control,tuneGrid = data.frame(size = 3, k = 1:2))
  # estimate variable importance
  importance <- varImp(model, scale=FALSE)
  # summarize importance
  print(importance)
  print(model)
  
  # Export results as rds
  fname <- paste0("tables/LVQ/",Baseline,"_",Target)
  saveRDS(importance, file=paste0(fname, "_importance.rds"))
  saveRDS(model, file=paste0(fname, "_model.rds"))
}

RunLVQ(Baseline = "NC", Target = "HCC", Data = datac, Metadata = metadata)
RunLVQ(Baseline = "NC", Target = "MM", Data = datac, Metadata = metadata)
RunLVQ(Baseline = "NC", Target = "MGUS", Data = datac, Metadata = metadata)
RunLVQ(Baseline = "HCC", Target = "MM", Data = datac, Metadata = metadata)
RunLVQ(Baseline = "MGUS", Target = "MM", Data = datac, Metadata = metadata)
RunLVQ(Baseline = "Cirr", Target = "HCC", Data = datac, Metadata = metadata)
RunLVQ(Baseline = "Cirr", Target = "NC", Data = datac, Metadata = metadata)
