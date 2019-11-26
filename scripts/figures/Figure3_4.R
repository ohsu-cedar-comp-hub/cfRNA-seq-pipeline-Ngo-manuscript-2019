# Script to generate figure 3 and 4 for cfRNA manuscript
library(reshape2)
library(tsne)
library(ggplot2)
library(MASS)
library(cowplot)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggpubr)

# Set the working directory
setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/figures/figure3/")

################################################
# Import Data 
################################################
data <- read.csv("../../tables/RPM_without_dates_updated.csv",sep=",", row.names = 1)
metadata <- read.csv("../../tables/PP_metadata_keep_FINAL_updated.csv",sep = ",", stringsAsFactors = FALSE)
colors <- read.delim("../../tables/colors.txt", stringsAsFactors = FALSE)

LVQ_HCCvsHD <- readRDS(file="../../tables/LVQ/HD_HCC_importance.rds")
LVQ_CirrvsHD <- readRDS(file="../../tables/LVQ/Cirr_HD_importance.rds")
LVQ_MMvsHD <- readRDS(file="../../tables/LVQ/HD_MM_importance.rds")
LVQ_MGUSvsHD <- readRDS(file="../../tables/LVQ/HD_MGUS_importance.rds")

################################################
# Format Data 
################################################
data$gene <- as.character(data$gene)

# Consolidate duplicate gene names
datam <- melt(data,id="gene")
datac <- dcast(datam,gene~variable,fun.aggregate = mean)
# Assign gene as rownames 
rownames(datac) <- datac[,1]
# Remove first column, gene names
datac <- datac[,-1]

################################################
# Boxplot for MG, MM and HD
################################################

# Import Arial type font
library(extrafont)
extrafont::font_import()
extrafont::fonts()
extrafont::loadfonts()

######################################################################
# Boxplot of top 10 LVQ genes from HD vs. MM pairwise comparison
######################################################################

# Select data
metadata <- metadata[match(colnames(datac),metadata$PP_ID),]
sample_sel <- metadata[metadata$Status=="HD"|metadata$Status=="MGUS"|metadata$Status=="MM",]$PP_ID
metadata_sub=metadata[metadata$Status=="HD"|metadata$Status=="MGUS"|metadata$Status=="MM",]

# Select top 10 LVQ genes
MM.vs.HD_top10lvq <- rownames(LVQ_MMvsHD$importance[order(LVQ_MMvsHD$importance$HD, decreasing = TRUE),])[1:10]

genes <- MM.vs.HD_top10lvq
datac_sel <- datac[rownames(datac) %in% genes,colnames(datac) %in% sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

# Transform data to log scale
genecount <- log(datac_sel+1,2)

# Format data
genecount$gene <- rownames(genecount)
genecountm <- melt(genecount,id="gene")
colnames(genecountm) <- c("gene","PP_ID","counts")
metadata_group <- metadata[,c(1,4)]
genecountm_group <- merge(genecountm,metadata_group,by="PP_ID")
genecountm_group$Status <- factor(genecountm_group$Status, levels = c("HD","MGUS","MM"))

my_comparisons <- list(c("HD","MGUS"),c("MGUS","MM"),c("HD","MM"))
colourHD <- Colors[Colors$Status=="HD",]$Colour
colourMGUS <- Colors[Colors$Status=="MGUS",]$Colour
colourMM <- Colors[Colors$Status=="MM",]$Colour

boxplot_MM <- ggplot(genecountm_group, aes(x=Status, y=counts,color=Status)) + 
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~gene, scale="free",nrow=2)+
  geom_jitter(width = 0.1, size = 0.4)+
  theme_bw()+
  theme(text = element_text(family = "Arial"),
        plot.margin = unit(c(0,0.7,0,0),"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Counts (log2(RPM+1))") +
  scale_colour_manual(values=c(colourHD, colourMGUS, colourMM))+
  stat_compare_means(method="t.test",comparisons = my_comparisons, label = "p.signif", size = 2) +
  scale_y_continuous(expand = c(0.1,0))

######################################################################
# Boxplot of top 10 LVQ genes from HD vs. HCC pairwise comparison
######################################################################

metadata <- metadata[match(colnames(datac),metadata$PP_ID),]
sample_sel <- metadata[metadata$Status=="HD"|metadata$Status=="Cirr"|metadata$Status=="HCC",]$PP_ID
metadata_sub <- metadata[metadata$Status=="HD"|metadata$Status=="Cirr"|metadata$Status=="HCC",]

# Select top 10 LVQ genes
HCC.vs.HD_top10lvq <- rownames(LVQ_HCCvsHD$importance[order(LVQ_HCCvsHD$importance$HD, decreasing = TRUE),])[1:10]

genes <- HCC.vs.HD_top10lvq
datac_sel <- datac[rownames(datac) %in% genes,colnames(datac) %in% sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

# Transform data to log scale
genecount <- log(datac_sel+1,2)

# Format data
genecount$gene <- rownames(genecount)
genecountm <- melt(genecount,id="gene")
colnames(genecountm) <- c("gene","PP_ID","counts")
metadata_group <- metadata[,c(1,4)]
genecountm_group <- merge(genecountm,metadata_group,by="PP_ID")
genecountm_group$Status <- factor(genecountm_group$Status, levels = c("HD","Cirr","HCC"))

my_comparisons <- list(c("HD","Cirr"),c("Cirr","HCC"),c("HD","HCC"))
colourCirr <- Colors[Colors$Status=="Cirr",]$Colour
colourHCC <- Colors[Colors$Status=="HCC",]$Colour

boxplot_HCC <- ggplot(genecountm_group, aes(x=Status, y=counts,color=Status)) + 
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~gene, nrow=2,scale="free")+
  geom_jitter(width = 0.1, size = 0.4)+ 
  theme_bw()+
  theme(text = element_text(family = "Arial"),
        plot.margin = unit(c(0,0.7,0,0),"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Counts (log2(RPM+1))") +
  scale_colour_manual(values=c(colourHD, colourCirr, colourHCC))+
  stat_compare_means(method="t.test",comparisons = my_comparisons, label = "p.signif", size = 2) +
  scale_y_continuous(expand = c(0.1,0))

###################################
# LDA Plots
###################################

### HD, MM, MGUS ###

sample_sel <- which(metadata$Status=='HD'|metadata$Status=='MGUS'|metadata$Status=='MM')
datac_sel <- datac[,sample_sel]
#make genename become legal name in R
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[sample_sel,]
row.names(metadata_sel) <- metadata_sel$PP_ID
metadata_sel <- metadata_sel[,c(4,2)]

# Select MG vs. HD top 10 LVQ genes
MG.vs.HD_top10lvq <- rownames(LVQ_MGUSvsHD$importance[order(LVQ_MGUSvsHD$importance[['HD']], decreasing = TRUE),])[1:10]

# Select combined gene set of top 10 most important gene from pairwise comparisons HD vs. MM and HD vs. MGUS
genes_sel <- c(MM.vs.HD_top10lvq,MG.vs.HD_top10lvq)
select_datac <- which( (row.names(datac)) %in% (genes_sel))

# Transform data to log scale and attach with metadata
data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) 
data_train <- data_train[,-2]
data_train$Status <- as.factor(data_train$Status)
table(data_train$Status)
names(data_train) <- make.names(names(data_train))

# Linear discriminant analysis
linear <- lda(Status~.,data_train)
#histogram
p <- predict(linear,data_train)
p

p1 <- predict(linear,data_train)$class
tab1 <- table(Predicted = p1, Actual = data_train$Status)
tab1
sum(diag(tab1))/sum(tab1)

Groups <- data_train$Status
my.lda <- data.frame(p$x[,1:2], Groups)
my.lda$Groups <- factor(my.lda$Groups, levels = c("HD","MGUS","MM"))

# Plot LDA 
LDA_MM <- ggplot(my.lda,aes(x=LD1,y=LD2,colour= Groups),size=5)+
  geom_point(size=3)+
  xlab("LD1")+
  ylab("LD2")+
  scale_color_manual(breaks = c("HD", "MGUS", "MM"),
                     values=c(colourHD, colourMGUS, colourMM)) +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        plot.margin = unit(c(0,0,0.7,0),"cm"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 5),
        legend.text = element_text(size = 8),
        legend.position="bottom",
        legend.direction="horizontal",
        aspect.ratio=1,
        legend.title = element_blank()) 

legend <- cowplot::get_legend(LDA_MM)

pdf("../figure3/Legend_LDA_MGUS_MM.pdf", 5, 2, useDingbats = FALSE)
grid.newpage()
grid.draw(legend)
dev.off()

### HD, Cirr, HCC ###

sample_sel <- which(metadata$Status=='HD'|metadata$Status=='Cirr'|metadata$Status=='HCC')
datac_sel <- datac[,sample_sel]
names(datac_sel) <- make.names(names(datac_sel))

metadata_sel <- metadata[sample_sel,]
row.names(metadata_sel) <- metadata_sel$PP_ID
metadata_sel <- metadata_sel[,c(4,2)]

# Select MG vs. HD top 10 LVQ genes
Cirr.vs.HD_top10lvq <- rownames(LVQ_CirrvsHD$importance[order(LVQ_CirrvsHD$importance[['HD']], decreasing = TRUE),])[1:10]

# Select combined gene set of top 10 most important gene from pairwise comparisons HD vs. HCC and HD vs. Cirr
genes_sel <- c(HCC.vs.HD_top10lvq,Cirr.vs.HD_top10lvq)
select_datac <- which((row.names(datac)) %in% (genes_sel))

# Transform data to log scale and attach with metadata
data_train <- cbind(metadata_sel,t(log(datac_sel[select_datac,]+1,2))) 
data_train <- data_train[,-2]
data_train$Status <- as.factor(data_train$Status)
table(data_train$Status)
names(data_train) <- make.names(names(data_train))

# Linear discriminant analysis
linear <- lda(Status~.,data_train)
#histogram
p <- predict(linear,data_train)
p

p1 <- predict(linear,data_train)$class
tab1 <- table(Predicted = p1, Actual = data_train$Status)
tab1
sum(diag(tab1))/sum(tab1)

Groups <- data_train$Status
my.lda <- data.frame(p$x[,1:2], Groups)
my.lda$Groups <- factor(my.lda$Groups, levels = c("HD","Cirr","HCC"))

#Plot LDA
LDA_HCC <- ggplot(my.lda,aes(x=LD1,y=LD2,colour= Groups),size=5)+
  geom_point(size=3)+
  xlab("LD1")+
  ylab("LD2")+
  scale_color_manual(breaks = c("HD", "Cirr", "HCC"),
                     values=c(colourHD, colourCirr, colourHCC)) +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        plot.margin = unit(c(0,0,0,0.7),"cm"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 5),
        legend.text = element_text(size = 8),
        legend.position="bottom",
        legend.direction="horizontal",
        aspect.ratio=1,
        legend.title = element_blank()) 

legend <- cowplot::get_legend(LDA_HCC)

pdf("../figure4/Legend_LDA_Cirr_HCC.pdf", 5, 2, useDingbats = FALSE)
grid.newpage()
grid.draw(legend)
dev.off()

# Plot as grid using cowplot
fig3 <- cowplot::plot_grid(boxplot_MM, LDA_MM, nrow = 1, rel_widths = c(1.3, 1))
fig4 <- cowplot::plot_grid(boxplot_HCC, LDA_HCC, nrow = 1, rel_widths = c(1.3, 1))

# Save plots
cowplot::save_plot("../figure3/Figure3_A_B.pdf", fig3, base_width = 8)
cowplot::save_plot("../figure4/Figure4_A_B.pdf", fig4, base_width = 8)
