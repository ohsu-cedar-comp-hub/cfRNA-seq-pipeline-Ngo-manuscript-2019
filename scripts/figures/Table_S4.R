# Generate Table S4

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

# Export top 10 LVQ genes for supplementary table S4
LuCa.vs.HD_top10lvq <- rownames(LVQ_LuCavsHD$importance[order(LVQ_LuCavsHD$importance$HD, decreasing = TRUE),])[1:10]
HCC.vs.HD_top10lvq <- rownames(LVQ_HCCvsHD$importance[order(LVQ_HCCvsHD$importance$HD, decreasing = TRUE),])[1:10]
MM.vs.HD_top10lvq <- rownames(LVQ_MMvsHD$importance[order(LVQ_MMvsHD$importance$HD, decreasing = TRUE),])[1:10]
MGUS.vs.HD_top10lvq <- rownames(LVQ_MGUSvsHD$importance[order(LVQ_MGUSvsHD$importance$HD, decreasing = TRUE),])[1:10]
Cirr.vs.HD_top10lvq <- rownames(LVQ_CirrvsHD$importance[order(LVQ_CirrvsHD$importance$HD, decreasing = TRUE),])[1:10]
LuCa.vs.MM_top10lvq <- rownames(LVQ_LuCavsMM$importance[order(LVQ_LuCavsMM$importance$MM, decreasing = TRUE),])[1:10]
HCC.vs.MM_top10lvq <- rownames(LVQ_HCCvsMM$importance[order(LVQ_HCCvsMM$importance$MM, decreasing = TRUE),])[1:10]
HCC.vs.LuCa_top10lvq <- rownames(LVQ_HCCvsLuCa$importance[order(LVQ_HCCvsLuCa$importance$LuCa, decreasing = TRUE),])[1:10]

write.table(as.data.frame(LuCa.vs.HD_top10lvq), file="../../tables/LVQ/LuCa_HD_top10_lvq.txt", 
            col.names = F, row.names = F, quote = F)
write.table(as.data.frame(HCC.vs.HD_top10lvq), file="../../tables/LVQ/HCC_HD_top10_lvq.txt", 
            col.names = F, row.names = F, quote = F)
write.table(as.data.frame(MM.vs.HD_top10lvq), file="../../tables/LVQ/MM_HD_top10_lvq.txt", 
            col.names = F, row.names = F, quote = F)
write.table(as.data.frame(LuCa.vs.MM_top10lvq), file="../../tables/LVQ/LuCa_MM_top10_lvq.txt", 
            col.names = F, row.names = F, quote = F)
write.table(as.data.frame(HCC.vs.MM_top10lvq), file="../../tables/LVQ/HCC_MM_top10_lvq.txt", 
            col.names = F, row.names = F, quote = F)
write.table(as.data.frame(Cirr.vs.HD_top10lvq), file="../../tables/LVQ/Cirr_HD_top10_lvq.txt", 
            col.names = F, row.names = F, quote = F)
write.table(as.data.frame(MGUS.vs.HD_top10lvq), file="../../tables/LVQ/MGUS_HD_top10_lvq.txt", 
            col.names = F, row.names = F, quote = F)
