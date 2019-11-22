# This is a script to re-generate nicer figures for the GO analysis
library(GO.db)
library(topGO)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(GenomicFeatures)
library(extrafont)
extrafont::font_import()
extrafont::fonts()
extrafont::loadfonts()

setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/figures/supplementary/GO/")

##-----------------------------------Functions--------------------------------------#
runGO <- function(geneList,xx=xx,otype,setName){
  tableDir <- "../../../tables/GO_tables/"
  
  # load annotation
  geneID2GO <- get(load("/Users/roskamsh/Desktop/CEDAR/anno/hg38.Ens_90.biomaRt.GO.external.geneID2GO.RData"))
  xx <- get(load("/Users/roskamsh/Desktop/cfRNA/PP_cohort_rough_work/GO_analysis/GO.db.Term.list.rda"))
  
  setLength       <- sum(as.numeric(levels(geneList))[geneList]) 
  fname           <- paste(sub("$", paste(setName, setLength, otype, "GO.txt", sep="_"), tableDir))
  GOData          <- new("topGOdata", ontology=otype, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO) ##run topGO
  resultFisher    <- runTest(GOData, algorithm = "classic", statistic = "fisher")## statistical test for topGO
  x               <- GenTable(GOData, classicFisher=resultFisher, topNodes=length(names(resultFisher@score)))## make go table for all terms
  x               <- data.frame(x)
  pVal            <- data.frame(pval=signif(resultFisher@score, 6)) ## get unrounded pvalue
  x$enrich        <- x$Significant/x$Expected ## calculate enrichment based on what you expect by chance
  x$p.unround     <- pVal[x$GO.ID,"pval"]## put unrounded pvalue in the table
  x$p.adj         <- signif(p.adjust(x$p.unround, method="BH"), 6)## calculate the adjusted pvalue with Benjamini & Hochberg correction
  x$log.p.adj     <- -log10(x$p.adj) ## convert adjusted p value to -log10 for plot magnitude
  x <- x[order(x$GO.ID),]
  write.table(x, file=fname, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE) ## save the table
  return(x)  
}

## function to make barplot of -log10 adjusted pvalues colored by enrichment

drawBarplot <- function(go, xx, ontology, setName, setSize, Baseline, Target){
  go <- go[!go$p.adj > 0.01,]
  if(nrow(go)>1){
    go$Term <- make.unique(paste(sapply(strsplit(as.character(substring(go$Term,1,50)), "\\,"), `[`, 1)))
    print(setName)
    print(setSize)
    title <- strsplit(setName, split='\\.')[[1]]
    go <- go[order(go$enrich, decreasing = TRUE),]
    go$Term <-factor(paste(go$Term), levels=rev(paste(go$Term))) 
    ptitle <- paste(Baseline, "vs", Target, title[[4]], title[[5]], title[[6]]) ## plot title
    pfname <- paste(setName,setSize,"pdf",sep=".")## name of png file
    if(nrow(go) < 20 ){
      toprange <- 1:nrow(go)
    }else{
      toprange <- 1:20
    }
    top <- go[toprange,]
    top$Term.full     <- sapply(top$GO.ID, FUN=function(n){Term(xx[[n]])})
    top <- top[order(top$enrich),]
    order <- as.vector(top$Term.full)
    top$Term.full <- factor(top$Term.full, levels=order)
    col <- colorRampPalette(c("white","navy"))(16)
    
    mar <- ifelse (Target == "LuCa", 1.8, 
                  ifelse (Target == "MM", 4.4, 0))
    
    p <- ggplot(top, aes(y=enrich, x=Term.full, fill=log.p.adj)) + ## ggplot barplot function
      geom_bar(stat="identity",colour="black") +
      ggtitle(ptitle) +
      xlab("") + ylab("enrichment") +
      scale_fill_gradient(low=col[2], high=col[15], name="-log10(FDR)", limits=c(0,ceiling(max(top$log.p.adj))))+
      coord_flip()+
      theme(text = element_text(family = "Arial"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(),
            plot.margin = unit(c(0,0,0,mar),"cm"),
            axis.text.x = element_text(vjust=1,color="black",size=5),
            axis.text.y = element_text(color="black",size=8),
            axis.title.x = element_text(size = 8),
            plot.title=element_blank(),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 8))
    
    return(p)
  }
}

GenerateGObarPlotUpGenes <- function(degFile, Baseline, Target, FC = 2, adjp = 0.01) {
  print("Loading differential expressed gene table")
  print(degFile)
  deg <- read.delim(file=degFile,header=TRUE,sep="\t", stringsAsFactors = FALSE)

  # load annotation
  geneID2GO <- get(load("/Users/roskamsh/Desktop/CEDAR/anno/hg38.Ens_90.biomaRt.GO.external.geneID2GO.RData"))
  xx <- get(load("/Users/roskamsh/Desktop/cfRNA/PP_cohort_rough_work/GO_analysis/GO.db.Term.list.rda"))
  
  print("get up genes and make geneList")
  up <- deg$padj < adjp & deg$log2FoldChange >= log2(FC)
  up <- unique(deg[up,]$GeneID)
  all <-unique(names(geneID2GO))
  up.geneList <-  factor(as.integer(all %in% up))
  names(up.geneList) <- all

  up.setsize <- sum(as.numeric(levels(up.geneList))[up.geneList])
  print("setsize for significant genes") 
  up.setsize
  
  adjplabel <- gsub("^0\\.","",adjp)
  
  print("make GO table for the up genes")
  #################################
  go.UP.BP <- runGO(geneList=up.geneList,xx=xx,otype="BP",setName=paste(Baseline, "vs", Target,"upFC",FC, "adjp", adjp, sep="."))
  
  print("make the png for the up genes")
  bp <- drawBarplot(go=go.UP.BP,xx=xx,ontology="BP", Baseline = Baseline, Target = Target,
                    setName=paste(Baseline, "vs", Target,"Up", "Regulated", "Genes",FC, "adjp", adjp, sep="."),setSize=up.setsize)
  
  return(bp)
}

LuCa_GO <- GenerateGObarPlotUpGenes(degFile = "../../../tables/DE/LuCa-vs-HD.diffexp.geneID.tsv", Baseline = "HD", Target = "LuCa")
MM_GO <- GenerateGObarPlotUpGenes(degFile = "../../../tables/DE/MM-vs-HD.diffexp.geneID.tsv", Baseline = "HD", Target = "MM")
HCC_GO <- GenerateGObarPlotUpGenes(degFile = "../../../tables/DE/HCC-vs-HD.diffexp.geneID.tsv", Baseline = "HD", Target = "HCC")

pdf("GO_grid.pdf", 8, 9, useDingbats = FALSE)
cowplot::plot_grid(LuCa_GO, MM_GO, HCC_GO, nrow = 3, rel_heights = c(1, 0.6, 1))
dev.off()
