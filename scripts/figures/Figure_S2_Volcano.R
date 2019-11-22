# This is a script for generating Figure-ready Volcano plots

setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/figures/supplementary/volcanoPlots")

# Import Arial type font
library(extrafont)
extrafont::font_import()
extrafont::fonts()
extrafont::loadfonts()

MakeVolcanoPlot <- function(degFile, Baseline, Target, adjp = 0.01, FC = 2) {
  downCol <- "#5159E8"
  upCol <- "#FF7E6E"
  ncCol <- "#CCCCCC"
  
  ## check if an rda file or tab sep
  deg <- read.delim(file=degFile)
  
  head(deg)
  dim(deg)
  
  ## select up regulated genes
  up <- deg$padj < adjp & deg$log2FoldChange > log2(FC)
  sum(up)
  
  ## select down regulated genes
  down <- deg$padj < adjp & deg$log2FoldChange < -log2(FC)
  sum(down)
  
  ## set labels for pdf
  comparison <- paste(Baseline, "vs", Target, sep=".")
  pdfFile <- paste(comparison,"VolcanoPlot.pdf",sep=".")
  print(pdfFile)
  
  ## calculate the -log10(adjp) for the plot
  deg$log10padj <- -log10(deg$padj)
  
  deg$Expression <- ifelse(down, 'down',
                           ifelse(up, 'up','NS'))
  deg$Expression <- factor(deg$Expression, levels=c("up","down","NS"))
  
  if (sum(up)==0 & sum(down)==0) {
    colours <- ncCol
  } else if (sum(up)==0) {
    colours <- c(downCol, ncCol)
  } else if (sum(down)==0) {
    colours <- c(upCol, ncCol)
  } else {
    colours <- c(upCol, downCol, ncCol)
  }
  
  p <- ggplot(data=deg, mapping=aes(x=log2FoldChange, y=log10padj, colour=Expression)) +
    geom_point(size = 0.8, alpha = 0.8) +
    geom_vline(xintercept = c(-log2(FC),log2(FC)), linetype="dashed", colour="gray45", size = 0.3) +
    geom_hline(yintercept = -log10(adjp), linetype="dashed", colour="gray45", size = 0.3) +
    ylab("-log10(FDR)") +
    xlab("log2(Fold Change)") +
    ggtitle(paste(Baseline, "vs", Target)) +
    scale_colour_manual(values=colours) +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5, face="plain", size = 10),
          aspect.ratio = 1,
          panel.background = element_blank(),
          axis.line = element_line(colour = "gray45"),
          legend.position = "none",
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 5))
  
  return(p)
}


LuCa_volcano <- MakeVolcanoPlot(degFile = "../../../tables/DE/LuCa-vs-HD.diffexp.geneID.tsv", Baseline = "HD", Target = "LuCa")
MM_volcano <- MakeVolcanoPlot(degFile = "../../../tables/DE/MM-vs-HD.diffexp.geneID.tsv", Baseline = "HD", Target = "MM")
HCC_volcano <- MakeVolcanoPlot(degFile = "../../../tables/DE/HCC-vs-HD.diffexp.geneID.tsv", Baseline = "HD", Target = "HCC")

pdf("volcanoPlots_A_C.pdf", 8.3, 2.6, useDingbats = FALSE)
cowplot::plot_grid(LuCa_volcano, MM_volcano, HCC_volcano, nrow = 1)
dev.off()
