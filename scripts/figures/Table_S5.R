# Generate supplemental tables for the permutation test (Table S5)
library(dplyr)
setwd("/Users/roskamsh/Desktop/cfRNA/manuscript/Re-run_updated_without_PP033/tables/permutation/supp_tables/")

GeneratePermutationTables <- function(Results, Assignment, Metadata, Baseline, Target) {
  numdiff <- read.csv(Results, row.names=1, header=TRUE)
  names(numdiff) <- c("NumDiffGenes", "Actual")
  assignment <- read.csv(Assignment, row.names=1, header=TRUE)
  
  # Generate subdata 
  md <- read.csv(file=Metadata, stringsAsFactors = FALSE)

  # Extract only the Types that we want in further analysis & only the PP_ID and Status informative columns
  md <- filter(md, Status == Baseline | Status == Target)
  
  # Keep only the PP_IDs of the types we have chosen in the metadata table above
  rownames(md) <- md$PP_ID
  md$PP_ID <- NULL
  
  # Obtain the number of genes that meet padj<0.01 for reference line in histogram
  md <- md[order(rownames(md)),]
  
  assignment[assignment==1] <- Baseline
  assignment[assignment==2] <- Target
  
  assignment <- as.data.frame(cbind(paste(md$Status), assignment))
  colnames(assignment)[1] <- "Correct Assignment"
  colnames(assignment) <- gsub("perm_", "Permutation ", colnames(assignment))
  
  write.table(t(assignment), paste0(Baseline, "-vs-", Target, "-permutation-assignments.txt"), sep = "\t", col.names = F, quote = F)
  
  numdiff$Actual <- NULL
  out_numdiff <- data.frame(paste("Permutation", seq(1, 500, 1)), numdiff$NumDiffGenes)
  colnames(out_numdiff) <- c("Permutation","Number of genes with padj < 0.01")
  
  write.table(out_numdiff, paste0(Baseline, "-vs-", Target, "-number-of-diff-genes.txt"), sep = "\t", row.names = F, quote = F)
}

GeneratePermutationTables(Results <- "../../../tables/permutation/LuCa-vs-HD.number.diff.genes.csv",
                          Assignment <- "../../../tables/permutation/LuCa-vs-HD.permutation.list.csv",
                          Metadata <- "../../../tables/PP_metadata_keep_FINAL_updated.csv",
                          Baseline <- "HD", Target <- "LuCa")
GeneratePermutationTables(Results = "../../../tables/permutation/MM-vs-HD.number.diff.genes.csv", 
                          Assignment = "../../../tables/permutation/MM-vs-HD.permutation.list.csv",
                          Metadata = "../../../tables/PP_metadata_keep_FINAL_updated.csv",
                          Baseline = "HD", Target = "MM")
GeneratePermutationTables(Results = "../../../tables/permutation/HCC-vs-HD.number.diff.genes.csv", 
                          Assignment = "../../../tables/permutation/HCC-vs-HD.permutation.list.csv",
                          Metadata = "../../../tables/PP_metadata_keep_FINAL_updated.csv",
                          Baseline = "HD", Target = "HCC")
