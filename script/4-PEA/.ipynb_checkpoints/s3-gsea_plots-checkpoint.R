# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


# EXPORT histogram of NES score / barplot of up and down regulated sets / gsea table / enrichment

# Load required libraries
#https://biostatsquid.com/fgsea-tutorial-gsea/
#https://www.genekitr.fun/plot-gsea-1

library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(data.table)

set.seed(123456)
setwd(paste0('/shared/projects/braf_mutated_melanoma','/mycovirus'))

rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F)

# SET PATHS
in_path <- "deseq_v2/"
out_path <- "PEA/results/"
bg_path <- "PEA/Background_genes/"

# -------------------------------- START  ------------------------------------------

matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}
create_dir <- function(dir_name) {
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
    print(paste("Directory", dir_name, "created."))
  } else {
    print(paste("Directory", dir_name, "already exists."))
  }
}

# Create directories
dir_names <- c("PEA/mainPathways", "PEA/imgs", "PEA/imgs/barplot", "PEA/imgs/enrichment", "PEA/imgs/GseaTable", "PEA/imgs/hist")
for (dir_name in dir_names) {
  create_dir(dir_name)
}

# --------------------------------------------------------------------------------------------

list.files(in_path)
list.files(bg_path)
recap <- c()
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)

for (file in list.files(out_path)){
  print(file)
  filename <- strsplit(file, ".", fixed = TRUE)[[1]]
  df <- read.csv(paste0(in_path, file), row.names = 1, sep = ',')[, -1]
  df$Symbol <- str_to_title(tolower(gsub("AFUA_", "Afu", df$Symbol)))
  my_genes <- df$Symbol
  numeric_cols <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue", "padj")
  for (col in numeric_cols) {
    df[[col]] <- gsub(",", ".", df[[col]])
    df[[col]] <- as.numeric(df[[col]])
  }
  gene_list <- df$log2FoldChange
  names(gene_list) = my_genes
  
  gene_list = sort(gene_list, decreasing = TRUE) # ca enleve en meme temps les na
  gene_list = gene_list[!duplicated(names(gene_list))]
  db <- prepare_gmt(gmt_files[1], df$Symbol, savefile = FALSE)
  
  GSEAres <-as.data.table(read.csv(paste0(out_path, file), sep = "\t", header = TRUE))  
  pval <- dim(subset(GSEAres, (pval < 0.05)))[1]
  padj <- dim(subset(GSEAres, (padj < 0.05)))[1]
  
  topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = 20), pathway]
  topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = 20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  collapsedPathways <- collapsePathways(GSEAres[order(padj)][padj < 0.05], db, gene_list)
  mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
  collapse <- length(mainPathways)
  
  write.csv(as.matrix(mainPathways), paste0("PEA/mainPathways/",file))
  best_term <- GSEAres[order(pval),][1,]$pathway
  leadingEdge <- length(strsplit(GSEAres[order(pval),][1,]$leadingEdge, " ")[[1]])
  path <- length(db[[GSEAres[order(pval),][1,]$pathway]])
  
  # -------------------------------------------------------------------------------------
  
  # HIST
  png(paste0(paste0("PEA/imgs/hist/",filename,".png")), width = 400, height = 300)
  hist(GSEAres$NES[GSEAres$pval<0.05], breaks = 20,main = file, col="#FF0000")
  dev.off()
  
  # GSEA TABLE
  pdf(file = paste0("PEA/imgs/GseaTable/",filename, '.pdf'), width = 20, height = 15)
  plotGseaTable(db[topPathways], stats = gene_list, fgseaRes = GSEAres, gseaParam = 0.5)
  dev.off()
  
  # BARPLOT
  setorder(GSEAres, NES)
  GSEAres_subset <- rbind(head(GSEAres, 20), tail(GSEAres, 20))
  gg<-ggplot(GSEAres_subset, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=filename) +
    theme_minimal()
  
  pdf(paste("PEA/imgs/barplot/", filename, ".pdf", sep = ""), width = 20, height = 15)
  print(gg)
  dev.off()
  
  # enrichment
  gg <- plotEnrichment(db[[head(GSEAres[order(pval), ], 1)$pathway]],
                       gene_list) + labs(title=paste0(head(GSEAres[order(pval), ], 1)$pathway,"\t",filename))
  pdf(paste("PEA/imgs/enrichment/", filename, ".pdf", sep = ""), width = 8, height = 6)
  print(gg)
  dev.off()
  
  my_list <- c(pval, padj, collapse, best_term, leadingEdge, path)
  recap <- cbind(recap,my_list)
}

# recapitulatif
recap <- as.data.frame(recap)
rownames(recap) <- c("significant term pval", "significant term padj", "NO term after collapse", "best term", "leadingEdge", "real pathway")
colnames(recap) <- list.files(out_path)
write.csv(recap, "PEA/gsea_recap.csv")



































