#!/usr/bin/env Rscript

# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


# EXPORT GSEAres

library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(data.table)

set.seed(123456)

rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F)

## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt --------------------------------------
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

# PATH TO CHANGE
in_path <- "deseq_v2/"
out_path <- "PEA/results/"
bg_path <- "PEA/Background_genes/"

ifelse(!dir.exists(file.path("PEA", "results")), dir.create(file.path("PEA", "results")), FALSE)


list.files(in_path)
list.files(bg_path)

sigs_set = c()
for (file in list.files(in_path)) {
  print(file)
  df <- read.csv(paste0(in_path, file), row.names = 1, sep = ',')[, -1]
  df$Symbol <- str_to_title(tolower(gsub("AFUA_", "Afu", df$Symbol)))
  my_genes <- df$Symbol
  gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
  numeric_cols <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue", "padj")
  for (col in numeric_cols) {
    df[[col]] <- gsub(",", ".", df[[col]])
    df[[col]] <- as.numeric(df[[col]])
  }
  gene_list <- df$log2FoldChange
  names(gene_list) = my_genes
  gene_list = sort(gene_list, decreasing = TRUE) # ca enleve en meme temps les na
  gene_list = gene_list[!duplicated(names(gene_list))]
  db <- prepare_gmt(gmt_files[1], my_genes, savefile = FALSE)
  GSEAres <- fgsea(pathways = db, # List of gene sets to check
                   stats = gene_list,
                   scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                   minSize = 2,
                   maxSize = 500,
                   nproc = 1) # for parallelisation
  
  x=length(rownames(subset(df, df$padj < 0.05 & abs(df$log2FoldChange) > 2)))
  print(x)
  sigs_set = append(sigs_set, x)
  fwrite(GSEAres, file=paste0(out_path, file), sep="\t", sep2=c("", " ", ""))
}

# Création du graphique à barres
par(mar = c(5, 5, 4, 2) + 6)
png("plots/significant_sets.png")
barplot(sigs_set, names.arg = "", ylab = "NO sign sets", main = "Barplot", las=1,ylim = c(0, 250), col=rainbow(length(sigs_set)))
text(x = barplot(sigs_set, plot = FALSE), y = -0.5, labels = list.files(in_path), srt = 40, adj = 1, xpd = TRUE, cex = 0.8)
grid(col = "gray", lty = "dotted")
dev.off()
