
# CONVERTIR gaf en gmt
# R
# OUAHBI Dina
# 28022024

# EXPORT GMT FILE 

library(dplyr)
setwd("/shared/ifbstor1/projects/mycovirus")
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(data.table)

names <- read.csv("go_term_names.csv")[,-1]

gaf_data <- read.delim("PEA/Background_genes/FungiDB-65_AfumigatusAf293_Curated_GO.gaf",
                       stringsAsFactors = FALSE,
                       header = FALSE,
                       comment.char = "!")

gmt_data <- gaf_data %>%
  dplyr::select(V2, V5) %>%
  group_by(V5) %>%
  summarize(gene_set = paste(unique(V2), collapse = "\t"))
gmt_data <- as.data.frame(gmt_data)
names(gmt_data)[names(gmt_data) == "V5"] <- "id"

merge_df <- merge(names, gmt_data, by = "id")
merge_df <- merge_df[, c('name', 'gene_set')]
merge_df$name <- gsub(",", ".", merge_df$name)

write.table(merge_df, "PEA/Background_genes/FungiDB-65_AfumigatusAf293_Curated_GO_v2.gmt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)