# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus

# REQUIREMENTS
library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)
library(vsn)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

library(gtsummary)
library(tidyverse)
library(rstatix)
library(ggpubr)

# ===================================================> PREPROCESS <=================================================== #

# READ AND PROCESS
setwd(paste0(getwd(),'/mycovirus'))

for (dir_path in c('pca', 'volcanoplots', 'deseq2', 'mycovirus_tables')) {
  if (!file.exists(dir_path)) {
    # If it doesn't exist, create the directory
    dir.create(dir_path)
    print(paste("Directory", dir_path, "created successfully."))
  } else {
    # If it already exists, print a message indicating so
    print(paste("Directory", dir_path, "already exists."))
  }
}


counts_data <- read.csv('results/counts_matrix.csv')
colData <- read.csv('results/colData.csv')
colData$mycovirus <- ifelse(colData$mycovirus == "", 'virusfree', colData$mycovirus)
select_cols <- c('group', 'condition', 'mycovirus', 'azoleR')
for (col in select_cols){
  colData[[col]] <- as.factor(colData[[col]])
}
counts_data <- counts_data[!duplicated(counts_data$X),]
rownames(counts_data) <- counts_data$X
counts_data$X <- NULL
rownames(colData) <- colData$sampleID
colData$sampleID <- NULL
colData = colData[-which(rownames(colData) == 'souche_af293'), ] # control
counts_data = counts_data[, -which(names(counts_data) %in% c('souche_af293'))]
counts_data <- round(counts_data)

# CHECK
all(colnames(counts_data) %in% rownames(colData)) # SAME PATIENTS?
all(colnames(counts_data) == rownames(colData)) # SAME ORDER?

# NAN
counts_data <- na.omit(counts_data)
col_missing_counts <- colSums(is.na(counts_data))
print(col_missing_counts)


# ===================================================> PCA <=================================================== #
df <- as.data.frame(t(counts_data)) # TRANSPOSE MY COUNT DATA

res.pca <- PCA(df, scale.unit = TRUE, ncp = 10, ind.sup = NULL, 
               quanti.sup = NULL, quali.sup = NULL, row.w = NULL, 
               col.w = NULL, graph = FALSE, axes = c(1,2))

fviz_eig(res.pca, addlabels = TRUE, barcolor = "#87F306", barfill = "#6A6564", ylim = c(0, 75)) # VARIANCE PROPORTIONS


# Color variables by cos
fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_var(res.pca, col.var = "cos2",
             col.circle = "Black",
             gradient.cols = c("blue", "yellow", "red"))



fviz_pca_ind(res.pca,
             col.ind = as.factor(colData$azoleR),
             #pointsize = "cos2",
             title = "azoleR",
             repel = TRUE # Avoid text overlapping (slow if many points)
)


# ===================================================> CORR <=================================================== #
pca_result <- prcomp(df, center = TRUE, scale = FALSE)
my_range <- 1:10
pc_df <-  c()
for (i in my_range) {
  grp <- pca_result$x[, i]
  pc_df <- cbind(pc_df, grp)
}
pc_df <- as.data.frame(pc_df)
colnames(pc_df) <- my_range

char_pc <- merge(pc_df, colData, 
      by = 'row.names', all = TRUE)
rownames(char_pc) <- char_pc$Row.names
char_pc$Row.names <- NULL


for (col in colnames(colData)) {
  tf <- paste0('pca/',col,'.rtf')
  char_pc %>% select(my_range, col)%>%
  tbl_summary(by = col,
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c("{median} ({p25}, {p75})"),
              missing_text = "NA"
  ) |>
    add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) |>
    add_n() |>
    add_q(method="BH") |>
    modify_header(label ~ "**Variable**") |>
    bold_labels() |>
    as_gt() |>
    gt::gtsave(filename = tf)
}





library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(pheatmap)


cuttof <- 1.5

deseq_by_comb <- function(design, comb, colData) {
  ref <- comb[1]
  sample <- colData[colData$temp %in% comb, ]
  counts_sample <- counts_data[rownames(sample)]
  dse <- DESeqDataSetFromMatrix(countData = counts_sample,
                                colData = sample,
                                design = ~temp)
  
  dse$temp = relevel(dse$temp, ref=ref)
  dds=DESeq(dse,   parallel = T)
  res = results(dds)
  title <- paste0(design,"_",ref,"_vs_",comb[2])
  write.csv2(res, paste0("deseq2/",title,".csv"), row.names = T)
  ##
  res_volcano <- as.data.frame(res[,c('log2FoldChange','pvalue','padj')])
  res_volcano$Gene = rownames(res_volcano) #indexer sur le nom
  res_volcano$color <- ifelse(res_volcano$pvalue<0.05,'p-value<0.05','NS')
  res_volcano$color <- ifelse(abs(res_volcano$log2FoldChange)>cuttof,paste0('abs(LFC)>',cuttof),res_volcano$color)
  res_volcano$color <- ifelse(res_volcano$pvalue<0.05 & abs(res_volcano$log2FoldChange)>cuttof,paste0('p-value <0.05 & abs(LFC)>', cuttof),res_volcano$color)
  res_volcano$color = as.factor(res_volcano$color)
  write.csv(table(res_volcano$color), paste0("mycovirus_tables/",title,".csv"))
  res_volcano = res_volcano[c(which(!is.na(res_volcano$pvalue))),]
  ##
  ggplot(res_volcano, aes(x = log2FoldChange, y = -log10(pvalue), colour=color)) + geom_point()+
    scale_colour_manual(values=c("#971A9B", "#424242", "#A72424", "#148811"))+
    theme_classic() +
    geom_text_repel(
      data = subset(res_volcano, (pvalue < 0.05 & abs(log2FoldChange)>cuttof)),
      aes(label = Gene),
      size = 4,
      colour="black",
      box.padding = unit(0.4, "lines"),
      point.padding = unit(0.2, "lines")
    )+
    theme_bw(base_size = 12) + theme(legend.title = element_blank(),legend.position = "bottom", legend.text = element_text(size=14)) +
    xlab('log2FoldChange')+ ylab('-log10 adjusted p-value')+
    ggtitle(paste0(design,"_",ref,"_vs_",comb[2]))
  ggsave(paste0("volcanoplots/",design,"_",ref,"_vs_",comb[2],".png"))
  
}

subset_group <- function(group){
  sub <- colData[colData$group %in% group, ]
  return(sub)
}

#---------------------------------------------------------------------------------------------------
#Analyse transcriptomique : groupe 

## BETWEEN GROUP
design <- "group"
colData$temp <- colData$group
deseq_by_comb(design,c('env', 'api'), colData)
deseq_by_comb(design,c('env', 'col'), colData)
deseq_by_comb(design,c('col', 'api'), colData)

## condition ALL GROUPS
design <- "condition"
colData$temp <- colData$condition
deseq_by_comb(design,c('virusfree', 'infected'), colData)

## azoleR ALL GROUPS
design <- "azoleR"
colData$temp <- colData$azoleR
deseq_by_comb(design,c('no', 'yes'), colData)

## PER GROUP
for (group in c('env', 'api', 'col')){
  design <- paste0("condition_",group)
  subset_gp <- subset_group(group)
  subset_gp$temp <- subset_gp$condition
  deseq_by_comb(design,c('virusfree', 'infected'), subset_gp)
}


#---------------------------------------------------------------------------------------------------
#Analyse transcriptomique : groupe 
colData <- read.csv('clustering/kmean_clustering.csv', sep=';')
rownames(colData) <- colData$X
colData$X <- NULL

# CHECK
all(colnames(counts_data) %in% rownames(colData)) # SAME PATIENTS?
all(colnames(counts_data) == rownames(colData)) # SAME ORDER?

## K2
design <- "X2"
colData$temp <- as.factor(colData$X2)
deseq_by_comb(design,c('1', '2'), colData)

## K3
design <- "X3"
colData$temp <- as.factor(colData$X3)
deseq_by_comb(design,c('1', '2'), colData)
deseq_by_comb(design,c('1', '3'), colData)
deseq_by_comb(design,c('2', '3'), colData)

## K4
design <- "X4"
colData$temp <- as.factor(colData$X4)
deseq_by_comb(design,c('1', '2'), colData)
deseq_by_comb(design,c('1', '3'), colData)
deseq_by_comb(design,c('1', '4'), colData)
deseq_by_comb(design,c('2', '3'), colData)
deseq_by_comb(design,c('2', '4'), colData)
deseq_by_comb(design,c('3', '4'), colData)




