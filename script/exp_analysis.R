

# ------------------------ AUTHOR  / OUAHBI DINA
# ------------------------ DATE  / 06-19-2023
# ------------------------ MAIN  / EXPRESSION GENETIQUE DIFFERENTIELLE ENTRE LES REPONDEUR ET NON A LA THERAPIE CIBLEE DANS LE CANCER CUTANEE

# my file list
list.files(path=".", pattern=NULL, all.files=FALSE,
           full.names=FALSE)

# REQUIREMENTS
library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(FactoMineR)
library(factoextra)
library("corrplot")

# READ DATAS
setwd("/shared/ifbstor1/projects/mycovirus/")
counts_data <- read.csv('results/counts_matrix.csv')
colData <- read.csv('results/colData.csv')

# indexs
rownames(counts_data) <- counts_data$X
counts_data$X <- NULL
rownames(colData) <- colData$sampleID
colData$sampleID <- NULL

# remove reference
counts_data = subset(counts_data, select = -c(souche_af293))
row_names_df_to_remove<-c("souche_af293")
colData<-colData[!(row.names(colData) %in% row_names_df_to_remove),]

# replace empty string in mycovirus 
colData[colData==""]<-"virusfree"

# stats ==========================================================================
set.seed(10)
str(colData)
table(colData$group)
table(colData$condition)
table(colData$mycovirus)
table(colData$azoleR)

# Prop of each category
prop.table(table(colData$group))

# barplot : frequency
par(mfrow=c(2,2))
barplot(prop.table(table(colData$group)), ylab = 'group')
barplot(prop.table(table(colData$condition)), ylab = 'condition')
barplot(prop.table(table(colData$mycovirus)), ylab = 'mycovirus')
barplot(prop.table(table(colData$azoleR)), ylab = 'azoleR')

# relationship two variables
table(colData$group, colData$condition)
table(colData$group, colData$mycovirus)
table(colData$group, colData$azoleR)

#contingency table
library(scales) # to show percentages on the y-axis
par(mfrow=c(1,3))
x = barplot(prop.table(table(colData$condition,colData$group)),
            col = rep(c("#00AFBB", "#FC4E07")),
            legend = TRUE,
            ylim = c(0, 0.5),
            yaxt = 'n', # remove y-axis
            ylab = 'Percent of stump')
# creating a y-axis with percentages
yticks = seq(0, 0.5, by = 0.05)
axis(2, at = yticks, lab = percent(yticks))
title("GROUP")
grid()

x = barplot(prop.table(table(colData$condition, colData$azoleR)),
            col = rep(c("#00AFBB","#FC4E07")),
            legend = TRUE,
            ylim = c(0, 1),
            yaxt = 'n', # remove y-axis
            ylab = 'Percent of stump')
# creating a y-axis with percentages
yticks = seq(0, 1, by = 0.05)
axis(2, at = yticks, lab = percent(yticks))
title('AZOLER')
grid()

x = barplot(prop.table(table(colData$condition, colData$mycovirus)),
            col = rep(c("#00AFBB", "#FC4E07")),
            legend = TRUE,
            ylim = c(0, 1),
            yaxt = 'n', # remove y-axis
            ylab = 'Percent of stump')
# creating a y-axis with percentages
yticks = seq(0, 1, by = 0.05)
axis(2, at = yticks, lab = percent(yticks))
title('MYCOVIRUS')
grid()



# PCA===========================================================================
df <- as.data.frame(t(counts_data))
res.pca <- PCA(df, graph = FALSE)

print(res.pca)

# variance per dim
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 80),barcolor = "black",barfill = "#FC4E07",)

# Color variables by groups
fviz_pca_var(res.pca, col.var = colData$condition, 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Condition")



fviz_pca_ind(res.pca, col.ind = colData$mycovirus, 
             repel = TRUE, # Avoid text overlapping (slow if many points)
             legend.title = "Group"
)

# Contributions of variables to PC1
par(mfrow=c(1,1))
fviz_contrib(res.pca, choice = "var", axes = 2, top = 200)

fviz_pca_var(PCA(df[,1:100], graph = FALSE), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE)

var <- get_pca_var(res.pca)
print(var)
contrib_df <- as.data.frame(var$contrib)
contrib_df_sorted <- contrib_df[order(-contrib_df$Dim.1), ]
tmp <- as.matrix(head(contrib_df, 10))
corrplot(tmp, method = c("number"), is.corr=FALSE, bg="#00AFBB",addgrid.col = T, tl.col = "black")


# CORRELATION TEST =============================================================
library(tidyverse)
library(rstatix)
library(ggpubr)
library(gtsummary)

pca_result <- prcomp(df, center = TRUE, scale = FALSE)
pc1_scores <- pca_result$x[, 1]

pc_df <- as.data.frame(pca_result$x[, 1])
colnames(pc_df) <- "pc1"
pc_df$pc2 <- pca_result$x[, 2]
pc_df$pc3 <- pca_result$x[, 3]
pc_df$pc4 <- pca_result$x[, 4]
pc_df$pc5 <- pca_result$x[, 5]

colData$group <- as.factor(colData$group)
colData$condition <- as.factor(colData$condition)
colData$mycovirus <- as.factor(colData$mycovirus)
colData$azoleR <- as.factor(colData$azoleR)

pc_df <- cbind(pc_df, colData)

# syndé les stat descriptive par groupe : VARIABLES CATEGORIELLE
pc_df%>%
  tbl_summary(by = group,
              #type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c("{median} ({p25}, {p75})"),
              missing_text = "NA"
  ) %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2))%>%
  add_n()%>%
  add_q(method="BH")%>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

# boxplot
# Visualisation : Boxplots avec p-values
pwc <- wilcox.test(pc_df$pc3~pc_df$condition)
print(pwc)
#pwc <- pwc %>% add_xy_position(x = "X1ECHEC...PSVr...2.4")
ggboxplot(pc_df %>% filter(!is.na(azoleR)), x = "group", y = "pc3", fill = 'group', notch = T) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))


#-------------------------------------------------------DESEQ2---------------------------------------------------------------------------------

# REQUIREMENTS
library(tximport)
library(DESeq2)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(plotly) #3D
library(EnhancedVolcano)
library(readxl)
library(FactoMineR)


all(colnames(counts_data) %in% rownames(colData)) # SAME PATIENTS?
all(colnames(counts_data) == rownames(colData)) # SAME ORDER?

counts_data <- round(counts_data)

# NANS
col_missing_counts <- colSums(is.na(counts_data))
print(col_missing_counts)

# SELECT INTERET COLS FROM COLData
samples = data.frame(sample = rownames(colData), 
                     azoleR = as.factor(colData$azoleR),
                     condition = as.factor(colData$condition))
rownames(samples) <- samples$sample

# ------------------------------  DDS OBJECTS -------------------------------------------

dseAzoler <- DESeqDataSetFromMatrix(countData = counts_data,
                                colData = samples,
                                design = ~ azoleR)

dseCondition <- DESeqDataSetFromMatrix(countData = counts_data,
                                 colData = samples,
                                 design = ~ condition)

# FILTER
dseAzoler <- dseAzoler[rowSums(counts(dseAzoler)) > 5, ] 
dseCondition <- dseCondition[rowSums(counts(dseCondition)) > 5, ] 

# Differential expression analysis
ddsAzoler=DESeq(dseAzoler,   parallel = T)
ddsCondition=DESeq(dseCondition, parallel = T)

# Results
resAzoler = results(ddsAzoler)
resCondition = results(ddsCondition)

# LCF SHRINKAGE
resultsNames(ddsAzoler)
resLFCAzoler <- lfcShrink(ddsAzoler, coef="azoleR_yes_vs_no", type="apeglm")
resultsNames(ddsCondition)
resLFCCondition <- lfcShrink(ddsCondition, coef="condition_virusfree_vs_infected", type="apeglm")

# Dispersion estimates
plotDispEsts(ddsAzoler, genecol = "black", fitcol = "red", finalcol = "green")
title(main="Azoler")

plotDispEsts(ddsCondition, genecol = "black", fitcol = "red", finalcol = "green")
title(main="Condition")

# MA PLOT
plotMA(resAzoler)
title(main=expression(azoleR_yes_vs_no))

plotMA(resCondition)
title(main=expression(condition_virusfree_vs_infected))

plotMA(resLFCAzoler)
title(main=expression(azoleR_yes_vs_no ~ "After shinkage"))

plotMA(resLFCCondition)
title(main=expression(condition_virusfree_vs_infected ~ "After shinkage"))

# VOLCANO PLOT
library(EnhancedVolcano)

# label genes with their name instead of their ensembl id
all(
  all(rownames(resAzoler) == rownames(resCondition)),
  all(rownames(resAzoler) == rownames(resLFCAzoler)),
  all(rownames(resAzoler) == rownames(resLFCCondition))
)


# Create a volcano plot
EnhancedVolcano(
  resLFCAzoler,
  lab = rownames(resLFCAzoler),  # Gene names
  x = 'log2FoldChange',    # Column with LFC values
  y = 'pvalue',            # Column with p-values
  drawConnectors = T,
  pCutoff = 1e-03,
  widthConnectors = 0.50,
  max.overlaps = Inf,
  endsConnectors = "first",
  lengthConnectors = unit(0, "npc"), labSize = 3
) + ggtitle(expression(azoleR_yes_vs_no ~ "After shinkage"))


#====================================================== groupe deux par deux
## SUBSET DATA : api vs col
table(colData$group)
colData1 <- colData[colData$group %in% c("col", "env"), ]
samples = data.frame(sample = rownames(colData1), 
                     group = as.factor(colData1$group))

counts_data1 <- counts_data[, names(counts_data) %in% rownames(colData1)]

dseGroup <- DESeqDataSetFromMatrix(countData = counts_data1,
                                    colData = samples,
                                    design = ~ group)

# FILTER
dseGroup <- dseGroup[rowSums(counts(dseGroup)) > 5, ] 
# Differential expression analysis
ddsGroup=DESeq(dseGroup,   parallel = T)
resGroup = results(ddsGroup)
resultsNames(ddsGroup)
resLFCGroup <- lfcShrink(ddsGroup, coef="group_env_vs_col", type="apeglm")
plotDispEsts(ddsGroup, genecol = "black", fitcol = "red", finalcol = "green")
title(main="group_env_vs_col")
plotMA(resGroup)
title(main=expression(group_env_vs_col))
plotMA(resLFCGroup)
title(main=expression(group_env_vs_col ~ "After shrinkage"))

# Create a volcano plot
EnhancedVolcano(
  resGroup,
  lab = rownames(resGroup),  # Gene names
  x = 'log2FoldChange',    # Column with LFC values
  y = 'pvalue',            # Column with p-values
  drawConnectors = T,
  widthConnectors = 0.50,
  pCutoff = 1e-03,
  max.overlaps = Inf,
  endsConnectors = "first",
  lengthConnectors = unit(0, "npc"), labSize = 4
) + ggtitle(expression(group_env_vs_col))



#====================================================== Comparer les mycovirus condition dans chaque groupe

table(colData$group)
colData1 <- colData[colData$group == "col", ]
samples = data.frame(sample = rownames(colData1), 
                     condition = as.factor(colData1$condition))
counts_data1 <- counts_data[, names(counts_data) %in% rownames(colData1)]

dseCondition <- DESeqDataSetFromMatrix(countData = counts_data1,
                                   colData = samples,
                                   design = ~ condition)
dseCondition <- dseCondition[rowSums(counts(dseCondition)) > 5, ] 
ddsCondition=DESeq(dseCondition,   parallel = T)
resCondition = results(ddsCondition)
resultsNames(ddsCondition)

# Create a volcano plot
EnhancedVolcano(
  resCondition,
  lab = rownames(resCondition),  # Gene names
  x = 'log2FoldChange',    # Column with LFC values
  y = 'padj',            # Column with p-values
  drawConnectors = T,
  widthConnectors = 0.50,
  #pCutoff = 1e-03,
  max.overlaps = Inf,
  endsConnectors = "first",
  lengthConnectors = unit(0, "npc"), labSize = 2
) + ggtitle(expression(condition_virusfree_vs_infected ~ "Col" ))

















ddsTxi <- DESeqDataSetFromMatrix(countData = counts_data,
                                 colData = samples,
                                 design =~ group)

#filter
ddsTC <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTC$group = relevel(ddsTC$group, ref="env") # SET REFERENCE
dds <- DESeq(ddsTC)
res <- results(dds)
res <- results(dds, pAdjustMethod = "BH")



# MA PLOT
plotMA(res, ylim=c(-2,2))

# NORMALISATION
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# ORDER DDS RESULTS
resOrdered <- res[order(res$pvalue),]


# ce qui nous importe le plus dans les res, c'est la p value et le fold change
# la p value pour reperer les genes sign
# le fold change c'est pour voir dans quel sensle gene est different 
summary(res)
# add title, label, new color to boxplot
boxplot(res$padj,
        main="condition",
        ylab="counts",
        xlab="P value adjuster",
        col="orange", notch=TRUE)


summary(res$padj)
summary(res$log2FoldChange)
hist(res$log2FoldChange) # si le foldchange est proche de 0 

# selectionner les genes diff avec un taux d'erreur de 5% et un logfold change >2
res[which(res$padj<0.05 & abs(res$log2FoldChange)>2),]

# --- Enregistrement ---
write.csv2(res, "deseq2/Res_allgenes_envApiColExcluded _Deseq2_RSEM_group.csv", row.names = T)

# ------------------------------  DDS / VOLCANO PLOT --------------------------------------------

res_volcano <- as.data.frame(res[,c('log2FoldChange','pvalue','padj')])
res_volcano$Gene = rownames(res_volcano) #indexer sur le nom
res_volcano$color <- ifelse(res_volcano$padj<0.05,'adjusted p-value<0.05','NS')
res_volcano$color <- ifelse(abs(res_volcano$log2FoldChange)>1,'abs(LFC)>1',res_volcano$color)
res_volcano$color <- ifelse(res_volcano$padj<0.05 & abs(res_volcano$log2FoldChange)>1,'ajusted p-value <0.05 & abs(LFC)>1',res_volcano$color)
res_volcano$color = as.factor(res_volcano$color)
table(res_volcano$color)
res_volcano = res_volcano[c(which(!is.na(res_volcano$padj))),]

# PLOTING
library(ggplot2); library(ggrepel)
ggplot(res_volcano, aes(x = log2FoldChange, y = -log10(padj), colour=color)) + geom_point()+
  scale_colour_manual(values=c("orange", "green", "red", "grey"))+
  theme_classic() +
  geom_text_repel(
    data = subset(res_volcano, (padj < 0.05 & abs(log2FoldChange)>1)),
    aes(label = Gene),
    size = 4,
    colour="black",
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.2, "lines")
  )+
  theme_bw(base_size = 12) + theme(legend.title = element_blank(),legend.position = "bottom", legend.text = element_text(size=14)) +
  xlab('log2FoldChange')+ ylab('-log10 adjusted p-value')+
  ggtitle("all vs env")



# ma plot is basicaly a scatter plot of the log fold change versus the mean of the normalise count. this plot tell us
# the genes that are differentialy expressed, and the genes in blues are sign defferently expressed genes and they have 
# adjusted p value less than 0.05
# you can also see some smaler triangles in the edges of the plot. this genes have higher FC, and the dirrection of triangles
# tell us the direction of FC.

# si presence de gene bleu dans le quadrant superieur droit ou inferieur droit


#------------------------------- HEATMAP

data = assay(ntd)
# selectionner des genes significative et tres differents (LFC > 3)
names_genes = res_volcano$Gene[which(res_volcano$padj<0.05 & abs(res_volcano$log2FoldChange)>10)]
data = data[c(which(rownames(data) %in% names_genes)),]

df <- data.frame(colData(dds)[,c("group")]) # select gene name and class
rownames(df) = colData(dds)[,1] ; names(df) = "group" #rename
ann_colors = list("group" = c('col' = "green", 'env' = "red", 'api' = 'blue', 'excluded' = 'yellow'))  #choose color between rep and non rep

data = data[, c(which(colnames(data) %in% rownames(df)[which(df$group=="col")]),
                which(colnames(data) %in% rownames(df)[which(df$group=="env")]),
                which(colnames(data) %in% rownames(df)[which(df$group=="api")]),
                which(colnames(data) %in% rownames(df)[which(df$group=="excluded")])
)]

pheatmap(data, cluster_rows=T, show_rownames=T,border_color="black",treeheight_col=0,
         color = colorRampPalette(rev(brewer.pal(n=10, name = "RdBu")))(100),
         cluster_cols=F, show_colnames = T, annotation_col=df, annotation_colors = ann_colors)





















































