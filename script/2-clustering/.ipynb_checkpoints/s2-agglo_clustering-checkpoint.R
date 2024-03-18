# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


# CHA
#https://www.datanovia.com/en/lessons/agglomerative-hierarchical-clustering/
library(factoextra)
library(FactoMineR)
library(factoextra)
library(stats)
library(dplyr)
library(gridExtra)


#Préparation des données (matrice numerique : row = obs / col = vars)
setwd(paste0(getwd(),'/mycovirus'))
counts_data <- read.csv('results/counts_matrix.csv')
rownames(counts_data) <- counts_data$X
counts_data$X <- NULL
counts_data = counts_data[, -which(names(counts_data) %in% c('souche_af293'))]

counts_data <- na.omit(counts_data)
col_missing_counts <- colSums(is.na(counts_data))
print(col_missing_counts)
df <- as.data.frame(t(counts_data))
df <- df[ , which(apply(df, 2, var) != 0)]
df_scale <- as.data.frame(scale(df)) 


#Calcul des informations de (dis)similarité entre chaque paire d'objets de l'ensemble de données.
res.dist <- dist(df_scale, method = "euclidean")
as.matrix(res.dist)[1:6, 1:6]

pcoa <- cmdscale(res.dist, k = 89, eig = TRUE)
data.pcoa <- pcoa$points %>%
  as.data.frame.matrix() %>%
  rename_with(~ gsub("V", "PCoA", .x, fixed = TRUE)) #' Evaluate number of kmeans cluster

#Utilisation de la fonction de liaison pour regrouper les objets dans une arborescence de cluster hiérarchique,
#en fonction des informations de distance générées à l'étape 1. Les objets/clusters proches les uns des autres 
#sont liés entre eux à l'aide de la fonction de liaison.
res.pca <- PCA(df_scale, scale.unit = TRUE, ncp = 10, ind.sup = NULL, 
               quanti.sup = NULL, quali.sup = NULL, row.w = NULL, 
               col.w = NULL, graph = FALSE, axes = c(1,2))

methods <- c('ward.D','ward.D2','single','complete','average','mcquitty','median','centroid')
k=2
hc_df <-  c()
dend_list <- list()
ind_list <- list()
for (method in methods) {
  res.hc <- hclust(d = res.dist, method = method)
  grp <- cutree(res.hc, k = k)
  hc_df <- cbind(hc_df, grp)
  res.coph <- cophenetic(res.hc) # hauteur = distance cophénétique : inversement proportionnelle a la similarité
  dend <- fviz_dend(res.hc, k=k, color_labels_by_k = T, cex = 0.5,
                    main = paste0("HC (",method,") Corr res.dist/res.coph : ",round(cor(res.dist, res.coph), 3))) +
    scale_color_manual(values = rainbow(3))
  ind <- fviz_pca_ind(res.pca,col.ind = as.factor(grp),pointsize = 1,title = method,repel = TRUE) +
    scale_color_manual(values = rainbow(3))
  
  dend_list[[length(dend_list) + 1]] <- dend
  ind_list[[length(ind_list) + 1]] <- ind
}
png("plots/clus_dend.png")
grid.arrange(grobs = dend_list, ncol = as.integer(length(dend_list)/2))
dev.off()

png("plots/clus_pca.png")
grid.arrange(grobs = ind_list, ncol = as.integer(length(ind_list)/2))
dev.off()
hc_df <- as.data.frame(hc_df)
colnames(hc_df) <- methods

write.csv2(hc_df, "clustering/aglomerative_clust_k2.csv", row.names = T)






