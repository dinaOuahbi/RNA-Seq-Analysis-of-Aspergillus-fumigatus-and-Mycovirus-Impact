# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


# CHA
#https://medium.com/@zullinira23/implementation-of-principal-component-analysis-pca-on-k-means-clustering-in-r-794f03ec15f
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cluster)
library(factoextra)
library(FactoMineR)
library(gridExtra)
library(devtools)

#Préparation des données (matrice numerique : row = obs / col = vars)
setwd(paste0(getwd(),'/mycovirus'))

# Check if the directory already exists

for (dir_path in c('plots', 'clustering')) {
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
rownames(counts_data) <- counts_data$X
counts_data$X <- NULL
counts_data = counts_data[, -which(names(counts_data) %in% c('souche_af293'))]

counts_data <- na.omit(counts_data)
col_missing_counts <- colSums(is.na(counts_data))
print(col_missing_counts)

df <- as.data.frame(t(counts_data))
which(apply(df, 2, var)==0)
df <- df[ , which(apply(df, 2, var) != 0)]

# PCA
res.pca <- PCA(df, scale.unit = TRUE, ncp = 10, graph = FALSE)
png("plots/screen_plot.png")
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100))
dev.off()

# SELECTIONNER LES 2 PREMIERS PC OPTIMALES
pca_df = prcomp(df, center = TRUE, scale = TRUE)
df_transform = as.data.frame(-pca_df$x[,1:2])

# AFFICHER LES CLUSTER OPTIMAL
plots_list <- list()
methods <- c('silhouette', 'wss', 'gap_stat')
for (method in methods) {
  plot <- fviz_nbclust(df_transform, kmeans, method = method, k.max = 10)
  plots_list[[length(plots_list) + 1]] <- plot
}
png("plots/k_optimal_cluster.png")
grid.arrange(grobs = plots_list, ncol = as.integer(length(plots_list)/2))
dev.off()

# AFFICHER LES CLUSTER DE SOUCHES / CLUSTER EXPORT
set.seed(0)
k_df <-  c()
plots_list <- list()
clusters <- c(2,3,4,5,6,8,10)
for (cluster in clusters) {
  res.kmeans <- kmeans(df_transform, centers = cluster, nstart = 50)
  k_df <- cbind(k_df, res.kmeans$cluster)
  plot <- fviz_cluster(res.kmeans, data = df_transform, main = paste0("kmean ",cluster, " clusters"),labelsize = 6, pointsize = 1,ellipse.alpha = 0)+
    scale_color_manual(values = rainbow(cluster))
  plots_list[[length(plots_list) + 1]] <- plot
}
png("plots/k_optimal_cluster_pca.png")
grid.arrange(grobs = plots_list, ncol = as.integer(length(plots_list)/2))
dev.off()

k_df <- as.data.frame(k_df)
colnames(k_df) <- clusters
table(k_df$`2`)
write.csv2(k_df, "clustering/kmean_clustering.csv", row.names = T)



