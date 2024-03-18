setwd("/shared/ifbstor1/projects/mycovirus/")
FILES <-list.files('genes_significatives/', pattern = "\\.csv$", full.names = F)


for (file in FILES) {
  file <- sub('\\.csv$', '', file) 
  df <- read.table(paste0('genes_significatives/',file,'.csv'), sep=',')
  colnames(df) <- df[1, ]
  df <- df[-1, ]
  names_genes <- df$Symbol
  genes_vect <- c(strsplit(names_genes, split = ' '))
  # export without double quotations
  writeLines(unlist(genes_vect), con = paste0('fungifun/',file,'.txt'))
}
