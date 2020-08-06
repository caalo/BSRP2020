# Week 4: Correlation Exercises with CRISPR Avana Data (in-class)
library(tidyverse)

# First, let's load in the CRISPR gene effect data
Achilles_gene_effect <- read_csv("Achilles_gene_effect.csv")

# Let's keep track of cell line metadata also
CCLE_metadata <- read_csv("sample_info.csv")

# What are the dimensions of this data? How many genes? How many cell lines?
dim(Achilles_gene_effect)
View(Achilles_gene_effect[1:5,1:5])

# Modify the column names to just get the gene (remove Entrez ID)
colnames(Achilles_gene_effect) <- sapply(colnames(Achilles_gene_effect), function(x){
  strsplit(x, split = " ", fixed = T)[[1]][1]
})

# Modify the Achilles_gene_effect data to rename cell lines to CCLE names
CCLE_metadata_ids <- CCLE_metadata %>%
  select(DepMap_ID, CCLE_Name)

CCLE_Achilles_merged <- merge(CCLE_metadata_ids, Achilles_gene_effect) %>%
  select(-DepMap_ID)

# What does one cell line's distribution look like?
### Select specific cell line (A549_LUNG), display a histogram
A549_Achilles <- CCLE_Achilles_merged %>%
  filter(CCLE_Name == "A549_LUNG") %>%
  gather(key = "gene", value = "dependency_score", -CCLE_Name)

ggplot(A549_Achilles, aes(x = dependency_score)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, color = "green") +
  geom_vline(xintercept = -1, color = "red") +
  labs(x = "Dependency score", y = "Number of genes", title = "Dependency score distribution in A549_LUNG")


# What does one gene's distribution look like?
### Select specific genes (KRAS, TP53, TCF7L1), display a histogram for each gene's distribution

Achilles_gathered <- CCLE_Achilles_merged %>%
  select(KRAS, TP53, TCF7L1) %>%
  gather(gene, dependency_score)
  
ggplot(data = subset(Achilles_gathered, gene == "KRAS"), aes(x = dependency_score)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = -1, color = "red") +
  geom_vline(xintercept = 0, color = "green") +
  labs(x = "Dependency score", y = "Number of cell lines", title = "Dependency score distribution for KRAS")

ggplot(data = subset(Achilles_gathered, gene == "TP53"), aes(x = dependency_score)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = -1, color = "red") +
  geom_vline(xintercept = 0, color = "green") +
  labs(x = "Dependency score", y = "Number of cell lines", title = "Dependency score distribution for TP53")

ggplot(data = subset(Achilles_gathered, gene == "TCF7L1"), aes(x = dependency_score)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = -1, color = "red") +
  geom_vline(xintercept = 0, color = "green") +
  labs(x = "Dependency score", y = "Number of cell lines", title = "Dependency score distribution for TCF7L1")


# How can we correlate these data together?
KRAS_correlations <- cor(CCLE_Achilles_merged$KRAS, CCLE_Achilles_merged[,-1])
View(t(KRAS_correlations))

BRCA1_correlations <- cor(CCLE_Achilles_merged$BRCA1, CCLE_Achilles_merged[,-1])
View(t(BRCA1_correlations))

# What does this correlation allow us to do?
# (Note: takes a long time! Go get a coffee or something)
gene_cor_matrix <- cor(CCLE_Achilles_merged[,-1], use = 'pairwise.complete') 
dim(gene_cor_matrix)

# How about this one? - note the use of t() to transpose the matrix (switch the rows and columns)
cell_cor_matrix <- cor(t(CCLE_Achilles_merged[,-1]), use = "pairwise.complete")
rownames(cell_cor_matrix) <- CCLE_Achilles_merged[,1] 
colnames(cell_cor_matrix) <- CCLE_Achilles_merged[,1]
dim(cell_cor_matrix)


# Let's make some heatmaps.
if(!require("heatmaply")){
  install.packages("heatmaply")
}

# Pick genes from a specific pathway


my_genes <- c("EGFR", "KRAS", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "TCF7L2")
heatmap(gene_cor_matrix[my_genes, my_genes])

# Show all cell line correlations
heatmap(cell_cor_matrix) # is this plot useful?

which_skin <- grep("SKIN", colnames(cell_cor_matrix)) # Which cell lines have the word "SKIN" in them?
heatmap(cell_cor_matrix[which_skin, which_skin]) # What does this plot show?