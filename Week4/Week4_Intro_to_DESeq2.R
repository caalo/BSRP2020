### An Intro to Differential Expression ###

# install.packages("BiocManager")
# BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)

# Load in the datasets
CCLE_metadata <- read_csv("sample_info.csv")
# Downloaded CCLE_RNAseq reads from here: https://ndownloader.figshare.com/files/22897973
CCLE_expression_raw <- read_csv("CCLE_RNAseq_reads.csv")
CCLE_mutation <- read_csv("CCLE_mutations.csv")

first_word_of_each <- function(x, pattern=" ", n=1){
  spl = strsplit(x, split = pattern, fixed = T)
  return(unlist(sapply(spl, `[[`, n)))
}

# Transforming count data to be DESeq2 compatible.
CCLE_expression <- CCLE_expression_raw %>%
  filter(X1 %in% CCLE_metadata$DepMap_ID) %>%
  column_to_rownames(var = "X1") %>%
  as.matrix() %>%
  round() %>%
  magrittr::set_rownames(CCLE_metadata$CCLE_Name[match(rownames(.), CCLE_metadata$DepMap_ID, nomatch = 0)]) %>%
  magrittr::set_colnames(first_word_of_each(colnames(.))) %>%
  t() # transpose to get samples on columns 

saveRDS(CCLE_expression, file = "CCLE_RNAseq_counts.rds")

# Here, we compare uterus and ovary cell lines to find the top differentially expressed genes between lineages
CCLE_metadata_lineage <- CCLE_metadata %>%
  filter(CCLE_metadata$CCLE_Name %in% colnames(CCLE_expression)) %>%
  filter(lineage %in% c("uterus", "ovary"))

CCLE_expression_lineage <- CCLE_expression[,match(colnames(CCLE_expression), CCLE_metadata_lineage$CCLE_Name, nomatch = 0)]

# Create a DESeq dataset (dds) from the integer counts matrix
dds_lineage <- DESeqDataSetFromMatrix(CCLE_expression_lineage, colData = CCLE_metadata_lineage, design = ~lineage)

dds_lineage_DE <- DESeq(dds_lineage)

res_lineage <- results(dds_lineage_DE)
res_lineage_df <- as.data.frame(res_lineage)

ggplot(data = res_lineage_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point()

plotDispEsts(dds_lineage_DE)



