### An Intro to Differential Expression ###

install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)

# Load in the datasets
CCLE_metadata <- read_csv("sample_info.csv")
CCLE_expression <- read_csv("CCLE_RNAseq_reads.csv")
CCLE_mutation <- read_csv("CCLE_mutations.csv")

# Create a DESeq dataset (dds)
CCLE_metadata_lineage <- CCLE_metadata %>%
  filter(CCLE_metadata$DepMap_ID %in% CCLE_expression$X1) %>%
  filter(lineage %in% c("uterus", "ovary"))

CCLE_expression_lineage <- CCLE_expression %>%
  filter(CCLE_expression$X1 %in% CCLE_metadata_lineage$DepMap_ID) %>%
  select(-X1) %>%
  as.matrix() %>%
  round()


dds <- DESeqDataSetFromMatrix(t(CCLE_expression_lineage), colData = CCLE_metadata_lineage, design = ~lineage)

dds <- DESeq(dds)

res <- results(dds)
res_df <- as.data.frame(res)

ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point()
