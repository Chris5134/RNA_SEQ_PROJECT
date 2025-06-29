#!/usr/bin/env Rscript
# scripts/sex_DE_analysis.R

# Load required libraries
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("DESeq2")
}
library(DESeq2)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# 1. Read in counts and metadata
counts <- as.matrix(
  read.csv("data/raw/counts.csv", row.names = 1, check.names = FALSE)
)
coldata <- read.csv(
  "data/raw/metadata.tsv",
  sep = "\t", row.names = 1, check.names = FALSE
)

# 2. Align sample names: replace counts columns with metadata rownames
colnames(counts) <- rownames(coldata)

# 3. Sanity check\stopifnot(all(colnames(counts) == rownames(coldata)))

# 4. Simplify genotype labels
coldata$genotype <- ifelse(
  coldata$genotype == "5xFADHEMI:CLUh2kbKIHO",
  "5xFAD_CLU",
  "5xFAD"
)

# 5. Build DESeq2 dataset and filter low counts
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ sex + genotype
)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]

# 6. Run DESeq2
dds <- DESeq(dds)

# 7. Overall genotype effect (5xFAD_CLU vs 5xFAD)
res_all <- results(dds, contrast = c("genotype", "5xFAD_CLU", "5xFAD"))
write.csv(as.data.frame(res_all), file = "results/tables/DE_all.csv")

# 8. Sex-specific analyses
for (s in c("male", "female")) {
  sel <- coldata$sex == s
  dds_s <- DESeqDataSetFromMatrix(
    countData = counts[, sel],
    colData = coldata[sel, ],
    design = ~ genotype
  )
  dds_s <- DESeq(dds_s)
  res_s <- results(dds_s, contrast = c("genotype", "5xFAD_CLU", "5xFAD"))
  write.csv(as.data.frame(res_s), file = sprintf("results/tables/DE_%s.csv", s))
}

# 9. PCA plot (optional)
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c("genotype", "sex"), returnData = TRUE)
p <- ggplot(pca_data, aes(PC1, PC2, color = genotype, shape = sex)) +
  geom_point(size = 3) + theme_bw()
ggsave("results/figures/PCA.png", p)

# 10. Volcano of overall results (optional)
df <- as.data.frame(res_all)
df$negLogP <- -log10(df$padj)
v <- ggplot(df, aes(log2FoldChange, negLogP)) +
  geom_point(alpha = 0.4) + theme_bw()
ggsave("results/figures/volcano_all.png", v)
