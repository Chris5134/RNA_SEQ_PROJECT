
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("DESeq2")
}
library(DESeq2)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

counts <- as.matrix(
  read.csv("data/raw/counts.csv", row.names = 1, check.names = FALSE)
)
coldata <- read.csv(
  "data/raw/metadata.tsv",
  sep = "\t", row.names = 1, check.names = FALSE
)

colnames(counts) <- rownames(coldata)


coldata$genotype <- ifelse(
  coldata$genotype == "5xFADHEMI:CLUh2kbKIHO",
  "5xFAD_CLU",
  "5xFAD"
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ sex + genotype
)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]

dds <- DESeq(dds)

res_all <- results(dds, contrast = c("genotype", "5xFAD_CLU", "5xFAD"))
write.csv(as.data.frame(res_all), file = "results/tables/DE_all.csv")

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

df <- as.data.frame(res_all)
df$negLogP <- -log10(df$padj)
v <- ggplot(df, aes(log2FoldChange, negLogP)) +
  geom_point(alpha = 0.4) + theme_bw()
ggsave("results/figures/volcano_all.png", v)
