#!/usr/bin/env Rscript
# scripts/sex_DE_analysis.R

# Library setup
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("DESeq2")
}
library(DESeq2)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
library(clusterProfiler)

if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db")
}
library(org.Mm.eg.db)

# Data import
counts <- as.matrix(
  read.csv("data/raw/counts.csv", row.names = 1, check.names = FALSE)
)
coldata <- read.csv(
  "data/raw/metadata.tsv", sep = "\t", row.names = 1, check.names = FALSE
)

# Align sample names
colnames(counts) <- rownames(coldata)
stopifnot(all(colnames(counts) == rownames(coldata)))

#Metadata cleanup
coldata$genotype <- factor(
  ifelse(coldata$genotype == "5xFADHEMI:CLUh2kbKIHO", "5xFAD_CLU", "5xFAD"),
  levels = c("5xFAD", "5xFAD_CLU")
)
coldata$sex <- factor(coldata$sex, levels = c("male", "female"))

#Build and filter DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ sex + genotype
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

#Run DESeq2
dds <- DESeq(dds)

#Extract results
res_all <- results(dds, contrast = c("genotype", "5xFAD_CLU", "5xFAD"))
write.csv(as.data.frame(res_all), "results/tables/DE_all.csv")

for (s in levels(coldata$sex)) {
  sel   <- coldata$sex == s
  dds_s <- DESeqDataSetFromMatrix(
    countData = counts[, sel],
    colData   = coldata[sel, ],
    design    = ~ genotype
  )
  dds_s <- DESeq(dds_s)
  res_s <- results(dds_s, contrast = c("genotype", "5xFAD_CLU", "5xFAD"))
  write.csv(as.data.frame(res_s), sprintf("results/tables/DE_%s.csv", s))
}

#Volcano plots with cutoffs
fc_cutoff <- log2(1.2)
p_cutoff  <- 0.05

make_volcano <- function(res, label) {
  df <- as.data.frame(res)
  df <- subset(df, !is.na(padj))
  df$negLogP <- -log10(df$padj)
  df$signif  <- abs(df$log2FoldChange) >= fc_cutoff & df$padj <= p_cutoff
  p <- ggplot(df, aes(log2FoldChange, negLogP)) +
       geom_point(aes(color = signif), alpha = 0.6) +
       scale_color_manual(values = c("grey70","red")) +
       theme_bw() +
       labs(title = paste("Volcano:", label), x = "log2 Fold Change", y = "-log10 adj. p-value") +
       theme(legend.position = "none")
  ggsave(sprintf("results/figures/volcano_%s.png", label), p, width = 7, height = 7)
  sig <- subset(df, signif)
  write.csv(sig, file = sprintf("results/tables/DE_%s_sig.csv", label), quote = FALSE)
}

make_volcano(res_all, "all")
res_male   <- read.csv("results/tables/DE_male.csv",   row.names = 1)
res_female <- read.csv("results/tables/DE_female.csv", row.names = 1)
make_volcano(res_male,   "male")
make_volcano(res_female, "female")

#GO enrichment for significant genes
# Overall
sig_all    <- read.csv("results/tables/DE_all_sig.csv",    row.names = 1)
genes_all  <- rownames(sig_all)
# Assume Ensembl IDs, strip version
genes_all  <- gsub("\\.\\d+$", "", genes_all)
entrez_all <- bitr(genes_all, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID

# Males
sig_male    <- read.csv("results/tables/DE_male_sig.csv",    row.names = 1)
genes_male  <- rownames(sig_male)
genes_male  <- gsub("\\.\\d+$", "", genes_male)
entrez_male <- bitr(genes_male, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID

# Females
sig_female    <- read.csv("results/tables/DE_female_sig.csv",    row.names = 1)
genes_female  <- rownames(sig_female)
genes_female  <- gsub("\\.\\d+$", "", genes_female)
entrez_female <- bitr(genes_female, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID

# Enrichment
ego_all    <- enrichGO(entrez_all,    OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
ego_male   <- enrichGO(entrez_male,   OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
ego_female <- enrichGO(entrez_female, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

# Plot top 10
dotplot(ego_all,    showCategory = 10) + ggtitle("GO BP – All samples")
ggsave("results/figures/GO_all.png", width = 8, height = 6)

dotplot(ego_male,   showCategory = 10) + ggtitle("GO BP – Male samples")
ggsave("results/figures/GO_male.png", width = 8, height = 6)

dotplot(ego_female, showCategory = 10) + ggtitle("GO BP – Female samples")
ggsave("results/figures/GO_female.png", width = 8, height = 6)
