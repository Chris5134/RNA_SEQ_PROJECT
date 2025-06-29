#!/usr/bin/env Rscript
# scripts/extract_metadata.R

# 1. 读入所有行
lines <- readLines("data/raw/GSE294133_series_matrix.txt")

# 2. 找到样本表头行，提取 sample IDs
header_line <- grep('^"ID_REF"', lines, value=TRUE)
samples <- strsplit(sub('^"ID_REF"\\t', '', header_line), '\\t')[[1]]
samples <- gsub('"', '', samples)

# 3. 找到 genotype 和 sex 那两行
gen_line <- grep('^!Sample_characteristics_ch1.*[Gg]enotype', lines, value=TRUE)[1]
sex_line <- grep('^!Sample_characteristics_ch1.*[Ss]ex', lines, value=TRUE)[1]

# 4. 解析出字段值
genotypes <- strsplit(sub('^!Sample_characteristics_ch1\\t', '', gen_line), '\\t')[[1]]
genotypes <- sub('^"genotype: ?', '', genotypes)
genotypes <- gsub('"', '', genotypes)

sexes <- strsplit(sub('^!Sample_characteristics_ch1\\t', '', sex_line), '\\t')[[1]]
sexes <- sub('^"[Ss]ex: ?', '', sexes)
sexes <- gsub('"', '', sexes)

# 5. 写出 TSV
meta <- data.frame(
  sample_id = samples,
  genotype  = genotypes,
  sex       = sexes,
  stringsAsFactors = FALSE
)
write.table(meta,
            file = "data/raw/metadata.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
cat("Wrote metadata for", nrow(meta), "samples to data/raw/metadata.tsv\n")
