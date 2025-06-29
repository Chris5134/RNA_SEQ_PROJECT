
lines <- readLines("data/raw/GSE294133_series_matrix.txt")

header_line <- grep('^"ID_REF"', lines, value=TRUE)
samples <- strsplit(sub('^"ID_REF"\\t', '', header_line), '\\t')[[1]]
samples <- gsub('"', '', samples)

gen_line <- grep('^!Sample_characteristics_ch1.*[Gg]enotype', lines, value=TRUE)[1]
sex_line <- grep('^!Sample_characteristics_ch1.*[Ss]ex', lines, value=TRUE)[1]

genotypes <- strsplit(sub('^!Sample_characteristics_ch1\\t', '', gen_line), '\\t')[[1]]
genotypes <- sub('^"genotype: ?', '', genotypes)
genotypes <- gsub('"', '', genotypes)

sexes <- strsplit(sub('^!Sample_characteristics_ch1\\t', '', sex_line), '\\t')[[1]]
sexes <- sub('^"[Ss]ex: ?', '', sexes)
sexes <- gsub('"', '', sexes)

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
