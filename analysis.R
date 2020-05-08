library(GEOquery)
library(dplyr)
library(tibble)
library(Seurat)

# Tab-delimited text file of log2 (incremented by 1 before log2) TPM RSEM expression values
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70240/suppl/GSE70240_Gmp.txt.gz", "data/gse70240.txt.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70236/suppl/GSE70236_Cmp.txt.gz", "data/gse70236.txt.gz")

gse70236_meta <- getGEO("GSE70236", GSEMatrix = TRUE)
gse70240_meta <- getGEO("GSE70240", GSEMatrix = TRUE)

gse70236_meta <- pData(phenoData(gse70236_meta[[1]]))
gse70240_meta <- pData(phenoData(gse70240_meta[[1]]))

# 确认两个gse基因完全相同
identical(gse70236$uid, gse70240$uid)

# 确认丰度列表和meta数据中sample一致

# 第一个字符大小，剩余字符小写 
only_first_upper <- function(x){
  paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))), sep = "")
}

sample1 <- strsplit(gse70236_meta$title, " ") %>% 
  purrr::map_chr(~ c(only_first_upper(.x[1]), .x[3]) %>% paste(collapse = "."))
identical(sample1, names(gse70236)[-1])
gse70236_meta$id <- sample1
gse70236_meta <- dplyr::select(gse70236_meta, id, everything())

sample2 <- strsplit(gse70240_meta$title, " ") %>% 
  purrr::map_chr(~ c(only_first_upper(.x[1]), .x[3]) %>% paste(collapse = "."))
identical(sample2, names(gse70240)[-1])
gse70240_meta$id <- sample2
gse70240_meta <- dplyr::select(gse70240_meta, id, everything())

# 合并两个gse
sample_abd <- dplyr::bind_cols(gse70236, gse70240[-1])
  
sample_meta <- dplyr::bind_rows(gse70236_meta, gse70240_meta)

readr::write_csv(sample_abd, "log2_abd.csv")
readr::write_csv(sample_meta, "sample_meta.csv")
