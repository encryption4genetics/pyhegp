library(dplyr)
library(qqman)
library(readr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
    write("Usage: Rscript jwas-manhattan.r GWAS-FILE", stderr())
    quit(status=1)
}
gwas_file = args[1]

gwas = mutate(read_tsv(gwas_file), antilog=10^(-modelfrequency))
manhattan(gwas, chr="chromosome", bp="position", snp="marker_ID", p="antilog")
