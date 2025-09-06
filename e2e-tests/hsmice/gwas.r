### pyhegp --- Homomorphic encryption of genotypes and phenotypes
### Copyright © 2025 Arun Isaac <arunisaac@systemreboot.net>
###
### This file is part of pyhegp.
###
### pyhegp is free software: you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### pyhegp is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with pyhegp. If not, see <https://www.gnu.org/licenses/>.

library(mixed.model.gwas)
library(dplyr)
library(qqman)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

mmgwas <- function(genotype, phenotype, phenotype_name) {
    ## Wrangle phenotype to vector.
    phenotype_vector = phenotype %>% pull(phenotype_name)
    sample_ids = phenotype %>% pull("sample-id")
    names(phenotype_vector) = sample_ids

    ## Wrangle genotype to matrix.
    genotype_matrix = genotype %>%
        select(all_of(sample_ids)) %>%
        as.matrix %>%
        t
    colnames(genotype_matrix) = genotype %>%
        transmute(snp.id=str_c(chromosome, ":", position)) %>%
        pull(snp.id)

    ## Compute GWAS.
    gwas = mixed.model.gwas(phenotype_vector, genotype_matrix)
    return(tibble(chromosome=as.integer(str_split_i(colnames(gwas$pval), ":", 1)),
                  position=as.integer(str_split_i(colnames(gwas$pval), ":", 2)),
                  snp.id=colnames(gwas$pval),
                  p=gwas$pval[1,]) %>%
           drop_na())
}

run_gwas_experiment <- function(genotype_file, phenotype_file, pvalues_file) {
    gwas = mmgwas(read_tsv(genotype_file),
                  read_tsv(phenotype_file) %>%
                  select("sample-id", "Biochem.ALP.resid"),
                  "Biochem.ALP.resid")
    write_tsv(gwas %>% select(chromosome, position, p),
              pvalues_file)
    manhattan(gwas, chr="chromosome", bp="position", snp="snp.id", p="p")
}

## Process command-line arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    write("Usage: Rscript hsmice.r GENOTYPE-FILE PHENOTYPE-FILE PVALUES-FILE", stderr())
    quit(status=1)
}
genotype_file = args[1]
phenotype_file = args[2]
pvalues_file = args[3]

## Run GWAS.
run_gwas_experiment(genotype_file, phenotype_file, pvalues_file)
