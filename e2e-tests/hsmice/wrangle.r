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

library(dplyr)
library(genio)
library(purrr)
library(readr)
library(tibble)
library(tidyr)

mean_sans_rm = partial(mean, na.rm=TRUE)

## Process command-line arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    write("Usage: Rscript hsmice.r INPUT_DIRECTORY OUTPUT_DIRECTORY", stderr())
    quit(status=1)
}
input_data_directory = args[1]
output_directory = args[2]

## Read phenotype data.
phenotype = read.table(file.path(input_data_directory, "HSmice.phe"),
                       header=TRUE) %>%
    select(id, sex, Anx.resid, BurrowedPelletWeight.resid, Context.resid,
           End.Weight.resid, Explore.resid, FN.preWeight.resid, Freeze.resid,
           Glucose.Weight.resid, OFT.Boli.resid, OFT.CenterTime.resid,
           OFT.Latency.resid, OFT.TotalActivity.resid, Start.Weight.resid,
           Weight.GrowthSlope.resid, Biochem.ALP.resid) %>%
    drop_na() %>%
    rename("sample-id"="id")
sample_ids = phenotype %>% pull("sample-id")

## Read genotype data, replace NAs by mean, and subset.
p = read_plink(file.path(input_data_directory, "HSmice.bed"))
genotype_matrix = t(apply(p$X, 1, function(x) ifelse(is.na(x), mean_sans_rm(x), x))) %>%
    data.frame() %>%
    select(all_of(sample_ids))
genotype = bind_cols(data.frame(chromosome=p$bim$chr,
                                position=p$bim$pos,
                                reference=p$bim$ref),
                     genotype_matrix)

## Write whole phenotype and the pieces.
midpoint = length(sample_ids) %/% 2
sample_ids1 = sample_ids[1:midpoint]
sample_ids2 = sample_ids[(midpoint+1):length(sample_ids)]
write_tsv(phenotype,
          file.path(output_directory, "phenotype.tsv"))
write_tsv(phenotype %>%
          inner_join(tibble("sample-id"=sample_ids1)),
          file.path(output_directory, "phenotype1.tsv"))
write_tsv(phenotype %>%
          inner_join(tibble("sample-id"=sample_ids2)),
          file.path(output_directory, "phenotype2.tsv"))

## Write whole genotype and the pieces.
write_tsv(genotype,
          file.path(output_directory, "genotype.tsv"))
write_tsv(genotype %>%
          select(chromosome, position, reference, all_of(sample_ids1)),
          file.path(output_directory, "genotype1.tsv"))
write_tsv(genotype %>%
          select(chromosome, position, reference, all_of(sample_ids2)),
          file.path(output_directory, "genotype2.tsv"))
