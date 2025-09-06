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

from collections import namedtuple
from functools import reduce
from pathlib import Path
import sys

import click
import numpy as np
import pandas as pd
from scipy.stats import special_ortho_group

from pyhegp.serialization import Summary, read_summary, write_summary, read_genotype, read_phenotype, write_genotype, write_phenotype, write_key, is_genotype_metadata_column

Stats = namedtuple("Stats", "n mean std")

def random_key(rng, n):
    return special_ortho_group.rvs(n, random_state=rng)

def standardize(matrix, mean, standard_deviation):
    m, _ = matrix.shape
    return ((matrix - np.tile(mean, (m, 1)))
            @ np.diag(1 / standard_deviation))

def unstandardize(matrix, mean, standard_deviation):
    m, _ = matrix.shape
    return ((matrix @ np.diag(standard_deviation))
            + np.tile(mean, (m, 1)))

def hegp_encrypt(plaintext, key):
    return key @ plaintext

def hegp_decrypt(ciphertext, key):
    return np.transpose(key) @ ciphertext

def genotype_summary(genotype):
    matrix = genotype.drop(columns=["chromosome", "position", "reference"]).to_numpy()
    return Summary(genotype.shape[0],
                   pd.DataFrame({"chromosome": genotype.chromosome,
                                 "position": genotype.position,
                                 "reference": genotype.reference,
                                 "mean": np.mean(matrix, axis=1),
                                 "std": np.std(matrix, axis=1)}))

def pool_stats(list_of_stats):
    sums = [stats.n*stats.mean for stats in list_of_stats]
    sums_of_squares = [(stats.n-1)*stats.std**2 + stats.n*stats.mean**2
                       for stats in list_of_stats]
    n = np.sum([stats.n for stats in list_of_stats])
    mean = np.sum(sums, axis=0) / n
    std = np.sqrt((np.sum(sums_of_squares, axis=0) - n*mean**2)
                  / (n - 1))
    return Stats(n, mean, std)

def pool_summaries(summaries):
    def pool_summaries2(summary1, summary2):
        # Drop any SNPs that are not in both summaries.
        data = pd.merge(summary1.data.rename(columns={"mean": "mean1",
                                                      "std": "std1"}),
                        summary2.data.rename(columns={"mean": "mean2",
                                                      "std": "std2"}),
                        how="inner",
                        on=("chromosome", "position", "reference"))
        pooled_stats = pool_stats([Stats(summary1.n,
                                         data.mean1.to_numpy(),
                                         data.std2.to_numpy()),
                                   Stats(summary2.n,
                                         data.mean2.to_numpy(),
                                         data.std2.to_numpy())])
        return Summary(pooled_stats.n,
                       pd.concat((data[["chromosome", "position", "reference"]],
                                  pd.DataFrame({"mean": pooled_stats.mean,
                                                "std": pooled_stats.std})),
                                 axis="columns"))
    pooled_summary = reduce(pool_summaries2, summaries)
    return Summary(pooled_summary.n,
                   pooled_summary.data.drop(columns=["reference"]))

def encrypt_genotype(genotype, key, summary):
    # Drop SNPs that have a zero standard deviation. Such SNPs have no
    # discriminatory power in the analysis and mess with our
    # standardization by causing a division by zero.
    summary = summary._replace(
        data=summary.data[~np.isclose(summary.data["std"], 0)])
    # Drop any SNPs that are not in both genotype and summary.
    common_genotype = pd.merge(genotype,
                               summary.data[["chromosome", "position"]],
                               on=("chromosome", "position"))
    sample_names = (common_genotype.drop(
        columns=["chromosome", "position", "reference"]).columns)
    genotype_matrix = common_genotype[sample_names].to_numpy().T
    encrypted_genotype_matrix = hegp_encrypt(standardize(
        genotype_matrix,
        summary.data["mean"].to_numpy(),
        summary.data["std"].to_numpy()),
                                             key)
    return pd.concat((common_genotype[["chromosome", "position"]],
                      pd.DataFrame(encrypted_genotype_matrix.T,
                                   columns=sample_names)),
                     axis="columns")

def encrypt_phenotype(phenotype, key):
    phenotype_matrix = phenotype.drop(columns=["sample-id"])
    sample_names = phenotype_matrix.columns
    return pd.concat((phenotype["sample-id"],
                      pd.DataFrame(hegp_encrypt(phenotype_matrix.to_numpy(),
                                                key),
                                   columns=sample_names)),
                     axis="columns")

def cat_genotype(genotypes):
    def cat2(df1, df2):
        return pd.merge(df1, df2,
                        how="inner",
                        on=list(filter(is_genotype_metadata_column,
                                       df1.columns)))
    match genotypes:
        # If there are no input data frames, return an empty data
        # frame with the chromosome and position columns.
        case []:
            return pd.DataFrame({"chromosome": pd.Series(dtype="str"),
                                 "position": pd.Series(dtype="int")})
        case _:
            return reduce(cat2, genotypes)

def cat_phenotype(phenotypes):
    match phenotypes:
        # If there are no input data frames, return an empty data
        # frame with the sample-id column.
        case []:
            return pd.DataFrame(data={"sample-id": []},
                                dtype="str")
        case _:
            return pd.concat(phenotypes)

@click.group()
def main():
    pass

@main.command("summary")
@click.argument("genotype-file", type=click.File("r"))
@click.option("--output", "-o", "summary_file",
              type=click.File("wb"),
              default="-",
              help="output file")
def summary_command(genotype_file, summary_file):
    write_summary(summary_file,
                  genotype_summary(read_genotype(genotype_file)))

@main.command("pool")
@click.option("--output", "-o", "pooled_summary_file",
              type=click.File("wb"),
              default="-",
              help="output file")
@click.argument("summary-files", type=click.File("rb"), nargs=-1)
def pool_command(pooled_summary_file, summary_files):
    summaries = [read_summary(file) for file in summary_files]
    pooled_summary = pool_summaries(summaries)
    max_snps = max(len(summary.data) for summary in summaries)
    if len(pooled_summary.data) < max_snps:
        dropped_snps = max_snps - len(pooled_summary.data)
        print(f"Dropped {dropped_snps} SNP(s)")
    write_summary(pooled_summary_file, pooled_summary)

@main.command("encrypt")
@click.argument("genotype-file", type=click.File("r"))
@click.argument("phenotype-file", type=click.File("r"), required=False)
@click.option("--summary", "-s", "summary_file", type=click.File("rb"),
              help="Summary statistics file")
@click.option("--key", "-k", "key_file", type=click.File("w"),
              help="Output key")
@click.option("--force", "-f", is_flag=True,
              help="Overwrite output files even if they exist")
def encrypt_command(genotype_file, phenotype_file, summary_file, key_file, force):
    def write_ciphertext(plaintext_path, writer):
        ciphertext_path = Path(plaintext_path + ".hegp")
        if ciphertext_path.exists() and not force:
            print(f"Output file {ciphertext_path} exists, cannot overwrite.")
            sys.exit(1)
        with ciphertext_path.open("w") as ciphertext_file:
            writer(ciphertext_file)

    genotype = read_genotype(genotype_file)
    if summary_file:
        summary = read_summary(summary_file)
    else:
        summary = genotype_summary(genotype)
    key = random_key(np.random.default_rng(),
                     len(genotype
                         .drop(columns=["chromosome", "position", "reference"])
                         .columns))
    if key_file:
        write_key(key_file, key)

    encrypted_genotype = encrypt_genotype(genotype, key, summary)
    if len(encrypted_genotype) < len(genotype):
        dropped_snps = len(genotype) - len(encrypted_genotype)
        print(f"Dropped {dropped_snps} SNP(s)")
    write_ciphertext(genotype_file.name,
                     lambda file: write_genotype(file, encrypted_genotype))

    if phenotype_file:
        write_ciphertext(phenotype_file.name,
                         lambda file:
                         write_phenotype(file, encrypt_phenotype(
                             read_phenotype(phenotype_file),
                             key)))

@main.command("cat-genotype")
@click.option("--output", "-o", "output_file",
              type=click.File("wb"),
              default="-",
              help="output file")
@click.argument("ciphertext-files", type=click.File("rb"), nargs=-1)
def cat_genotype_command(output_file, ciphertext_files):
    write_genotype(output_file,
                   cat_genotype([read_genotype(file)
                                 for file in ciphertext_files]))

@main.command("cat-phenotype")
@click.option("--output", "-o", "output_file",
              type=click.File("wb"),
              default="-",
              help="output file")
@click.argument("ciphertext-files", type=click.File("rb"), nargs=-1)
def cat_phenotype_command(output_file, ciphertext_files):
    write_phenotype(output_file,
                    cat_phenotype([read_phenotype(file)
                                   for file in ciphertext_files]))

if __name__ == "__main__":
    main()
