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

import click
import numpy as np
from scipy.stats import special_ortho_group

from pyhegp.serialization import Summary, read_summary, write_summary, read_genotype

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

def hegp_encrypt(plaintext, mean, standard_deviation, key):
    return key @ standardize(plaintext, mean, standard_deviation)

def hegp_decrypt(ciphertext, mean, standard_deviation, key):
    return unstandardize(np.transpose(key) @ ciphertext,
                         mean, standard_deviation)

def pool_stats(list_of_stats):
    sums = [stats.n*stats.mean for stats in list_of_stats]
    sums_of_squares = [(stats.n-1)*stats.std**2 + stats.n*stats.mean**2
                       for stats in list_of_stats]
    n = np.sum([stats.n for stats in list_of_stats])
    mean = np.sum(sums, axis=0) / n
    std = np.sqrt((np.sum(sums_of_squares, axis=0) - n*mean**2)
                  / (n - 1))
    return Stats(n, mean, std)

@click.group()
def main():
    pass

@main.command()
@click.argument("genotype-file", type=click.File("r"))
@click.option("--output", "-o", "summary_file",
              type=click.File("wb"),
              default="-",
              help="output file")
def summary(genotype_file, summary_file):
    genotype = read_genotype(genotype_file)
    write_summary(summary_file,
                  Summary(genotype.shape[0],
                          np.mean(genotype, axis=0),
                          np.std(genotype, axis=0)))

@main.command()
@click.option("--output", "-o", "pooled_summary_file",
              type=click.File("wb"),
              default="-",
              help="output file")
@click.argument("summary-files", type=click.File("rb"), nargs=-1)
def pool(pooled_summary_file, summary_files):
    summaries = [read_summary(file) for file in summary_files]
    pooled_stats = pool_stats([Stats(summary.n, summary.mean, summary.std)
                               for summary in summaries])
    write_summary(pooled_summary_file,
                  Summary(pooled_stats.n,
                          pooled_stats.mean,
                          pooled_stats.std))

@main.command()
@click.argument("genotype-file", type=click.File("r"))
@click.option("--summary", "-s", "summary_file", type=click.File("rb"),
              help="Summary statistics file",
              required=True)
@click.option("--key", "-k", "key_path", type=click.Path(),
              help="Output key",
              required=True)
@click.option("--output", "-o", "ciphertext_path", type=click.Path(),
              default="-",
              help="Output ciphertext")
def encrypt(genotype_file, summary_file, key_path, ciphertext_path):
    genotype = read_genotype(genotype_file)
    summary = read_summary(summary_file)
    rng = np.random.default_rng()
    key = random_key(rng, len(genotype))
    encrypted_genotype = hegp_encrypt(genotype, summary.mean,
                                      summary.std, key)
    np.savetxt(key_path, key, delimiter=",", fmt="%f")
    np.savetxt(ciphertext_path, encrypted_genotype, delimiter=",", fmt="%f")

if __name__ == "__main__":
    main()
