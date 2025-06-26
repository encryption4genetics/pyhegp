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

import click
import numpy as np
from scipy.stats import special_ortho_group

def random_key(rng, n):
    return special_ortho_group.rvs(n, random_state=rng)

def standardize(genotype_matrix, maf):
    m, _ = genotype_matrix.shape
    return ((genotype_matrix - np.tile(maf, (m, 1)))
            @ np.diag(1 / np.sqrt(2 * maf * (1 - maf))))

def hegp_encrypt(plaintext, maf, key):
    return key @ plaintext
    # FIXME: Add standardization.
    # return key @ standardize(plaintext, maf)

def hegp_decrypt(ciphertext, key):
    return np.transpose(key) @ ciphertext

def read_genotype(genotype_file):
    return np.loadtxt(genotype_file, delimiter=",")
    # snps = genotype_file.readline().split(",")
    # return np.loadtxt(genotype_file, delimiter=",", skiprows=1, usecols=range(1, 1+len(snps)))

@click.group()
def main():
    pass

@main.command()
@click.argument("genotype-file", type=click.File("r"))
@click.argument("maf-file", type=click.File("r"))
@click.argument("key-path", type=click.Path())
@click.argument("ciphertext-path", type=click.Path())
def encrypt(genotype_file, maf_file, key_path, ciphertext_path):
    genotype = read_genotype(genotype_file)
    maf = np.loadtxt(maf_file)
    rng = np.random.default_rng()
    key = random_key(rng, len(genotype))
    encrypted_genotype = hegp_encrypt(genotype, maf, key)
    np.savetxt(key_path, key, delimiter=",", fmt="%f")
    np.savetxt(ciphertext_path, encrypted_genotype, delimiter=",", fmt="%f")

@main.command()
@click.argument("key-file", type=click.File("r"))
@click.argument("ciphertext-file", type=click.File("r"))
@click.argument("plaintext-path", type=click.Path())
def decrypt(key_file, ciphertext_file, plaintext_path):
    key = np.loadtxt(key_file, delimiter=",")
    ciphertext = np.loadtxt(ciphertext_file, delimiter=",")
    genotype = hegp_decrypt(ciphertext, key)
    np.savetxt(plaintext_path, genotype, delimiter=",", fmt="%f")

if __name__ == "__main__":
    main()
