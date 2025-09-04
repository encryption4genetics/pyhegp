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

from itertools import pairwise
import math
from pathlib import Path
import shutil

from click.testing import CliRunner
from hypothesis import given, strategies as st
from hypothesis.extra.numpy import arrays, array_shapes
import numpy as np
import pandas as pd
from pytest import approx, mark

from pyhegp.pyhegp import Stats, main, hegp_encrypt, hegp_decrypt, random_key, pool_stats, standardize, unstandardize, genotype_summary, encrypt_genotype, encrypt_phenotype, cat_genotype, cat_phenotype
from pyhegp.serialization import Summary, read_summary, read_genotype, is_genotype_metadata_column, is_phenotype_metadata_column
from pyhegp.utils import negate

from helpers.strategies import *

@given(st.lists(st.lists(arrays("float64",
                                st.shared(array_shapes(min_dims=1, max_dims=1),
                                          key="pool-vector-length"),
                                elements=st.floats(min_value=-100, max_value=100)),
                         min_size=2),
                min_size=1))
def test_pool_stats(pools):
    combined_pool = sum(pools, [])
    pooled_stats = pool_stats([Stats(len(pool),
                                     np.mean(pool, axis=0),
                                     np.std(pool, axis=0, ddof=1))
                               for pool in pools])
    assert (pooled_stats.n == len(combined_pool)
            and pooled_stats.mean == approx(np.mean(combined_pool, axis=0),
                                            rel=1e-6)
            and pooled_stats.std == approx(np.std(combined_pool, axis=0, ddof=1),
                                           rel=1e-6))

def test_encrypt_command(tmp_path):
    shutil.copy("test-data/encrypt-test-genotype.tsv", tmp_path)
    ciphertext = tmp_path / "encrypt-test-genotype.tsv.hegp"
    result = CliRunner().invoke(main, ["encrypt",
                                       "-s", "test-data/encrypt-test-summary",
                                       str(tmp_path / "encrypt-test-genotype.tsv")])
    assert result.exit_code == 0
    assert ciphertext.exists()
    assert "Dropped 1 SNP(s)" in result.output
    with ciphertext.open("rb") as genotype_file:
        encrypted_genotype = read_genotype(genotype_file)
    # TODO: Properly compare encrypted genotype data frame with
    # expected output once it is possible to specify the key.
    assert len(encrypted_genotype) == 3

def no_column_zero_standard_deviation(matrix):
    return not np.any(np.isclose(np.std(matrix, axis=0), 0))

@given(st.one_of(
    arrays("int32",
           array_shapes(min_dims=2, max_dims=2, min_side=2),
           elements=st.integers(min_value=0, max_value=2)),
    # The array above is the only realistic input, but we test more
    # kinds of inputs for good measure.
    arrays("int32",
           array_shapes(min_dims=2, max_dims=2, min_side=2),
           elements=st.integers(min_value=0, max_value=100)),
    arrays("float64",
           array_shapes(min_dims=2, max_dims=2, min_side=2),
           elements=st.floats(min_value=0, max_value=100))))
def test_hegp_encryption_decryption_are_inverses(plaintext):
    rng = np.random.default_rng()
    key = random_key(rng, len(plaintext))
    assert hegp_decrypt(hegp_encrypt(plaintext, key), key) == approx(plaintext)

@given(arrays("float64",
              array_shapes(min_dims=2, max_dims=2),
              elements=st.floats(min_value=0, max_value=100))
       # Reject matrices with zero standard deviation columns since
       # they trigger a division by zero.
       .filter(no_column_zero_standard_deviation))
def test_standardize_unstandardize_are_inverses(matrix):
    mean = np.mean(matrix, axis=0)
    standard_deviation = np.std(matrix, axis=0)
    assert unstandardize(standardize(matrix, mean, standard_deviation),
                         mean, standard_deviation) == approx(matrix)

def square_matrices(order, elements=None):
    def generate(draw):
        n = draw(order)
        return draw(arrays("float64", (n, n), elements=elements))
    return st.composite(generate)

def is_singular(matrix):
    # We want to avoid nearly singular matrices as well. Hence, we set
    # a looser absolute tolerance.
    return math.isclose(np.linalg.det(matrix), 0, abs_tol=1e-6)

@given(square_matrices(st.shared(st.integers(min_value=2, max_value=7),
                                 key="n"),
                       elements=st.floats(min_value=0, max_value=10))()
       .filter(negate(is_singular)),
       arrays("float64",
              st.shared(st.integers(min_value=2, max_value=7),
                        key="n"),
              elements=st.floats(min_value=0, max_value=10)))
def test_conservation_of_solutions(genotype, phenotype):
    rng = np.random.default_rng()
    key = random_key(rng, len(genotype))
    assert (approx(np.linalg.solve(genotype, phenotype),
                   abs=1e-6, rel=1e-6)
            == np.linalg.solve(hegp_encrypt(genotype, key),
                               hegp_encrypt(phenotype, key)))

@mark.xfail
@given(genotype_frames(st.shared(st.integers(min_value=2, max_value=10),
                                 key="number-of-samples"),
                       reference_present=st.just(True)),
       keys(st.shared(st.integers(min_value=2, max_value=10),
                      key="number-of-samples")))
def test_encrypt_genotype_does_not_produce_na(genotype, key):
    assert not encrypt_genotype(genotype,
                                key,
                                genotype_summary(genotype)).isna().any(axis=None)

@given(phenotype_frames(st.shared(st.integers(min_value=2, max_value=10),
                                  key="number-of-samples")),
       keys(st.shared(st.integers(min_value=2, max_value=10),
                      key="number-of-samples")))
def test_encrypt_phenotype_does_not_produce_na(phenotype, key):
    assert not encrypt_phenotype(phenotype, key).isna().any(axis=None)

def test_pool_command(tmp_path):
    columns = ["chromosome", "position", "reference", "mean", "std"]
    complete_summary = tmp_path / "complete-summary"
    result = CliRunner().invoke(main, ["pool",
                                       "-o", complete_summary,
                                       "test-data/pool-test-summary1",
                                       "test-data/pool-test-summary2"],
                                catch_exceptions=True)
    assert result.exit_code == 0
    assert complete_summary.exists()
    assert "Dropped 2 SNP(s)" in result.output
    with complete_summary.open("rb") as summary_file:
        pooled_summary = read_summary(summary_file)
    with open("test-data/pool-test-complete-summary", "rb") as summary_file:
        expected_pooled_summary = read_summary(summary_file)
    pd.testing.assert_frame_equal(pooled_summary.data,
                                  expected_pooled_summary.data)
    assert pooled_summary.n == expected_pooled_summary.n

def split_data_frame(draw, df, metadata_columns):
    metadata = df[metadata_columns]
    data_columns = [column
                    for column in df.columns
                    if column not in metadata_columns]
    data = df[data_columns]
    split_points = sorted(draw(st.lists(st.integers(min_value=0,
                                                    max_value=len(data_columns)),
                                        min_size=0,
                                        ## Something reasonably small.
                                        max_size=len(data_columns))))
    return [pd.concat((metadata, data[data_columns[start:end]]),
                      axis="columns")
            for start, end
            in pairwise([0] + split_points + [len(data_columns)])]

@st.composite
def catenable_genotype_frames(draw):
    genotype = draw(genotype_frames())
    return ([genotype]
            + split_data_frame(draw,
                               genotype,
                               list(filter(is_genotype_metadata_column,
                                           genotype.columns))))

@given(catenable_genotype_frames())
def test_cat_genotype(genotypes):
    complete_genotype, *split_genotypes = genotypes
    pd.testing.assert_frame_equal(complete_genotype,
                                  cat_genotype(split_genotypes))

@st.composite
def catenable_phenotype_frames(draw):
    phenotype = draw(phenotype_frames())
    return ([phenotype]
            + split_data_frame(draw,
                               phenotype,
                               list(filter(is_phenotype_metadata_column,
                                           phenotype.columns))))

@given(catenable_phenotype_frames())
def test_cat_phenotype(phenotypes):
    complete_phenotype, *split_phenotypes = phenotypes
    pd.testing.assert_frame_equal(complete_phenotype,
                                  cat_phenotype(split_phenotypes))

def test_simple_workflow(tmp_path):
    shutil.copy(f"test-data/genotype.tsv", tmp_path)
    ciphertext = tmp_path / "genotype.tsv.hegp"
    result = CliRunner().invoke(main,
                                ["encrypt", str(tmp_path / "genotype.tsv")])
    assert result.exit_code == 0
    assert ciphertext.exists()

def test_joint_workflow(tmp_path):
    runner = CliRunner()
    for i in range(4):
        shutil.copy(f"test-data/genotype{i}.tsv", tmp_path)
        summary = tmp_path / f"summary{i}"
        result = runner.invoke(
            main, ["summary", str(tmp_path / f"genotype{i}.tsv"),
                   "-o", summary])
        assert result.exit_code == 0
        assert summary.exists()
    complete_summary = tmp_path / "complete-summary"
    result = runner.invoke(
        main, ["pool",
               "-o", complete_summary,
               *(str(tmp_path / f"summary{i}") for i in range(4))])
    assert result.exit_code == 0
    assert complete_summary.exists()
    for i in range(4):
        ciphertext = tmp_path / f"genotype{i}.tsv.hegp"
        result = runner.invoke(
            main, ["encrypt",
                   "-s", complete_summary,
                   str(tmp_path / f"genotype{i}.tsv")])
        assert result.exit_code == 0
        assert ciphertext.exists()
    complete_ciphertext = tmp_path / "complete-genotype.tsv.hegp"
    result = runner.invoke(
        main, ["cat-genotype",
               "-o", complete_ciphertext,
               *(str(tmp_path / f"genotype{i}.tsv.hegp") for i in range(4))])
    assert result.exit_code == 0
    assert complete_ciphertext.exists()
