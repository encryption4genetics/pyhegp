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

from hypothesis import given, settings, strategies as st
from hypothesis.extra.numpy import arrays, array_shapes
import numpy as np
from pytest import approx

from pyhegp.pyhegp import Stats, hegp_encrypt, hegp_decrypt, random_key, pool_stats, standardize, unstandardize

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
            and pooled_stats.mean == approx(np.mean(combined_pool, axis=0))
            and pooled_stats.std == approx(np.std(combined_pool, axis=0, ddof=1)))

def no_column_zero_standard_deviation(matrix):
    return not np.any(np.isclose(np.std(matrix, axis=0), 0))

@given(st.one_of(
    arrays("int32",
           array_shapes(min_dims=2, max_dims=2),
           elements=st.integers(min_value=0, max_value=2)),
    # The array above is the only realistic input, but we test more
    # kinds of inputs for good measure.
    arrays("int32",
           array_shapes(min_dims=2, max_dims=2),
           elements=st.integers(min_value=0, max_value=100)),
    arrays("float64",
           array_shapes(min_dims=2, max_dims=2),
           elements=st.floats(min_value=0, max_value=100)))
       # Reject matrices with zero standard deviation columns since
       # they trigger a division by zero.
       .filter(no_column_zero_standard_deviation))
def test_hegp_encryption_decryption_are_inverses(plaintext):
    mean = np.mean(plaintext, axis=0)
    standard_deviation = np.std(plaintext, axis=0)
    rng = np.random.default_rng()
    key = random_key(rng, len(plaintext))
    assert hegp_decrypt(hegp_encrypt(plaintext, mean, standard_deviation, key),
                        mean, standard_deviation, key) == approx(plaintext)

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
