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

from pyhegp.pyhegp import hegp_encrypt, hegp_decrypt, random_key

@given(st.one_of(
    arrays("int32",
           array_shapes(min_dims=2, max_dims=2, min_side=2, max_side=100),
           elements=st.integers(min_value=0, max_value=2)),
    # The array above is the only realistic input, but we test more
    # kinds of inputs for good measure.
    arrays("int32",
           array_shapes(min_dims=2, max_dims=2, min_side=2, max_side=100),
           elements=st.integers(min_value=0, max_value=100)),
    arrays("float64",
           array_shapes(min_dims=2, max_dims=2, min_side=2, max_side=100),
           elements=st.floats(min_value=0, max_value=100)))
)
@settings(deadline=None)
def test_hegp_encryption_decryption_are_inverses(plaintext):
    rng = np.random.default_rng()
    key = random_key(rng, len(plaintext))
    # FIXME: We don't use maf at the moment.
    maf = None
    assert hegp_decrypt(hegp_encrypt(plaintext, maf, key), key) == approx(plaintext)
