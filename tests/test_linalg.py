### pyhegp --- Homomorphic encryption of genotypes and phenotypes
### Copyright © 2026 Arun Isaac <arunisaac@systemreboot.net>
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

from hypothesis import given, strategies as st
from hypothesis.extra.numpy import arrays
import numpy as np
from pytest import approx

from pyhegp.linalg import BlockDiagonalMatrix

@st.composite
def block_diagonal_matrices(draw, max_block_size=10, max_number_of_blocks=None):
    return BlockDiagonalMatrix(
        [draw(arrays("float64", (n, n),
                     elements=st.floats(min_value=-10,
                                        max_value=10,
                                        allow_nan=False,
                                        allow_infinity=False)))
         for n in draw(st.lists(st.integers(min_value=1,
                                            max_value=max_block_size),
                                min_size=1,
                                max_size=max_number_of_blocks))])

@given(block_diagonal_matrices(max_number_of_blocks=10))
def test_block_diagonal_matrix_transpose(block_diagonal_matrix):
    assert (np.transpose(block_diagonal_matrix).__array__()
            == approx(np.transpose(block_diagonal_matrix.__array__())))

@st.composite
def block_diagonal_matrix_product_multiplicands(draw):
    block_diagonal_matrix = draw(block_diagonal_matrices())
    return (block_diagonal_matrix,
            draw(arrays("float64",
                        # Either a vector or a matrix
                        (len(block_diagonal_matrix),
                         *draw(st.one_of(
                             st.just(()),
                             st.builds(lambda x: (x,),
                                       st.integers(min_value=1,
                                                   max_value=100))))),
                        elements=st.floats(min_value=-10,
                                           max_value=10,
                                           allow_nan=False,
                                           allow_infinity=False))))

@given(block_diagonal_matrix_product_multiplicands())
def test_block_diagonal_matrix_product(multiplicands):
    block_diagonal_matrix, multiplier = multiplicands
    assert ((block_diagonal_matrix @ multiplier)
            == approx(block_diagonal_matrix.__array__() @ multiplier))
