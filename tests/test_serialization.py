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

import tempfile

from hypothesis import given, strategies as st
from hypothesis.extra.numpy import arrays, array_shapes
from pytest import approx

from pyhegp.serialization import Summary, read_summary, write_summary, read_summary_headers

@given(st.integers(),
       arrays("float64",
              st.shared(array_shapes(max_dims=1), key="number-of-snps"),
              elements=st.floats()),
       arrays("float64",
              st.shared(array_shapes(max_dims=1), key="number-of-snps"),
              elements=st.floats()))
def test_read_write_summary_are_inverses(n, mean, std):
    with tempfile.TemporaryFile() as file:
        write_summary(file, Summary(n, mean, std))
        file.seek(0)
        summary = read_summary(file)
        assert ((summary.n == n) and
                (summary.mean == approx(mean, nan_ok=True)) and
                (summary.std == approx(std, nan_ok=True)))

@st.composite
def properties_and_whitespace(draw):
    n = draw(st.integers(min_value=0, max_value=10))
    return (draw(st.dictionaries(st.text(st.characters(codec="ascii",
                                                       exclude_categories=("Cc", "Zs")),
                                         min_size=1)
                                 .filter(lambda key: not key.startswith("#")),
                                 st.text(st.characters(codec="ascii",
                                                       exclude_categories=("Cc",)),
                                         min_size=1),
                                 min_size=n, max_size=n)),
            draw(st.integers(min_value=0, max_value=10)),
            draw(st.lists(st.integers(min_value=0, max_value=10),
                          min_size=n, max_size=n)))

@given(properties_and_whitespace())
def test_read_summary_headers_variable_whitespace(properties_and_whitespace):
    properties, header_whitespace, key_value_whitespace = properties_and_whitespace
    with tempfile.TemporaryFile() as file:
        file.write(b"#")
        file.write(b" " * header_whitespace)
        file.write(b"pyhegp summary file version 1\n")
        for (key, value), key_value_whitespace in zip(properties.items(), key_value_whitespace):
            file.write(b"#")
            file.write(b" " * key_value_whitespace)
            file.write(f"{key} {value}\n".encode("ascii"))
        file.seek(0)
        assert properties == read_summary_headers(file)
