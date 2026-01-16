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
import pandas as pd
from pytest import approx

from pyhegp.serialization import read_summary, write_summary, read_summary_headers, read_genotype, write_genotype, read_phenotype, write_phenotype, read_key, write_key

from helpers.strategies import summaries, genotype_frames, phenotype_frames, keys

@given(summaries())
def test_read_write_summary_are_inverses(summary):
    with tempfile.TemporaryFile() as file:
        write_summary(file, summary)
        file.seek(0)
        recovered_summary = read_summary(file)
        pd.testing.assert_frame_equal(summary.data,
                                      recovered_summary.data)
        assert summary.n == recovered_summary.n

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

@given(genotype_frames())
def test_read_write_genotype_are_inverses(genotype):
    with tempfile.TemporaryFile() as file:
        write_genotype(file, genotype)
        file.seek(0)
        pd.testing.assert_frame_equal(genotype, read_genotype(file))

@given(phenotype_frames())
def test_read_write_phenotype_are_inverses(phenotype):
    with tempfile.TemporaryFile() as file:
        write_phenotype(file, phenotype)
        file.seek(0)
        pd.testing.assert_frame_equal(phenotype, read_phenotype(file))

@given(keys(st.integers(min_value=2, max_value=10)))
def test_read_write_key_are_inverses(key):
    with tempfile.TemporaryFile() as file:
        write_key(file, key)
        file.seek(0)
        assert key == approx(read_key(file), nan_ok=True)
