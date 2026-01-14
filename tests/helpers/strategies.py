### pyhegp --- Homomorphic encryption of genotypes and phenotypes
### Copyright © 2025–2026 Arun Isaac <arunisaac@systemreboot.net>
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

from hypothesis import strategies as st
from hypothesis.extra.pandas import column, columns, data_frames, range_indexes
import pandas as pd
from scipy.stats import special_ortho_group
from typing import assert_never

from pyhegp.serialization import Summary, is_genotype_metadata_column, is_phenotype_metadata_column
from pyhegp.utils import negate

tabless_printable_ascii_text = st.text(
    # Exclude control characters and tab.
    st.characters(codec="ascii",
                  exclude_categories=("Cc",),
                  exclude_characters=("\t",)),
    min_size=1)

chromosomes = tabless_printable_ascii_text
positions = st.integers(min_value=0,
                        max_value=10*10**9)
references = st.text(st.characters(codec="ascii",
                                   categories=(),
                                   include_characters=("A", "G", "C", "T")),
                     min_size=1)

sample_names = (tabless_printable_ascii_text
                .filter(negate(is_genotype_metadata_column)))

def genotype_metadata(draw, number_of_snps, reference_present):
    match list(zip(*draw(st.lists(st.tuples(chromosomes, positions, references)
                                  if reference_present
                                  else st.tuples(chromosomes, positions),
                                  min_size=number_of_snps,
                                  max_size=number_of_snps,
                                  unique=True)))):
        case []:
            return pd.DataFrame({"chromosome": pd.Series(dtype="str"),
                                 "position": pd.Series(dtype="int")}
                                | ({"reference": pd.Series(dtype="str")}
                                   if reference_present else {}))
        case chromosomes_lst, positions_lst, *references_lst:
            return pd.DataFrame({"chromosome": pd.Series(chromosomes_lst, dtype="str"),
                                 "position": pd.Series(positions_lst, dtype="int")}
                                | ({"reference": pd.Series(*references_lst, dtype="str")}
                                   if reference_present else {}))
        case _ as unreachable:
            assert_never(unreachable)

@st.composite
def summaries(draw):
    stats = draw(data_frames(
        columns=columns(["mean", "std"],
                        dtype="float64",
                        elements=st.floats(allow_nan=False))))
    return Summary(draw(st.integers()),
                   pd.concat((genotype_metadata(draw,
                                                len(stats),
                                                draw(st.booleans())),
                              stats),
                             axis="columns"))

@st.composite
def genotype_frames(draw,
                    number_of_samples=st.integers(min_value=0,
                                                  max_value=10),
                    reference_present=st.booleans()):
    _number_of_samples = draw(number_of_samples)
    dosages = draw(data_frames(
        columns=columns(draw(st.lists(sample_names,
                                      min_size=_number_of_samples,
                                      max_size=_number_of_samples,
                                      unique=True)),
                        dtype="float64",
                        elements=st.floats(min_value=0,
                                           max_value=100,
                                           allow_nan=False))))
    return pd.concat((genotype_metadata(draw,
                                        len(dosages),
                                        draw(reference_present)),
                      dosages),
                     axis="columns")

phenotype_names = st.lists(tabless_printable_ascii_text
                           .filter(negate(is_phenotype_metadata_column)),
                           unique=True)

@st.composite
def phenotype_frames(draw,
                     number_of_samples=st.integers(min_value=0,
                                                   max_value=10),
                     intercept_present=st.booleans()):
    _number_of_samples = draw(number_of_samples)
    return draw(data_frames(
        columns=([column(name="sample-id",
                         dtype="str",
                         elements=tabless_printable_ascii_text,
                         unique=True)]
                 + ([column(name="intercept",
                            dtype="float64",
                            elements=st.floats(min_value=-1,
                                               max_value=1,
                                               allow_nan=False))]
                    if draw(intercept_present)
                    else [])
                 + columns(draw(phenotype_names),
                           dtype="float64",
                           elements=st.floats(min_value=-1000,
                                              max_value=1000,
                                              allow_nan=False))),
        index=range_indexes(min_size=_number_of_samples,
                            max_size=_number_of_samples)))

@st.composite
def keys(draw, size):
    return (special_ortho_group(draw(size),
                                seed=draw(st.integers(min_value=0,
                                                      max_value=2**32-1)))
            .rvs())
