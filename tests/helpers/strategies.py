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

from hypothesis import strategies as st
from hypothesis.extra.pandas import column, columns, data_frames
from scipy.stats import special_ortho_group

from pyhegp.serialization import Summary, is_genotype_metadata_column, is_phenotype_metadata_column
from pyhegp.utils import negate

tabless_printable_ascii_text = st.text(
    # Exclude control characters and tab.
    st.characters(codec="ascii",
                  exclude_categories=("Cc",),
                  exclude_characters=("\t",)),
    min_size=1)

chromosome_column = column(name="chromosome",
                           dtype="str",
                           elements=tabless_printable_ascii_text)

position_column = column(name="position",
                         dtype="int")

reference_column = column(name="reference",
                          dtype="str",
                          elements=st.text(
                              st.characters(codec="ascii",
                                            categories=(),
                                            include_characters=("A", "G", "C", "T")),
                              min_size=1))

sample_names = (tabless_printable_ascii_text
                .filter(negate(is_genotype_metadata_column)))

@st.composite
def summaries(draw):
    return Summary(draw(st.integers()),
                   draw(data_frames(
                       columns=([chromosome_column, position_column]
                                + ([reference_column] if draw(st.booleans()) else [])
                                + columns(["mean", "std"],
                                          dtype="float64",
                                          elements=st.floats(allow_nan=False))))))

@st.composite
def genotype_frames(draw,
                    number_of_samples=st.integers(min_value=0,
                                                  max_value=10),
                    reference_present=st.booleans()):
    _number_of_samples = draw(number_of_samples)
    genotype = draw(data_frames(
        columns=([chromosome_column, position_column]
                 + ([reference_column]
                    if draw(reference_present)
                    else [])
                 + columns(draw(st.lists(sample_names,
                                         min_size=_number_of_samples,
                                         max_size=_number_of_samples,
                                         unique=True)),
                           dtype="float64",
                           elements=st.floats(allow_nan=False)))))
    return genotype.drop_duplicates(subset=list(
        filter(is_genotype_metadata_column,
               genotype.columns)),
                                    ignore_index=True)

phenotype_names = st.lists(tabless_printable_ascii_text
                           .filter(negate(is_phenotype_metadata_column)),
                           unique=True)

@st.composite
def phenotype_frames(draw):
    phenotype = draw(data_frames(
        columns=([column(name="sample-id",
                         dtype="str",
                         elements=tabless_printable_ascii_text)]
                 + columns(draw(phenotype_names),
                           dtype="float64",
                           elements=st.floats(allow_nan=False)))))
    return phenotype.drop_duplicates(subset=list(
        filter(is_phenotype_metadata_column,
               phenotype.columns)),
                                     ignore_index=True)

@st.composite
def keys(draw, size):
    return (special_ortho_group(draw(size),
                                seed=draw(st.integers(min_value=0,
                                                      max_value=2**32-1)))
            .rvs())
