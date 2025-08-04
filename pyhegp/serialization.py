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
import csv
from itertools import takewhile

import pandas as pd

SUMMARY_HEADER = b"# pyhegp summary file version 1\n"

Summary = namedtuple("Summary", "n data")

def peek(file):
    c = file.read(1)
    file.seek(-1, 1)
    return c

def header_lines(file):
    while peek(file) == b"#":
        yield file.readline()

def read_summary_headers(file):
    assert (file.readline().decode("ascii").lstrip("#").lstrip()
            == SUMMARY_HEADER.decode("ascii").lstrip("#").lstrip())
    return dict(line.decode("ascii").rstrip("\n").lstrip("#").lstrip().split(" ", maxsplit=1)
                for line in header_lines(file))

def read_summary(file):
    headers = read_summary_headers(file)
    return Summary(int(headers["number-of-samples"]),
                   pd.read_csv(file,
                               sep="\t",
                               header=0,
                               dtype={"chromosome": "str",
                                      "position": "int",
                                      "reference": "str",
                                      "mean": "float",
                                      "standard-deviation": "float"},
                               na_filter=False)
                   .rename(columns={"standard-deviation": "std"}))

def write_summary(file, summary):
    file.write(SUMMARY_HEADER)
    file.write(f"# number-of-samples {summary.n}\n".encode("ascii"))
    (summary.data
     .rename(columns={"std": "standard-deviation"})
     .to_csv(file,
             sep="\t",
             float_format="%.8g",
             index=False))

def read_genotype(file):
    df = pd.read_csv(file,
                     quoting=csv.QUOTE_NONE,
                     sep="\t",
                     na_filter=False)
    sample_columns = [column
                      for column in df.columns
                      if column not in ["chromosome", "position", "reference"]]
    df.chromosome = df.chromosome.astype("str")
    df.position = df.position.astype("int")
    if "reference" in df:
        df.reference = df.reference.astype("str")
    df[sample_columns] = df[sample_columns].astype("float")
    return df

def write_genotype(file, genotype):
    (genotype
     .to_csv(file,
             quoting=csv.QUOTE_NONE,
             sep="\t",
             float_format="%.8g",
             index=False))
