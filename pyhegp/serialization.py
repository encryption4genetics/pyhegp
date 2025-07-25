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
from itertools import takewhile

import numpy as np

SUMMARY_HEADER = b"# pyhegp summary file version 1\n"

Summary = namedtuple("Summary", "n mean std")

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
                   *np.loadtxt(file, ndmin=2))

def write_summary(file, summary):
    file.write(SUMMARY_HEADER)
    file.write(f"# number-of-samples {summary.n}\n".encode("ascii"))
    np.savetxt(file,
               np.row_stack((summary.mean, summary.std)),
               fmt="%.8g")

def read_genotype(genotype_file):
    return np.loadtxt(genotype_file, delimiter=",")

def write_genotype(file, genotype):
    np.savetxt(file, genotype, delimiter=",", fmt="%f")
