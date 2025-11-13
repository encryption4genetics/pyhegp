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

import sys

import pandas as pd

data_file = sys.argv[1]
query_expression = sys.argv[2]

if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], sep="\t")
    qtl = df.query(query_expression)
    # Assert that the QTL is on chromosome 4.
    assert (qtl.chromosome == 4).all()
    # Assert that the QTL is within 2 Mb of the expected position.
    assert ((qtl.position - 137715608).abs() < 2*10**6).all()
