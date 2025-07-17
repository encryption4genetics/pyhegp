pyhegp is a Python library and CLI utility implementing homomorphic encryption of genotypes and phenotypes as described in [Private Genomes and Public SNPs: Homomorphic Encryption of Genotypes and Phenotypes for Shared Quantitative Genetics](https://academic.oup.com/genetics/article/215/2/359/5930450).

# Install development version

In a new directory, create a python virtual environment and activate it.
```
mkdir pyhegp
cd pyhegp
python3 -m venv .venv
source .venv/bin/activate
```
Install the development version of pyhegp into the virtual environment.
```
pip install git+https://github.com/encryption4genetics/pyhegp
```

# Run tests

Run the test suite using
```
python3 -m pytest
```

# License

pyhegp is free software released under the terms of the [GNU General Public License](https://www.gnu.org/licenses/gpl.html), either version 3 of the License, or (at your option) any later version.
