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

import numpy as np

from itertools import accumulate, pairwise
from scipy.linalg import block_diag

class BlockDiagonalMatrix:
    def __init__(self, _blocks):
        self.blocks = _blocks
        self.shape = (sum(len(block) for block in self.blocks),) * 2

    def __len__(self):
        return self.shape[0]

    def __repr__(self):
        return f"BlockDiagonalMatrix{self.blocks}"

    def __array_function__(self, func, types, args, kwargs):
        if ((func is np.transpose)
            and (all(issubclass(type, BlockDiagonalMatrix)
                     for type in types))):
            return BlockDiagonalMatrix([np.transpose(block)
                                        for block in self.blocks])
        else:
            return NotImplemented

    def __array__(self):
        return block_diag(*self.blocks)

    def __matmul__(self, multiplier):
        return np.concatenate(
            [block @ multiplier[start:stop, ...]
             for (start, stop), block
             in zip(pairwise(accumulate((len(block) for block in self.blocks),
                                        initial=0)),
                    self.blocks)])

    def savetxt(self, file, *args, **kwargs):
        return np.savetxt(file, self.to_ndarray(), *args, **kwargs)
