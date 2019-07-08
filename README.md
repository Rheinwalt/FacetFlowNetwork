# FacetFlowNetwork

Flow accumulation on TINs of point-cloud data 

## Install

The core of this module is written in C (see main.c). It is
wrapped by Cython into a Python module, i.e., you will need Python,
Cython and a C compiler.

    git clone https://github.com/Rheinwalt/FacetFlowNetwork.git
    cd FacetFlowNetwork
    sudo python setup.py install

## Usage

This module has on class called *ffn*, and all functionality is
structured inside this class by class methods.

~~~~~~~~~~~~~~~~~ {.python .numberLines}
from FacetFlowNetwork import ffn
help(ffn)
~~~~~~~~~~~~~~~~~


## Tutorial

~~~~~~~~~~~~~~~~~ {.python .numberLines}
import numpy as np
from matplotlib import pyplot as pl
from FacetFlowNetwork import ffn

def GaussianHill(n = 1000):
    """
    Create a Gaussian hill sampling point cloud with n points.
    """
    x = np.random.random(n)
    y = np.random.random(n)
    z = np.exp(-x*x-y*y)
    return (x, y, z)

# construct FFN from a Gaussian hill
x, y, z = GaussianHill()
G = ffn(x, y, z)

# visualize the specific catchment area (SCA) for each facet of the FFN
pl.title('Gaussian hill FFN SCA estimate')
pl.tripcolor(x, y, G.tri, facecolors = G.sca(), vmax = 1)
cb = pl.colorbar()
cb.set_label('SCA [m]')
pl.xlabel('x [m]')
pl.ylabel('y [m]')
pl.gca().set_aspect('equal')
pl.savefig('Gauss.pdf')

# store FFN to disk
G.save('Gauss.hdf')

# load FFN from disk
F = ffn(fname = 'Gauss.hdf')

# compare the differences, should be zero since both are identical
dsca = G.sca() - F.sca()
print('Sum of deviations: %f' % np.sum(dsca))
del F

# export FFN SCA to LAS file
rgbc = G.sca()
pnts = G.facet_centroids()
rgbc[rgbc > 1] = 1
G.export('Gauss.las', rgbc, pnts)

# alternatively, average FFN SCA to the original point cloud and
# export that to a LAS file
rgbc = G.fatp(G.sca())
rgbc[rgbc > 1] = 1
G.export('Gauss_fatp.las', rgbc)
~~~~~~~~~~~~~~~~~

## Bugs

The C routines inside this module might crash for large point clouds
if the stack size on your system is too small. In that case it helps to
have an unlimited stack for your session:

   ulimit -s unlimited

