# FacetFlowNetwork

Flow accumulation on TINs of point-cloud data by facet flow networks (FFNs).

## Install

The core of this module is written in C (see main.c). It is
wrapped by Cython into a Python module, i.e., you will need Python,
Cython and a C compiler.

    git clone https://github.com/Rheinwalt/FacetFlowNetwork.git
    cd FacetFlowNetwork
    sudo python setup.py install

## Usage

This module has one class called *ffn*, and all functionality is
structured inside this class by class methods.

~~~~~~~~~~~~~~~~~ {.python .numberLines}
from FacetFlowNetwork import ffn
help(ffn)
~~~~~~~~~~~~~~~~~

## Tutorial

First tutorial on synthetic point-cloud data of a Gaussian hill surface.
It covers the generation of 1000 points on a 1 by 1 meter region of interest,
FFN construction, specific catchment area (SCA) estimation, visualization,
input / output of FFNs, and export to LAS files. We recommend
[displaz](https://github.com/c42f/displaz "a hackable lidar viewer") as a
LAS file viewer.

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
~~~~~~~~~~~~~~~~~

![Gaussian hill FFN SCA](tutorial/Gauss.png?raw=true "Gaussian hill FFN SCA")

~~~~~~~~~~~~~~~~~ {.python .numberLines}

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

## Publication

This software is associated to *A network-based flow accumulation algorithm for point clouds: Facet-Flow Networks (FFN)*,
Journal of Geophysical  Research: Earth Surface, 2019 ([10.1029/2018JF004827](https://doi.org/10.1029/2018JF004827)). 
