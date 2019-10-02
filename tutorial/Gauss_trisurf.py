import numpy as np
from matplotlib import pyplot as pl
from mpl_toolkits.mplot3d import proj3d
import matplotlib.tri as mtri

def GaussianHill(n = 1000):
    """
    Create a Gaussian hill sampling point cloud with n points.
    """
    x = np.random.random(n)-0.5
    y = np.random.random(n)-0.5
    x *= 3
    y *= 3
    z = np.exp(-x*x-y*y)
    return (x, y, z)

x, y, z = GaussianHill()
tri = mtri.Triangulation(x, y)

fg, ax = pl.subplots(1, 1, subplot_kw = dict(projection = '3d'))
pl.title('Gaussian hill TIN')
ax.plot_trisurf(x, y, z, triangles = tri.triangles,
                cmap = pl.cm.magma)
pl.show()
