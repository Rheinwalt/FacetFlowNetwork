import numpy as np
from matplotlib import pyplot as pl
from matplotlib.colors import LogNorm
import matplotlib.tri as mtri
from FacetFlowNetwork import ffn

# load FFN from disk
G = ffn(fname = 'Alluvial_fan.hdf')

# LAS export with average SCA at points
var = np.log10(G.fatp(G.sca()))
G.export('Alluvial_fan_sca.las', var)

# LAS export with only the upper quartile SCA
pts = np.transpose((G.x, G.y, G.z))
sel = var > np.percentile(var, 75)
pts = pts[sel]
var = var[sel]
G.export('Alluvial_fan_sca_p75.las', var, pts)

# detailed map view of SCA around xc, yc
xc, yc = 451708, 4060410
ww = 50
sca = G.sca()
pts = G.facet_centroids()
xf, yf = pts[:, 0], pts[:, 1]
sel = (xc-ww <= xf)*(xf < xc+ww)*(yc-ww <= yf)*(yf < yc+ww)
tri, sca = G.tri[sel], sca[sel]

pl.figure(1, (8, 6))
pl.title('Alluvial fan FFN SCA estimate')
pl.tripcolor(G.x-xc+ww, G.y-yc+ww, tri, facecolors = sca,
             cmap = pl.cm.magma_r,
             norm = LogNorm(vmin = 0.1, vmax = 1e5))
cb = pl.colorbar()
cb.set_label('SCA [m]')
pl.xlim((0, ww+ww))
pl.ylim((0, ww+ww))
pl.xlabel('x - %i [m]' % (xc-ww))
pl.ylabel('y - %i [m]' % (yc-ww))
pl.gca().set_aspect('equal')
pl.savefig('Alluvial_fan_sca.pdf')

pl.close('all')
pl.figure(1, (8, 6))
pl.title('Alluvial fan FFN SCA estimate')
pl.tripcolor(G.x-xc+ww, G.y-yc+ww, tri, facecolors = sca,
             cmap = pl.cm.magma_r,
             norm = LogNorm(vmin = 0.1, vmax = 1e5))
cb = pl.colorbar()
cb.set_label('SCA [m]')
pl.xlim((0, ww+ww))
pl.ylim((0, ww+ww))
pl.xlabel('x - %i [m]' % (xc-ww))
pl.ylabel('y - %i [m]' % (yc-ww))
pl.gca().set_aspect('equal')
pl.savefig('Alluvial_fan_sca.png')

# detailed map view of SCA with gouraud shading
sca = G.fatp(G.sca())
x, y = G.x, G.y
sel = (xc-ww <= x)*(x < xc+ww)*(yc-ww <= y)*(y < yc+ww)
x, y, sca = x[sel]-xc+ww, y[sel]-yc+ww, sca[sel]
tri = mtri.Triangulation(x, y)

pl.close('all')
pl.figure(1, (8, 6))
pl.title('Alluvial fan FFN SCA estimate')
pl.tripcolor(tri, sca, shading = 'gouraud',
             cmap = pl.cm.magma_r,
             norm = LogNorm(vmin = 0.1, vmax = 1e5))
cb = pl.colorbar()
cb.set_label('SCA [m]')
pl.xlim((0, ww+ww))
pl.ylim((0, ww+ww))
pl.xlabel('x - %i [m]' % (xc-ww))
pl.ylabel('y - %i [m]' % (yc-ww))
pl.gca().set_aspect('equal')
pl.savefig('Alluvial_fan_sca_gouraud.pdf')

pl.close('all')
pl.figure(1, (8, 6))
pl.title('Alluvial fan FFN SCA estimate')
pl.tripcolor(tri, sca, shading = 'gouraud',
             cmap = pl.cm.magma_r,
             norm = LogNorm(vmin = 0.1, vmax = 1e5))
cb = pl.colorbar()
cb.set_label('SCA [m]')
pl.xlim((0, ww+ww))
pl.ylim((0, ww+ww))
pl.xlabel('x - %i [m]' % (xc-ww))
pl.ylabel('y - %i [m]' % (yc-ww))
pl.gca().set_aspect('equal')
pl.savefig('Alluvial_fan_sca_gouraud.png')

