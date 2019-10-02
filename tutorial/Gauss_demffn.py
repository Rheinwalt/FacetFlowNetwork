import numpy as np
from matplotlib import pyplot as pl
from FacetFlowNetwork import demffn

pw = 0.01
xr = np.arange(-1.5, 1.5, pw)
yr = np.arange(-1.4, 1.4, pw)

x, y = np.meshgrid(xr, yr)
r = np.sqrt(x*x + y*y)

dem = np.exp(-r*r)
tsca = r / 2.0

f = demffn(dem, pw)
sca = f.fatp(f.sca())
sca.shape = dem.shape

pl.figure(1, (8, 6))
im = pl.imshow(sca, origin = 'lower',
        interpolation = 'none', extent = [0,1.5,0,1.4],
        cmap = pl.cm.viridis_r,
        vmin = 0, vmax = tsca.max())
cb = pl.colorbar()
cb.set_label('SCA [m]')
pl.gca().set_aspect('equal')
pl.savefig('Gauss_demffn_sca.png')

pl.close('all')

pl.figure(1, (8, 6))
im = pl.imshow(100*(sca-tsca)/tsca, origin = 'lower',
        interpolation = 'none', extent = [0,1.5,0,1.4],
        cmap = pl.cm.bwr,
        vmin = -25, vmax = 25)
cb = pl.colorbar()
cb.set_label(r'Relative SCA error [%]')
pl.gca().set_aspect('equal')
pl.savefig('Gauss_demffn_sca_error.png')
