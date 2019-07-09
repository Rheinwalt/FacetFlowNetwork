import numpy as np
from matplotlib import pyplot as pl
from FacetFlowNetwork import ffn

def thinning(x, y, z, r):
    from scipy.spatial import cKDTree as kdtree

    b = np.ones(len(x), dtype = 'bool')   
    tree = kdtree(np.transpose((x, y)))
    lsts = tree.query_ball_tree(tree, r = r)
    for l in lsts:
        if len(l) > 1:
            la = np.array(l)
            dz = np.abs(z[la] - np.median(z[la]))
            b[la] = False
            b[la[np.argmin(dz)]] = True

    print('keeping %.3f %%' % (100.*len(b[b])/len(b)))
    return (x[b], y[b], z[b])

def data(fname):
    from laspy.file import File

    print('loading %s ..' % fname)
    f = File(fname, mode = 'r')
    x = f.x
    y = f.y
    z = f.z
    f.close()
    print('we have %.2e points' % len(x))

    print('add low noise ..')
    x += np.random.random(len(x)) / 1000.0
    y += np.random.random(len(x)) / 1000.0
    z += np.random.random(len(x)) / 1000.0

    print('remove points that are closer than 10cm ..')
    x, y, z = thinning(x, y, z, 0.1)
    print('we have %.2e points' % len(x))    
    return (x, y, z)

x, y, z = data('Alluvial_fan.laz')
G = ffn(x, y, z)
G.save('Alluvial_fan.hdf')
