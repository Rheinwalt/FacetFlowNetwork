import sys
import numpy as np
from libc.stdlib cimport free

# for LAS file RGB colors
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as pl

__version__ = '0.1'
__author__ = 'Aljoscha Rheinwalt'
__author_email__ = 'aljoscha.rheinwalt@uni-potsdam.de'

cdef extern from "main.c":
    void network(unsigned int *net, double *w, double *a, double *d, double *p, double *t,
            const unsigned int *tri, const double *x, const double *y, const double *z,
            const unsigned int m, const unsigned int n)
    void tunnel(unsigned int *net, double *w, double *a,
            const unsigned int *tri, const double *x, const double *y, const double *z,
            const unsigned int m, const unsigned int n, const unsigned int ml)
    void linkthroughput(double *ltp,
            const unsigned int *net, const double *w, const double *a,
            const unsigned int m)
    unsigned int * upstreamnetwork(const unsigned int *net,
            const unsigned int m, unsigned int *l)

class ffn:
    """
    Facet flow network class
    """
    def __init__(self, x = None, y = None, z = None,
            tri = None, fname = None,
            tunneling = True, tunnelmaxlen = 100):
        cdef double[:] xv, yv, zv
        cdef unsigned int[:,:] tv
        cdef unsigned int[:,:] nv
        cdef unsigned int m, n
        cdef double[:,:] wv, av
        cdef double[:] dv, pv, sv
        if fname is None:
            if tri is None:
                print('Delaunay triangulation ..')
                from scipy.spatial import Delaunay
                tri = Delaunay(np.transpose((x-np.mean(x), y-np.mean(y))),
                        qhull_options = 'Qt').simplices.astype('uint32')
            if x is None:
                sys.exit('error: fname and x both not specified!')

            self.x = np.array(x)
            self.y = np.array(y)
            self.z = np.array(z)
            self.tri = tri
            self.tunneled = False
            self.ltp = None
            self.m = tri.shape[0]
            self.n = len(x)
            self.net = np.zeros((self.m, 2), dtype = 'uint32')
            self.w = np.zeros((self.m, 2), dtype = 'float64')
            self.a = np.zeros((self.m, 2), dtype = 'float64')
            self.d = np.zeros(self.m, dtype = 'float64')
            self.aspect = np.zeros(self.m, dtype = 'float64')
            self.slope = np.zeros(self.m, dtype = 'float64')
            m = self.m
            n = self.n
            nv = self.net
            wv = self.w
            av = self.a
            dv = self.d
            pv = self.aspect
            sv = self.slope
            xv = self.x
            yv = self.y
            zv = self.z
            tv = self.tri
            print('FFN construction ..')
            network(&nv[0,0], &wv[0,0], &av[0,0], &dv[0], &pv[0], &sv[0],
                    &tv[0,0], &xv[0], &yv[0], &zv[0], m, n)
            self.aspect *= 180 / np.pi
            self.slope *= 180 / np.pi
            if tunneling:
                print('resolving sinks by tunneling ..')
                tunnel(&nv[0,0], &wv[0,0], &av[0,0],
                       &tv[0,0], &xv[0], &yv[0], &zv[0],
                       m, n, tunnelmaxlen)
                self.tunneled = True

        else:
            from h5py import File

            f = File(fname, 'r')
            self.x = f['x'][:]
            self.y = f['y'][:]
            self.z = f['z'][:]
            self.tri = f['tri'][:]
            self.net = f['net'][:]
            self.w = f['w'][:]
            self.a = f['a'][:]
            self.d = f['d'][:]
            self.aspect = f['aspect'][:]
            self.slope = f['slope'][:]
            flags = f['flags'][:]
            self.tunneled = bool(flags[0])
            self.ltp = None
            self.m = self.tri.shape[0]
            self.n = self.x.shape[0]

    def __str__(self):
        return 'ffn(%.1e points, %.1e facets)' % (self.n, self.m)

    def tunneling(self, tunnelmaxlen = 100):
        """
        Tunnel flow out of sinks
        """
        cdef double[:] xv, yv, zv
        cdef unsigned int[:,:] tv
        cdef unsigned int[:,:] nv
        cdef unsigned int m, n
        cdef double[:,:] wv, av
        if not self.tunneled:
            m = self.m
            n = self.n
            nv = self.net
            wv = self.w
            av = self.a
            xv = self.x
            yv = self.y
            zv = self.z
            tv = self.tri
            tunnel(&nv[0,0], &wv[0,0], &av[0,0],
                   &tv[0,0], &xv[0], &yv[0], &zv[0],
                   m, n, tunnelmaxlen)
            self.tunneled = True

    def upstream(self):
        """
        Returns the FFN in reverse, links point upstream
        """
        cdef unsigned int[:] rv
        cdef unsigned int[:,:] nv
        cdef unsigned int m, l
        cdef unsigned int *p
        cdef unsigned int[:] pv

        nv = self.net
        m = self.m
        p = upstreamnetwork(&nv[0,0], m, &l)
        pv = <unsigned int[:l]> p
        r = np.asarray(pv).copy()
        free(p)
        return r

    def sca(self):
        """
        Specific catchment area
        """
        cdef double[:,:] lv, wv, av
        cdef unsigned int[:,:] nv
        cdef unsigned int m
        if self.ltp is None:
            self.ltp = np.zeros((self.m, 2), dtype = 'float64')
            lv = self.ltp
            nv = self.net
            wv = self.w
            av = self.a
            m = self.m
            linkthroughput(&lv[0,0], &nv[0,0], &wv[0,0], &av[0,0], m)
        a = self.ltp.sum(axis = 1)
        return a / self.d

    def tda(self):
        """
        Total drainage area
        """
        cdef double[:,:] lv, wv, av
        cdef unsigned int[:,:] nv
        cdef unsigned int m
        if self.ltp is None:
            self.ltp = np.zeros((self.m, 2), dtype = 'float64')
            lv = self.ltp
            nv = self.net
            wv = self.w
            av = self.a
            m = self.m
            linkthroughput(&lv[0,0], &nv[0,0], &wv[0,0], &av[0,0], m)
        a = self.ltp.sum(axis = 1)
        return a

    def fatp(self, var):
        """
        Facet values averaged to point values
        """
        cdef unsigned int i, k
        cdef double[::] rv
        cdef double[:] vv
        cdef unsigned int[::] cv
        cdef unsigned int[:, :] tv
        assert len(var) == self.m, '%i facets and %i values' % (self.m, len(var))
        r = np.zeros(self.n, dtype = 'float64')
        c = np.zeros(self.n, dtype = 'uint32')
        rv = r
        cv = c
        tv = self.tri
        vv = var
        for i in range(self.m):
            for k in range(3):
                rv[tv[i, k]] += vv[i]
                cv[tv[i, k]] += 1
        return r / c

    def fmtp(self, var):
        """
        Facet values maximum to point values
        """
        cdef unsigned int i, k
        cdef double[::] rv
        cdef double[:] vv
        cdef unsigned int[:, :] tv
        assert len(var) == self.m, '%i facets and %i values' % (self.m, len(var))
        r = np.zeros(self.n, dtype = 'float64')
        rv = r
        tv = self.tri
        vv = var
        for i in range(self.m):
            for k in range(3):
                if vv[i] > rv[tv[i, k]]:
                    rv[tv[i, k]] = vv[i]
        return r

    def save(self, fname, compr = 'lzf'):
        """
        Save FFN to a compressed HDF file
        """
        import os
        from h5py import File

        fn, fe = os.path.splitext(fname)
        f = File('%s.hdf' % fn, 'w')
        f.create_dataset('x', data = self.x, compression = compr)
        f.create_dataset('y', data = self.y, compression = compr)
        f.create_dataset('z', data = self.z, compression = compr)
        f.create_dataset('tri', data = self.tri, compression = compr)
        f.create_dataset('net', data = self.net, compression = compr)
        f.create_dataset('w', data = self.w, compression = compr)
        f.create_dataset('a', data = self.a, compression = compr)
        f.create_dataset('d', data = self.d, compression = compr)
        f.create_dataset('aspect', data = self.aspect, compression = compr)
        f.create_dataset('slope', data = self.slope, compression = compr)
        flags = np.zeros(2)
        if self.tunneled:
            flags[0] = 1
        f.create_dataset('flags', data = flags, compression = compr)
        f.close()

    def export(self, fname, var, pnts = None, cmap = pl.cm.magma_r):
        """
        Export point cloud to LAS file with RGB colors according to var
        """
        import laspy
        import os
        from subprocess import call

        fn, fe = os.path.splitext(fname)
        v = var - np.min(var)
        v /= v.max()
        rgb = cmap(v)
        rgb = rgb[:, :3]
        rgb *= 65535
        rgb = rgb.astype('uint')
        header = laspy.header.Header()
        header.data_format_id = 2
        f = laspy.file.File('%s.las' % fn, mode = 'w', header = header)
        f.header.scale = [0.001, 0.001, 0.001]
        if pnts is None:
            f.header.offset = [self.x.min(), self.y.min(), self.z.min()]
            f.x = self.x
            f.y = self.y
            f.z = self.z
        else:
            f.header.offset = [pnts[:,0].min(), pnts[:,1].min(), pnts[:,2].min()]
            f.x = pnts[:, 0]
            f.y = pnts[:, 1]
            f.z = pnts[:, 2]
        f.set_red(rgb[:, 0])
        f.set_green(rgb[:, 1])
        f.set_blue(rgb[:, 2])
        f.close()
        try:
            r = call(['laszip', '%s.las' % fn])
            if not r:
                r = call(['rm', '%s.las' % fn])
        except:
            pass

