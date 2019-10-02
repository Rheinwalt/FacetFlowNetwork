from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules = [Extension('FacetFlowNetwork',
                        ['FacetFlowNetwork.pyx'],
                        include_dirs = [numpy.get_include()],
                        extra_compile_args = ['-fopenmp'],
                        extra_link_args = ['-lgomp'],
                        )]

setup(name = 'FacetFlowNetwork',
      author = 'Aljoscha Rheinwalt',
      author_email = 'aljoscha.rheinwalt@uni-potsdam.de',
      ext_modules = cythonize(ext_modules,
      compiler_directives={'language_level': 3}))

