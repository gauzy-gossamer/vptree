from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

ext = Extension("vptree", ["vptree.pyx"],
    include_dirs = [numpy.get_include()])

setup(
    name = "vptree",
    ext_modules = cythonize('*.pyx'),
    include_dirs = [numpy.get_include()]
)
#cythonize('*.pyx')
