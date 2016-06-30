from setuptools import setup, Extension
import os

from Cython.Build import cythonize
import numpy as np

# Supplies __version__
from ch_frb_io import __version__



REQUIRES = ['numpy', 'scipy', 'cython', 'h5py', 'bitshuffle',]

# Generate test data.
from ch_frb_io.tests.data import generate
generate()

COMPILE_FLAGS = ['-O3', '-ffast-math', '-march=native']

EXTENSIONS = [
        ]


extensions = cythonize(EXTENSIONS,
        include_path=[np.get_include()],
        )

setup(
    name = 'ch_frb_io',
    version = __version__,
    packages = ['ch_frb_io', 'ch_frb_io.tests'],
    ext_modules = extensions,
    scripts=['bin/decompress-chfrb-data'],
    install_requires=REQUIRES,
    package_data = {'ch_frb_io.tests' : ['data/*']},

    author = "Kiyoshi Masui, Kendrick Smith",
    author_email = "kiyo@physics.ubc.ca",
    description = "CHIME FRB IO",
    url = "https://github.com/CHIMEFRB/ch_frb_io"
)
