from distutils.core import setup, Extension
from Cython.Build import cythonize


setup(
	ext_modules = cythonize(Extension(
		name = "Matrix",
		sources = ["Matrix.pyx", "vMatrix.cpp"],
		extra_compile_args=["-std=c++11"],
		extra_link_args=["-std=c++11"])))