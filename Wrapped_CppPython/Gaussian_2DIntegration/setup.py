from distutils.core import setup, Extension

'''
Setup file to include c++ files as modules for use in python.

> python setup.py build_ext --inplace
   --inplace: This tells distutils to put the extension lib in the current directory.
              Otherwise, it will put it inside a build hierarchy, and you'd have to move it to use it.
'''

GInt2D = Extension('_GaussInt2D', # name of the extension
                    sources=['SwigGaussInt2D.i', 'GaussianIntegration.cpp'], # a list of source filenames. May be C,C++,Objective-C,SWIG
                    swig_opts=["-c++"],
                  )

setup (name = 'GaussInt2D',
       author      = 'Tony Alberti',
       description = '2D_Gaussian_Integration',
       ext_modules = [GInt2D],
       )
