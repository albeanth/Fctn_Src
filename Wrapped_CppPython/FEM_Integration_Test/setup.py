from distutils.core import setup, Extension

'''
Setup file to include c++ files as modules for use in python.
> python setup.py build_ext --inplace
   --inplace: This tells distutils to put the extension lib in the current directory.
              Otherwise, it will put it inside a build hierarchy, and you'd have to move it to use it.
'''

FEMInt = Extension('_FEMInt', # name of the extension
                    sources=['SwigFEMInt.i', 'FEM_Integration.cpp'], # a list of source filenames. May be C,C++,Objective-C,SWIG
                    swig_opts=["-c++"],
                  )

setup (name = 'FEMInt',
       ext_modules = [FEMInt],
       author      = 'Tony Alberti',
       description = 'FEM_Integration_Test',
       )
