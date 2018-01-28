from distutils.core import setup, Extension

'''
Setup file to include c++ files as modules for use in python.
> python setup.py build_ext --inplace
   --inplace: This tells distutils to put the extension lib in the current directory.
              Otherwise, it will put it inside a build hierarchy, and you'd have to move it to use it.
'''

GInt = Extension('_GaussInt_XT', # name of the extension
                    sources=['SwigGaussInt.i', 'Mesh_InterfaceXT.cpp'], # a list of source filenames. May be C,C++,Objective-C,SWIG
                    swig_opts=["-c++"],
                    extra_compile_args=["-fopenmp"],
                    extra_link_args=["-fopenmp"],
                  )

setup (name = 'GaussInt_XT',
       ext_modules = [GInt],
       author      = 'Tony Alberti',
       description = 'Gaussian_Integration',
       )
