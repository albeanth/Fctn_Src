from distutils.core import setup, Extension

'''
Setup file to include c++ files as modules for use in python.
> python setup.py build_ext --inplace
   --inplace: This tells distutils to put the extension lib in the current directory.
              Otherwise, it will put it inside a build hierarchy, and you'd have to move it to use it.
'''

InterpTest = Extension(name = '_InterpTest', # name of the library (.so file) made
                    sources=['Interp_Test.i', 'Interp_Test.cpp'], # a list of source filenames. May be C,C++,Objective-C,SWIG
                    swig_opts=['-c++',],
                    include_dirs=['/usr/local/include'],
                    library_dirs=['/usr/local/lib'],
                    libraries=['gsl', 'gslcblas', 'm'],
                    runtime_library_dirs=['/usr/local/lib'],
                    # extra_compile_args=['-fopenmp'],
                    )

setup(name = 'InterpTest',
      ext_modules = [InterpTest],
      author      = 'Tony_Alberti',
      description = 'Interpolation_Test',
      )
