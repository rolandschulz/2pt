from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension("auto_corr3_c",
              ["auto_corr3_c.pyx"],
              library_dirs=['/ccs/home/z8g/software/gromacs/lib'],
              libraries=["m","gmx"]) # Unix-like specific
    
    ]

setup(
    name = "Demos",
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules,
    include_dirs=['/ccs/home/z8g/software/gromacs/include/gromacs',
                  '/sw/analysis-x64/scipy/0.9.0rc5/centos5.5_gcc4.1.2/install_dir_numpy/lib/python2.7/site-packages/numpy/core/include'],
    )
