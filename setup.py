from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension("auto_corr3_c",
              ["auto_corr3_c.pyx"],
              library_dirs=['/ccs/home/z8g/download/gromacs-4.0.5/lib'],
              libraries=["m","gmx"]) # Unix-like specific
    
    ]

setup(
    name = "Demos",
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules,
    include_dirs=['/ccs/home/z8g/download/gromacs-4.0.5/include/',
                  '/sw/yona/python/2.7.1/centos5.5_gnu4.1.2/2.7_install_dir/include/python2.7/',
                  '/sw/yona/python/2.7.1/centos5.5_gnu4.1.2/2.7_install_dir/lib/python2.7/site-packages/numpy/core/include/']
    )
