from ctypes import c_float,cdll,c_double,c_int,byref

from os import environ

mNumPy=1
try:
    import numpy as N
    #from numpy.ctypeslib import ndpointer
    auto_mode=mNumPy                 
except:
    auto_mode=0

if "GROMPYDOUBLE" in environ:
    isdouble=True
    print("Loading grompy with double precision library")
    c_real = c_double
    if auto_mode&mNumPy:
        N_real = N.float64
    libmd=cdll.LoadLibrary("libmd_d.so")
    libgmx=cdll.LoadLibrary("libgmx_d.so")
else:
    isdouble=False
    print("Loading grompy with single precision library")
    c_real = c_float
    if auto_mode&mNumPy:
        N_real = N.float32
    libmd=cdll.LoadLibrary("libmd.so")
    libgmx=cdll.LoadLibrary("libgmx.so")


#sets it to fs - would need to be fixed in libgmx.
#No oother exported function but parse_args allows to set the time
libgmx.default_time() 
                         
rvec=c_real*3
matrix=c_real*3*3

class GMXctypesError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

 
