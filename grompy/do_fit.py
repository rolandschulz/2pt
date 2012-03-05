from ctypes import c_int,POINTER
from grompy import libgmx,matrix
from grompy import c_real,auto_mode,mNumPy

libgmx.calc_similar_ind.restype = c_real
if mNumPy:
    import numpy as N
    from grompy import N_real
    from numpy.ctypeslib import ndpointer
        

def rmsdev(natoms,indices,masses,atomsvec1,atomsvec2):
    return libgmx.calc_similar_ind(c_int(0),natoms,indices,masses,atomsvec1,atomsvec2)

def rhodev(natoms,indices,masses,atomsvec1,atomsvec2):
    return libgmx.calc_similar_ind(c_int(1),natoms,indices,masses,atomsvec1,atomsvec2)


def calc_fit_R(weights,atomsvec1,atomsvec2,natoms=None,dimensions=3,mtx=None):
    if mNumPy:
        if natoms==None:
            if not len(weights)==len(atomsvec1)==len(atomsvec2):
                raise TypeError("Array dimension don't match")
            natoms=len(weights)
        if mtx==None: mtx = N.empty((3,3),N_real)
    else:
        if mtx==None: mtx = matrix()

    libgmx.calc_fit_R(dimensions,natoms,weights,atomsvec1,atomsvec2,mtx)
    return mtx
if mNumPy: libgmx.calc_fit_R.argtypes=[c_int,c_int,ndpointer(dtype=N_real),ndpointer(ndim=2,dtype=N_real),ndpointer(ndim=2,dtype=N_real),ndpointer(ndim=2,dtype=N_real)]


def do_fit(weights,atomsvec1,atomsvec2,natoms=None,dimensions=3):
    if mNumPy:
        if natoms==None:
            if not len(weights)==len(atomsvec1)==len(atomsvec2):
                raise TypeError("Array dimension don't match")
            natoms=len(weights)
    libgmx.do_fit_ndim(dimensions,natoms,weights,atomsvec1,atomsvec2)
if mNumPy: libgmx.do_fit_ndim.argtypes=[c_int,c_int,ndpointer(dtype=N_real),ndpointer(ndim=2,dtype=N_real),ndpointer(ndim=2,dtype=N_real)]
    
def reset_x(ncmatoms,cmindices,atomsvec,masses,nrtoreset,dimensions=3,rindices=None):
    if not rindices:
        rindices=POINTER(c_int)()
    libgmx.reset_x_ndim(dimensions,ncmatoms,cmindices,nrtoreset,rindices,atomsvec,masses)
    
