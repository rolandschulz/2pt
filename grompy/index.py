from ctypes import *
from grompy import libgmx,auto_mode,mNumPy

if auto_mode&mNumPy:
    import numpy as N

def rd_index(statfile,ngroups,isizep=None,indexp=None,grpnamesp=None):
    filename=c_char_p(statfile)
    if not isizep:
        isize=(c_int*ngroups)()
        isizep=pointer(isize)
    if not indexp:
        index=(POINTER(c_int)*ngroups)()
        indexp=pointer(index)
    if not grpnamesp:  
        grpnames=(c_char_p*ngroups)()
        grpnamesp=pointer(grpnames)
    
    libgmx.rd_index(filename,ngroups,isizep,indexp,grpnamesp)
    if auto_mode&mNumPy:
        nindex = [0]*ngroups
        for i in range(ngroups):
            tmp = (c_int*isize[i]).from_address(addressof(index[i].contents))
            nindex[i] = N.frombuffer(tmp,dtype=N.int32) #count=3*self.f.natoms
        index = nindex  
    return isize,index,grpnames
  
def get_index(atoms,ngroups,filename=c_char_p(),isizep=None,indexp=None,grpnamesp=None):
    if type(filename)==str:
        filename=c_char_p(filename)

    if not isizep:
        isize=(c_int*ngroups)()
        isizep=pointer(isize)
    if not indexp:
        index=(POINTER(c_int)*ngroups)()
        indexp=pointer(index)
    if not grpnamesp:  
        grpnames=(c_char_p*ngroups)()
        grpnamesp=pointer(grpnames)

    libgmx.get_index(byref(atoms),filename,ngroups,isizep,indexp,grpnamesp)
    if auto_mode&mNumPy:
        nindex = [0]*ngroups
        for i in range(ngroups):
            tmp = (c_int*isize[i]).from_address(addressof(index[i].contents))
            nindex[i] = N.frombuffer(tmp,dtype=N.int32) #count=3*self.f.natoms
        index = nindex
    return isize,index,grpnames

class index:
    pass
        
