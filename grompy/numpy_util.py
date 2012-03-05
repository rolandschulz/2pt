from ctypes import *
from grompy import *
from grompy.types import t_topology,t_atom,t_atoms

if auto_mode==mNumPy:
    import numpy as N
    def maprvec(i,natoms):
        """convert POINTER(rvec) to numpy 2d arry"""
        #using addess, because the straig forward from_buffer is new in 2.6
        x = (rvec*natoms).from_address(addressof(i.contents)) 
        x = N.frombuffer(x,dtype=N_real) #count=3*self.f.natoms
        x.shape = (natoms,3)
        return x

    def ctype_to_numpy(t):   #TODO: add auto alignment
        _ctypes_to_numpy = {
            c_real: N_real, 
            c_char : N.int8,
            c_wchar : N.int16,
            c_byte : N.int8,
            c_ubyte : N.uint8,
            c_short : N.int16,
            c_ushort : N.uint16,
            c_int : N.int32,
            c_uint : N.uint32,
            c_long : N.dtype("i"+str(sizeof(c_long))), 
            c_ulong : N.dtype("u"+str(sizeof(c_long))),
            c_float : N.float32,
            c_double : N.float64
            }
        ty = [(_x,_ctypes_to_numpy[_y]) for _x,_y in t]
        return ty


    def mapatom(a,natoms):
        a = (t_atom*natoms).from_address(addressof(a.contents))
        ty = ctype_to_numpy(t_atom._fields_)
        ty +=[('',N.dtype('V3'))]  #add alignment bytes, would be nice to have a function finding places for required alignment from offset and size and adds it automatically
        return N.frombuffer(a,ty)  #why no aligned option?
