from ctypes import c_char_p,POINTER,c_int,byref,addressof
from copy   import copy, deepcopy
from grompy import libgmx,rvec,matrix
from grompy import c_real, auto_mode, mNumPy
from grompy.types  import TRX_READ_X,t_trxframe

if auto_mode==mNumPy:
	from grompy.numpy_util import maprvec
	
   
def open_xtc(filename,mode):
	fileptr=c_char_p(filename)
	modeptr=c_char_p(mode)
	return libgmx.open_xtc(fileptr,modeptr)
	
def close_xtc(handle):
	libgmx.close_xtc(handle)

def read_first_xtc(handle):
	natoms = c_int()
	step = c_int()
	time = c_real()
	prec = c_real()
	bOK = c_int()
	box =  matrix()
	xpointer = POINTER(rvec)()
	ret=libgmx.read_first_xtc(handle,byref(natoms),byref(step),byref(time),box,byref(xpointer),byref(prec),byref(bOK))
#	if ret!=1:
#		raise GMXctypesError,"read_first_xtc did not return 1"
	return natoms.value,step.value,time.value,box,xpointer,prec.value,bOK.value,ret
	
def read_next_xtc(handle, natoms, xpointer):
	step = c_int()
	time = c_real()
	prec = c_real()
	bOK = c_int()
	box =  matrix()
	ret=libgmx.read_next_xtc(handle,natoms,byref(step),byref(time),box,xpointer,byref(prec),byref(bOK))
#	if ret!=1:
#		raise GMXctypesError,"read_next_xtc did not return 1"

	return step.value,time.value,box,prec.value,bOK.value,ret

def write_xtc(handle,natoms,step,time,box,xpointer,prec):
	ret=libgmx.write_xtc(handle,natoms,step,c_real(time),box,xpointer,c_real(prec))
	return ret
#	if ret!=1:
#		raise GMXctypesError,"write_xtc did not return 1"

class frame:
    pass

class xtcfile:
    def __init__(self,filename):
        self.fh = open_xtc(filename,"r")
    def __iter__(self):
        n = frame()
        self.natoms,n.step,n.time,n.box,n.x,n.prec,bOK,ret = read_first_xtc(self.fh)
        while bOK and n.prec!=0: #HACK: prec should normally never be 0 and is 0 at EOF
          yield n
          l = n
          n = frame()
          n.x=l.x   #reuse memory for coordinates
          del l.x   #make sure the old coordinates are not used by clients (because they get overwritten)
          n.step,n.time,n.box,n.prec,bOK,ret = read_next_xtc(self.fh,self.natoms,n.x)

            
class gmxframe(frame,object):
    def __init__(self, o=None):  #before Python 2.6 copy.copy is broken for Structure, 
                                 #so we need an unphytonic copy constructor
        self.fr = t_trxframe()
        if o!=None:
          attr = zip(*t_trxframe._fields_)[0]
          for i in attr:            #not all needed but coping all is safe
              setattr(self.fr,i,getattr(o.fr,i))
          if auto_mode==mNumPy:
              if self.fr.bX: self.x = o.x 
              if self.fr.bV: self.v = o.v 
              if self.fr.bF: self.f = o.f
        self.__invalid = False

    def invalidate(self,x):
        """removes array attributes - needed after update - e.g. next frame"""
        self.__invalid = x
        if auto_mode==mNumPy and x:
            if self.fr.bX: del self.x
            if self.fr.bV: del self.v
            if self.fr.bF: del self.f
    invalid=property(None,invalidate)

    def __getattr__(self,name):
        if not self.__invalid or not name in ["x","v","f","atoms"]:#don't allow access to arrays after update
          return getattr(self.fr,name)
        else:
          raise AttributeError(name)

class gmxfile:
    def __init__(self,filename,flags=TRX_READ_X):
        self.status = c_int()
        self.f = gmxframe()
        self.ret = libgmx.read_first_frame(byref(self.status),filename,byref(self.f.fr),flags)
	if not self.ret:
		raise IOError("Failed to open %s. File not found or not a valid GROMACS trajectory."%(filename,))
        self.natoms = self.f.natoms
        if auto_mode==mNumPy:
            if self.f.bX: self.f.x = maprvec(self.f.fr.x,self.natoms)
            if self.f.bV: self.f.v = maprvec(self.f.fr.v,self.natoms)
            if self.f.bF: self.f.f = maprvec(self.f.fr.f,self.natoms)

    #for efficienty we allocate large frame arrays only once and reuse between iterations. We invalidate the last iteration object to prohibit access (because it would return the the new data not the original data)
    def __iter__(self):
        while self.ret:
            c = gmxframe(self.f)  #copy of frame object (to make sure user is not changing our copy)
            yield self.f
            self.f.invalid = True #old copy is now invalid (shared date is being overwritten)
            self.f = c  #point to copy
            self.ret = libgmx.read_next_frame(self.status,byref(self.f.fr))


