from __future__ import print_function, division
import sys
import numpy as N
import numpy.linalg as LA
import scipy.signal
from optparse import OptionParser
from grompy import N_real
from grompy.fio import gmxfile
from grompy.index import rd_index
import grompy.do_fit
from auto_corr3_c import calc_fit_R,compute_omega
#from grompy.tpxio import tpxfile
from grompy import libgmx
from grompy.types  import TRX_NEED_X, TRX_NEED_V, TRX_NEED_F
from math import sqrt
from matrix import *
from itertools import izip

#hardcoded constants - should be improved that those are read from input files
mass=N.array([15.99940,1.00790,1.00790],dtype=N_real).reshape(3,1)  #this is for water. Has to be changed if other molecules should be computed
dt = 0.001  #timestep in ps
ds = 4      #sampling step
nP = 5000   #number sampling points
bCom=True
bNormalize=True

#compute center of mass and rotation autocorrelation for water


"""Computes angular velocity of molecules.
   Input: coordinates, velocities, coordinates along prinzipal axis, principal moments of intertia"""
def compute_omega_py(frm_x,frm_v,ic,I,mass):
    assert(frm_x.shape==frm_v.shape)
    R = N.empty((3,3),N_real)
    w = N.empty((len(frm_x),3))
    for m in range(len(frm_x)):
        #print(m)
        x = frm_x[m]
        v = frm_v[m]

        #calc_fit_R(3,3,mass,ic,x,R)
        #print(mass,ic,x,R)
        #libgmx.calc_fit_R(3,3,mass.ctypes.data,ic.ctypes.data,x.ctypes.data,R.ctypes.data)
        R=rotate_matrix(ic,x)
        #libgmx.calc_fit_R(3,3,N.ones(3,dtype=N_real),ic,x,R)
        #continue
        
        if __debug__:
            #R2=rotate_matrix(ic,x)
            #N.testing.assert_almost_equal(R,R2)

            evec=LA.eigh(inertia_tensor(x,mass))[1]
            M=N.dot(R,evec)   #transformation from R to evec transformation should be reflection at coordinate system planes
            #N.testing.assert_almost_equal(N.dot(N.diag(N.diag(M)),M),N.eye(3),decimal=4) #thus should be diagonal
            N.testing.assert_almost_equal(abs(N.diag(M)),1,decimal=4) #and elements either 1 or -1 (this means it is closely diagonal too)

        L=N.sum(N.cross(x,v)*mass,axis=0)
        w[m]=N.dot(R,L.reshape(3,1)).reshape(3)/I

        if __debug__:
            w2 = cross(x[0],v[0])/N.dot(x[0],x[0]) #*(N.dot(x[i],x[i])*mass[i]); hydrgens to inaccurate
            w2 = N.dot(w2,R.T)
            N.testing.assert_almost_equal(w2[::2],w[m,::2],decimal=0)  # y axis goes through Oxygen

def norm_corr(corr_inp):
   corr = corr_inp.reshape(len(corr_inp),-1)
   #corr = N.sum(corr,axis=1)
   corr = N.mean(corr,axis=1)
   corr/=N.arange(len(corr)*2,len(corr),-1)  #correct for different number of summands (ends are padded in numpy)
   if bNormalize:
       corr/=corr[0]
   return corr


def main():
   parser = OptionParser()
   #parser.add_option("-s",metavar="FILE",help="Reference Structure file",default="topol.tpr")
   parser.add_option("-n",metavar="FILE",help="Index file")
   parser.add_option("-f",metavar="FILE",help="Trajectory with coordinates, velocies, and forces")
   parser.add_option("-x",metavar="FILE",help="Trajectory with coordinates")
   parser.add_option("-v",metavar="FILE",help="Trajectory with velocities")
   parser.add_option("-k",metavar="FILE",help="Trajectory with forces")
   parser.add_option("-o",metavar="FILE",help="Output",default="corr.dat")
   parser.add_option("-r",metavar="FILE",help="Output",default="rot.dat")
   (options, args) = parser.parse_args()

   #open files
   #sf = tpxfile(options.s)

   #read index
   if options.n:
       isize,idx,iname = rd_index(options.n,1)
       print("Using group %s with %s atoms"%(iname[0],isize[0]))
       idx=idx[0]

   if not options.f and not(options.x and options.v and options.k):
       print("Requires either one trajectory containing all input (-f) or three (-x, -v, -f)")
       sys.exit(1);
   if options.f and (options.x or options.v or options.k):
       print("-f cannot be used with x,v, or k")
       sys.exit(1);

   if options.f:
       tf = gmxfile(options.f,flags=TRX_NEED_X|TRX_NEED_V|TRX_NEED_F)
       natoms = tf.natoms
   else:
       xf = gmxfile(options.x,flags=TRX_NEED_X)
       vf = gmxfile(options.v,flags=TRX_NEED_X) #DCD doesn't have V/F so we expect all 3 files have X field
       kf = gmxfile(options.k,flags=TRX_NEED_X)
       print("Warning: NAMD units are expected for input files!")
       natoms = xf.natoms
       tf = izip(xf,vf,kf)
       if natoms!=vf.natoms or natoms!=kf.natoms:
           print("Number of atoms doesn't match between files!")
           sys.exit(1);

   if (natoms%3!=0):
       print("This program only supports water. Number of atoms has to be multiple of 3.");
       sys.exit(1);

   temp_v=[]
   temp_x=[]
   temp_f=[]
   trn_v=N.empty((nP,natoms/3,3,3),dtype=N_real)
   trn_x=N.empty((nP,natoms/3,3,3),dtype=N_real)

   for i,frm in enumerate(tf):
       if not options.f:
           frm_t = frm
           class frm:  #convert tuple into attributes
               x = frm_t[0].x
               v = frm_t[1].x*20.45482706
               f = frm_t[2].x*(100/4.184)
       if options.n:
           frm_x = frm.x[idx]
           frm_v = frm.v[idx]
           frm_f = frm.f[idx]
       else:
           frm_x = frm.x
           frm_v = frm.v
           frm_f = frm.f

       frm_v=frm_v.reshape(-1,3,3)
       frm_x=frm_x.reshape(-1,3,3)
       frm_f=frm_f.reshape(-1,3,3)

       temp_v += [frm_v.copy()]
       temp_x += [frm_x.copy()]
       temp_f += [frm_f.copy()]
       temp_v = temp_v[-3:]
       temp_x = temp_x[-3:]
       temp_f = temp_f[-3:]

       if i%ds!=2: continue  #assumes ds>=3, 
       if __debug__:
           a=temp_x[1]
           b=temp_x[0]+temp_v[1]*dt
           #print(temp_v[1]*dt)
           #print(temp_x[0])
           #idx=N.abs(a-b).max(axis=1).max(axis=1)>.0005
           #print(a.shape,temp_v[1].shape)
           #print(idx)
           #print(a[idx])
           #print(b[idx])
           N.testing.assert_almost_equal(a,b,decimal=2,verbose=True)  #low accuracy dominated by pressure coupling

       i=i//ds
       trn_v[i] = .5*(temp_v[1]+temp_v[2])-dt/16/mass*(temp_f[2]-temp_f[0])
       #trn_v[i] = .5*(temp_v[1]+temp_v[2])
       ##trn_v2 = 2*(trn_v[1:-1]+trn_v[2:])-3/2/dt*(trn_x[2:]-trn_x[:-2])
       trn_x[i] = temp_x[1]

       if i==nP-1: break

       if __debug__:
           a=(temp_v[0]+temp_v[1])/2#trn_v2[n]
           b=trn_v[i]
           #idx=N.abs((a-b))>.4+.2*N.abs(a+b)   #.6 .4  old values: for 2fs: a few,, 1fs: none
           idx=N.abs((a-b))>.6+.4*N.abs(a+b)   #.6 .4  old values: for 2fs: a few,, 1fs: none
           idx2=N.abs((a-b))<.1+.05*N.abs(a+b) #.1 .05 : 75% - 1fs: 91%        
           assert(N.sum(idx)/N.prod(idx2.shape)<1/400)
           #assert(N.sum(idx2)/N.prod(idx2.shape)>.9)
           assert(N.sum(idx2)/N.prod(idx2.shape)>.70)

   assert(i==nP-1)





   N.set_string_function(N.array_repr,False)

   x=trn_x[0][0]  #1st frame, 1st molecule
   x=x-N.sum(x*mass,axis=0)/N.sum(mass)
   I=inertia_tensor(x,mass)
   I,evec=LA.eigh(I)
   print()
   print("I",I)
   ic = N.array(N.dot(x,evec),dtype=N_real) #initial water coordinates (oriented according to I); orientation x: H1->H2, y: O->(H1+H2)/2, z: perpendicular to plane


   print(trn_x.shape,trn_v.shape)




   trn_w=N.empty((len(trn_v),natoms/3,3),dtype=N_real)

   #COM
   trn_x -= N.sum(trn_x*mass,axis=2).reshape(nP,natoms/3,1,3)/N.sum(mass)  #reshape doesn't change shape but reinserts the summed axis
   trn_v_com = trn_v - N.sum(trn_v*mass,axis=2).reshape(nP,natoms/3,1,3)/N.sum(mass)

   for n in range(len(trn_v)):
       trn_w[n] = compute_omega(trn_x[n], trn_v_com[n], ic, I, mass.reshape(3))
       #trn_w[n] = compute_omega_py(trn_x[n], trn_v_com[n], ic, I, mass)



   trn_w=trn_w[:nP]
   corr = N.apply_along_axis(lambda a: scipy.signal.fftconvolve(a,a[::-1],'same'),axis=0,arr=trn_w)
   corr = corr[len(corr)/2:] #remove negative lag
   corr=norm_corr(corr)
   rf=file(options.r,"w")
   map(lambda x:print(x,file=rf),corr)



   #if len(trn_v)%2==1:  #remove last if uneven number of trn-data points
   #
   trn_v=trn_v[:nP]  #TODO: make configurable, switch here for trn/rot



   print(trn_v.shape)
   if bCom:
       #trn_v = trn_v.reshape((-1,natoms/3,3,3)) #group as waters
       trn_v = N.mean(trn_v*mass, axis=2) #average over waters, now shape is: frame, mol, coor
   print(trn_v.shape)
   corr = N.apply_along_axis(lambda a: scipy.signal.fftconvolve(a,a[::-1],'same'),axis=0,arr=trn_v)
   #corr = N.apply_along_axis(lambda a: N.correlate(a,a,'same'),axis=0,arr=trn_v) #slower identical alternative
   print(corr.shape)
   corr = corr[len(corr)/2:] #remove negative lag
   if bCom:
       if not bNormalize:
           corr/=N.sum(mass) #we multiplied trn with mass and because correlation is trn*trn the correlation was mass^2
   else:
       corr = corr.reshape(-1,natoms/3,3,3) #frames, molecule, atoms, xyz; so that we can multiply with mass
       print(corr.shape)
       corr *= mass
   corr=norm_corr(corr)
   print(corr.shape)
   of=file(options.o,"w")
   map(lambda x:print(x,file=of),corr)

if __name__ == "__main__":
    main()
