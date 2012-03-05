import numpy as np
cimport numpy as np
cimport cython
#from libc.math cimport sqrtf
import numpy.linalg as LA
from matrix import cross, inertia_tensor
ctypedef float real

cdef extern from "types/simple.h":
    pass

cdef extern from "math.h":
    float sqrtf(float x)

cdef unsigned int XX=0,YY=1,ZZ=2

#def cprod(np.ndarray[np.float32_t,ndim=1] a not None,np.ndarray[np.float32_t,ndim=1] b not None,np.ndarray[np.float32_t,ndim=1] c not None):
#    c_cprod(<real*>a.data,<real*>b.data,<real*>c.data)

cdef extern from "vec.h":
    ctypedef float *rvec
    ctypedef float *matrix
    void cprod (rvec a, rvec b, rvec c)
    void clear_rvec (rvec a)
    void mvmul(matrix a,rvec src,rvec dest) 

cdef extern from "do_fit.h":
    ctypedef float rvec
    ctypedef float *matrix
    void c_calc_fit_R "calc_fit_R" (int ndim,int natoms,rvec *w_rls,rvec *xp,rvec *x, matrix R)
    
def calc_fit_R(ndim, natoms, np.ndarray[np.float32_t,ndim=2] w_rls, np.ndarray[np.float32_t,ndim=2] xp, np.ndarray[np.float32_t,ndim=2] x, np.ndarray[np.float32_t,ndim=2] R):
    c_calc_fit_R(ndim, natoms, <float*> w_rls.data, <rvec*> xp.data, <rvec*> x.data, <rvec*> R.data)

N_real = np.float32
"""Computes angular velocity of molecules.
   Input: coordinates, velocities, coordinates along prinzipal axis, principal moments of intertia"""
def compute_omega(np.ndarray[np.float32_t,ndim=3] frm_x,
                  np.ndarray[np.float32_t,ndim=3] frm_v,
                  np.ndarray[np.float32_t,ndim=2] ic,
                  np.ndarray[np.float32_t,ndim=1] I,
                  np.ndarray[np.float32_t,ndim=1] mass):

    cdef np.ndarray[np.float32_t,ndim=2] R = np.empty((3,3),N_real)
    cdef np.ndarray[np.float32_t,ndim=2] w = np.empty((len(frm_x),3),N_real)
    cdef np.ndarray[np.float32_t,ndim=2] x,v
    cdef np.ndarray[np.float32_t,ndim=1] L = np.empty(3,N_real)
    cdef np.ndarray[np.float32_t,ndim=1] t = np.empty(3,N_real)
    cdef unsigned int m,n,i
    for m in range(len(frm_x)):
        c_calc_fit_R(3,3,&mass[0],<rvec*>&ic[0,0],<rvec*>&frm_x[m,0,0],<rvec*>&R[0,0])
        
        if __debug__:
            x = frm_x[m]
            v = frm_v[m]
            evec=LA.eigh(inertia_tensor(x,mass.reshape(3,1)))[1]
            M=np.dot(R,evec)   #transformation from R to evec transformation should be reflection at coordinate system planes
            np.testing.assert_almost_equal(abs(np.diag(M)),1,decimal=4) #and elements either 1 or -1 (this means it is closely diagonal too)

        clear_rvec(&L[0])
        for n in range(3): #over atoms
            cprod(&frm_x[m,n,0],&frm_v[m,n,0],&t[0])  #t=frm_x x frm_v
            for i in range(3): L[i]+=t[i]*mass[n]
            #L+=t*mass[n]   #implicit over coordinates
        #L=np.sum(np.cross(x,v)*mass.reshape(3,1),axis=0)

        mvmul(<rvec*>&R[0,0],&L[0],&w[m,0]) #w=RL
        for i in range(3): w[m,i]/=sqrtf(I[i]) #not actually w (should be devided by I not sqrt(I) but we need w*sqrt(I)
        #w[m]=np.dot(R,L.reshape(3,1)).reshape(3)/I

        if __debug__:
            w2 = cross(x[0],v[0])/np.dot(x[0],x[0]) #*(np.dot(x[i],x[i])*mass[i]); hydrgens to inaccurate
            w2 = np.dot(w2,R.T)
            np.testing.assert_almost_equal(w2[::2],w[m,::2],decimal=0)  # y axis goes through Oxygen
    return w

    
