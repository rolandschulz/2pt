import numpy as N
import numpy.linalg as LA
from grompy import N_real

#numpy is actually slower for norm and cross if it is for just one vector
def norm(v):
    x,y,z=v
    return sqrt(x*x+y*y+z*z)

def cross(v1, v2):
    assert(v1.shape==(3,) and v2.shape==(3,))
    x0, y0, z0 = v1
    x1, y1, z1 = v2
    return N.array((y0*z1-y1*z0, z0*x1-z1*x0, x0*y1-x1*y0),dtype=N_real)

""" angle between the two input vector """
def angle(x1,x2):
    return N.arccos(N.dot(x1,x2)/norm(x1)/norm(x2))/N.pi*180

""" inertia tensor of list of atoms """
def inertia_tensor(lv,mass):
    t=N.zeros((3,3),dtype=N_real)
    for i in range(len(lv)):
        x,y,z=lv[i]
        t+=mass[i]*[[y*y+z*z,0,0],[-y*x,x*x+z*z,0],[-x*z,-y*z,x*x+y*y]]
    return t                      

""" get the rotation matrix which rotates M2 to (inverse of) M
    both matrices are supposed to be 3x3 matrices and have to be identical besides a rotation
    M should be the pseudo inverse of the coordinates described by M
    Only works in special cases so far!  Left here because much faster and probably can be fixed
    """
def rotation_pinv(M,M2):
    x=N.dot(M,M2)
    x2=cross(x[:,0],x[:,1]) #the 3rd column is undefined because of singularity but we now we can compute it using cross product because it is a rotation
    x2=x2/norm(x2)
    x[2]=x2
    return x

""" computes rotation matrix for 3*N coordinates relative to their COM """
def rotate_matrix(x,c):
    t=N.tensordot(x,c,(0,0))
    u,s,v=LA.svd(t)
    if __debug__:
        R0=N.dot(u,v)
        #R is different whether rotating around COG or COM
        #t=N.tensordot(x-N.mean(x,axis=0),c-N.mean(c,axis=0),(0,0))
        #u,s,v=LA.svd(t)
        #R2=N.dot(u,v)
        #N.testing.assert_almost_equal(R,R2)

        #checking that det and cross method to gurantee det(R)=1 are equaivalen
        R2=N.dot(N.dot(u,N.diag([1,1,LA.det(R0)])),v)
    u[2]=cross(u[0],u[1])
    v[2]=cross(v[0],v[1])
    R=N.dot(u,v)
    if __debug__:
        N.testing.assert_almost_equal(R,R2,decimal=5)
    return R
