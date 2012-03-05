from ctypes import c_char_p, c_int, c_char,c_ubyte, c_ushort,Structure, Union, POINTER
from grompy import c_real, rvec, matrix
#typedef int atom_id;
atom_id = c_int
#typedef int t_functype;
t_functype = c_int
#typedef atom_id t_iatom;
t_iatom = atom_id

DIM = 3
MAXFORCEPARAM = 12
NR_RBDIHS = 6

#enum {
(
  F_BONDS,
  F_G96BONDS,
  F_MORSE,
  F_CUBICBONDS,
  F_CONNBONDS,
  F_HARMONIC,
  F_FENEBONDS,
  F_TABBONDS,
  F_TABBONDSNC,
  F_ANGLES,
  F_G96ANGLES,
  F_CROSS_BOND_BONDS,
  F_CROSS_BOND_ANGLES,
  F_UREY_BRADLEY,
  F_QUARTIC_ANGLES,
  F_TABANGLES,
  F_PDIHS,
  F_RBDIHS,
  F_FOURDIHS,
  F_IDIHS,
  F_PIDIHS,
  F_TABDIHS,
  F_LJ14,
  F_COUL14,
  F_LJC14_Q,
  F_LJC_PAIRS_NB,
  F_LJ,
  F_BHAM,
  F_LJ_LR,
  F_BHAM_LR,
  F_DISPCORR,
  F_COUL_SR,
  F_COUL_LR,
  F_RF_EXCL,
  F_COUL_RECIP,
  F_DPD,
  F_POLARIZATION,
  F_WATER_POL,
  F_THOLE_POL,
  F_POSRES,
  F_DISRES,
  F_DISRESVIOL,
  F_ORIRES,
  F_ORIRESDEV,
  F_ANGRES,
  F_ANGRESZ,
  F_DIHRES,
  F_DIHRESVIOL,
  F_CONSTR,
  F_CONSTRNC,
  F_SETTLE,
  F_VSITE2,
  F_VSITE3,
  F_VSITE3FD,
  F_VSITE3FAD,
  F_VSITE3OUT,
  F_VSITE4FD,
  F_VSITE4FDN,
  F_VSITEN,
  F_COM_PULL,
  F_EQM,
  F_EPOT,
  F_EKIN,
  F_ETOT,
  F_ECONSERVED,
  F_TEMP,
  F_PRES,
  F_DVDL,
  F_DKDL,
  F_DGDL_CON,
  F_NRE) = list(map(c_int, range(71)))
#/* This number is for the total number of energies	*/
#};

#/* For trxframe.flags, used in trxframe read routines.
 #* When a READ flag is set, the field will be read when present,
 #* but a frame might be returned which does not contain the field.
 #* When a NEED flag is set, frames not containing the field will be skipped.
 #*/
TRX_READ_X  =  (1<<0)
TRX_NEED_X  =  (1<<1)
TRX_READ_V  =  (1<<2)
TRX_NEED_V  =  (1<<3)
TRX_READ_F  =  (1<<4)
TRX_NEED_F  =  (1<<5)
#/* Useful for reading natoms from a trajectory without skipping */
TRX_DONT_SKIP = (1<<6)




class t_block(Structure):
    _fields_ = [("nr", c_int),
	    ("index", POINTER(atom_id)),
	    ("nalloc_index", c_int)]
#  int nr;			/* The number of blocks			*/
#  atom_id *index;		/* Array of indices (dim: nr+1) 	*/
#  int nalloc_index;             /* The allocation size for index        */


class t_blocka(Structure):
    _fields_ = [("nr", c_int),
	    ("index", POINTER(atom_id)),
	    ("nra", c_int),
	    ("a", POINTER(atom_id)),
	    ("nalloc_index", c_int),
	    ("nalloc_a", c_int)]
#  int nr;			/* The number of blocks			*/
#  atom_id *index;		/* Array of indices in a (dim: nr+1)	*/
#  int nra;			/* The number of atoms 			*/
#  atom_id *a;			/* Array of atom numbers in each group 	*/
#				/* (dim: nra)				*/
#				/* Block i (0<=i<nr) runs from		*/
#				/* index[i] to index[i+1]-1. There will */
#				/* allways be an extra entry in index	*/
#				/* to terminate the table		*/
#  int nalloc_index;             /* The allocation size for index        */
#  int nalloc_a;                 /* The allocation size for a            */

class t_symbuf(Structure):
    pass
  
t_symbuf._fields_ = [("bufsize", c_int),
	    ("buf", POINTER(c_char_p)),
	    ("next", POINTER(t_symbuf))]
#  int bufsize;
#  char **buf;
#  struct symbuf *next;


class t_symtab(Structure):
    _fields_ = [("nr", c_int),
	    ("symbuf", POINTER(t_symbuf))]
#  int      nr;
#  t_symbuf *symbuf;


class t_iparams(Union):
    class bham(Structure):
        _fields_ = [("a", c_real),
	      ("b", c_real),
	      ("c", c_real)]

    class harmonic(Structure):
        _fields_ = [("rA", c_real),
	      ("krA", c_real),
	      ("rb", c_real),
	      ("krB", c_real)]

    class cubic(Structure):
        _fields_ = [("b0", c_real),
	      ("kb", c_real),
	      ("kcub", c_real)]

    class fene(Structure):
        _fields_ = [("bm", c_real),
	      ("kb", c_real)]

    class cross_bb(Structure):
        _fields_ = [("r1e", c_real),
	      ("r2e", c_real),
	      ("krr", c_real)]

    class cross_ba(Structure):
        _fields_ = [("r1e", c_real),
	      ("r2e", c_real),
	      ("r3e", c_real),
	      ("krt", c_real)]

    class u_b(Structure):
        _fields_ = [("theta", c_real),
	      ("ktheta", c_real),
	      ("r13", c_real),
	      ("kUB", c_real)]

    class qangle(Structure):
        _fields_ = [("theta", c_real),
	      ("c", c_real * 5)]

    class polarize(Structure):
        _fields_ = [("alpha", c_real)]

    class wpol(Structure):
        _fields_ = [("al_x", c_real),
	      ("al_y", c_real),
	      ("al_z", c_real),
	      ("rOH", c_real),
	      ("rHH", c_real),
	      ("rOD", c_real)]

    class thole(Structure):
        _fields_ = [("a", c_real),
	      ("alpha1", c_real),
	      ("alpha2", c_real),
	      ("rfac", c_real)]

    class lj(Structure):
        _fields_ = [("c6", c_real),
	      ("c12", c_real)]

    class lj14(Structure):
        _fields_ = [("c6A", c_real),
	      ("c12A", c_real),
	      ("c6B", c_real),
	      ("c12B", c_real)]

    class ljc14(Structure):
        _fields_ = [("fqq", c_real),
	      ("qi", c_real),
	      ("qj", c_real),
	      ("c6", c_real),
	      ("c12", c_real)]

    class ljcnb(Structure):
        _fields_ = [("qi", c_real),
	      ("qj", c_real),
	      ("c6", c_real),
	      ("c12", c_real)]

    class pdihs(Structure):
        _fields_ = [("phiA", c_real),
	      ("cpA", c_real),
	      ("mult", c_int),
	      ("phiB", c_real),
	      ("cpB", c_real)]
  
    class constr(Structure):
        _fields_ = [("dA", c_real),
	      ("dB", c_real)]

    class settle(Structure):
        _fields_ = [("doh", c_real), 	      
	      ("dhh", c_real)]

    class morse(Structure):
        _fields_ = [("b0", c_real),
	      ("cb", c_real),
	      ("beta", c_real)]
  
    class posres(Structure):
        _fields_ = [("pos0A", c_real * DIM),
	      ("fcA", c_real * DIM),
	      ("pos0B", c_real * DIM),
	      ("fcB", c_real * DIM)]

    class rbdihs(Structure):
        _fields_ = [("rbcA", c_real * NR_RBDIHS),
	      ("rbcB", c_real * NR_RBDIHS)]

    class vsite(Structure):
        _fields_ = [("a", c_real),
	      ("b", c_real),
	      ("c", c_real),
	      ("d", c_real),
	      ("e", c_real),
	      ("f", c_real)]

    class vsiten(Structure):
        _fields_ = [("n", c_int),
	      ("a", c_real)]

    class disres(Structure):
        _fields_ = [("low", c_real),
	      ("up1", c_real),
	      ("up2", c_real),
	      ("kfac", c_real),
	      ("type", c_int),
	      ("label", c_int)]

    class dihres(Structure):
        _fields_ = [("phi", c_real),
	      ("dphi", c_real),
	      ("kfac", c_real),
	      ("label", c_int),
	      ("power", c_int)]

    class orires(Structure):
        _fields_ = [("ex", c_int),
	      ("power", c_int),
	      ("label", c_int),
	      ("c", c_real),
	      ("obs", c_real),
	      ("kfac", c_real)]

    class tab(Structure):
        _fields_ = [("table", c_int),
	      ("kA", c_real),
	      ("kB", c_real)]

    class generic(Structure):
        _fields_ = [("buf", c_real * MAXFORCEPARAM)]


#  /* Some parameters have A and B values for free energy calculations.
#   * The B values are not used for regular simulations of course.
#   * Free Energy for nonbondeds can be computed by changing the atom type.
#   * The harmonic type is used for all harmonic potentials:
#   * bonds, angles and improper dihedrals
#   */
#  struct {real a,b,c;	                                   } bham;
#  struct {real rA,krA,rB,krB;           	           } harmonic; 
#  /* No free energy supported for cubic bonds, FENE, WPOL or cross terms */ 
#  struct {real b0,kb,kcub;                                 } cubic;
#  struct {real bm,kb;                                      } fene;
#  struct {real r1e,r2e,krr;                                } cross_bb;
#  struct {real r1e,r2e,r3e,krt;                            } cross_ba;
#  struct {real theta,ktheta,r13,kUB;                       } u_b;
#  struct {real theta,c[5];                                 } qangle; 
#  struct {real alpha;                                      } polarize;
#  struct {real al_x,al_y,al_z,rOH,rHH,rOD;                 } wpol;
#  struct {real a,alpha1,alpha2,rfac;                       } thole;
#  struct {real c6,c12;				           } lj;
#  struct {real c6A,c12A,c6B,c12B;		           } lj14;
#  struct {real fqq,qi,qj,c6,c12;	                   } ljc14;
#  struct {real qi,qj,c6,c12;		                   } ljcnb;
#  /* Proper dihedrals can not have different multiplicity when
#   * doing free energy calculations, because the potential would not
#   * be periodic anymore.
#   */ 
#  struct {real phiA,cpA;int mult;real phiB,cpB;            } pdihs;
#  struct {real dA,dB;		        	           } constr;
#  /* Settle can not be used for Free energy calculations of water bond geometry.
#   * Use shake (or lincs) instead if you have to change the water bonds.
#   */
#  struct {real doh,dhh;                                   } settle;
#  /* No free energy supported for morse bonds */ 
#  struct {real b0,cb,beta;                        	  } morse;
#  struct {real pos0A[DIM],fcA[DIM],pos0B[DIM],fcB[DIM];   } posres;
#  struct {real rbcA[NR_RBDIHS], rbcB[NR_RBDIHS];          } rbdihs;
#  struct {real a,b,c,d,e,f;                               } vsite;   
#  struct {int  n; real a;                                 } vsiten;   
#  struct {real low,up1,up2,kfac;int type,label;           } disres; 
#  struct {real phi,dphi,kfac;int label,power;             } dihres;  
#  struct {int  ex,power,label; real c,obs,kfac;           } orires;
#  struct {int  table;real kA;real kB;                     } tab;
#  struct {real buf[MAXFORCEPARAM];	  	          } generic; /* Conversion */





class t_ilist(Structure):
    _fields_ = [("nr", c_int),
	    ("nr_nonperturbed", c_int),
	    ("iatoms", POINTER(t_iatom)),
	    ("nalloc", c_int)]
#  int nr;
#  int nr_nonperturbed;
#  t_iatom *iatoms;
#  int nalloc;




class t_idef(Structure):
    _fields_ = [("ntypes", c_int),
	    ("atnr", c_int),
	    ("t_functype", POINTER(t_functype)),
	    ("t_iparams", POINTER(t_iparams)),
	    ("fudgeQQ", c_real),
	    ("iparams_posres", POINTER(t_iparams)),
	    ("iparams_posres_nalloc", c_int),
	    ("il", t_ilist * F_NRE.value),
	    ("ilsort", c_int)
	    ]

#  int ntypes;
#  int atnr;
#  t_functype *functype;
#  t_iparams  *iparams;
#  real fudgeQQ;
#  t_iparams  *iparams_posres;
#  int iparams_posres_nalloc;
#  t_ilist il[F_NRE];
#  int ilsort;


class t_atom(Structure):
    _fields_ = [("m", c_real),
	    ("q", c_real),
	    ("mB", c_real),
	    ("qB", c_real),
	    ("type", c_ushort),
	    ("typeB", c_ushort),
	    ("ptype", c_int),
	    ("resnr", c_int),
	    ("atomnumber", c_int),
	    ("chain", c_ubyte)]
#  real 		m,q;		/* Mass and charge			*/
#  real 		mB,qB;		/* Mass and charge for Free Energy calc */
#  unsigned short type;		/* Atom type				*/
#  unsigned short typeB;		/* Atom type for Free Energy calc	*/
#  int           ptype;		/* Particle type			*/
#  int 		resnr;		/* Residue number			*/
#  int           atomnumber;     /* Atomic Number or NOTSET              */
#  unsigned char chain;          /* chain identifier                     */

class t_pdbinfo(Structure):
    _fields_ = [("type", c_int),
	    ("atomnr", c_int),
	    ("altloc", c_char),
	    ("atomnm", c_char * 6),
	    ("pdbresnr", c_char * 6),
	    ("occup", c_real),
	    ("bfac", c_real),
	    ("bAnisotropic", c_int),
	    ("uij", c_int * 6)
	    ]
#  int  type;                    /* PDB record name                      */
#  int  atomnr;                  /* PDB atom number                      */
#  char altloc;                  /* Alternate location indicator         */
#  char atomnm[6];               /* True atom name including spaces      */
#  char pdbresnr[6];             /* PDB res number                       */
#  real occup;                   /* Occupancy                            */
#  real bfac;                    /* B-factor                             */
#  bool bAnisotropic;            /* (an)isotropic switch                 */
#  int  uij[6];                  /* Anisotropic B-factor                 */

class t_grps(Structure):
    _fields_ = [("nr", c_int),
	    ("nm_ind", POINTER(c_int))
	    ]

#  int  nr;			/* Number of different groups		*/
#  int  *nm_ind;                 /* Index in the group names             */

class t_atoms(Structure):
    _fields_ = [("nr", c_int),
	    ("atom", POINTER(t_atom)),
	    ("atomname", POINTER(POINTER(c_char_p))),
	    ("atomtype", POINTER(POINTER(c_char_p))),
	    ("atomtypeB", POINTER(POINTER(c_char_p))),
	    ("nres", c_int),
	    ("resname", POINTER(POINTER(c_char_p))),
	    ("pdbinfo", c_char_p)
	    ]
#  int           nr;             /* Nr of atoms                          */
#  t_atom	*atom;		/* Array of atoms (dim: nr)		*/
#				/* The following entries will not 	*/
#				/* allways be used (nres==0)	 	*/
#  char		***atomname;	/* Array of pointers to atom name	*/
#				/* use: (*(atomname[i]))		*/
#  char		***atomtype;	/* Array of pointers to atom types	*/
#				/* use: (*(atomtype[i]))		*/
#  char		***atomtypeB;	/* Array of pointers to B atom types	*/
#				/* use: (*(atomtypeB[i]))		*/
#  int		nres;		/* Nr of residue names			*/
#  char		***resname; 	/* Array of pointers to residue names 	*/
#				/* use: (*(resname[i]))	       	*/
#  t_pdbinfo     *pdbinfo;       /* PDB Information, such as aniso. Bfac */


class t_atomtypes(Structure):
    _fields_ = [("nr", c_int),
	    ("radius", POINTER(c_real)),
	    ("vol", POINTER(c_real)),
	    ("surftens", POINTER(c_real)),
	    ("atomnumber", POINTER(c_int))
	    ]
#  int           nr;              /* number of atomtypes                     */
#  real         *radius;         /* GBSA radius for each atomtype        */
#  real         *vol;            /* GBSA efective volume for each atomtype   */
#  real         *surftens;       /* implicit solvent surftens for each atomtype */
#  int          *atomnumber;     /* Atomic number, used for QM/MM */


class t_topology(Structure):
    _fields_ = [("name", POINTER(c_char_p)),
	      ("idef", t_idef),
	      ("atoms", t_atoms),
	      ("atomtypes", t_atomtypes),
	      ("cgs", t_block),
	      ("mols", t_block),
	      ("excls", t_blocka),
	      ("symtab", t_symtab)]

#  char  	**name;		/* Name of the topology	       	        */
#  t_idef	idef;		/* The interaction function definition	*/
#  t_atoms	atoms;		/* The atoms		       	        */
#  t_atomtypes   atomtypes;      /* Atomtype properties                  */
#  t_block       cgs;            /* The charge groups                    */
#  t_block       mols;           /* The molecules                        */
#  t_blocka      excls;          /* The exclusions                       */
#  t_symtab	symtab;		/* The symbol table			*/


#/* The bools indicate whether a field was read from the trajectory.
 #* Do not try to use a pointer when its bool is FALSE, as memory might
 #* not be allocated.
 #*/ 
 
class t_trxframe(Structure):
  _fields_ = [(   "flags", c_int),           # flags for read_first/next_frame 
              (  "not_ok", c_int),           # integrity flags (see statutil.h 
              ( "bDouble", c_int),           # Double precision?               
              (  "natoms", c_int),           # number of atoms (atoms, x, v, f)
              (      "t0", c_real),          # time of the first frame, needed 
                                           # for skipping frames with -dt    
              (     "tpf", c_real),          # time of the previous frame, not 
                                           # the read, but real file frames  
              (    "tppf", c_real),          # time of two frames ago          
                                           # tpf and tppf are needed to      
                                           # correct rounding errors for -e  
              (  "bTitle", c_int),           # 
              (   "title", c_char_p),        # title of the frame         
              (   "bStep", c_int),
              (    "step", c_int),           # MD step number                  
              (   "bTime", c_int),
              (    "time", c_real),          # time of the frame               
              ( "bLambda", c_int),
              (  "Lambda", c_real),          # free energy perturbation lambda 
              (  "bAtoms", c_int),
              (   "atoms", POINTER(t_atoms)),# atoms struct (natoms)        
              (   "bPrec", c_int),
              (    "prec", c_real),          # precision of x, fraction of 1 nm
              (      "bX", c_int),
              (       "x", POINTER(rvec)),   # coordinates (natoms)            
              (      "bV", c_int),
              (       "v", POINTER(rvec)),   # velocities (natoms)             
              (      "bF", c_int),
              (       "f", POINTER(rvec)),   # forces (natoms)                 
              (    "bBox", c_int),
              (     "box", matrix),          # the 3 box vectors               
              (    "bPBC", c_int),
              (    "ePBC", c_int)]           # the type of pbc                 