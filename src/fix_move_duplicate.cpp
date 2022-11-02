
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_move_duplicate.h"
#include "compute.h"
#include "domain.h"
#include "modify.h"
#include "group.h"
#include "atom.h"
#include "input.h"
#include "update.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include "comm.h"




using namespace LAMMPS_NS;
using namespace FixConst;

enum{COMPUTE};

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

FixMoveDupl::FixMoveDupl(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), displace(NULL),locdisplace(NULL)
{

int i;

create_attribute=1;
restart_pbc = 1;
nevery=1;

not_found_value=1;

int expand = 0;
char **earg;

nvalues=1;
nvalues = utils::expand_args(FLERR,nvalues,&arg[3],1,earg,lmp);

if (earg != &arg[3]) expand = 1;
arg = earg;

  which = new int[nvalues];
  argindex = new int[nvalues];
  ids = new char*[nvalues];
  value2index = new int[nvalues];

  printf("NVALUES1 %d \n",nvalues);

  for (int i = 0; i < nvalues; i++) {
    if (strncmp(arg[i],"c_",2) == 0 ) {
          if (arg[i][0] == 'c') which[i] = COMPUTE;

          int n = strlen(arg[i]);
          char *suffix = new char[n];
          strcpy(suffix,&arg[i][2]);

// The C library function char *strchr(const char *str, int c) searches for the first occurrence of the character c (an unsigned char) in the string pointed to by the argument str.

// The C library function int atoi(const char *str) converts the string argument str to an integer (type int).

         char *ptr = strchr(suffix,'[');
         if (ptr) {
            if (suffix[strlen(suffix)-1] != ']')
                       error->all(FLERR,"Illegal fix move/dupl command");
            argindex[i] = atoi(ptr+1);

            *ptr = '\0';
         } else argindex[i] = 0;

         n = strlen(suffix) + 1;
         ids[i] = new char[n];
         strcpy(ids[i],suffix);
         delete [] suffix; 
         not_found_value=0;   
    }else error->all(FLERR,"Illegal move/dupl command");
  }

if(not_found_value)
    error->all(FLERR,"Illegal fix move/dupl command");

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }

for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix move/dupl does not exist");
      if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
        error->all(FLERR,"Fix ave/atom compute does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_peratom_cols)
        error->all(FLERR,"Fix ave/atom compute array is accessed out-of-range");
        }
    }

  // this fix produces either a per-atom vector or array

  peratom_flag = 1;
  size_peratom_cols = nvalues;
  nmax = atom->natoms;


  buff_size=(nmax+1)*nvalues;



// fpl = fopen("dupl.log","w");

/*
  // perform initial allocation of atom-based array
  // register with Atom class

  xoriginal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // xoriginal = initial unwrapped positions of atoms

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  time_origin = update->ntimestep;*/
}

/* ---------------------------------------------------------------------- */

FixMoveDupl::~FixMoveDupl()
{

  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;

  memory->destroy(locdisplace);
  memory->destroy(displace);
}

/* ---------------------------------------------------------------------- */

void FixMoveDupl::allocate(){

 memory->create(locdisplace,nmax+1,nvalues,"fixmovedupl:locdisplace");
 memory->create(displace,nmax+1,nvalues,"fixmovedupl:displace");

}

/* ---------------------------------------------------------------------- */

int FixMoveDupl::setmask()
{

  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMoveDupl::init()
{


  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix move/dupl does not exist");
      value2index[m] = icompute;
        }
    }

  allocate();

  for (int i = 0; i < nmax+1  ; i++)
      for (int j = 0; j < nvalues  ; j++)
	locdisplace[i][j] = 0.0;


}

/* ---------------------------------------------------------------------- */

void FixMoveDupl:: collect_displacement(){

  int n,j,i,k;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int natoms=atom->natoms;
  tagint *tag = atom->tag;

 int me = comm->me;


/*
   clear invoked flag of all computes
   called everywhere that computes are used, before computes are invoked
   invoked flag used to avoid re-invoking same compute multiple times
   and to flag computes that store invocation times as having been invoked
 */

  modify->clearstep_compute();

  for (i=0;i<natoms;i++){
       locdisplace[i][0]=0.0;
       locdisplace[i][1]=0.0;
       locdisplace[i][2]=0.0;
       locdisplace[i][3]=0.0;
  }


  for (int m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];


    if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= INVOKED_PERATOM;
    	       }

        int jm1 = j - 1;
        double **compute_array = compute->array_atom;
        for (i = 0; i < nlocal; i++)
          if (! (mask[i] & groupbit)){
              k=tag[i];
              locdisplace[k][m] = compute_array[i][jm1];
//printf("atomo %d %d %d %d %f  \n",k,m,i,jm1,compute_array[i][jm1]);
               }
          }
    }

// printf ("QUA FIX DIOCANE %d %f %f \n",buff_size,displace[0][0],locdisplace[0][0]);


//MPI_Reduce(&locdisplace[0][0],&displace[0][0],buff_size,MPI_DOUBLE,MPI_SUM,0,world);
//if(me == 0){for (i=1; i<  atom->natoms+1;i++)printf ("QUA FIX DIOCANE %d %f %f %f \n",i,displace[i][0],displace[i][1],displace[i][2]);}
//MPI_Bcast(&displace[0][0],buff_size,MPI_DOUBLE,0,world);

// if(me == 0){for (i=1; i<  atom->natoms+1;i++)printf ("BROADCAST %d %f %f %f \n",i,displace[i][0],displace[i][1],displace[i][2]);}

  MPI_Allreduce(&locdisplace[0][0],&displace[0][0],buff_size,MPI_DOUBLE,MPI_SUM,world);


//if(me == 0){for (i=1; i<  atom->natoms+1;i++)fprintf (fpl,"BROADCAST %d %f %f %f \n",i,displace[i][0],displace[i][1],displace[i][2]);}


}

/* ---------------------------------------------------------------------- */

void FixMoveDupl::post_integrate()
{



  int i,j,k;
  double xold[3];

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *type = atom->type;
  imageint *image = atom->image;
  int ndiff;


//  fprintf(fpl,"TIMESTEP %d\n\n\n\n\n",update->ntimestep);

//printf("GIIOVANNONA\n");

  collect_displacement();

//printf("GIIOVANNONA\n");

    for (i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
          j=tag[i];
	  k=j-1;
            xold[0]=x[i][0];
            xold[1]=x[i][1];
            xold[2]=x[i][2];

          x[i][0] = x[i][0] + displace[k][0]; 
          x[i][1] = x[i][1] + displace[k][1]; 
          x[i][2] = x[i][2] + displace[k][2];

// fprintf(fpl,"prima %d %d %f %f %f %f %f %f %f %f %f \n",i,j,xold[0],xold[1],xold[2],displace[k][0],displace[k][1],displace[k][2],x[i][0],x[i][1],x[i][2]);

      domain->remap(x[i],image[i]);
//	  domain->remap_near(x[i],xold);

//fprintf(fpl,"dopo %f %f %f  \n",x[i][0],x[i][1],x[i][2]);
        }
    }

}



/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixMoveDupl::set_arrays(int i)
{

/*
  double **x = atom->x;
  imageint *image = atom->image;

  int *mask = atom->mask;

  // particle not in group

  if (!(mask[i] & groupbit)) {
    xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
    return;
  }

  // current time still equal fix creation time

  if (update->ntimestep == time_origin) {
    domain->unmap(x[i],image[i],xoriginal[i]);
    return;
  }
*/
  tagint *tag = atom->tag;
  int j=tag[i];

    displace[j][0] = locdisplace[j][0] = 0.0;
    displace[j][1] = locdisplace[j][1] = 0.0;
    displace[j][2] = locdisplace[j][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixMoveDupl::reset_dt()
{
  error->all(FLERR,"Resetting timestep size is not allowed with fix move/duplicate");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixMoveDupl::memory_usage()
{
  double bytes = 0.0;
  int maxatom = atom->nmax;
  bytes = 3 * maxatom * 3 * sizeof(double);  // displace, localdisplace, compute_array
  return bytes;
}
