/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.


This compute calculates the interactions among all atom types 
   Author: Di Pasquale Nicodemo
   University of Leicester, March 2020
   email: nicodemo.dipasquale@gmail.com    
   
   The documentation for this pair potential can be browsed at the following link:
   https://demonico85.github.io/cleaving/

 ----------------------------------------------------------------------*/


#include <mpi.h>
#include <string.h>
#include "compute_paircleav.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

ComputePairCl::ComputePairCl(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  pstyle(NULL),  pair(NULL)
{
  if (narg < 4 || narg > 7) error->all(FLERR,"Illegal compute pair command");

  array_flag = 1;
  scalar_flag = 1;
  vector_flag = 1;
  extscalar   = 1; //global scalar is extensive
  extvector   = 1; //same for the vector
  extarray    = 1; //same for the array
  peratom_flag = 0; //no per atom sub
  peflag = 1;  // this allows to have eflag != 0 in the potential computation. Timestep 0 energy is always calculated, for ts > 0 it depends from this flag
  timeflag = 1;
  norm = 0;


  int n = strlen(arg[3]) + 1;
  if (lmp->suffix) n += strlen(lmp->suffix) + 1;
  pstyle = new char[n];
  strcpy(pstyle,arg[3]);

  

  if (narg > 4) {
    if (strcmp(arg[4],"norm") == 0) norm = 1;
    else 
        error->all(FLERR,"Illegal compute pair command (norm)");
  

    nwalls = utils::numeric(FLERR,arg[5],false,lmp);


    if(strcmp(arg[6],"x") != 0 && strcmp(arg[6],"y") != 0 && strcmp(arg[6],"z") != 0)
     error->all(FLERR,"Illegal compute pair command (direction)");

// For Now working only with NVT/NVE, because the size of the box is taken 
// at the beginning at the simulation
// We need the surface perpendicular to the direction arg[6]

    if(strcmp(arg[6],"x") == 0){
        edge1 = domain->yprd;
        edge2 = domain->zprd;
        }
    else if(strcmp(arg[6],"y") == 0){
        edge1 = domain->xprd;
        edge2 = domain->zprd;
        }
    else if(strcmp(arg[6],"z") == 0){
        edge1 = domain->xprd;
        edge2 = domain->yprd;
        }

    }
  // check if pair style with and without suffix exists

  pair = force->pair_match(pstyle,1);
  if (!pair && lmp->suffix) {
    strcat(pstyle,"/");
    strcat(pstyle,lmp->suffix);
    pair = force->pair_match(pstyle,1);
  }

  if (!pair)
    error->all(FLERR,"Unrecognized pair style in compute pair command");
  
  npair = pair->nextra;
  
  




  ntypes = atom->ntypes;
  
  int dmy;
  
  dmy=(npair-1)/ntypes/ntypes;
  
  printf("QUAAAAA %d %d %d \n",ntypes,dmy,npair);
  
    if(dmy == 4)ntypes=2*ntypes;
  
  
  size_array_rows = ntypes+1;
  size_array_cols = ntypes+1;


  allocate();


}

/* ---------------------------------------------------------------------- */

ComputePairCl::~ComputePairCl()
{
  delete [] pstyle;
  delete [] one;
  delete [] vector;
  delete [] array;
}

/* ---------------------------------------------------------------------- */

void ComputePairCl::init()
{

  double areaxy;

  // recheck for pair style in case it has been deleted

  pair = force->pair_match(pstyle,1);
  if (!pair)
    error->all(FLERR,"Pair style not found");

  // Calculation of the area to be used for normalization

    if(norm){
        areaxy= nwalls*edge1*edge2; // 
        iareaxy = 1./areaxy;
        }
    else
        iareaxy = 0.0;
}

/* ---------------------------------------------------------------------- */

void ComputePairCl::allocate()
{

int i;

  if (npair) {
    size_vector = npair;
    extvector = 1;
    one = new double[npair];
    vector = new double[npair];
  } else one = vector = NULL;


  if (size_array_rows && size_array_cols) {
    array = new double*[size_array_rows];
    for (i = 0; i < size_array_rows; i++)
        array[i] = new double[size_array_rows];
    for (i = 0; i < size_array_rows; i++){
        array[i][0] = 0.0;
        array[0][i] = 0.0;
        }   
  } else array = NULL;


}

/* ---------------------------------------------------------------------- */

void ComputePairCl::compute_vector(){

//retrieve(){


  for (int i = 0; i < npair; i++)
    one[i] = pair->pvector[i];
  MPI_Allreduce(one,vector,npair,MPI_DOUBLE,MPI_SUM,world);


}

/* ---------------------------------------------------------------------- */

double ComputePairCl::compute_scalar(){

int i,j,n;

  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
      error->all(FLERR,"Energy was not tallied on needed timestep");

   compute_vector();

   scalar = 0.0;

    for (i = 1; i < size_array_rows; i++){
        for (j = 1; j < size_array_cols; j++){
            n=j+(i-1)*ntypes;
            scalar += vector[n];
            }
        }

    if(norm) scalar *= iareaxy;



    return scalar;

}


/* ---------------------------------------------------------------------- */

void ComputePairCl::compute_array()
{

int i,j,n;

  invoked_array = update->ntimestep;
  if (update->eflag_global != invoked_array){
//    printf ("%d %d\n\n",invoked_array,update->eflag_global);
    error->all(FLERR,"Energy was not tallied on needed timestep");}


    compute_vector();

    for (i = 1; i < size_array_rows; i++){
        for (j = 1; j < size_array_cols; j++){
            n=i+(j-1)*ntypes;
            array[i][j] = vector[n];
            }
        }


    for (i = 1; i < size_array_rows; i++){
        for (j = i+1; j < size_array_cols; j++){
            array[i][j] += array[j][i];
            array[j][i] = 0.0;
        }
    }

    if(norm) {
        for (i = 1; i < size_array_rows; i++){
            for (j = i; j < size_array_cols; j++){
             array[i][j] *= iareaxy;
                } 
            }  
        }


}




































