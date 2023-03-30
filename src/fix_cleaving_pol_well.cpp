/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.


This fix calculates the well cleaving interactions

   Author: Di Pasquale Nicodemo
   University of Leicester, March 2020
   email: nicodemo.dipasquale@gmail.com    
   
   The documentation for this pair potential can be browsed at the following link:
   https://demonico85.github.io/cleaving/

 ----------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_cleaving_pol_well.h"
#include "comm.h"
#include "domain.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixForcePolWell::FixForcePolWell(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  int iarg, me, strl;
  
//  Fix(lmp, narg, arg)

  dynamic_group_allow = 0;
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 0;

  poswell  = NULL;
  namefile = NULL;
  fpw      = NULL;
  ewell    = NULL;
  ewell_all= NULL;

  ntypes = atom->ntypes;
  size_vector = ntypes+1;

// Specific for this subroutine 

  if (narg < 7) error->all(FLERR,"Illegal fix wellforce command");

  dw     = utils::numeric(FLERR,arg[3],false,lmp);
  rw     = utils::numeric(FLERR,arg[4],false,lmp);
  expon  = utils::numeric(FLERR,arg[5],false,lmp);
  lambda = utils::numeric(FLERR,arg[6],false,lmp);


/*
             _________________________________   ___________________   
            |                                                       |                  
narg          4       5         6          7         8      9            
type         dw        rw     expon     lambda      file name_file        
pair_coeff   1.0      1.0      3            0.25      
*/

  iarg=0;
  me = comm->me;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wellforce command");
        if (me == 0) {  

          strl = strlen (arg[iarg+1])+1;
          namefile = new char[strl];
//          memory->create(namefile,strl,"fix:namefile");
          strcpy(namefile,arg[iarg+1]);
          fpw = fopen(namefile,"r");
          if (fpw == NULL) {
            char str[128];
            sprintf(str,"Cannot open fix wellforce file %s",arg[iarg+1]);
            error->one(FLERR,str);
                }
            }
        }
      iarg += 1;
    }


}


/* ---------------------------------------------------------------------- */

FixForcePolWell::~FixForcePolWell()
{

  memory->destroy(poswell);
  memory->destroy(ewell);
  memory->destroy(ewell_all);

//  delete [] namefile;
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void FixForcePolWell::allocate()
{


  memory->create(ewell,size_vector,"fix:ewall");
  memory->create(ewell_all,size_vector,"fix:ewall_all");
  memory->create(poswell,nwells,3,"fix:wells");

}

/* ---------------------------------------------------------------------- */

int FixForcePolWell::setmask()
{

  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}



/* ----------------------------------------------------------------------
   init variables
------------------------------------------------------------------------- */

void FixForcePolWell::init()
{
  // check variables

double areaxy;



cleavwork = 0.0;
localwork = 0.0;

xedge = domain->boxhi[0] - domain->boxlo[0];
yedge = domain->boxhi[1] - domain->boxlo[1];
zedge = domain->boxhi[2] - domain->boxlo[2];

xhalf = xedge*0.5;
yhalf = yedge*0.5;
zhalf = zedge*0.5;
 
areaxy= 2.0*xedge*yedge; // there are two walls
iareaxy = 1./areaxy;


//  if(dw <= 0.0)
//    error->all(FLERR,"Illegal Wall fix addforce command (dw)");
  if(rw <= 0.0)
    error->all(FLERR,"Illegal Wall fix addforce command (rw)");
  if(expon <= 0.0)
    error->all(FLERR,"Illegal Wall fix addforce command (expon)");

  // Calculate init quantities

  rwsq  = rw*rw;
  irwsq  = 1./rwsq;
  M = expon-1.0;
  prefactor = -2.0*expon*dw*lambda*irwsq;

    ReadFWells();


}


/* ---------------------------------------------------------------------- */

void FixForcePolWell::post_force(int vflag)
{
  int i,j,itype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,forcelj;
  double Phi,fx,fy,fz;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  cleavwork = 0.0;
  localwork = 0.0;
  eflag = 0;

for(i=0; i<size_vector; i++) ewell[i]=0.0;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

      Phi=0.0;
      fx=0.0;
      fy=0.0;
      fz=0.0;

        for (j = 0; j < nwells; j++) {

          delx = xtmp - poswell[j][0];
          if (delx >   xhalf ) delx = delx - xedge;
          if (delx <= -xhalf ) delx = delx + xedge;


          dely = ytmp - poswell[j][1];
          if (dely >   yhalf ) dely = dely - yedge;
          if (dely <= -yhalf ) dely = dely + yedge;

          delz = ztmp - poswell[j][2];
         if (delz >   zhalf ) delz = delz - zedge;
         if (delz <= -zhalf ) delz = delz + zedge;

          rsq = delx*delx + dely*dely + delz*delz;

            if ( rsq < rwsq ) {

                fpair = rsq*irwsq-1.0;
                forcelj = prefactor*pow(fpair,M);
                Phi += dw*pow(fpair,expon);

                fx += forcelj * delx;
                fy += forcelj * dely;
                fz += forcelj * delz;

                ewell[itype] += Phi;


                }
            }
                f[i][0]   += fx;
                f[i][1]   += fy;
                f[i][2]   += fz; 
                localwork += Phi;  
        }
    }


}

/* ---------------------------------------------------------------------- */

void FixForcePolWell::setup(int vflag)
{
  post_force(vflag);
}


/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixForcePolWell::compute_scalar()
{
  // only sum across procs one time

  MPI_Allreduce(&localwork, &cleavwork,1,MPI_DOUBLE,MPI_SUM,world);
  
  cleavwork=cleavwork*iareaxy;

  return cleavwork;
}

/* ----------------------------------------------------------------------
   FCC111 walls
------------------------------------------------------------------------- */


void FixForcePolWell::ReadFWells()
{
    

  int me,cnt;

  me = comm->me;

    if(poswell != NULL)
               memory->destroy(poswell);

    if (me == 0) {
      fpw=fopen(namefile, "r");
        if (fpw == NULL) {
          char str[128];
          sprintf(str,"Cannot open data file %s",namefile);
          error->one(FLERR,str);
        }

        fscanf(fpw, "%d", &nwells);
    }

    MPI_Bcast(&nwells, 1, MPI_INT, 0, world);
    allocate();

    if( me == 0){
        cnt=0;
         while(cnt < nwells){
            for(int i=0;i<3;i++) fscanf (fpw, "%lf", &poswell[cnt][i]);
            cnt+=1;
        }
        fclose(fpw);
    }

    fpw = NULL;
    MPI_Bcast(&cnt, 1, MPI_INT, 0, world);

    if(nwells != cnt) error->all(FLERR,"Check input file wells, different number of elements");

    MPI_Bcast(&poswell[0][0], nwells*3, MPI_DOUBLE, 0, world);



}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixForcePolWell::compute_vector(int i)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewell,ewell_all,size_vector,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewell_all[i];
}

