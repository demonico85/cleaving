/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.


This fix calculates the wall cleaving interactions

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
#include "fix_cleaving_wall.h"
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

FixForceLJWall::FixForceLJWall(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
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
  nevery      = 1;
 
// Init some Variables

 lposwall = NULL;
 hposwall = NULL;
 namefile = NULL;
 fpw      = NULL;
 ewall    = NULL;
 ewall_all= NULL;
/*

narg          4       5         6            7         8   
type       epsilon  sigma  dist_W_from cent   delta   rcut    
*/


  ntypes = atom->ntypes;
  size_vector = ntypes+1;
  allocate();

  iarg=0;
  me = comm->me;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wallforce command");
        if (me == 0) {  

          strl = strlen (arg[iarg+1])+1;
          namefile = new char[strl];
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


  if (narg < 9) error->all(FLERR,"Illegal fix addforce command");

  epsilon       = utils::numeric(FLERR,arg[3],false,lmp);
  sigma         = utils::numeric(FLERR,arg[4],false,lmp);
  zw            = utils::numeric(FLERR,arg[5],false,lmp);
  delta         = utils::numeric(FLERR,arg[6],false,lmp);
  rw            = utils::numeric(FLERR,arg[7],false,lmp);
  cleav_wall    = utils::numeric(FLERR,arg[8],false,lmp);





}


/* ---------------------------------------------------------------------- */

FixForceLJWall::~FixForceLJWall()
{

  memory->destroy(lposwall);
  memory->destroy(hposwall);
  memory->destroy(ewall);
  memory->destroy(ewall_all);

}

/* ---------------------------------------------------------------------- */

void FixForceLJWall::allocate()
{


  memory->create(ewall,size_vector,"fix:ewall");
  memory->create(ewall_all,size_vector,"fix:ewall_all");

}

/* ---------------------------------------------------------------------- */

int FixForceLJWall::setmask()
{

  int mask = 0;
  mask |= POST_FORCE;
//  mask |= THERMO_ENERGY;
  return mask;
}



/* ----------------------------------------------------------------------
   init variables
------------------------------------------------------------------------- */

void FixForceLJWall::init()
{
  // check variables

double prec;
double areaxy;


prec = 1.e-10;
cleavwork = 0.0;
localwork = 0.0;


   xprd = domain->xprd;
   yprd = domain->yprd;
   zprd = domain->zprd;
 
areaxy= 2.0*xprd*yprd; // there are two walls
iareaxy = 1./areaxy;


  if(epsilon < 0.0)
    error->all(FLERR,"Illegal Wall fix addforce command (eps)");
  if(sigma < 0.0)
    error->all(FLERR,"Illegal Wall fix addforce command (sigma)");
  if (strcmp(update->unit_style,"lj") == 0) {
    if(epsilon < 1.0-prec || epsilon > 1.0+prec)
        error->all(FLERR,"Illegal Wall fix addforce command (eps and lj units)");
    if(epsilon < 1.0-prec || epsilon > 1.0+prec)
        error->all(FLERR,"Illegal Wall fix addforce command (sigma and lj units)");   
    }
  if(delta <= 0.0)
    error->all(FLERR,"Illegal Wall fix addforce command (delta)");
  invdelta = 1.0/delta;
  if(zw < 0.0)
    error->all(FLERR,"Illegal Wall fix addforce command (zw)");
  if(rw < 0.0)
    error->all(FLERR,"Illegal Wall fix addforce command (rw)");
   
   rwsq=rw*rw;



  // Calculate init quantities

  ReadFWalls();

  lj1 = 48.0 * epsilon * pow(sigma,12.0);
  lj2 = 24.0 * epsilon * pow(sigma, 6.0);
  lj3 =  4.0 * epsilon * pow(sigma,12.0);
  lj4 =  4.0 * epsilon * pow(sigma, 6.0);
  
  i2delta = 1.0/delta/2.0;
  i4delta = 1.0/delta/4.0;
  pconst = (2.0-delta)/4.0;


}

/* ---------------------------------------------------------------------- */

void FixForceLJWall::setup(int vflag)
{
  post_force(vflag);
}


/* ---------------------------------------------------------------------- */

void FixForceLJWall::post_force(int vflag)
{




  
  int i,j,itype,ptype;
  double xtmp,ytmp,ztmp,ldelx,hdelx,ldely,hdely,fpair;
  double ldelz,hdelz;
// double zsq1, zsq2;
  double lrsq,hrsq,r2inv,r6inv,forcelj;
  double phi1,phi2,Phi,fx1,fy1,fz1,fx2,fy2,fz2,modmw;
  double invphis,udiff,usum,C1,C2;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;



  cleavwork = 0.0;
  localwork = 0.0;
  eflag = 0;

  for(i=0; i<size_vector; i++) ewall[i]=0.0;


  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
        itype = type[i];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];

        Phi=0.0;
        phi1=0.0;
        phi2=0.0;
        fx1=0.0;
        fy1=0.0;
        fz1=0.0;
        fx2=0.0;
        fy2=0.0;
        fz2=0.0;
    
//        zsq1 = fabs(ztmp-cleav_wall-zw);
//        zsq2 = fabs(ztmp-cleav_wall+zw);

        for (j = 0; j < natomsinwalls[0]; j++) {

          ldelx = xtmp - lposwall[j][0];
          ldely = ytmp - lposwall[j][1];
          ldelz = ztmp - lposwall[j][2];

          domain->minimum_image(ldelx,ldely,ldelz);

          lrsq = ldelx*ldelx + ldely*ldely + ldelz*ldelz;


      if ( lrsq < rwsq ) {
        r2inv = 1.0/lrsq;
        r6inv = r2inv*r2inv*r2inv;
        phi1 += r6inv*(lj3*r6inv - lj4)+epsilon;
        fpair = r6inv * (lj1*r6inv - lj2);
        forcelj = fpair*r2inv;

        fx1 += forcelj * ldelx;
        fy1 += forcelj * ldely;
        fz1 += forcelj * ldelz;
        }
    }

    for (j = 0; j < natomsinwalls[1]; j++) {

          hdelx = xtmp - hposwall[j][0];
          hdely = ytmp - hposwall[j][1];
          hdelz = ztmp - hposwall[j][2];

          domain->minimum_image(hdelx,hdely,hdelz);

          hrsq = hdelx*hdelx + hdely*hdely + hdelz*hdelz;

      if ( hrsq < rwsq ) {
        r2inv = 1.0/hrsq;
        r6inv = r2inv*r2inv*r2inv;
        phi2 += r6inv*(lj3*r6inv - lj4)+epsilon;
        fpair = r6inv * (lj1*r6inv - lj2);
        forcelj = fpair*r2inv;
        
        fx2 += forcelj * hdelx;
        fy2 += forcelj * hdely;
        fz2 += forcelj * hdelz;
        }
    }

    ptype=-10;

    if(phi1 > 0.0 && phi2 > 0.0){
    ptype = WallP(phi1,phi2);

    if(ptype == 0){
        Phi=phi1;
        f[i][0]   += fx1;
        f[i][1]   += fy1;
        f[i][2]   += fz1; 
        localwork -= fz1;  
        ewall[itype]+=Phi;
        }
    else if(ptype == 1){
        Phi=phi2;
        f[i][0]   += fx2;
        f[i][1]   += fy2;
        f[i][2]   += fz2; 
        localwork += fz2;  
        ewall[itype]+=Phi;
        }
    else if (ptype == 2){
		udiff = phi1-phi2;
		usum  = phi1+phi2;
        invphis    = udiff/usum;

        Phi        = pconst*usum - i4delta*udiff*udiff/usum;      
		C1         = pconst + i4delta*invphis*invphis;
		C2         = i2delta*invphis;
        modmw      = C1*(fz2-fz1) + C2*(fz1+fz2);

        f[i][0]   += C2*(fx2-fx1) + C1*(fx1+fx2);
        f[i][1]   += C2*(fy2-fy1) + C1*(fy1+fy2);;
        f[i][2]   += C2*(fz2-fz1) + C1*(fz1+fz2);;

        ewall[itype] +=Phi;
        localwork += modmw;
            }
        }
    }
}

}



/* ---------------------------------------------------------------------- */

int FixForceLJWall::WallP(double x, double y)
{

double mindist,hidist;

  mindist=y-delta*(x+y); 
  hidist=y+delta*(x+y);


  if(x <= mindist) return 0;
  else if(x >= hidist) return 1;
  else return 2;

}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixForceLJWall::compute_scalar()
{
  // only sum across procs one time

  MPI_Allreduce(&localwork, &cleavwork,1,MPI_DOUBLE,MPI_SUM,world);
  
  cleavwork=cleavwork*iareaxy;

  return cleavwork;
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixForceLJWall::compute_vector(int i)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,size_vector,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[i];
}


/* ----------------------------------------------------------------------
   FCC111 walls
------------------------------------------------------------------------- */


void FixForceLJWall::ReadFWalls()
{
    

  int i,me,cnt;
  char str[128];
  
  me = comm->me;

    if(lposwall != NULL)memory->destroy(lposwall);
    if(hposwall != NULL)memory->destroy(hposwall);

    if (me == 0) {
      fpw=fopen(namefile, "r");
        if (fpw == NULL) {
          sprintf(str,"Cannot open data file %s",namefile);
          error->one(FLERR,str);
        }

        if(fscanf(fpw, "%d", &natomsinwalls[0]) != 1){ 
          sprintf(str,"Cannot read number of walls in %s",namefile);
          error->one(FLERR,str);}
    }

    MPI_Bcast(&natomsinwalls[0], 1, MPI_INT, 0, world);
    
    memory->create(lposwall,natomsinwalls[0],3,"fix:lwalls");

    if( me == 0){
        cnt=0;
         while(cnt < natomsinwalls[0]){
            for(int i=0;i<3;i++) {
            	if(fscanf (fpw, "%lf", &lposwall[cnt][i]) != 1){
            	sprintf(str,"Cannot read position of wall %d in %s",i,namefile);
         	error->one(FLERR,str);}
         		}
            cnt+=1;
        }
        if (fscanf(fpw, "%d", &natomsinwalls[1]) != 1){ 
          sprintf(str,"Cannot read number of wells in %s",namefile);
          error->one(FLERR,str);}
    }

    MPI_Bcast(&cnt, 1, MPI_INT, 0, world);
    if(natomsinwalls[0] != cnt) error->all(FLERR,"Check input file walls, different number of elements");

    MPI_Bcast(&lposwall[0][0], natomsinwalls[0]*3, MPI_DOUBLE, 0, world);


    MPI_Bcast(&natomsinwalls[1], 1, MPI_INT, 0, world);
    memory->create(hposwall,natomsinwalls[1],3,"fix:hwalls");




    if( me == 0){
        cnt=0;
        while(cnt < natomsinwalls[1]){
            for(int i=0;i<3;i++) {
            	if(fscanf (fpw, "%lf", &hposwall[cnt][i]) != 1){
            	sprintf(str,"Cannot read position of wall %d in %s",i,namefile);
         	error->one(FLERR,str);}
         		}
            cnt+=1;
        }
        fclose(fpw);
    }

    MPI_Bcast(&cnt, 1, MPI_INT, 0, world);
    if(natomsinwalls[1] != cnt) error->all(FLERR,"Check input file walls, different number of elements");

    MPI_Bcast(&hposwall[0][0], natomsinwalls[1]*3, MPI_DOUBLE, 0, world);


// Considering flat surface 


for (i=0; i<natomsinwalls[0]; i++){
    lposwall[i][2] = lposwall[i][2] + cleav_wall-zw; 

    }
for(i=0; i<natomsinwalls[1]; i++){
    hposwall[i][2] = hposwall[i][2] + cleav_wall+zw;
    }



}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixForceLJWall::memory_usage()
{
  double bytes = 0.0;
  bytes = natomsinwalls[0] * natomsinwalls[1] * 3 * sizeof(double);  // lposwall, hposwall
  return bytes;
}

