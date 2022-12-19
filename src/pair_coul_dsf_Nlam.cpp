/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ORNL)
   References: Fennell and Gezelter, JCP 124, 234104 (2006)

   Contributing author: Di Pasquale Nicodemo
    University of Leicester, March 2020

% Force single not yet changed to include lambda

------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_coul_dsf_Nlam.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairCoulNLDSF::PairCoulNLDSF(LAMMPS *lmp) : Pair(lmp) {

int n=atom->ntypes;

    nextra = n*n+1; // the zero-th position is not considered
    pvector = new double[nextra];
    natoms=atom->natoms+1; //array does start from zero

   xprd = domain->xprd;
   yprd = domain->yprd;
   zprd = domain->zprd;
   xy = domain->xy;
   yz = domain->yz;
   xz = domain->xz;

}

/* ---------------------------------------------------------------------- */

PairCoulNLDSF::~PairCoulNLDSF()
{
  if (copymode) return;

  delete [] pvector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(powlambda);
    memory->destroy(powDlambda);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulNLDSF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,m;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double r,rsq,r2inv,forcecoul,factor_coul,flam,fDlam;
  double prefactor,erfcc,erfcd,t;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int scaling,k1,k2;

  ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int n=atom->ntypes;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;


  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  tagint *molecule = atom->molecule;
  int *ivector = atom->ivector[index];

  global_boundary();

  for(i=0; i<nextra ; i++){
    pvector[i] = 0.0;
    }



  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
 

    if (eflag) {
      double e_self = -(e_shift/2.0 + alpha/MY_PIS) * qtmp*qtmp*qqrd2e;
      ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_coulsq) {
        r2inv = 1.0/rsq;

        r = sqrt(rsq);
        prefactor = qqrd2e*qtmp*q[j]/r;
        erfcd = exp(-alpha*alpha*rsq);
        t = 1.0 / (1.0 + EWALD_P*alpha*r);
        erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
        forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd +
                                 r*f_shift) * r;
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;


        k1=molecule[i];
        k2=molecule[j];


        scaling=0;
        if(k1 != k2)
             scaling = find_scaling(ivector[i],ivector[j],i,j,x[j]);

          flam = 1.0;
          fDlam = 1.0;

          if(scaling){
              flam = powlambda[itype][jtype];
              fDlam =  powDlambda[itype][jtype];
            }


        fpair = flam * forcecoul * r2inv;



        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          m=jtype+(itype-1)*n;
          ecoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
          if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
	   if(scaling){
             pvector[m] += fDlam * ecoul;
             ecoul *= flam ;
            }

        } else ecoul = 0.0;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             0.0,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}
/* ----------------------------------------------------------------------
 *  *   Global Boundary
 *   * ------------------------------------------------------------------------- */

void PairCoulNLDSF::global_boundary(){

int i,k;
int ibox[3];
int nlocal = atom->nlocal;
int loc_box[natoms];
double xx[3];
double **x = atom->x;

imageint *image = atom->image;
imageint locim;
tagint *tag = atom->tag;

     for(i=0; i<natoms; i++){
        loc_box[i]=0;
    }


     for(i=0; i<nlocal; i++){
        k=tag[i];

        ibox[0] = (image[i] & IMGMASK) - IMGMAX;
        ibox[1] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ibox[2] = (image[i] >> IMG2BITS) - IMGMAX;


      if (domain->triclinic == 0) {
        xx[0] = x[i][0] + xprd*ibox[0];
        xx[1] = x[i][1] + yprd*ibox[1];
        xx[2] = x[i][2] + zprd*ibox[2];
      } else {
        xx[0] = x[i][0] + xprd*ibox[0] + ibox[1]*xy + ibox[2]*xz;
        xx[1] = x[i][1] + yprd*ibox[1] + ibox[2]*yz;
        xx[2] = x[i][2] + zprd*ibox[2];
      }


       locim = ((imageint) IMGMAX << IMG2BITS) |
           ((imageint) IMGMAX << IMGBITS) | IMGMAX;


        domain->remap(xx,locim);

        ibox[0] = ( locim & IMGMASK) - IMGMAX;
        ibox[1] = ( locim >> IMGBITS & IMGMASK) - IMGMAX;
        ibox[2] = ( locim >> IMG2BITS) - IMGMAX;

        loc_box[k]=ibox[ind_dir];


        }


 MPI_Allreduce(&loc_box[0],&gbox[0],natoms,MPI_INT,MPI_SUM,world);


}

/* ----------------------------------------------------------------------
  find scaling
------------------------------------------------------------------------- */

int PairCoulNLDSF::find_scaling(int imvec, int jmvec, int i, int j,double *xj){


tagint *tag = atom->tag;

int nlocal = atom->nlocal;
int k1,k2,ipbc,jpbc,ghost;
double lamda[3],xx[3];

ghost=1;
k1=tag[i];
k2=tag[j];

 if(j < nlocal)ghost=0;

ipbc = gbox[k1];
jpbc = gbox[k2];



if(ipbc != 0 && jpbc == 0) {
	if(imvec != jmvec)
	   return 1;
    else
	   return 0;
    }
else if (ipbc == 0 && jpbc != 0){
        if(imvec != jmvec)
           return 1;
        else
           return 0;
    }
else if (ipbc == 0 && jpbc == 0){
     if(imvec != jmvec){
       if(ghost){
// Passing by reference the value is modified we don't want it
    	   xx[0] = xj[0];
    	   xx[1] = xj[1];
    	   xx[2] = xj[2];
           domain->x2lamda(xx,lamda);
           if( lamda[ind_dir] < 0 || lamda[ind_dir] > 1.0 )
              return 1;
           else
              return 0;
        }
        else
           return 0;
        }
     else
       return 0;
    }
else if(ipbc !=0 && jpbc != 0 ){
    if(imvec != jmvec)
        return 1;
    else
        return 0;
    }
else
  return -1001;

}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulNLDSF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(lam,n+1,n+1,"pair:lambda");
  memory->create(gbox,natoms,"pair:gbox");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(powlambda,n+1,n+1,"pair:powlambda");
  memory->create(powDlambda,n+1,n+1,"pair:powDlambda");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulNLDSF::settings(int narg, char **arg)
{

  if (narg != 6) error->all(FLERR,"Illegal pair_style command");

  alpha = utils::numeric(FLERR,arg[0],false,lmp);
  cut_coul = utils::numeric(FLERR,arg[1],false,lmp);

 lambda = utils::numeric(FLERR,arg[2],false,lmp);

  int n = strlen(arg[3]) + 1;
  idflag = new char[n];
  strcpy(idflag,arg[3]);


  int flag, ncols;

  index = atom->find_custom( idflag, flag, ncols);

  if(flag) error->all(FLERR,"Custom property must be integer for dsf pair_style command");

  if(strcmp(arg[4],"x") != 0 && strcmp(arg[4],"y") != 0 && strcmp(arg[4],"z") != 0)
     error->all(FLERR,"Illegal LJ-BJ wells_command (direction)");
// For Now working only with NVT/NVE, because the size of the box is taken
// at the beginning at the simulation


  if(strcmp(arg[4],"x") == 0){
     ind_dir=0;
      }
  else if(strcmp(arg[4],"y") == 0){
     ind_dir=1;
      }
  else if(strcmp(arg[4],"z") == 0){
     ind_dir=2;
     }

    npow =  utils::numeric(FLERR,arg[5],false,lmp);

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulNLDSF::coeff(int narg, char **arg)
{

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double l = lambda;
  l = utils::numeric(FLERR,arg[2],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      lam[i][j] = l;
      powlambda[i][j] = pow(l,npow);
      powDlambda[i][j] = npow*pow(l,npow-1.0);
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulNLDSF::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/dsf requires atom attribute q");

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;
  double erfcc = erfc(alpha*cut_coul);
  double erfcd = exp(-alpha*alpha*cut_coul*cut_coul);
  f_shift = -(erfcc/cut_coulsq + 2.0/MY_PIS*alpha*erfcd/cut_coul);
  e_shift = erfcc/cut_coul - f_shift*cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulNLDSF::init_one(int i, int j)
{
  lam[j][i] = lam[i][j];
  powlambda[j][i] = powlambda[i][j];
  powDlambda[j][i] = powDlambda[i][j];

  return cut_coul;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulNLDSF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulNLDSF::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulNLDSF::write_restart_settings(FILE *fp)
{
  fwrite(&alpha,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&npow,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulNLDSF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&alpha,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&npow,sizeof(double),1,fp);
  }
  MPI_Bcast(&alpha,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&npow,1,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

double PairCoulNLDSF::single(int i, int j, int itype, int jtype, double rsq,
                           double factor_coul, double factor_lj,
                           double &fforce)
{
  double r2inv,r,erfcc,erfcd,prefactor,t;
  double forcecoul,phicoul;

  r2inv = 1.0/rsq;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    prefactor = factor_coul * force->qqrd2e * atom->q[i]*atom->q[j]/r;
    erfcd = exp(-alpha*alpha*rsq);
    t = 1.0 / (1.0 + EWALD_P*alpha*r);
    erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;

    forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS*erfcd +
      r*f_shift) * r;

    phicoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
    eng += phicoul;
  } else forcecoul = 0.0;

  fforce = forcecoul * r2inv;

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairCoulNLDSF::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  return NULL;
}
