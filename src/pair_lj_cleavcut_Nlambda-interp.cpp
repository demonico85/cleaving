/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------


This subroutine is written from the original pair_lj_cut in LAMMPS adapted for
cleaving calculations

   Contributing author: Di Pasquale Nicodemo
   University of Leicester, March 2020
   email: nicodemo.dipasquale@gmail.com    
   
   The documentation for this pair potential can be browsed at the following link:
   https://demonico85.github.io/cleaving/
    

------------------------------------------------------------------------- */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_lj_cleavcut_Nlambda-interp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace MathConst;



/* ---------------------------------------------------------------------- */

PairLJCleavCutNLInt::PairLJCleavCutNLInt(LAMMPS *lmp) : Pair(lmp) , idflag(NULL)
{

  writedata = 1;
  pallocation = 0;
  int n=atom->ntypes;
  natoms=atom->natoms; //array does start from zero
  natoms=natoms+1; 
  if(n > 0){
    nextra = n*n+1; // the zero-th position is not considered
    pvector = new double[nextra];
    pallocation = 1;
    }

   xprd = domain->xprd;
   yprd = domain->yprd;
   zprd = domain->zprd;
   xy = domain->xy;
   yz = domain->yz;
   xz = domain->xz;
 
}

/* ---------------------------------------------------------------------- */

PairLJCleavCutNLInt::~PairLJCleavCutNLInt()
{

    delete [] pvector;
    delete [] idflag;

  if (allocated) {
    memory->destroy(gbox);
    memory->destroy(radm);
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(lam);
    memory->destroy(powlambda);
    memory->destroy(powDlambda);
    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    memory->destroy(lju0);    
    memory->destroy(lju1);    
    memory->destroy(lju2); 
    memory->destroy(lju3);    
    memory->destroy(lju4);     
  }

}

/* ---------------------------------------------------------------------- */

void PairLJCleavCutNLInt::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,m,k1,k2;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,flam,fDlam;
  double rsq,r2inv,r6inv,forcelj,factor_lj,rfourth;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int scaling;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int n = atom->ntypes;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  tagint *molecule = atom->molecule;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
 
  global_boundary(); 
 
  int *ivector = atom->ivector[index];

  for(i=0; i<nextra ; i++){
    pvector[i] = 0.0;
    }


  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < radm[itype][jtype]){
    
          rfourth=rsq*rsq;
          
          forcelj = -(lju3[itype][jtype]*sqrt(rsq) + lju4[itype][jtype]*rsq);
          
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

          fpair = flam * factor_lj * forcelj;
 
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
            evdwl = lju0[itype][jtype]+lju1[itype][jtype]*rsq+lju2[itype][jtype]*rfourth;
            evdwl *= factor_lj;
            if(scaling){
              pvector[m] += fDlam * evdwl;
              evdwl *= flam;
              }
          }

          if (evflag) ev_tally(i,j,nlocal,newton_pair,
                               evdwl,0.0,fpair,delx,dely,delz);
      }
      else if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);

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



        fpair =  flam * factor_lj * forcelj * r2inv;
 
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
          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
          offset[itype][jtype];
          evdwl *= factor_lj;
          if(scaling){
            pvector[m] += fDlam * evdwl;
            evdwl *= flam;
            }
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
//      else if (rsq < )
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 *  *   Global Boundary
 *   * ------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::global_boundary(){

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

int PairLJCleavCutNLInt::find_scaling(int imvec, int jmvec, int i, int j,double *xj){


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

void PairLJCleavCutNLInt::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  if(!pallocation) {nextra = n*n+1; // the zero-th position is not considered
    pvector = new double[nextra]; 
    pallocation = 1;}

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(gbox,natoms,"pair:gbox");
  memory->create(lam,n+1,n+1,"pair:lambda");
  memory->create(powlambda,n+1,n+1,"pair:powlambda");
  memory->create(powDlambda,n+1,n+1,"pair:powDlambda");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(radm,n+1,n+1,"pair:radm");
  memory->create(lju0,n+1,n+1,"pair:lju0");
  memory->create(lju1,n+1,n+1,"pair:lju1");
  memory->create(lju2,n+1,n+1,"pair:lju2");  
  memory->create(lju3,n+1,n+1,"pair:lju3");
  memory->create(lju4,n+1,n+1,"pair:lju4");  
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset"); 

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::settings(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
  lambda = utils::numeric(FLERR,arg[1],false,lmp);


  int n = strlen(arg[2]) + 1;
  idflag = new char[n];
  strcpy(idflag,arg[2]);

  int flag= -1, ncols;
  index = atom->find_custom( idflag, flag, ncols);



  if(flag != 0) error->all(FLERR,"Custom property must be integer for lj pair_style command");

  if(strcmp(arg[3],"x") != 0 && strcmp(arg[3],"y") != 0 && strcmp(arg[3],"z") != 0)
     error->all(FLERR,"Illegal LJ-BJ wells_command (direction)");

// For Now working only with NVT/NVE, because the size of the box is taken 
// at the beginning at the simulation

    if(strcmp(arg[3],"x") == 0){
        ind_dir=0;
    }
    else if(strcmp(arg[3],"y") == 0){
        ind_dir=1;
    }
    else if(strcmp(arg[3],"z") == 0){
        ind_dir=2;
    }
    

    npow =  utils::numeric(FLERR,arg[4],false,lmp);

    rspace = utils::numeric(FLERR,arg[5],false,lmp);


  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++){
        if (setflag[i][j]) cut[i][j] = cut_global;
        lam[i][j] = lambda;
    }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

/* IMPORTANT:
      The order of the argument is essential. You cannot put the cut_off for the pairs
      if you don't insert lambda as well, otherwise the cut-off will be read as lambda. 
      
      If you want to insert the value of lambda *only* then you can avoid to include the
      cut-off since it is not read  
*/

void PairLJCleavCutNLInt::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients A");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);

  double l = lambda;
  if (narg > 4 )   l = utils::numeric(FLERR,arg[4],false,lmp);

  double cut_one = cut_global;
  if (narg == 6) cut_one = utils::numeric(FLERR,arg[5],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      lam[i][j] = l;
      radm[i][j] = sigma[i][j]*rspace;
      powlambda[i][j] = pow(l,npow);
      powDlambda[i][j] = npow*pow(l,npow-1.0);
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients B");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::init_style()
{
  // request regular or rRESPA neighbor lists

  int irequest;

    if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
        error->all(FLERR,"There is no respa implementation yet");
    }

    irequest = neighbor->request(this,instance_me);


}

/* ----------------------------------------------------------------------
   Coefficients of the interpolation
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::calculate_int_coefficients(int i, int j)
{

    double r2inv,r6inv,philj,dphilj,ddphilj;

 
  
  r2inv = 1.0/radm[i][j]/radm[i][j];
  r6inv = r2inv*r2inv*r2inv;
  
  philj = r6inv*4*epsilon[i][j]*( pow(sigma[i][j],12.0)*r6inv-pow(sigma[i][j],6.0));
  dphilj = -24.0*epsilon[i][j]*r6inv*(2.0*pow(sigma[i][j],12.0)*r6inv-pow(sigma[i][j],6.0))/radm[i][j];
  ddphilj = 24.0*epsilon[i][j]*r6inv*r2inv*(26.0*pow(sigma[i][j],12.0)*r6inv-7*pow(sigma[i][j],6.0));

  
//a3=rc*cljp-1/3*rc^2*cljpp;
//a4=rc^2/4*cljpp-1/2*rc*cljp;
//a0 = clj-a3 -a4;

  
    lju1[i][j] = radm[i][j]*dphilj - 1.0/3.0*radm[i][j]*radm[i][j]*ddphilj;
    lju2[i][j] = 0.25*radm[i][j]*radm[i][j]*ddphilj - 0.5*radm[i][j]*dphilj;  
    lju0[i][j] = philj - lju1[i][j] - lju2[i][j];
    lju3[i][j] = 3.0*lju1[i][j]/pow(radm[i][j],3.0);
    lju4[i][j] = 4.0*lju2[i][j]/pow(radm[i][j],4.0);
    

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCleavCutNLInt::init_one(int i, int j)
{



  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j]   = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j]     = mix_distance(cut[i][i],cut[j][j]); 
    radm[i][j]    = mix_distance(radm[i][i],radm[j][j]); 
  }


  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  
  
  calculate_int_coefficients(i, j);
 
    
  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = powlambda[i][j] * 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lju0[j][i] = lju0[i][j];
  lju1[j][i] = lju1[i][j];
  lju2[j][i] = lju2[i][j];
  lju3[j][i] = lju3[i][j];
  lju4[j][i] = lju4[i][j];
  
  
  radm[j][i] = radm[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  lam[j][i] = lam[i][j];
  offset[j][i] = offset[i][j];
  cutsq[i][j] = cut[i][j]*cut[i][j];
  cutsq[j][i] =  cutsq[i][j];

  powlambda[j][i] = powlambda[i][j];
  powDlambda[j][i] = powDlambda[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&lam[i][j],sizeof(double),1,fp);
        fwrite(&powlambda[i][j],sizeof(double),1,fp);
        fwrite(&powDlambda[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error); 
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&lam[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&powlambda[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&powDlambda[i][j],sizeof(double),1,fp,nullptr,error);
              }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lam[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&powlambda[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&powDlambda[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&rspace,sizeof(double),1,fp);  
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
  fwrite(&ind_dir,sizeof(int),1,fp);
  fwrite(&npow,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&rspace,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&tail_flag,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&ind_dir,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&npow,sizeof(double),1,fp,nullptr,error);  
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&rspace,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ind_dir,1,MPI_INT,0,world);
  MPI_Bcast(&npow,1,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,epsilon[i][i],sigma[i][i],lam[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCleavCutNLInt::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j],lam[i][j] );
}

/* ---------------------------------------------------------------------- */

/*
    Scaling not considered. All pair potentials here are scaled by lambda
    
*/

double PairLJCleavCutNLInt::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r2inv,r6inv,forcelj,philj,rfourth;


  if (rsq < radm[itype][jtype]){
      rfourth=rsq*rsq;
      forcelj = lju3[itype][jtype] + lju4[itype][jtype]*rsq;
      fforce = npow*pow(lam[itype][jtype],npow-1.0) * factor_lj * forcelj;
      philj =  lju0[itype][jtype]+lju1[itype][jtype]*rsq+lju2[itype][jtype]*rfourth;
      philj = pow(lam[itype][jtype],npow)*philj;
  }
  else {
  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  fforce = npow*pow(lam[itype][jtype],npow-1.0)*factor_lj*forcelj*r2inv;

  philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
    offset[itype][jtype];
  philj = pow(lam[itype][jtype],npow)*philj;
    }
  return factor_lj*philj;
}


/* ---------------------------------------------------------------------- */

void *PairLJCleavCutNLInt::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return NULL;
}
