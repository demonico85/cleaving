/* ----------------------------------------------------------------------
   Contributing author: Nicodemo Di Pasquale (nicodemo.dipasquale@gmail.com) 
------------------------------------------------------------------------- */
/*
      __________________________________________________________________
      | A1        #  a2             |            a3      #         A4   |                               
      |           #                 |                    #              |
      |           #                 |                    #              |
      |           #                 |                    #              |
      |           #                 |                    #              |
      |           #                 |                    #              |
      |           #                 |                    #              |
      |           #                 |                    #              |
      |  a7       #   A5            |            A6      #        a8    |
      -------------------------------------------------------------------

Capital case: phase 1
Lower case: phase 2

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_lj_BGcleav.h"
#include "domain.h"
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

using namespace LAMMPS_NS;
using namespace MathConst;

// Definition coefficients 
// Davidchack, R. L. and Laird, B. B. Journal of Chemical Physics,  118, 7651:7657, (2003)

#define  C1  0.016132
#define  C2  3136.6
#define  C3 -68.069
#define  C4 -0.083312
#define  C5  0.74689


/* ---------------------------------------------------------------------- */

PairLJBGcleav::PairLJBGcleav(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  respa_enable  = 0;
  manybody_flag = 0;
  writedata     = 1; 
  restartinfo   = 0;
  one_coeff     = 1;
  pallocation   = 0;
  //eflag_global  = 1;

  int n=atom->ntypes;
  if(n > 0){
    nextra = n*n+1; // the zero-th position is not considered
    pvector = new double[nextra];
    pallocation = 1;
    }
    
    natoms=atom->natoms+1; //array does start from zero
    
    xprd = domain->xprd;
    yprd = domain->yprd;
    zprd = domain->zprd;
    xy = domain->xy;
    yz = domain->yz;
    xz = domain->xz;


}


/* ---------------------------------------------------------------------- */

PairLJBGcleav::~PairLJBGcleav()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(cutsq_in);
    memory->destroy(cut_in);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(lj5);
    memory->destroy(lj6);
    memory->destroy(lj7);
    memory->destroy(lj8);
    memory->destroy(lj9);
    memory->destroy(lj10);
    memory->destroy(lamcoeff);
    memory->destroy(c1);
    memory->destroy(c5);
    delete [] pvector;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJBGcleav::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

 if(!pallocation) {nextra = n*n+1; // the zero-th position is not considered
    pvector = new double[nextra]; 
    pallocation = 1;}

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cutsq_in,n+1,n+1,"pair:cutsq_in");
  memory->create(cut_in,n+1,n+1,"pair:cut2");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(lj5,n+1,n+1,"pair:lj5");
  memory->create(lj6,n+1,n+1,"pair:lj6");
  memory->create(lj7,n+1,n+1,"pair:lj7");
  memory->create(lj8,n+1,n+1,"pair:lj8");
  memory->create(lj9,n+1,n+1,"pair:lj9");
  memory->create(lj10,n+1,n+1,"pair:lj10");
  memory->create(c1,n+1,n+1,"pair:c1");
  memory->create(c5,n+1,n+1,"pair:c5");
  memory->create(gbox,natoms,"pair:gbox");
  memory->create(lamcoeff,n+1,n+1,"pair:lamcoeff");
  for(int i=0; i<n+1; i++ ){
     for(int j=i; j<n+1; j++){
        lamcoeff[i][j]=-10.0;
        lamcoeff[j][i]=-10.0;
       }
    }
}


/* ---------------------------------------------------------------------- */

void PairLJBGcleav::compute(int eflag, int vflag)
{

  int i,j,ii,jj,inum,jnum,itype,jtype,m,scaling;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj,flam;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int n = atom->ntypes;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  for(i=0; i<nextra ; i++){
    pvector[i] = 0.0;
    }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

    
    global_boundary();


//fpl =NULL;
//fpl=fopen("forces.log", "a");
/*
if(update->ntimestep == 1)if(comm->me == 0){fprintf(fpl," TIMESTEP %d  \n",update->ntimestep);}
if(update->ntimestep < 10){
//  for(i=0; i<n+1 ; i++){
//    for(j=0; j<n+1; j++)fprintf(fpl,"%d %d \n %f %f %f %f %f \n %f %f %f %f %f \n %f %f %f %f\n\n",i,j,lj1[i][j],lj2[i][j],lj3[i][j],lj4[i][j],lj5[i][j],lj6[i][j],lj7[i][j],lj8[i][j],lj9[i][j],lj10[i][j],cutsq[i][j],cutsq2[i][j],cut[i][j],cut2[i][j],lamcoeff[i][j]);
//    }}


for (i=0;i<nlocal;i++)
    fprintf(fpl,"%d %f %f %f\n",i,f[i][0],f[i][1],f[i][2]);
}


MPI_Barrier(world); */

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


      if (rsq <= cutsq_in[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);


          scaling=0;
        flam = 1.0;

          if(itype != jtype)
               scaling = find_scaling(itype,jtype,i,j,x[j]);
          
          if(scaling)
            flam = lamcoeff[itype][jtype];
//if(itype != jtype) fprintf(fpl,"A %d %d %d \n %f %f\n",itype,jtype,scaling,ztmp,x[j][2]); 

        fpair = flam*factor_lj*forcelj*r2inv;

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
            evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) + c1[itype][jtype];
            evdwl *= factor_lj;
            if(scaling)pvector[m] += evdwl;
            evdwl *= flam; 
          }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);

        }
      else if (rsq < cutsq[itype][jtype]){

        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj5[itype][jtype]*r6inv + lj6[itype][jtype]);


        flam = 1.0;

          if(itype != jtype)
               scaling = find_scaling(itype,jtype,i,j,x[j]);
          
        if(scaling)
            flam = lamcoeff[itype][jtype];

//if(itype != jtype) fprintf(fpl,"B %d %d %d \n %f %f\n",itype,jtype,scaling,ztmp,x[j][2]);
        fpair = flam*factor_lj*(forcelj*r2inv - lj7[itype][jtype]);

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
            evdwl = r6inv*(lj8[itype][jtype]*r6inv + lj9[itype][jtype]) + lj10[itype][jtype]*rsq + c5[itype][jtype];
            evdwl *= factor_lj;
            if(scaling)pvector[m] += evdwl;
            evdwl *= flam; 
            }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);

      }
    }
  }
//fclose (fpl);

  if (vflag_fdotr) virial_fdotr_compute();


}

/* ----------------------------------------------------------------------
 *  *   Global Boundary
 *   * ------------------------------------------------------------------------- */

void PairLJBGcleav::global_boundary(){

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
 *   find scaling
 *   ------------------------------------------------------------------------- */

int PairLJBGcleav::find_scaling(int imvec, int jmvec, int i, int j,double *xj){


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
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJBGcleav::write_restart(FILE *fp)
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
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJBGcleav::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void  PairLJBGcleav::write_restart_settings(FILE *fp)
{

  fwrite(&cut_global_in,sizeof(double),1,fp);
  fwrite(&cut_global_out,sizeof(int),1,fp);

}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void  PairLJBGcleav::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global_in,sizeof(double),1,fp);
    fread(&cut_global_out,sizeof(double),1,fp);
  }
  MPI_Bcast(&cut_global_in,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global_out,1,MPI_INT,0,world);
}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJBGcleav::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");


  cut_global_in  = force->numeric(FLERR,arg[0]);
  cut_global_out = force->numeric(FLERR,arg[1]);
  lambda      = force->numeric(FLERR,arg[2]);
  cleavwall = force->numeric(FLERR,arg[3]);  

  // reset cutoffs that have been explicitly set

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

    
    
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
      if (setflag[i][j]){
           cut[i][j] = cut_global_out;
           cut_in[i][j] = cut_global_in;
      }
  }

    
    


}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
/*

come funziona setflag:

se nell file .in sono specificati i,j allora setta le mixed sigma,eps,cut tutti 
ai valori indicati nel file per tutte le coppie i,j. 
Altrimenti per quella coppia setflag[i][j] = 0 e vengono calcolate (dopo) le mixed interactions


narg          1     2        3       4         5        6       7 
type          I     J      epsilon  sigma     cut1     cut2      
pair_coeff    1     2       1.0      1.0      1.1      0.25      1.1

*/

void PairLJBGcleav::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (narg == 6)
    error->all(FLERR,"Incorrect args for pair coefficients (you need two cut-offs)");


  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  double lm = lambda;
  if (narg > 4)lm = force->numeric(FLERR,arg[4]);

  double cutinner = cut_global_in;
   if (narg > 5) cutinner = force->numeric(FLERR,arg[5]);
  double cutouter = cut_global_out;
  if (narg == 7) cutouter = force->numeric(FLERR,arg[6]);


  if(cutinner > cutouter)
         error->all(FLERR,"First cut_off must be smaller than second one in the definition");

//  printf("COEFF %s %s %d %f %f\n",arg[0],arg[1],narg,lm,ftr);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_in[i][j] = cutinner;
      cut[i][j] = cutouter;
      setflag[i][j] = 1;
      lamcoeff[i][j] = lm;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");


}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJBGcleav::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJBGcleav::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

/*

lj1 = coefficient Force r^-12 first cut-off
lj2 = coefficient Force r^-6 first cut-off
lj3 = coefficient Potential r^-12 first cut-off
lj4 = coefficient Potential r^-6 first cut-off
lj5 = coefficient Force r^-12 second cut-off
lj6 = coefficient Force r^-6 second cut-off
lj7 = coefficient Force r^-2 second cut-off
lj8 = coefficient Potential r^-12 second cut-off
lj9 = coefficient Potential r^-6 second cut-off
lj10 = coefficient Force r^-2 second cut-off

*/

double PairLJBGcleav::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
    cut_in[i][j] = mix_distance(cut_in[i][i],cut_in[j][j]);
  }


  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj5[i][j] = 12.0 * C2 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj6[i][j] =  6.0 * C3 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj7[i][j] =  2.0 * C4 * epsilon[i][j] * 1./pow(sigma[i][j],2.0);
  lj8[i][j] =  C2 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj9[i][j] =  C3 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj10[i][j] = C4 * epsilon[i][j] * 1./pow(sigma[i][j],2.0);
    cutsq_in[i][j] = cut_in[i][j]*cut_in[i][j];


  lj1[j][i]  = lj1[i][j];
  lj2[j][i]  = lj2[i][j];
  lj3[j][i]  = lj3[i][j];
  lj4[j][i]  = lj4[i][j];
  lj5[j][i]  = lj5[i][j];
  lj6[j][i]  = lj6[i][j];
  lj7[j][i]  = lj7[i][j];
  lj8[j][i]  = lj8[i][j];
  lj9[j][i]  = lj9[i][j];
  lj10[j][i] = lj10[i][j];


    cutsq_in[j][i] = cutsq_in[i][j];
  lamcoeff[j][i] = lamcoeff[i][j];

  c1[i][j] = C1*epsilon[i][j];
  c1[j][i] = c1[i][j];
  c5[i][j] = C5*epsilon[i][j];
  c5[j][i] = c5[i][j];


  return cut[i][j];
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJBGcleav::init_style()
{

  // request regular neighbor lists

  int irequest;

    if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
        error->all(FLERR,"There is no respa implementation yet");
    }

    irequest = neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
 *   
 *   Single subroutine 
 *
 * ------------------------------------------------------------------------- */


double PairLJBGcleav::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r2inv,r6inv,forcelj,philj;


  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;

  if( rsq < cut_global_in*cut_global_in){
  	philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) + c1[itype][jtype];
	forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  }
  else if(rsq < cut_global_out*cut_global_out){
        philj   = r6inv*(lj8[itype][jtype]*r6inv + lj9[itype][jtype]) + lj10[itype][jtype]*rsq + c5[itype][jtype];
        forcelj = r6inv * (lj5[itype][jtype]*r6inv + lj6[itype][jtype]);
	}

  fforce = lamcoeff[itype][jtype]*factor_lj*forcelj*r2inv;

  return lamcoeff[itype][jtype]*factor_lj*philj;
}






















