// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Amalie Frischknecht and Ahmed Ismail (SNL)
   simpler force assignment added by Rolf Isele-Holder (Aachen University)
------------------------------------------------------------------------- */

#include "pair_lj_cut_tip4p_short_cleav.h"

#include "compute_chunk_atom.h"
#include "compute_com_chunk.h"
#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "modify.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCutTIP4PShortCleav::PairLJCutTIP4PShortCleav(LAMMPS *lmp) :
  PairLJCutCoulLong(lmp)
{
  tip4pflag = 1;
  ewaldflag = pppmflag = 1;  // for clarity, though inherited from parent class

  single_enable = 0;
  respa_enable = 0;
  writedata = 1;

  nmax = 0;
  hneigh = nullptr;
  newsite = nullptr;

  // TIP4P cannot compute virial as F dot r
  // due to finding bonded H atoms which are not near O atom

  no_virial_fdotr_compute = 1;
  
  dubtypes=ntypes*2;


  pallocation = 0;
  natoms=natoms+1;
  if(ntypes > 0){
    nextra = dubtypes*dubtypes+1; // the zero-th position is not considered
    pvector = new double[nextra];
    pallocation = 1;
    }

   xprd = domain->xprd;
   yprd = domain->yprd;
   zprd = domain->zprd;
   xy = domain->xy;
   yz = domain->yz;
   xz = domain->xz;

  nchunk = 0;  
}

/* ---------------------------------------------------------------------- */

PairLJCutTIP4PShortCleav::~PairLJCutTIP4PShortCleav()
{
  memory->destroy(hneigh);
  memory->destroy(newsite);
  
  if (allocated) {  
    memory->destroy(lam);
    memory->destroy(powlambda);
    memory->destroy(powDlambda);
    memory->destroy(gbox);
    memory->destroy(lamC);
    memory->destroy(powlambdaC);
    memory->destroy(powDlambdaC);
    memory->destroy(giflag);
    }
  
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(lam,n+1,n+1,"pair:lambda");
  memory->create(powlambda,n+1,n+1,"pair:powlambda");
  memory->create(powDlambda,n+1,n+1,"pair:powDlambda");
  memory->create(lamC,n+1,n+1,"pair:lambdaC");
  memory->create(powlambdaC,n+1,n+1,"pair:powlambdaC");
  memory->create(powDlambdaC,n+1,n+1,"pair:powDlambdaC");
  memory->create(giflag,natoms,"pair:giflag");  
  
  
}

/* ---------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable,key;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double fraction,table;
  double r,r2inv,r6inv,forcecoul,forcelj,cforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double fO[3],fH[3],fd[3],v[6];
  double *x1,*x2,*xH1,*xH2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  // reallocate hneigh & newsite if necessary
  // initialize hneigh[0] to -1 on steps when reneighboring occurred
  // initialize hneigh[2] to 0 every step

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(hneigh);
    memory->create(hneigh,nmax,3,"pair:hneigh");
    memory->destroy(newsite);
    memory->create(newsite,nmax,3,"pair:newsite");
  }
  if (neighbor->ago == 0)
    for (i = 0; i < nall; i++) hneigh[i][0] = -1;
  for (i = 0; i < nall; i++) hneigh[i][2] = 0;

  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  tagint *tag = atom->tag;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double cut_coulsqplus = (cut_coul+2.0*qdist) * (cut_coul+2.0*qdist);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  tagint *molecule = atom->molecule;
  int t1, t2, k1, k2, scaling,m;  
  double flam,fDlam;
  global_boundary();

  for(i=0; i<nextra ; i++){
    pvector[i] = 0.0;
    }

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    // if atom I = water O, set x1 = offset charge site
    // else x1 = x of atom I

    if (itype == typeO) {
      if (hneigh[i][0] < 0) {
        iH1 = atom->map(tag[i] + 1);
        iH2 = atom->map(tag[i] + 2);
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        // set iH1,iH2 to closest image to O
        iH1 = domain->closest_image(i,iH1);
        iH2 = domain->closest_image(i,iH2);
        compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
        hneigh[i][0] = iH1;
        hneigh[i][1] = iH2;
        hneigh[i][2] = 1;

      } else {
        iH1 = hneigh[i][0];
        iH2 = hneigh[i][1];
        if (hneigh[i][2] == 0) {
          hneigh[i][2] = 1;
          compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
        }
      }
      x1 = newsite[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // LJ interaction based on true rsq

      if (rsq < cut_ljsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        
        
        scaling=0;
         k1=molecule[i];
         k2=molecule[j];
         t1=tag[i];
         t2=tag[j];
         if (switchlj){
         	if(k1 != k2)
                	scaling = find_scaling(giflag[t1],giflag[t2],i,j,x[j]);
			}
          flam = 1.0;
          fDlam = 1.0;

          if(scaling){
              flam = powlambda[itype][jtype];
              fDlam =  powDlambda[itype][jtype];
            }        
        
        
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        forcelj *= flam * factor_lj * r2inv;

        f[i][0] += delx*forcelj;
        f[i][1] += dely*forcelj;
        f[i][2] += delz*forcelj;
        f[j][0] -= delx*forcelj;
        f[j][1] -= dely*forcelj;
        f[j][2] -= delz*forcelj;

        if (eflag) {
          m=jtype+(itype-1)*dubtypes;
          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
            offset[itype][jtype];
          evdwl *= factor_lj;
		  if(scaling){
               pvector[m] += fDlam * evdwl;
               evdwl *= flam;
            }
        } else evdwl = 0.0;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,forcelj,delx,dely,delz);
      }

      // adjust rsq and delxyz for off-site O charge(s) if necessary
      // but only if they are within reach

      if (rsq < cut_coulsqplus) {
        if (itype == typeO || jtype == typeO) {

          // if atom J = water O, set x2 = offset charge site
          // else x2 = x of atom J

          if (jtype == typeO) {
            if (hneigh[j][0] < 0) {
              jH1 = atom->map(tag[j] + 1);
              jH2 = atom->map(tag[j] + 2);
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              // set jH1,jH2 to closest image to O
              jH1 = domain->closest_image(j,jH1);
              jH2 = domain->closest_image(j,jH2);
              compute_newsite(x[j],x[jH1],x[jH2],newsite[j]);
              hneigh[j][0] = jH1;
              hneigh[j][1] = jH2;
              hneigh[j][2] = 1;

            } else {
              jH1 = hneigh[j][0];
              jH2 = hneigh[j][1];
              if (hneigh[j][2] == 0) {
                hneigh[j][2] = 1;
                compute_newsite(x[j],x[jH1],x[jH2],newsite[j]);
              }
            }
            x2 = newsite[j];
          } else x2 = x[j];

          delx = x1[0] - x2[0];
          dely = x1[1] - x2[1];
          delz = x1[2] - x2[2];
          rsq = delx*delx + dely*dely + delz*delz;
        }

        // Coulombic interaction based on modified rsq

        if (rsq < cut_coulsq) {
          r2inv = 1 / rsq;
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) {
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }

/* This is a repetition of what we wrote for the LJ part, but it is needed 
 * if coul_cut_off > lj_cut_off then we are not calculating the scaling for some pairs*/

         scaling=0;
         k1=molecule[i];
         k2=molecule[j];
         t1=tag[i];
         t2=tag[j];
         if (switchcoul){
         	if(k1 != k2)
                	scaling = find_scaling(giflag[t1],giflag[t2],i,j,x[j]);
	      }
	      
          flam = 1.0;
          fDlam = 1.0;

// It can be avoided if we enforce coul_cut_off <= lj_cut_off

          if(scaling){
              flam = powlambdaC[itype][jtype];
              fDlam =  powDlambdaC[itype][jtype];
            }



          cforce = flam*forcecoul * r2inv;

          // if i,j are not O atoms, force is applied directly
          // if i or j are O atoms, force is on fictitious atom & partitioned
          // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
          // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
          // preserves total force and torque on water molecule
          // virial = sum(r x F) where each water's atoms are near xi and xj
          // vlist stores 2,4,6 atoms whose forces contribute to virial

          n = 0;
          key = 0;

          if (itype != typeO) {
            f[i][0] += delx * cforce;
            f[i][1] += dely * cforce;
            f[i][2] += delz * cforce;

            if (vflag) {
              v[0] = x[i][0] * delx * cforce;
              v[1] = x[i][1] * dely * cforce;
              v[2] = x[i][2] * delz * cforce;
              v[3] = x[i][0] * dely * cforce;
              v[4] = x[i][0] * delz * cforce;
              v[5] = x[i][1] * delz * cforce;
            }
            vlist[n++] = i;

          } else {
            key++;

            fd[0] = delx*cforce;
            fd[1] = dely*cforce;
            fd[2] = delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[i][0] += fO[0];
            f[i][1] += fO[1];
            f[i][2] += fO[2];

            f[iH1][0] += fH[0];
            f[iH1][1] += fH[1];
            f[iH1][2] += fH[2];

            f[iH2][0] += fH[0];
            f[iH2][1] += fH[1];
            f[iH2][2] += fH[2];

            if (vflag) {
              xH1 = x[iH1];
              xH2 = x[iH2];
              v[0] = x[i][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
              v[1] = x[i][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
              v[2] = x[i][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
              v[3] = x[i][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
              v[4] = x[i][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
              v[5] = x[i][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
            }
            vlist[n++] = i;
            vlist[n++] = iH1;
            vlist[n++] = iH2;
          }

          if (jtype != typeO) {
            f[j][0] -= delx * cforce;
            f[j][1] -= dely * cforce;
            f[j][2] -= delz * cforce;

            if (vflag) {
              v[0] -= x[j][0] * delx * cforce;
              v[1] -= x[j][1] * dely * cforce;
              v[2] -= x[j][2] * delz * cforce;
              v[3] -= x[j][0] * dely * cforce;
              v[4] -= x[j][0] * delz * cforce;
              v[5] -= x[j][1] * delz * cforce;
            }
            vlist[n++] = j;

          } else {
            key += 2;

            fd[0] = -delx*cforce;
            fd[1] = -dely*cforce;
            fd[2] = -delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[j][0] += fO[0];
            f[j][1] += fO[1];
            f[j][2] += fO[2];

            f[jH1][0] += fH[0];
            f[jH1][1] += fH[1];
            f[jH1][2] += fH[2];

            f[jH2][0] += fH[0];
            f[jH2][1] += fH[1];
            f[jH2][2] += fH[2];

            if (vflag) {
              xH1 = x[jH1];
              xH2 = x[jH2];
              v[0] += x[j][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
              v[1] += x[j][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
              v[2] += x[j][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
              v[3] += x[j][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
              v[4] += x[j][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
              v[5] += x[j][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
            }
            vlist[n++] = j;
            vlist[n++] = jH1;
            vlist[n++] = jH2;
          }

          if (eflag) {
 			m=(jtype+ntypes)+(ntypes+itype-1)*dubtypes;
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
            if(scaling){
//fprintf(fp,"Scaled %d %d %f %f \n",tag[i],tag[j],ecoul,pvector[m]);
               pvector[m] += fDlam * ecoul;
               ecoul *= flam;
//fprintf(fp,"Scaled %f %f %f \n",fDlam,ecoul,pvector[m]);
            }
          } else ecoul = 0.0;

          if (evflag) ev_tally_tip4p(key,vlist,v,ecoul,alpha);
        }
      }
    }
  }
}



/* ----------------------------------------------------------------------
   Global Boundary
------------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::global_boundary(){

int i,k,m;
int ibox[3];
int nlocal = atom->nlocal;
int loc_box[natoms],loc_iflag[natoms];
double xx[3];
double **x = atom->x;

imageint *image = atom->image;
imageint locim;
tagint *tag = atom->tag;

ccom->compute_array();
//nchunk = cchunk->nchunk;

double **com = ccom->array;
int *ichunk = cchunk->ichunk;

     for(i=0; i<natoms; i++){
        loc_box[i]=0;
        loc_iflag[i]=0;
    }


     for(i=0; i<nlocal; i++){
         m = ichunk[i]-1;
         k=tag[i];

         
         if(com[m][ind_dir] > lboxhalf )loc_iflag[k]=1;

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

//fprintf(fp,"GB %d %d %f %d %d \n",k,i,x[i][2],loc_iflag[k],loc_box[k]);

        }


 MPI_Allreduce(&loc_box[0],&gbox[0],natoms,MPI_INT,MPI_SUM,world);
 MPI_Allreduce(&loc_iflag[0],&giflag[0],natoms,MPI_INT,MPI_SUM,world);

}


/* ----------------------------------------------------------------------
    Find scaling
------------------------------------------------------------------------- */

int PairLJCutTIP4PShortCleav::find_scaling(int imvec, int jmvec, int i, int j,double *xj){


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
   global settings
------------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::settings(int narg, char **arg)
{
  if (narg < 10 || narg > 14) error->all(FLERR,"Illegal pair_style command");

  typeO = utils::inumeric(FLERR,arg[0],false,lmp);
  typeH = utils::inumeric(FLERR,arg[1],false,lmp);
  typeB = utils::inumeric(FLERR,arg[2],false,lmp);
  typeA = utils::inumeric(FLERR,arg[3],false,lmp);
  qdist = utils::numeric(FLERR,arg[4],false,lmp);

  lambda = utils::numeric(FLERR,arg[5],false,lmp);

  if(strcmp(arg[6],"x") != 0 && strcmp(arg[6],"y") != 0 && strcmp(arg[6],"z") != 0)
     error->all(FLERR,"Illegal LJ-BJ wells_command (direction)");

    if(strcmp(arg[6],"x") == 0){
        ind_dir=0;
        lboxhalf=xprd/2.0;
    }
    else if(strcmp(arg[6],"y") == 0){
        ind_dir=1;
        lboxhalf=yprd/2.0;
    }
    else if(strcmp(arg[6],"z") == 0){
        ind_dir=2;
        lboxhalf=zprd/2.0;
    }


    npow =  utils::numeric(FLERR,arg[7],false,lmp);
    idchunk = utils::strdup(arg[8]);
    idcom = utils::strdup(arg[9]);
    
    

  

cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
if (!cchunk)
  error->all(FLERR,"Chunk/atom compute {} does not exist or is not chunk/atom style", idchunk);

ccom = dynamic_cast<ComputeCOMChunk *>(modify->get_compute_by_id(idcom));
if (!ccom)
  error->all(FLERR,"Com/chunk compute {} does not exist or is not com/chunk style", idcom);

// check that idchunk is consistent with ccom->idchunk

if (ccom && (strcmp(idchunk,ccom->idchunk) != 0))
         error->all(FLERR,"Fix spring/chunk chunk ID {} not the same as compute com/chunk chunk ID {}", idchunk, ccom->idchunk);




  cut_lj_global = utils::numeric(FLERR,arg[10],false,lmp);
  cut_coul = utils::numeric(FLERR,arg[11],false,lmp);


  switchcoul = 1; // default switch: on 
  switchlj   = 1;

  if (narg > 12){ 
    switchcoul = utils::numeric(FLERR,arg[12],false,lmp);  
    switchlj = utils::numeric(FLERR,arg[13],false,lmp);
	}

  cut_coulsq = cut_coul * cut_coul;

//printf("ALLOCATED %d\n",allocated); 

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {cut_lj[i][j] = cut_lj_global;
        	lam[i][j] = lambda;
        	lamC[i][j]= lambda;
     }
  }
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style lj/cut/tip4p/long requires atom IDs");
  if (!force->newton_pair)
    error->all(FLERR,
               "Pair style lj/cut/tip4p/long requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,
               "Pair style lj/cut/tip4p/long requires atom attribute q");
  if (force->bond == nullptr)
    error->all(FLERR,"Must use a bond style with TIP4P potential");
  if (force->angle == nullptr)
    error->all(FLERR,"Must use an angle style with TIP4P potential");

  PairLJCutCoulLong::init_style();

  // set alpha parameter

  const double theta = force->angle->equilibrium_angle(typeA);
  const double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (cos(0.5*theta) * blen);

  const double mincut = cut_coul + qdist + blen + neighbor->skin;
  if (comm->get_comm_cutoff() < mincut) {
    if (comm->me == 0)
      error->warning(FLERR, "Increasing communication cutoff to {:.8} for TIP4P pair style",
                     mincut);
    comm->cutghostuser = mincut;
  }
}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);

  double l = lambda;
  if (narg > 4 )   l = utils::numeric(FLERR,arg[4],false,lmp);

  lambdaC = lambda;
  if (narg == 6)  lambdaC = utils::numeric(FLERR,arg[5],false,lmp);

  double cut_lj_one = cut_lj_global;
  if (narg == 7) cut_lj_one = utils::numeric(FLERR,arg[6],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
    
//printf("BBBB %d %d %f %d  \n",i,j);    
    
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_lj[i][j] = cut_lj_one;
      lam[i][j] = l;
      lamC[i][j] = lambdaC;
      powlambda[i][j] = pow(l,npow);
      powDlambda[i][j] = npow*pow(l,npow-1.0);
      powlambdaC[i][j] = pow(lambdaC,npow);
      powDlambdaC[i][j] = npow*pow(lambdaC,npow-1.0);
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutTIP4PShortCleav::init_one(int i, int j)
{
  double cut = PairLJCutCoulLong::init_one(i,j);
  

  lam[j][i] = lam[i][j];
  powlambda[j][i] = powlambda[i][j];
  powDlambda[j][i] = powDlambda[i][j];

  lamC[j][i] = lamC[i][j];
  powlambdaC[j][i] = powlambdaC[i][j];
  powDlambdaC[j][i] = powDlambdaC[i][j];

  // check that LJ epsilon = 0.0 for water H
  // set LJ cutoff to 0.0 for any interaction involving water H
  // so LJ term isn't calculated in compute()

  if ((i == typeH && epsilon[i][i] != 0.0) ||
      (j == typeH && epsilon[j][j] != 0.0))
    error->all(FLERR,"Water H epsilon must be 0.0 for "
               "pair style lj/cut/tip4p/short/cleav");

  if (i == typeH || j == typeH)
    cut_ljsq[j][i] = cut_ljsq[i][j] = 0.0;
  


  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::write_restart_settings(FILE *fp)
{
  fwrite(&typeO,sizeof(int),1,fp);
  fwrite(&typeH,sizeof(int),1,fp);
  fwrite(&typeB,sizeof(int),1,fp);
  fwrite(&typeA,sizeof(int),1,fp);
  fwrite(&qdist,sizeof(double),1,fp);

  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&typeO,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&typeH,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&typeB,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&typeA,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&qdist,sizeof(double),1,fp,nullptr,error);

    utils::sfread(FLERR,&cut_lj_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&ncoultablebits,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&tabinner,sizeof(double),1,fp,nullptr,error);
  }

  MPI_Bcast(&typeO,1,MPI_INT,0,world);
  MPI_Bcast(&typeH,1,MPI_INT,0,world);
  MPI_Bcast(&typeB,1,MPI_INT,0,world);
  MPI_Bcast(&typeA,1,MPI_INT,0,world);
  MPI_Bcast(&qdist,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void PairLJCutTIP4PShortCleav::compute_newsite(double *xO, double *xH1,
                                         double *xH2, double *xM)
{
  double delx1 = xH1[0] - xO[0];
  double dely1 = xH1[1] - xO[1];
  double delz1 = xH1[2] - xO[2];

  double delx2 = xH2[0] - xO[0];
  double dely2 = xH2[1] - xO[1];
  double delz2 = xH2[2] - xO[2];

  xM[0] = xO[0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = xO[1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = xO[2] + alpha * 0.5 * (delz1 + delz2);
}

/* ---------------------------------------------------------------------- */

void *PairLJCutTIP4PShortCleav::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"qdist") == 0) return (void *) &qdist;
  if (strcmp(str,"typeO") == 0) return (void *) &typeO;
  if (strcmp(str,"typeH") == 0) return (void *) &typeH;
  if (strcmp(str,"typeA") == 0) return (void *) &typeA;
  if (strcmp(str,"typeB") == 0) return (void *) &typeB;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of hneigh
------------------------------------------------------------------------- */

double PairLJCutTIP4PShortCleav::memory_usage()
{
  double bytes = (double)maxeatom * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  bytes += (double)2 * nmax * sizeof(double);
  return bytes;
}
