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

#include "pair_lj_cut_tip4p_long_cleav_types_splitlj.h"

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

PairLJCutTIP4PLongCleavTypesSplitLJ::PairLJCutTIP4PLongCleavTypesSplitLJ(LAMMPS *lmp) :
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
  
  ntypes=atom->ntypes;
  natoms=atom->natoms; //array does start from zero



  if(natoms == 0){
      error->all(FLERR,"Pair style must be declared after box generation");
  }
  

  pallocation = 0;
  natoms=natoms+1;
  dubtypes=ntypes*2;
      
  if(ntypes > 0){
    nextra = dubtypes*dubtypes+1; // the zero-th position is not considered
    pvector = new double[nextra];
    pallocation = 1;
    }


}

/* ---------------------------------------------------------------------- */

PairLJCutTIP4PLongCleavTypesSplitLJ::~PairLJCutTIP4PLongCleavTypesSplitLJ()
{
  memory->destroy(hneigh);
  memory->destroy(newsite);
  
  
  if (allocated) {  
    memory->destroy(lam);
    memory->destroy(powlambdaLJb);
    memory->destroy(powDlambdaLJb);
    memory->destroy(powlambdaLJa);
    memory->destroy(powDlambdaLJa);
    memory->destroy(scalingC);
    memory->destroy(lamC);
    memory->destroy(powlambdaC);
    memory->destroy(powDlambdaC);
    memory->destroy(scalingLJ);
    memory->destroy(SCoffset);
    memory->destroy(SCDoffset);
    memory->destroy(unSCoffset);
    }
    
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongCleavTypesSplitLJ::allocate()
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
  memory->create(powlambdaLJa,n+1,n+1,"pair:powlambdaLJa");
  memory->create(powDlambdaLJa,n+1,n+1,"pair:powDlambdaLJa");
  memory->create(powlambdaLJb,n+1,n+1,"pair:powlambdaLJb");
  memory->create(powDlambdaLJb,n+1,n+1,"pair:powDlambdaLJb");
  memory->create(lamC,n+1,n+1,"pair:lambdaC");
  memory->create(powlambdaC,n+1,n+1,"pair:powlambdaC");
  memory->create(powDlambdaC,n+1,n+1,"pair:powDlambdaC");
  memory->create(scalingC,n+1,n+1,"pair:powlambdaC");
  memory->create(scalingLJ,n+1,n+1,"pair:powDlambdaC");  
  memory->create(SCoffset,n+1,n+1,"pair:SCoffset");
  memory->create(SCDoffset,n+1,n+1,"pair:SCoffset");
  memory->create(unSCoffset,n+1,n+1,"pair:unSCoffset");  

  ntypes=n;
  dubtypes=n*2;
  if(!pallocation){
    nextra = dubtypes*dubtypes+1; // the zero-th position is not considered
    pvector = new double[nextra];
    pallocation = 1;
    }

}


/* ---------------------------------------------------------------------- */

void PairLJCutTIP4PLongCleavTypesSplitLJ::compute(int eflag, int vflag)
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
  

  int t1, t2, scaling,m;  
  double flamA,fDlamA,flamB,fDlamB;
    

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
        
        
          flamA  = 1.0;
          fDlamA = 1.0;
          flamB  = 1.0;
          fDlamB = 1.0;
          
          scaling = scalingLJ[itype][jtype];

          if(scaling){
              flamA  = powlambdaLJa[itype][jtype];
              fDlamA =  powDlambdaLJa[itype][jtype];
              flamB  = powlambdaLJb[itype][jtype];
              fDlamB =  powDlambdaLJb[itype][jtype];              
            }     
            

    
        forcelj = r6inv * (flamA*lj1[itype][jtype]*r6inv - flamB*lj2[itype][jtype]);
        forcelj *= factor_lj * r2inv;
                
//        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
//        forcelj *= flamA * factor_lj * r2inv;
        
// if (comm->me == 0)   printf("%d %d %d %d %f %f \n LJ %f \n ",i,j,itype,jtype, flamA,flamB,forcelj);         

        f[i][0] += delx*forcelj;
        f[i][1] += dely*forcelj;
        f[i][2] += delz*forcelj;
        f[j][0] -= delx*forcelj;
        f[j][1] -= dely*forcelj;
        f[j][2] -= delz*forcelj;

 

 
        if (eflag) {        
	      if(scaling){
	        m=jtype+(itype-1)*dubtypes;  
            pvector[m] += r6inv*(fDlamA * lj3[itype][jtype]*r6inv- fDlamB * lj4[itype][jtype]) -
            	SCDoffset[itype][jtype];
            evdwl = r6inv*(flamA * lj3[itype][jtype] * r6inv- flamB * lj4[itype][jtype]) -
            	SCoffset[itype][jtype];
                        }
            else{                 
        		evdwl = r6inv*( lj3[itype][jtype] * r6inv - lj4[itype][jtype]) -
            		unSCoffset[itype][jtype];
            	evdwl *= factor_lj;	          
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

        
          flamA = 1.0;
          fDlamA = 1.0;
          
          scaling = scalingC[itype][jtype];

          if(scaling){
              flamA = powlambdaC[itype][jtype];
              fDlamA =  powDlambdaC[itype][jtype];
            }     
            


          cforce = flamA*forcecoul * r2inv;
 // if (comm->me == 0)  printf("%d %d %d %d %f %f \n COUL %f %f %f\n ",i,j,itype,jtype,  powlambda[itype][jtype],powDlambda[itype][jtype],cut_coulsq,cut_coulsqplus,cforce); 

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
               pvector[m] += fDlamA * ecoul;
               ecoul *= flamA;
//                                      if (comm->me == 0)   printf("%d %d %d %d SC %d \n C %f %f %f %f \n ",i,j,itype,jtype,scaling,ecoul,flam,fDlam,pvector[m] );  
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
   global settings
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongCleavTypesSplitLJ::settings(int narg, char **arg)
{
  if (narg < 8 || narg > 11)  error->all(FLERR,"Illegal pair_style command");
  

  typeO = utils::inumeric(FLERR,arg[0],false,lmp);
  typeH = utils::inumeric(FLERR,arg[1],false,lmp);
  typeB = utils::inumeric(FLERR,arg[2],false,lmp);
  typeA = utils::inumeric(FLERR,arg[3],false,lmp);
  qdist = utils::numeric(FLERR,arg[4],false,lmp);

  lambda = utils::numeric(FLERR,arg[5],false,lmp);

    npowlja =  utils::numeric(FLERR,arg[6],false,lmp);
    npowljb =  utils::numeric(FLERR,arg[7],false,lmp);
    npowc   =  utils::numeric(FLERR,arg[8],false,lmp);



  cut_lj_global = utils::numeric(FLERR,arg[9],false,lmp);
  cut_coul = cut_lj_global;
  if (narg > 10)  cut_coul = utils::numeric(FLERR,arg[10],false,lmp);


  cut_coulsq = cut_coul * cut_coul;

//printf("ALLOCATED %d %d %d \n",npowlja ,npowljb,npowc); 

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

void PairLJCutTIP4PLongCleavTypesSplitLJ::init_style()
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
    printf("\n\n CHECKCOULOMB %f %f %f %f \n\n",cut_coul,qdist,blen, neighbor->skin);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongCleavTypesSplitLJ::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 9)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);

   int switchC = 0, switchLJ = 0;
  
  if (narg > 4 ){
  switchC = utils::inumeric(FLERR,arg[4],false,lmp);  
  switchLJ = utils::inumeric(FLERR,arg[5],false,lmp);    
 	} 
  double l = lambda;
  if (narg > 6) l = utils::numeric(FLERR,arg[6],false,lmp);
 
  
  lambdaC = lambda;
  if (narg > 7)  lambdaC = utils::numeric(FLERR,arg[7],false,lmp);
 
 
  double cut_lj_one = cut_lj_global;
  if (narg == 9) cut_lj_one = utils::numeric(FLERR,arg[8],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
   
    
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_lj[i][j] = cut_lj_one;
      lam[i][j] = l;
      lamC[i][j] = lambdaC;
      powlambdaLJb[i][j] = pow(l,npowljb);
      powDlambdaLJb[i][j] = npowljb*pow(l,npowljb-1.0);
      powlambdaLJa[i][j] = pow(l,npowlja);
      powDlambdaLJa[i][j] = npowlja*pow(l,npowlja-1.0);
      powlambdaC[i][j] = pow(lambdaC,npowc);
      powDlambdaC[i][j] = npowc*pow(lambdaC,npowc-1.0);
      scalingC[i][j] = switchC;
      scalingLJ[i][j] = switchLJ;      
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutTIP4PLongCleavTypesSplitLJ::init_one(int i, int j)
{  double cut = PairLJCutCoulLong::init_one(i,j);


  lam[j][i] = lam[i][j];
  
  powlambdaLJa[j][i] = powlambdaLJa[i][j];
  powDlambdaLJa[j][i] = powDlambdaLJa[i][j];
  powlambdaLJb[j][i] = powlambdaLJb[i][j];
  powDlambdaLJb[j][i] = powDlambdaLJb[i][j];

  lamC[j][i] = lamC[i][j];
  powlambdaC[j][i] = powlambdaC[i][j];
  powDlambdaC[j][i] = powDlambdaC[i][j];


      scalingC[j][i] = scalingC[i][j] ;
      scalingLJ[j][i] = scalingLJ[i][j] ;     

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    SCoffset[i][j] = 4.0 * epsilon[i][j] * (powlambdaLJa[j][i] * pow(ratio, 12.0) - powlambdaLJb[j][i] * pow(ratio, 6.0));
    SCDoffset[i][j] = 4.0 * epsilon[i][j] * (powDlambdaLJa[j][i] * pow(ratio, 12.0) - powDlambdaLJb[j][i] * pow(ratio, 6.0));
    unSCoffset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));    
  } else{
	unSCoffset[i][j] = 0.0;  
    SCoffset[i][j] = 0.0;
    SCDoffset[i][j] = 0.0;
    	}
    	
//    	printf("%d %d %f %f %f %f %f %f %f\n",i,j,powlambdaLJa[i][j],powDlambdaLJa[i][j],powlambdaLJb[i][j],powDlambdaLJb[i][j],unSCoffset[i][j] ,SCoffset[i][j],SCDoffset[i][j] );
    
      SCoffset[j][i] = SCoffset[i][j];
      SCDoffset[j][i] = SCoffset[i][j];
      unSCoffset[j][i] = unSCoffset[i][j];

  // check that LJ epsilon = 0.0 for water H
  // set LJ cutoff to 0.0 for any interaction involving water H
  // so LJ term isn't calculated in compute()

  if ((i == typeH && epsilon[i][i] != 0.0) ||
      (j == typeH && epsilon[j][j] != 0.0))
    error->all(FLERR,"Water H epsilon must be 0.0 for "
               "pair style lj/cut/tip4p/long/cleav");

  if (i == typeH || j == typeH)
    cut_ljsq[j][i] = cut_ljsq[i][j] = 0.0;

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongCleavTypesSplitLJ::write_restart_settings(FILE *fp)
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

void PairLJCutTIP4PLongCleavTypesSplitLJ::read_restart_settings(FILE *fp)
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

void PairLJCutTIP4PLongCleavTypesSplitLJ::compute_newsite(double *xO, double *xH1,
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

void *PairLJCutTIP4PLongCleavTypesSplitLJ::extract(const char *str, int &dim)
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

double PairLJCutTIP4PLongCleavTypesSplitLJ::memory_usage()
{
  double bytes = (double)maxeatom * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  bytes += (double)2 * nmax * sizeof(double);
  return bytes;
}
