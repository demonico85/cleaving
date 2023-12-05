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
------------------------------------------------------------------------- */

#include "pppm_tip4p_dummy.h"

#include "atom.h"
#include "domain.h"
#include "force.h"
#include "error.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

#define OFFSET 16384

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMTIP4Pdummy::PPPMTIP4Pdummy(LAMMPS *lmp) : PPPM(lmp)
{
  triclinic_support = 1;
  tip4pflag = 1;
}

/* ---------------------------------------------------------------------- */

void PPPMTIP4Pdummy::init()
{
  // TIP4P PPPM requires newton on, b/c it computes forces on ghost atoms

  if (force->newton == 0)
    error->all(FLERR,"Kspace style pppm/tip4p requires newton on");

  PPPM::init();
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

void PPPMTIP4Pdummy::particle_map()
{
  int nx,ny,nz,iH1,iH2;
  double *xi,xM[3];

  int *type = atom->type;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);
      xi = xM;
    } else xi = x[i];

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((xi[0]-boxlo[0])*delxinv+shift) - OFFSET;
    ny = static_cast<int> ((xi[1]-boxlo[1])*delyinv+shift) - OFFSET;
    nz = static_cast<int> ((xi[2]-boxlo[2])*delzinv+shift) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
        nz+nlower < nzlo_out || nz+nupper > nzhi_out) flag++;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Out of range atoms - cannot compute PPPM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMTIP4Pdummy::make_rho()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz,iH1,iH2;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  double *xi,xM[3];

  // clear 3d density array

  FFT_SCALAR *vec = &density_brick[nzlo_out][nylo_out][nxlo_out];
  for (i = 0; i < ngrid; i++) vec[i] = ZEROF;

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  int *type = atom->type;
  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);
      xi = xM;
    } else xi = x[i];

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (xi[0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (xi[1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (xi[2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          density_brick[mz][my][mx] += x0*rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMTIP4Pdummy::fieldforce_ik()
{
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void PPPMTIP4Pdummy::fieldforce_ad()
{
 
}


/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
------------------------------------------------------------------------- */

void PPPMTIP4Pdummy::fieldforce_peratom()
{
}


void PPPMTIP4Pdummy::poisson()
{
}


void PPPMTIP4Pdummy::poisson_ad()
{
}

void PPPMTIP4Pdummy::poisson_ik()
{
}



/* ----------------------------------------------------------------------
  Fix handling of TIP4P dipole compared to PPPMDisp::slabcorr
------------------------------------------------------------------------- */

#define SMALL 0.00001

void PPPMTIP4Pdummy::slabcorr()
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double zprd_slab = domain->zprd*slab_volfactor;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double *xi, xM[3]; int iH1, iH2;  //for TIP4P virtual site

  // sum local contributions to get global dipole moment
  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);
      xi = xM;
    } else xi = x[i];
    dipole += q[i]*xi[2];
  }

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
    for (int i = 0; i < nlocal; i++)
      dipole_r2 += q[i]*x[i][2]*x[i][2];

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all -
    qsum*dipole_r2 - qsum*qsum*zprd_slab*zprd_slab/12.0)/volume;
  const double qscale = force->qqrd2e * scale;

  if (eflag_global) energy_1 += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd_slab*zprd_slab/12.0);
  }

  // add on force corrections

  double ffact = qscale * (-4.0*MY_PI/volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    double fzi_corr = ffact * q[i]*(dipole_all - qsum*x[i][2]);
    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);
      f[i][2] += fzi_corr*(1 - alpha);
      f[iH1][2] += 0.5*alpha*fzi_corr;
      f[iH2][2] += 0.5*alpha*fzi_corr;
    }
    else f[i][2] += fzi_corr;
  }
}

/* ----------------------------------------------------------------------
   find 2 H atoms bonded to O atom i
   compute position xM of fictitious charge site for O atom
   also return local indices iH1,iH2 of H atoms
------------------------------------------------------------------------- */

void PPPMTIP4Pdummy::find_M(int i, int &iH1, int &iH2, double *xM)
{
  double **x = atom->x;

  iH1 = atom->map(atom->tag[i] + 1);
  iH2 = atom->map(atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1) error->one(FLERR,"TIP4P hydrogen is missing");
  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one(FLERR,"TIP4P hydrogen has incorrect atom type");

  if (triclinic) {

    // need to use custom code to find the closest image for triclinic,
    // since local atoms are in lambda coordinates, but ghosts are not.

    int *sametag = atom->sametag;
    double xo[3],xh1[3],xh2[3],xm[3];
    const int nlocal = atom->nlocal;

    for (int ii = 0; ii < 3; ++ii) {
      xo[ii] = x[i][ii];
      xh1[ii] = x[iH1][ii];
      xh2[ii] = x[iH2][ii];
    }

    if (i < nlocal) domain->lamda2x(x[i],xo);
    if (iH1 < nlocal) domain->lamda2x(x[iH1],xh1);
    if (iH2 < nlocal) domain->lamda2x(x[iH2],xh2);

    double delx = xo[0] - xh1[0];
    double dely = xo[1] - xh1[1];
    double delz = xo[2] - xh1[2];
    double rsqmin = delx*delx + dely*dely + delz*delz;
    double rsq;
    int closest = iH1;

    // no need to run lamda2x() here -> ghost atoms

    while (sametag[iH1] >= 0) {
      iH1 = sametag[iH1];
      delx = xo[0] - x[iH1][0];
      dely = xo[1] - x[iH1][1];
      delz = xo[2] - x[iH1][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rsqmin) {
        rsqmin = rsq;
        closest = iH1;
        xh1[0] = x[iH1][0];
        xh1[1] = x[iH1][1];
        xh1[2] = x[iH1][2];
      }
    }
    iH1 = closest;

    closest = iH2;
    delx = xo[0] - xh2[0];
    dely = xo[1] - xh2[1];
    delz = xo[2] - xh2[2];
    rsqmin = delx*delx + dely*dely + delz*delz;

    while (sametag[iH2] >= 0) {
      iH2 = sametag[iH2];
      delx = xo[0] - x[iH2][0];
      dely = xo[1] - x[iH2][1];
      delz = xo[2] - x[iH2][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rsqmin) {
        rsqmin = rsq;
        closest = iH2;
        xh2[0] = x[iH2][0];
        xh2[1] = x[iH2][1];
        xh2[2] = x[iH2][2];
      }
    }
    iH2 = closest;

    // finally compute M in real coordinates ...

    double delx1 = xh1[0] - xo[0];
    double dely1 = xh1[1] - xo[1];
    double delz1 = xh1[2] - xo[2];

    double delx2 = xh2[0] - xo[0];
    double dely2 = xh2[1] - xo[1];
    double delz2 = xh2[2] - xo[2];

    xm[0] = xo[0] + alpha * 0.5 * (delx1 + delx2);
    xm[1] = xo[1] + alpha * 0.5 * (dely1 + dely2);
    xm[2] = xo[2] + alpha * 0.5 * (delz1 + delz2);

    // ... and convert M to lamda space for PPPM

    domain->x2lamda(xm,xM);

  } else {

    // set iH1,iH2 to index of closest image to O

    iH1 = domain->closest_image(i,iH1);
    iH2 = domain->closest_image(i,iH2);

    double delx1 = x[iH1][0] - x[i][0];
    double dely1 = x[iH1][1] - x[i][1];
    double delz1 = x[iH1][2] - x[i][2];

    double delx2 = x[iH2][0] - x[i][0];
    double dely2 = x[iH2][1] - x[i][1];
    double delz2 = x[iH2][2] - x[i][2];

    xM[0] = x[i][0] + alpha * 0.5 * (delx1 + delx2);
    xM[1] = x[i][1] + alpha * 0.5 * (dely1 + dely2);
    xM[2] = x[i][2] + alpha * 0.5 * (delz1 + delz2);
  }
}
