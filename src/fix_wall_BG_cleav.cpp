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

#include "fix_wall_BG_cleav.h"

#include "atom.h"
#include "error.h"
#include "math_special.h"

#include <cmath>


// Definition coefficients 
// Davidchack, R. L. and Laird, B. B. Journal of Chemical Physics,  118, 7651:7657, (2003)

#define  C1  0.016132
#define  C2  3136.6
#define  C3 -68.069
#define  C4 -0.083312
#define  C5  0.74689

using namespace LAMMPS_NS;
using MathSpecial::powint;

/* ---------------------------------------------------------------------- */

FixWallBGCleav::FixWallBGCleav(LAMMPS *lmp, int narg, char **arg) : FixWall(lmp, narg, arg)
{
  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

void FixWallBGCleav::precompute(int m)
{
  coeff1[m] = 48.0 * epsilon[m] * powint(sigma[m], 12);
  coeff2[m] = 24.0 * epsilon[m] * powint(sigma[m], 6);
  coeff3[m] = 4.0 * epsilon[m] * powint(sigma[m], 12);
  coeff4[m] = 4.0 * epsilon[m] * powint(sigma[m], 6);
  coeff5[m] = 12.0 * C2 * epsilon[m] * powint(sigma[m],12);
  coeff6[m] =  6.0 * C3 * epsilon[m] * powint(sigma[m],6);
  coeff7[m] =  2.0 * C4 * epsilon[m] * 1./powint(sigma[m],2);
  coeff8[m] =  C2 * epsilon[m] * powint(sigma[m],12);
  coeff9[m] =  C3 * epsilon[m] * powint(sigma[m],6);
  coeff10[m] = C4 * epsilon[m] * 1./powint(sigma[m],2);
  c1[m] = C1*epsilon[m];  
  c5[m] = C5*epsilon[m];  
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallBGCleav::wall_particle(int m, int which, double coord)
{
  double delta, rinv, r2inv, r6inv, fwall;
  double vn, rsq;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0)
        delta = x[i][dim] - coord;
      else
        delta = coord - x[i][dim];
      if (delta >= cutoffout[m]) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }
      rsq   = delta * delta;
      rinv  = 1.0 / delta;
      r2inv = rinv * rinv;
      r6inv = r2inv * r2inv * r2inv;
      
      if(delta < cutoff[m]){
          fwall = side * r6inv * (coeff1[m] * r6inv - coeff2[m]) * rinv;
      	  f[i][dim] -= fwall;
          ewall[0] += r6inv * (coeff3[m] * r6inv - coeff4[m]) + c1[m];;
          ewall[m + 1] += fwall;
          }
      else{
          fwall = side * (r6inv * (coeff5[m] * r6inv + coeff6[m]) * r2inv - coeff7[m]);
      	  f[i][dim] -= fwall;
          ewall[0] += r6inv * (coeff8[m] * r6inv + coeff9[m]) + coeff10[m] * rsq + c5[m];          
          ewall[m + 1] += fwall;   
      }


      if (evflag) {
        if (side < 0)
          vn = -fwall * delta;
        else
          vn = fwall * delta;
        v_tally(dim, i, vn);
      }
    }

  if (onflag) error->one(FLERR, "Particle on or inside fix wall surface");
}
