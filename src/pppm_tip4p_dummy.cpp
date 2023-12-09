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

}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMTIP4Pdummy::make_rho()
{

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
}

/* ----------------------------------------------------------------------
   find 2 H atoms bonded to O atom i
   compute position xM of fictitious charge site for O atom
   also return local indices iH1,iH2 of H atoms
------------------------------------------------------------------------- */

void PPPMTIP4Pdummy::find_M(int i, int &iH1, int &iH2, double *xM)
{
}
