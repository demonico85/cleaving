/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.


This compute calculates the distance (Delta x, Delta y, Delta z) crossed by atoms in the specified group
  at each  time step
   Author: Di Pasquale Nicodemo
   University of Leicester, March 2020
   email: nicodemo.dipasquale@gmail.com    
   
   The documentation for this pair potential can be browsed at the following link:
   https://demonico85.github.io/cleaving/

 ----------------------------------------------------------------------*/


#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(displace/atom_cleav,ComputeDisplaceAtomCleav);
// clang-format on
#else

#ifndef LMP_COMPUTE_DISPLACE_ATOM_CLEAV_H
#define LMP_COMPUTE_DISPLACE_ATOM_CLEAV_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDisplaceAtomCleav : public Compute {
 public:
  ComputeDisplaceAtomCleav(class LAMMPS *, int, char **);
  ~ComputeDisplaceAtomCleav() override;
  void init() override;
  void compute_peratom() override;
  void set_arrays(int) override;
  void refresh() override;
  double memory_usage() override;

 private:
  int nmax;
  double **displace;
  char *id_fix;
  class FixStoreAtom *fix;

  int refreshflag, ivar, nvmax;    // refresh option is enabled
  char *rvar;                      // for incremental dumps
  double *varatom;
};

}    // namespace LAMMPS_NS

#endif
#endif
