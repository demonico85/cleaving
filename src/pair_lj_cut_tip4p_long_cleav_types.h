/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/tip4p/long/cleav/types,PairLJCutTIP4PLongCleavTypes);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_CLEAV_TYPES_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_CLEAV_TYPES_H

#include "pair_lj_cut_coul_long.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PLongCleavTypes : public PairLJCutCoulLong {
 public:
  PairLJCutTIP4PLongCleavTypes(class LAMMPS *);
  ~PairLJCutTIP4PLongCleavTypes() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void init_style() override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart_settings(FILE *fp) override;
  void read_restart_settings(FILE *fp) override;
  void *extract(const char *, int &) override;
  double memory_usage() override;

 protected:
 
 FILE *fp;
 
  int natoms,pallocation,ind_dir,npow,ntypes;
  int dubtypes;
  int **scalingC, **scalingLJ;
  double xprd,yprd,zprd,xy,yz,xz,lboxhalf, lambda,lambdaC;
  double **powlambda, **powDlambda,**lam;
  double **powlambdaC, **powDlambdaC,**lamC;
 
 
  int typeH, typeO;    // atom types of TIP4P water H and O atoms
  int typeA, typeB;    // angle and bond types of TIP4P water
  double alpha;        // geometric constraint parameter for TIP4P

  int nmax;            // info on off-oxygen charge sites
  int **hneigh;        // 0,1 = indices of 2 H associated with O
                       // 2 = 0 if site loc not yet computed, 1 if yes
  double **newsite;    // locations of charge sites

  void compute_newsite(double *, double *, double *, double *);
  
  
  void allocate() override;  
  
};

}    // namespace LAMMPS_NS

#endif
#endif
