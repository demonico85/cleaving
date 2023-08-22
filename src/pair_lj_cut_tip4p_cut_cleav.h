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
PairStyle(lj/cut/tip4p/cut/cleav,PairLJCutTIP4PCutCleav);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_CUT_H_CLEAV
#define LMP_PAIR_LJ_CUT_TIP4P_CUT_H_CLEAV

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PCutCleav : public Pair {
 public:
  PairLJCutTIP4PCutCleav(class LAMMPS *);
  ~PairLJCutTIP4PCutCleav() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void *extract(const char *, int &) override;
  double memory_usage() override;

 protected:
//FILE *fp;
  int    natoms,pallocation,index,ind_dir,npow,ntypes;
  int *gbox,*giflag;
    int nchunk;
  double xprd,yprd,zprd,xy,yz,xz, lambda,lambdaC;
  double cut_lj_global, cut_coul_global,lboxhalf;
  double cut_coul, cut_coulsq;
  double cut_coulsqplus;    // extended value for cut_coulsq
  double **cut_lj, **cut_ljsq;
  double **epsilon, **sigma;
  double **powlambda, **powDlambda,**lam;
  double **powlambdaC, **powDlambdaC,**lamC;
  double **lj1, **lj2, **lj3, **lj4, **offset;

  char *idchunk, *idcom;
  class ComputeChunkAtom *cchunk;
  class ComputeCOMChunk *ccom;

  int typeH, typeO;    // atom types of TIP4P water H and O atoms
  int typeA, typeB;    // angle and bond types of TIP4P water
  double alpha;        // geometric constraint parameter for TIP4P
  double qdist;

  int nmax;            // info on off-oxygen charge sites
  int **hneigh;        // 0,1 = indices of 2 H associated with O
                       // 2 = 0 if site loc not yet computed, 1 if yes
  double **newsite;    // locations of charge sites

  void allocate();
  void compute_newsite(double *, double *, double *, double *);
  void global_boundary();
  int find_scaling(int,int,int,int,double *);
  
};
}    // namespace LAMMPS_NS

#endif
#endif
