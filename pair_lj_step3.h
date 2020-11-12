/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/cleavs3,PairLJStep3)

#else

#ifndef LMP_PAIR_CLEAV_STEP3
#define LMP_PAIR_CLEAV_STEP3

#include "pair.h"

namespace LAMMPS_NS {

class PairLJStep3 : public Pair {
 public:
  PairLJStep3(class LAMMPS *);
  virtual ~PairLJStep3();
  virtual void compute(int, int);

  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 protected:

  int pallocation;
  double cut_global_out, cut_global_in, lambda;
  double delta, Dfac;
  double **cut,**cut_in, **cutsq_in;
  double **epsilon,**sigma, **lamcoeff, **Dlambda, **Dfactorlam;
  double **lj1,**lj2,**lj3,**lj4, **lj5, **lj6, **lj7, **lj8, **lj9, **lj10,**c1,**c5;
//  void set_combinations();
  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
