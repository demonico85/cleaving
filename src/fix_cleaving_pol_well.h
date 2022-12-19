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

#ifdef FIX_CLASS

FixStyle(wellPforce,FixForcePolWell)

#else

#ifndef LMP_FIX_CLEAV_PWELL
#define LMP_FIX_CLEAV_PWELL

#include "fix.h"

namespace LAMMPS_NS {

class FixForcePolWell : public Fix {
 public:
  FixForcePolWell(class LAMMPS *, int, char **);
  ~FixForcePolWell();

  int setmask();
  void init();
  void post_force(int );
  void setup(int );
  double compute_scalar ();
  double compute_vector(int);

 protected:

  int nwells, ntypes, eflag;
  double dw,expon,lambda;
  double rw,rwsq,irwsq, prefactor,M;
  double iareaxy,cleaving_wall;
  double cut_global, zw, a0,zwsq,invdelta;
  double *ewell,*ewell_all;
  double **poswell;
  double xhalf, yhalf, zhalf,xedge,yedge,zedge;
  double cleavwork, localwork;

  char *namefile;
  FILE *fpw;
  void ReadFWells();
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
