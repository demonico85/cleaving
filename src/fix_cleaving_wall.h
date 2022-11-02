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

FixStyle(wallforce,FixForceLJWall)

#else

#ifndef LMP_FIX_CLEAV_WALL
#define LMP_FIX_CLEAV_WALL

#include "fix.h"

namespace LAMMPS_NS {

class FixForceLJWall : public Fix {
 public:
  FixForceLJWall(class LAMMPS *, int, char **);
  ~FixForceLJWall();

  int setmask();
  void init();
  void post_force(int );
  void setup(int);
  double compute_scalar ();
  double compute_vector(int);
  double memory_usage();

 protected:

  int natomsinwalls[2], ntypes, eflag;
  double sigma, epsilon, delta,rw,iareaxy,cleavpos[2];
  double cut_global, zw, a0,zwsq,invdelta,rwsq,cleav_wall;
  double i2delta, i4delta, pconst;
  double **lposwall, **hposwall;
  double *ewall,*ewall_all;
  double xprd,yprd,zprd;
  double lj1,lj2,lj3,lj4; 
  double cleavwork, localwork;
  char *namefile;

  FILE *fpw;
  int WallP(double, double);
  void ReadFWalls();
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
