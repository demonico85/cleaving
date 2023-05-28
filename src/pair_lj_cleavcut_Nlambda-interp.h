/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------


This subroutine is written from the original pair_lj_cut in LAMMPS adapted for
cleaving calculations

   Contributing author: Di Pasquale Nicodemo
   University of Leicester, March 2020
   email: nicodemo.dipasquale@gmail.com    
   
   The documentation for this pair potential can be browsed at the following link:
   https://demonico85.github.io/cleaving/
    

------------------------------------------------------------------------- */


#ifdef PAIR_CLASS

PairStyle(lj/cleavcutNlint,PairLJCleavCutNLInt)

#else

#ifndef LMP_PAIR_LJ_CLVCUT_NL_INT_H
#define LMP_PAIR_LJ_CLVCUT_NL_INT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCleavCutNLInt : public Pair {
 public:
  PairLJCleavCutNLInt(class LAMMPS *);
  virtual ~PairLJCleavCutNLInt();
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
  void *extract(const char *, int &);

 protected:
  int natoms,pallocation,index,ind_dir,npow;
  int *gbox;
  double rspace;
//  double u0,u1,u2;
  double cut_global,lambda,cleavwall;
  double xprd,yprd,zprd,xy,yz,xz;
  double **cut;
  double **epsilon,**sigma, **lam;
  double **lj1,**lj2,**lj3,**lj4,**offset,**radm;
  double **powlambda, **powDlambda;
  double **lju0, **lju1, **lju2, **lju3, **lju4;  
  char *idflag;


  virtual void allocate();
  void global_boundary();
  void calculate_int_coefficients(int, int);
  int find_scaling(int,int,int ,int,double *);

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
