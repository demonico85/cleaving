/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ORNL)
   References: Fennell and Gezelter, JCP 124, 234104 (2006)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------


This subroutine is written from the original pair_lj_cut in LAMMPS adapted for
cleaving calculations

   Contributing author: Di Pasquale Nicodemo
   University of Leicester, March 2020
   email: nicodemo.dipasquale@gmail.com    
   
   The documentation for this pair potential can be browsed at the following link:
   https://demonico85.github.io/cleaving/
    
% Known issues:
   Force single not yet changed to include lambda

------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(coul/dsfNl,PairCoulNLDSF)


#else

#ifndef LMP_PAIR_COUL_NL_DSF_H
#define LMP_PAIR_COUL_NL_DSF_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulNLDSF : public Pair {
 public:
  PairCoulNLDSF(class LAMMPS *);
  ~PairCoulNLDSF();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:

  int *gbox;
  int ind_dir,index,natoms,npow;
  double **lam;
  double xprd,yprd,zprd,xy,yz,xz;
  double cut_coul,cut_coulsq;
  double alpha,lambda;
  double f_shift,e_shift;
  double **powlambda, **powDlambda;

  char *idflag;

  int find_scaling(int,int,int ,int,double *);

  void global_boundary();


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

E: Pair style coul/dsf requires atom attribute q

The atom style defined does not have this attribute.

*/
