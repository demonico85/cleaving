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

FixStyle(move/dupl,FixMoveDupl)

#else

#ifndef LMP_FIX_MOVE_DUPL
#define LMP_FIX_MOVE_DUPL

#include "fix.h"


namespace LAMMPS_NS {

class FixMoveDupl : public Fix {
 public:
  FixMoveDupl(class LAMMPS *, int, char **);
  ~FixMoveDupl();

  int setmask();
  void post_integrate();
  void init();

  void set_arrays(int);
  double memory_usage();
  void reset_dt();

 protected:

  int nvalues, nrestart,buff_size,nmax;
  int  icompute,not_found_value;
  int time_origin;
  int *which,*argindex,*value2index;
  double **displace, **locdisplace;
  char **ids;  

//  FILE *fpl;

 void allocate();
 void collect_displacement();
};

}

#endif
#endif
