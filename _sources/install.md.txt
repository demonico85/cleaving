# Installation

## Requirements

The package requires a version of LAMMPS $\geq$ Sep2022, which can be downloaded from [here](https://www.lammps.org/download.html). Older versions are **not** compatible.

## Compilation

1. Compile Lammps with `make` following the [instructions](https://docs.lammps.org/Build_make.html)
2. Copy all files in `/cleaving/src` directory to the lammps `src` folder (*e.g.* `/home/user/lammps-22Dec2022/src`)
3. Run again `make mpi` to patch Lammps with the cleaving package

### Known issues

## Testing
