# Installation

## Requirements

The package requires a version of LAMMPS $\geq$ Mar2023, which can be downloaded from [here](https://www.lammps.org/download.html). Older versions are **not** compatible.

## Compilation

1. Compile Lammps with `make` following the [instructions](https://docs.lammps.org/Build_make.html)
2. Copy all files in `/cleaving/src` directory to the lammps `src` folder (*e.g.* `/home/user/lammps-28Mar2023/src`)
3. Run again `make mpi` to patch Lammps with the cleaving package

### Known issues

## Testing

In order to test the installation, the folder `./test` contains a script that will run a sample of cleaving calculation. 

1. After compiling LAMMPS with the cleaving sources, enter the ddirectory `./test`
2. Run the script `./cleav_test.sh` with the option  `-l` which contains the name with the full path of the LAMMPS executable to be tested
3. An optional argument `-n` can be given to select the number of the MPI processes to be used, if not provided, the default value is 4
4. The whole test should take around 2 minutes with 4 cores. Please, note that the result obtained are *not* going to be physically correct. This script only tests if the pairs/fix/compute of the cleaving package are running without errors. In order to obtain physically meaningful results, please refers to the examples provided in [lj_SV](./example_SV_wwlls.md) and [lj_SL](./example_SL_walls.md).

