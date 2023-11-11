# Installation

## Requirements

The package requires a version of LAMMPS $\geq$ Mar2023, which can be downloaded from [here](https://www.lammps.org/download.html). Older versions are **not** compatible.

## Compilation

1. Compile Lammps with `make` following the [instructions](https://docs.lammps.org/Build_make.html)
2. Copy all files in `/cleaving/src` directory to the lammps `src` folder (*e.g.* `/home/user/lammps-28Mar2023/src`)
3. Run again `make mpi` to patch Lammps with the cleaving package

### Known issues

The list of known issues can be browsed online [here](https://github.com/demonico85/cleaving/issues).

## Testing

The installation can be tested by heading over to `/test`, which contains a script that runs a short sample of cleaving calculation.

1. LAMMPS should be patched with the cleaving files and compiled as detailed [above](#compilation).
2. Run the script `./cleav_test.sh -l LAMMPS_BINARY`, where `LAMMPS_BINARY`, which defaults to `lmp_mpi`, is the path of the LAMMPS executable to be tested.
3. An optional argument `-n` can be given to select the number of the MPI processes to be used. The default value is 4.
4. The whole test should take a few minutes on 4 cores. Please, note that the results obtained are *not* going to be physically correct. This script only tests if the pairs/fix/compute of the cleaving package run without errors. In order to obtain physically meaningful results, please refer to the [examples](example.md).

## Testing calculations

The results of the cleaving code can be tested with a short test in the directory `/test`, which allows to run a sample calculation for the case of solid-liquid and solid-vacuum interface with Broughton and Gilmer LJ potential.

1. LAMMPS should be patched with the cleaving files and compiled as detailed [above](#compilation).
2. Run the script `./cleav_regression.sh -l LAMMPS_BINARY`, where `LAMMPS_BINARY`, which defaults to `lmp_mpi`, is the path of the LAMMPS executable to be tested.
3. An optional argument `-n` can be given to select the number of the MPI processes to be used. The default value is 4.
4. The whole test should take a few minutes on 4 cores. The results that should be obtained are
   * Solid-Vacuum  Surface Free Energy is:  5.7299
   * Solid-Liquid Surface Free Energy is:  0.0957595

Please, note that the results obtained are *not* going to be physically correct. This script only tests if the cleaving code gives consistent results and can be used to check if any modifications affects the known behaviour. In order to obtain physically meaningful results, please refer to the [examples](example.md).

