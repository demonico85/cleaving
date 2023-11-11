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

The installation can be tested by heading over to `/test`. We provide two different tests: the first one checks that the installation (and therefore the main LAMMPS cleaving commands) works by running very short simulations, while the second one runs a full cleaving calculation yielding the solid-liquid and solid-vacuum surface free energy of the Broughton and Gilmer LJ potential.

The tests can be run by executing `./cleav_run.sh` (for the installation check) or `./cleav_full.sh` (for the scientific test) following these instructions:

1. LAMMPS should be patched with the cleaving files and compiled as detailed [above](#compilation).
2. Run the chosen script with a `-l LAMMPS_BINARY` argument, where `LAMMPS_BINARY`, which defaults to `lmp_mpi`, is the path of the LAMMPS executable to be tested.
3. An optional argument `-n` can be given to select the number of the MPI processes to be used, which defaults to 4.
4. Both tests should take a few minutes on 4 cores. 

Please, note that the results obtained in either case are *not* going to be physically correct or meaningful (for which we refer to the [examples](example.md)). For `cleav_run.sh`, the script only tests if the pairs/fix/compute of the cleaving package run without errors. By contrast, `cleav_full.sh` tests if the cleaving code gives consistent results, and can be used to check if any modifications affects the behaviour in an unwanted way. The results that should be obtained are
* Solid-Vacuum Surface Free Energy: 5.7299
* Solid-Liquid Surface Free Energy: 0.0957595
