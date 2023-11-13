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

The installation can be tested by heading over to `/test`. We provide two different tests: the first one checks that the installation (and therefore the main LAMMPS cleaving commands) works by performing one-time-step simulations, while the second one runs slightly longer simulations to perform "interaction" checks for each cleaving step for the solid-liquid and solid-vacuum Broughton and Gilmer LJ potential.

The tests can be run by executing `./cleav_run.sh` (for the installation check) or `./cleav_int.sh` (for the interaction test):

1. LAMMPS should be patched with the cleaving files and compiled as detailed [above](#compilation).
2. Run the chosen script with a `-l LAMMPS_BINARY` argument, where `LAMMPS_BINARY`, which defaults to `lmp_mpi`, is the path of the LAMMPS executable to be tested.
3. An optional argument `-n` can be given to select the number of the MPI processes to be used, which defaults to 4.
4. Either test should run in less than a minute on 4 cores. 

Please, note that the results obtained in either case are *not* going to be physically correct or meaningful (for which we refer to the [examples](example.md)). For `cleav_run.sh`, the script only tests if the pairs/fix/compute of the cleaving package run without errors. By contrast, `cleav_int.sh` tests if the cleaving code gives consistent results, and should be used to check if any modifications affects the behaviour in an unwanted way. The final output (which is printed at the end of the procedure, but also stored in the `test_iRT1/results.dat` and `test_iRT2/results.dat` files) should be:

    #####################################
    Results from iRT1 (solid-vacuum LJ):
    Wells interactions step1: -0.00965234
    Switch-off interactions work step3: -0.66751
    Wells interactions step4: -2.81862
    #####################################

    #####################################
    Results from iRT2 (solid-liquid LJ):
    Walls interactions step1: -0.00162756
    Walls interactions step2: -0.231522
    Switch interaction step3: -0.859394
    Walls interactions step4: -0.47436
    #####################################
