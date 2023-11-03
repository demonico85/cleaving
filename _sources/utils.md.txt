# Utils

Here, we report a list of the subroutines included in the directory `utils/` along with a brief comment on what is their purpose.

````{note}
These programs are useful for cleaving pre/post processing calculations (e.g., extracting and averaging useful quantities). These small programs are given without warranty
````

* `TzerogatherEn.sh`: Extract the value of $\gamma$ as difference in energy for a zero temperature LJ crystal at different strains
* `Tzerogathermin.sh`: Extract the value of $\gamma$ from step3 for a zero temperature LJ crystal at different strains
* `allgatheresults.sh`: Collect work in all the four steps, for different run of cleaving, at different temperature, with different crystal orientation and different cleaving potential in a LJ crystal
* `athenawork.sh`: Collect work in all the four steps for a LJ crystal, version to be used in a queue system to run the postprocessing as a batch job
* `blockAv.py`: Script to calculate block average and some statistics of a list of measurements
* `c3cryst.f90`: Rewrite the file inters3.* obtained in step3 in a format similar to step1/2/4 to be read using the same subroutines used for these steps. Valid for zero temperature LJ crystals
* `c3cryst_wells.f90`: Same as `c3cryst.f90`, to be used with wells cleaving potential
* `calcintegrals.f90`: Trapezoidal rule to calculate work integrals 
* `calcworkstep3.f90`: Same as `c3cryst.f90`, valid for solid-liquid interfaces in a LJ crystal
* `calcworkstep3_nodupl.f90`: Same as `calcworkstep3.f90`, do not consider the duplicate atoms in step3
* `check_sim.sh`: Check if the simulation in a queue system are still running or if they have stopped before the end of the calculation completed
* `cleaving.sh`: Automated script to run the entire cleaving model for a solid-liquid interface. It needs `preparation.sh`
* `compile.sh`: Compile the fortran subroutines in the `utils/` directory
* `dump2.data.sh`: Convert dump configuration to data file. It needs `dump2data.f90`
* `dump2.data.sh`: Convert dump configuration to data file. It is called by `dump2data.f90`
* `fixedlayer.f90`:  Create a block of FCC crystal with certain orientation
* `gatherwork.sh`: Check if cleaving calculation is running, if not, it collects the work file at each step
* `multiple_sim.sh`: Submit cleaving calculations in different condition to a queue system
* `new_step3IN.f90`: Create input file for step3 in a LJ crystal
* `preparation.sh`: Script called by `cleaving.sh` to initialize the different directories used for cleaving
* `pressureprofile.sh`: Script to calculate the ensemble average of the pressure (stress) profile in a box
* `remove.sh`: Clean directory
* `restartmann.sh`: restart cleaving calculation on a queue system (written for mannitol system)
* `runcleav.sh`: Submit a new cleaving calculation or restart it to a queue system
* `runcleavcry.sh`: Same as `cleaving.sh` but for a solid-vacuum LJ crystal at finite temperature
* `runzerocleavcry.sh`: Same as `runcleavcry.sh` but for a system at zero temperature
* `s3zeroinpljcry.sh`: Set up step3 for solid-vacuum interface for a LJ crystal at zero temperature
* `s4inpljcry.sh`: Set up step4 for solid-vacuum interface for a LJ crystal 
* `selwall.sh`: Calculate work separately for the two pairs of walls in step4
* `singleres.sh`: Calculate work in all the steps for a single configuration of the inputs (e.g., orientation, temperature) for walls cleaving potential
* `singleres_wells.sh`: Same as `singleres.sh` but for wells cleaving potential
* `step3IN.f90`: Create input file for step3 in a LJ crystal
* `step4IN.f90`: Create input file for step4 in a LJ crystal
* `submit.sh`: Submit a range of cleaving calculation as batch jobs
* `work.sh`: Calculate block average of the work at each step of cleaving for a LJ system
* `work_new_s4.sh`: Same as `work.sh` but modified to give the work for each pairs of walls separately
* `work_s3.sh`: Same as work, modified with to use a different convention for the fixes name identifying the cleaving work
* `workmannitol.sh`: Same as `work.sh` but modified for mannitol system