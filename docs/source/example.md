# Examples


## Solid-Vacuum interface of a Lennard-Jones crystal with wells

In this example we will set up the cleaving calculation for the calculation of the SFE of a Lennard-Jones crystal in contact with vacuum

The input files for the whole calculations are already given in the directory `/examples/lj_SV` but in this tutorial we will go through the writing of such files from scratch

### Step 1 

1. Copy the file `bulk.in` from the directory `/examples/lj_S`, we will build the cleaving calculations by modifying this file. 

2.  Create the data file with the initial system. In our case we will use a Lennard-Jones crystall in fcc configuration along the direction (111) at the (reduced) temperature of 0.1. Copy the file `fcc111-T01.lmp` from the directory `/examples/lj_systems/`

3. Create the input file for the wells. The exact format is given in the description of the appropriate fix. Copy the file `fcc    fcc111-T01-wells.lmp` from the directory `/examples/lj_systems/`

4. Create a file for the variation of the strength of the wells. This file contains a sequence of increasing consecutive numbers in the interval $[0,1]$ (extremes included). You can copy the file `/examples/lj_SV/step1/lambda_wells.dat` (if you call it differently change its name in point 6)

5. The walls in the system are introduced using the new fix `wellsPforce`: 
    
   ```
   fix f2 all wellPforce ${dw} ${rw} ${expp} ${lambda} file fcc111-T01-wells.lmp 
   ```
where the explanation of the different parameters is given in the description of the new fix.

6. In the wells version of the cleaving model, the wells are introduced by switching the parameter lambda from 0 (no interactions between atoms and wells) and 1 (full interactions between atoms and wells). Lammps  allows the creation of an input file which can perform several runs in a row by changing the parameter between the different runs. The code to include in the `bulk.in` file is

```
variable lam file lambda_wells.dat
variable  cnt   equal  1
variable a0 equal exp(1/3*ln(4/1.05604))
variable  rw    equal sqrt(2)*${a0}/4.0*1.2


label here
variable    lambda equal ${lam}

fix f2 all wellPforce 6.0 ${rw} ${expp} ${lambda} file fcc111-T01-wells.lmp 

print       "Well depth ${lambda}"
run   5000

fix  f5 all ave/time 100 5 500  c_thermo_temp c_thermo_pe f_f2 v_lambda file  out/ave.F.${cnt}.out

run  15000

unfix f5
unfix f2 

write_data data/Fstep1.${cnt}.data nocoeff
variable    cntt equal ${cnt}+1
variable    cnt  equal ${cntt}
next lam
jump SELF here
```