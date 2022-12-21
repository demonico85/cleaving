# Examples


## Solid-Vacuum interface of a Lennard-Jones crystal with wells

In this example we will set up the cleaving calculation for the calculation of the SFE of a Lennard-Jones crystal in contact with vacuum

The input files for the whole calculations are already given in the directory `/examples/lj_SV` but in this tutorial we will go through the writing of such files from scratch

### Step 1 

1. Copy the file `bulk.in` from the directory `/examples/lj_S`, we will build the cleaving calculations by modifying this file. 

2.  Create the data file with the initial system. In our case we will use a Lennard-Jones crystall in fcc configuration along the direction (111) at the (reduced) temperature of 0.1. Copy the file `fcc111-T01.lmp` from the directory `/examples/lj_systems/`

3. Create the input file for the wells. The exact format is given in the description of the appropriate fix. Copy the file `fcc    fcc111-T01-wells.lmp` from the directory `/examples/lj_systems/`


   E.g.:
```
0.0
0.001
0.002
0.05
0.1
...
0.998
0.999
1.0
```

    Note: 
   * The file must start at 0 and end at 1. 
   * There is no internal control in the code that the boundaries are correct

4. Create a file for the variation of the strength of the wells. This file contains a sequence of increasing consecutive numbers in the interval $[0,1]$ (extremes included). You can copy the file `/examples/lj_SV/step1/lambda_wells.dat` (if you call it differently change its name in point 6)

5. The walls in the system are introduced using the new fix `wellsPforce`: 
    
   ```
   fix f2 all wellPforce ${dw} ${rw} ${expp} ${lambda} file fcc111-T01-wells.lmp 
   ```
    where the explanation of the different parameters is given in the description of the new fix.

6. In the wells version of the cleaving model, the wells are introduced by switching the parameter lambda from 0 (no interactions between atoms and wells) and 1 (full interactions between atoms and wells). Lammps  allows the creation of an input file which can perform several runs in a row by changing the parameter between the different runs. The code to include in the `bulk.in` file is

```
variable lam file lambda_wells.dat
variable i   equal  1
variable a0 equal exp(1/3*ln(4/1.05604))
variable rw    equal sqrt(2)*${a0}/4.0*1.2


label here
variable    lambda equal ${lam}

fix f2 all wellPforce 6.0 ${rw} ${expp} ${lambda} file fcc111-T01-wells.lmp 

print       "Well depth ${lambda}"
run   5000

fix  f5 all ave/time 100 5 500  c_thermo_temp c_thermo_pe f_f2 v_lambda file  out/ave.F.${i}.out

run  15000

unfix f5
unfix f2 

write_data data/Fstep1.${i}.data nocoeff
variable    ii equal ${i}+1
variable    i  equal ${ii}
next lam
jump SELF here
```

   To keep the main directory clean from all the output files generated during the run, we include such files in two subdirs `out` and `data` which must be created before running the simulation, otherwise Lammps will throw an error.
   We refer to the [Lammps documentation](https://docs.lammps.org/jump.html) for the use of the use of the `jump` command to create a loop. Each iteration of the loops produces the following files:
   * `Fstep1.${cnt}.data`: data file containing the last configuration of the i-th iteration
   * `ave.F.${cnt}.out`: File which contains a summary of  the properties of the system, including the work (f_f2) 

### Step 2

In the second step we are switching off the interactions among the two sides of the cleaving wall. 

1. Copy the last `Fstep1.${i}.data` file from the Step1 folder

2. The switching off is implemented directly in the definition of the pair interactions. We therefore need to change the pair interaction to the new defined type

```
pair_style lj/BGcleavwellspbc ${cutoff1} ${cutoff2} z
```

   all the parameters (cutoff-1, cutoff2, epsilon, sigma) remain identical to those used in the Step1. Note that there is a third parameter in this pair_style, the direction normal to the cleaving plane (z in this case). 

4. Add the command for the walls in the lammps script file: `fix f2 all wellPforce 6.0 ${rw} ${expp} 1.0 file fcc111-T01-wells.lmp`. In this step the strenght of the walls remains constant, therefore we replace 1 to the variable ${lambda}

3. Create a file for the switching off of the interactions across teh cleaving plane. This is obtained by specifying in a file the amount of the coordinate (z in our case) the box system will be moved away from its periodic images. The file needs to always start with zero, which is the first point of the switching

   E.g.:
```
0.0
0.001
0.002
0.05
0.1
...
2.0
2.2
2.4
2.8
3.2
...
```

   Note: in this case there is no upper boundary. The limit mnust ensure that there are no more interactions between the box and its periodic images. Usually a value of 2.6 is enough.

4. The actual switching off is obtained through another loop which increases the size of the box. Since the atoms are kept in place by the cleaving potential, increasing the size of the box creates a vacuum space. When the vacuum space is larger than cutoff2 than the box and its image across the cleaving wall do not interact anymore.

```

variable initzhi equal $(zhi)
variable initzlo equal $(zlo)

label here
variable dhi equal ${zc}
variable newzhi equal ${dhi}+${initzhi}
variable newzlo equal ${initzlo}

print       "Forward Interaction ${cnt}: ${zc} "
print       " ${initzhi} ${dhi} zhi "

change_box all z final ${newzlo}  ${newzhi}

print       "B ${initzhi} ${dhi} zhi "

run   ${eqnts}

fix   f6 all ave/time ${Nevery} ${Nrepeat} ${Nfreq}  c_1[*] file dat/inters3.${i}.dat mode vector
fix   f7 all ave/time ${Nevery} ${Nrepeat} ${Nfreq} c_thermo_temp c_thermo_pe v_totW v_lambda f_f2 v_dhi file  out/ave.F.${i}.out

run   ${nts}

unfix f6
unfix f7

write_data  data/Fstep3.${i}.data  nocoeff
variable    ii equal ${i}+1
variable    i  equal ${ii}
next zc
jump SELF here
```

   We refer to the [Lammps documentation](https://docs.lammps.org/change_box.html) for the exact explanation of the `change_box` command


### Step 3

In the last step we are removing the wells leaving free the newly created interfaces. 

1. Copy the last dat file from the Step2 directory

2. Create the new directories `/out` and `/data`

3. Create the input file for the wells. The exact format is given in the description of the appropriate fix. Copy the file `fcc    fcc111-T01-wells.lmp` from the directory `/examples/lj_systems/`


   E.g.:
```
1.0
0.99
0.98
...
0.02
0.01
0.0
```

   Note:
   * The file must start at 1 and end at 0
   * There is no internal control in the code on the boundaries
   * It does not need to be the reverse of the file used in Step1


4. The loop in Step3 in analogous to the loop in Step1 run backwards


```
label here
variable    lambda equal ${lam}

fix f2 all wellPforce ${dw} ${rw} ${expp} ${lambda} file fcc111-T01-wells.lmp
variable totW   equal "f_f2"

print       "Well depth ${lambda}"

run   ${eqnts}

fix  f6 all ave/time  ${Nevery} ${Nrepeat} ${Nfreq}  c_thermo_temp c_thermo_pe v_totW v_lambda f_f2   file out/ave.F.$i.out

run  ${nts}

unfix f2
unfix f6
write_data data/Fstep4.$i.data nocoeff
variable    ii equal $i+1
variable    i  equal ${ii}
next lam
jump SELF here
```



#### Calculation of the SFE 

The SFE is obtained by summing the work performed in the Step1, Step2, Step3

1. The files `.out` generated in Step1 contains the quantity `f_f2` which is the work performed. An average of that quantity for each lambda gives the variation of the energy in Step1. The integration of the results quantity over lambda gives the total work in Step1

2. The files `.dat` generated in Step2 contains the interactions _switched-off_ during the Step2. By averaging these values for each value of zw we obtain the variation of the energy in Step2. The integration of the results over zw gives the total work in Step2


3. Step3 is analogous to Step1



