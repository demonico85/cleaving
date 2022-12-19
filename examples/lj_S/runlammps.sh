#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --job-name=ljbulk
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --account=a12-vfl
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ndp8@le.ac.uk
#SBATCH --export=PATH

module purge

module load gcc/6.3.0/1 openmpi/3.0.1/01 

lmp="./lmp_mpi"


    rm -r ./dump/ 2> /dev/null
    rm -r ./restart/ 2> /dev/null

    mkdir dump
    mkdir restart


mpirun -np 28 $lmp  < bulk.in 

exit $?
