#!/bin/bash -l

#$ -P Gold
#$ -A Brunel_allocation

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=15:00:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=1G

# Set the working directory to somewhere in your scratch space.  
#$ -cwd 

# Set the name of the job.
#$ -N ljcrys 

# Request 40 cores.
#$ -pe mpi 40 
export $OMP_NUM_THREADS=1


module purge
module load gcc-libs
module load compilers/intel/2018
module load mpi/intel/2018
module load gerun



# Parameters to be changed

dirscripts="/home/mmm1133/Scratch/LJ_sol_liq/utils/"
lmp="lmp_mpi_03Nov22"
lmpath="/home/mmm1133/Scratch/LJ_sol_liq/"
lmp="$lmpath/$lmp"

trap "exit 15" TERM
export TOP_PID=$$

#step="step4"

if [ ! -z "$step" ];
  then
    echo "Restarting from $step"
    $dirscripts/cleaving.sh -n $NSLOTS -l $lmp -s $step
else
    echo "Starting from step1"
    $dirscripts/cleaving.sh -n $NSLOTS -l $lmp
fi
