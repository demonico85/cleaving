#! /bin/bash

step="step1"

cd ./area_A0
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A0 ~/bin/inputs/runcleavcry.sh -s $step

cd ../area_A1
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A1 ~/bin/inputs/runcleavcry.sh -s $step



cd ../area_A2
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A2 ~/bin/inputs/runcleavcry.sh -s $step



cd ../area_A3
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A3 ~/bin/inputs/runcleavcry.sh -s $step

cd ../area_A4
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A4 ~/bin/inputs/runcleavcry.sh -s $step

cd ../area_A5
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A5 ~/bin/inputs/runcleavcry.sh -s $step

cd ../area_A6
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A6 ~/bin/inputs/runcleavcry.sh -s $step

cd ../area_A7
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A7 ~/bin/inputs/runcleavcry.sh -s $step

cd ../area_A8
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A8 ~/bin/inputs/runcleavcry.sh -s $step


cd ../area_A9
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A9 ~/bin/inputs/runcleavcry.sh -s $step


cd ../area_A10
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A10 ~/bin/inputs/runcleavcry.sh -s $step



cd ../area_A11
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A11 ~/bin/inputs/runcleavcry.sh -s $step


cd ../area_A12
echo $(pwd)
rm slurm-* 2> /dev/null

sbatch -J A12 ~/bin/inputs/runcleavcry.sh -s $step







