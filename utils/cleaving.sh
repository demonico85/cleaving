#! /bin/bash

currdir=$(pwd)
inpdir=/home/nico/lammps/inputs_8Jul2020
nproc=20

for or in 111 #100 110 #111 
do
   for T in T1 # T2 T3
   do
     for ty in walls # wells
       do
	 cd $currdir
	 fold=$(echo "fcc"$or"_"$T"_"$ty)
         cd $fold
	 $inpdir/runcleaving.sh -n $nproc
      done
   done
done
