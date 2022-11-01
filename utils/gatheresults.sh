#! /bin/bash


currdir=$(pwd)
inpdir=/home/nico/lammps/inputs_8Jul2020
nbins=5
replace=0


if [ -d results ] && [ $replace -eq 1 ];
then
   rm -r results
fi

if [ ! -d results ];
  then
    mkdir results
fi

for or in 110 #111 #110 #100  
do
   for T in T1 # T2 T3
   do
     for ty in walls # wells
       do
         cd $currdir
         fold=$(echo "fcc"$or"_"$T"_"$ty)
         cd $fold
         for st in step4 #step1 step2 step4 
           do
             cd $currdir/$fold
             cd $st/out
             $inpdir/work.sh F $nbins
	     if [ ! -f F-work.dat ];
	       then
		       echo "Results non produced in $fold for $st"
		       echo "Chek it better"
	       else
		       cp F-work.dat $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st".dat")
	       fi
	   done
	   cd $currdir/$fold
	   st=step3
             cd $st/dat
	     $inpdir/s3work
	     echo $inpdir/s3work
	     if [ ! -x $inpdir/s3work ];
               then
                echo "ERROR: s3work not found or not executable"
		echo "Exiting..."
		exit -1
             fi
             $inpdir/work.sh F $nbins
             if [ ! -f F-work.dat ];
               then
                       echo "Results non produced in $fold for $st"
                       echo "Chek it better"
               else
                       cp F-work.dat $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st".dat")
               fi	   
      done
   done
done

