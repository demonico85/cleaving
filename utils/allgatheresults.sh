#! /bin/bash


currdir=$(pwd)
inpdir=./utils
nbins=5
replace=0
in=1
fin=10

if [ -d results ] && [ $replace -eq 1 ];
then
   rm -r results
fi

if [ ! -d results ];
  then
    mkdir results
fi

for i in `seq $in $fin`
do
for or in 111 110 100  
do
   for T in T1 T2 T3
   do
     for ty in walls #wells
       do
         cd $currdir
         fold=$(echo "fcc"$or"_"$T"_"$ty"_"$i)
         echo "Trying dir: $fold"
         if [ ! -d $fold ];
           then
             echo "dir not found. Skipping..."
             continue
         fi
         cd $fold
         for st in step1 step2 step4 
           do
	     echo $st
             cd $currdir/$fold
             cd $st/out
             $inpdir/work.sh F $nbins > work.log
	     if [ ! -f F-work.dat ];
	       then
		       echo "Results non produced in $fold for $st"
		       echo "Chek it better"
	       else
		       cp F-work.dat $currdir/results/$(echo $i"_fcc"$or"_"$T"_"$ty"_"$st".dat")
	       fi
	   done
	   cd $currdir/$fold
	   st=step3
	   echo $st
             cd $st/dat
	     echo $inpdir/compiled/s3work
	     if [ ! -x $inpdir/compiled/s3work ];
               then
                echo "ERROR: s3work not found or not executable"
		echo "Exiting..."
		exit -1
             fi
             $inpdir/compiled/s3work
             $inpdir/work.sh F $nbins > work.log
             if [ ! -f F-work.dat ];
               then
                       echo "Results non produced in $fold for $st"
                       echo "Chek it better"
               else
                       cp F-work.dat $currdir/results/$(echo $i"_fcc"$or"_"$T"_"$ty"_"$st".dat")
               fi	   
      done
   done
done
done
