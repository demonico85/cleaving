#! /bin/bash


insim=1
totsim=5
currdir=$(pwd)
inpdir="/home/mmm1133/Scratch/LJ_sol_liq/utils"
nbins=5
replace=0

#extp=("walls" "wells")
extp=("walls")

#temp=("T1" "T2" "T3")
temp=("T3")

#orientation=("111" "100" "110")
orientation=("100")

if [ -d results ] && [ $replace -eq 1 ];
then
   rm -r results
fi

if [ ! -d results ];
  then
    mkdir results
fi

for or in "${orientation[@]}" 
do
   for T in "${temp[@]}" 
   do
     for ty in "${extp[@]}"
       do
         for i in `seq $insim $totsim `
           do
             fold=$(echo "fcc"$or"_"$T"_"$ty"_"$i)
             if [ ! -d $fold ];
                then
                  continue
             fi
             $inpdir/singleres.sh $or $T $ty $i
         done	   
      done
   done
done

