#! /bin/bash


currdir=$(pwd)
inpdir="/home/mmm1133/interface/inputs_8Jul2020"
insim=1
totsim=1
nproc=20
replace=1


if [ ! -d cases ];
  then
    echo "ERROR: No cases dir found"
    echo "Exiting..."
    exit
fi

for or in 111 # 100 110 
do
   for T in T1 #T2 T3
   do
     for ty in walls # wells
       do
         cd $currdir
         fold=$(echo "fcc"$or"_"$T"_"$ty)
         for i in `seq $insim $totsim `
           do
		cd $currdir
		newfold=$(echo $fold"_"$i)
                if [ $replace -eq 1 ];
                  then
                   echo "Deleting old $fold v $i"
                   rm -r $newfold
                fi
                if [ ! -d $newfold ];
                  then
                    echo "No $newfold"
                    echo "creating new simulations"
 
		    cp -r ./cases/$fold ./$newfold
		    cd $newfold
                    qsub -N $( echo "fcc"$or"_"$T"_"$ty"_"$i) $inpdir/runcleaving.sh
                else
                   echo "Skipping: $newfold"
                fi
         done
      done
   done
done

