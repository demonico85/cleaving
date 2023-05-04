#! /bin/bash

function checkcommands () {
    if ! command -v $1 &> /dev/null
      then
        echo "$1 could not be found"
        exit
    fi
}

currdir=$(pwd)
inpdir="/home/mmm1133/Scratch/LJ_sol_liq/utils/"
insim=2  # starting index
totsim=5 # total number of sims
nproc=8
replace=1

step="step4"


let endsim=insim+totsim-1


if [ ! -z $step ];
  then
    echo "WARNING: you are running the simulation starting with $step"
    echo "If that is the intendend behaviour ignore this warning"
    echo "Otherwise kill the script and qdel the jobs"
    echo ""
    replace=0
fi


checkcommands  "mpirun"

if [ ! -d cases ];
  then
    echo "ERROR: No cases dir found"
    echo "Exiting..."
    exit
fi

for or in 100 #100 110 111 
do
   for T in T3  
   do
     for ty in  walls
       do
         if [ $ty == "walls" ]
           then
             prfx="wa"
         else 
             prfx="we"
         fi
         cd $currdir
         fold=$(echo "fcc"$or"_"$T"_"$ty)
         for i in `seq $insim $endsim `
           do
		cd $currdir
		newfold=$(echo $fold"_"$i)
                if [ ! -d $newfold ] && [ ! -z $step ];
                  then
                     echo "$fold not found. Skipping"
                     continue
                fi

                if [ $replace -eq 1 ];
                  then
                   echo "Deleting old $fold v.$i"
                   rm -r $newfold 2> /dev/null
                fi
                if [ ! -d $newfold ];
                  then
                    echo "No $newfold"
                    echo "creating new simulations"
 
                    if [ ! -d ./cases/$fold ];
                       then
                         echo "No case $fold found in the dir"
                         echo "Skipping: $fold"
                         continue
                    fi
		         cp -r ./cases/$fold ./$newfold
                elif [  -d $newfold ] && [ -z "$step" ];
                  then
                   echo "Skipping: $newfold"
                   continue
                fi

                cd $newfold
                qsub -N $( echo $prfx$or$T"_"$i) -v step=$step $inpdir/runcleaving.sh

         done
      done
   done
done

