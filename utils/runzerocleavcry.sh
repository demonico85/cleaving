#! /bin/bash
#SBATCH --time=100:00:00
#SBATCH --job-name=zw06T3
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --account=a12-vfl
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ndp8@le.ac.uk
#SBATCH --export=PATH

# ---------------------- SOME FUNCTION DEFINITION ------------------------------

function preparedir {


    rm $1/ave.F.* 2> /dev/null
    rm $1/ave.B.* 2> /dev/null
    rm $1/force.dump 2> /dev/null
    rm -r $1/dump 2> /dev/null
    mkdir $1/dump

    rm -r $1/group 2> /dev/null
    mkdir $1/group

    rm -r $1/restart 2> /dev/null
    mkdir $1/restart


    rm -r $1/data 2> /dev/null
    mkdir $1/data

    rm -r $1/out 2> /dev/null
    mkdir $1/out

    rm $1/forces.log 2> /dev/null
}

####################################################################
####################################################################

#      -------- Load some modules for the simulations ------------
#        module purge
#        module load gcc/6.3.0/1 openmpi/3.0.1/01


#   PBS_O_WORKDIR=$SLURM_SUBMIT_DIR
#   PBS_NUM_PPN=$SLURM_NTASKS

PBS_NUM_PPN=28
PBS_O_WORKDIR=$(pwd)


lmp="lmp_mpi_06Feb"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $PBS_O_WORKDIR
echo "Entering woring dir..."
echo $(pwd)
echo

lmpath=$(echo "/home/n/ndp8/shuttleworth/Tzero/$lmp")
if [ ! -e  $lmpath ];
  then
    echo "ERROR: Lammps executable not found in $lmpath"
    exit
fi

temp="T0"
orient="no"
startstep="step3"
runsim=1
now=$(date +"%m_%d_%Y_%s")
log=$(echo $now.log)
cp $lmpath $PBS_O_WORKDIR


while getopts s:o:t: option
do
 case "${option}"
 in
 o) orient=${OPTARG}
    if [ "${#OPTARG}" -eq 0 ];
      then
        echo "ERROR: -o needs an argument:"
        echo "<111> or <110> or <100>"
        echo "exiting..."
        exit 2
    fi
    if [ "$orient" != "111" ] && [ "$orient" != "110" ] && [ "$orient" != "100" ];
      then
        echo "ERROR: the input should be 111, 110, 100"
        echo
        exit 4
    fi
  ;;
   \?)
      echo "ERROR: Invalid option: -$OPTARG" 
      echo "Exiting ..."
      exit 4
      ;;
 esac
done


if [ $orient == "no" ];
  then
  echo "ERROR: you must specify orientation with option -o"
  echo "Exiting..."
  exit -2
fi



lmp=$(echo "../"$lmp)

suffarea=$(echo $PBS_O_WORKDIR | rev | cut -d/ -f1 | cut -d_ -f1 | rev )

# Some inputs that can be changed

zwf=0.62

startstep="step3"



step="step3"
        rm -r $PBS_O_WORKDIR/$step 2> /dev/null
        cp -r ../wallscry/$step .
        cp -r ../systems/$(echo $suffarea"_fcc"$orient"-"$temp"-walls.lmp") ./$step/.

if [ ! -e  ../systems/$(echo $suffarea"_fcc"$orient"-"$temp"-walls.lmp") ]
 then
	echo "systems does not exist"
	exit
fi


#for step in step1 step3 step4
#      do
#        rm -r $PBS_O_WORKDIR/$step 2> /dev/null
#        cp -r ../wallscry/$step .
#        if [ $step == "step1" ];
#          then
#            cp -r ../systems/$(echo $suffarea"_fcc111-T1.lmp") ./$step/.
#	fi
#	cp -r ../systems/$(echo $suffarea"_fcc111-T1-walls.lmp") ./$step/.
#done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 3

step="step3"

echo $step
echo
cd $PBS_O_WORKDIR
dir=$PBS_O_WORKDIR/$step
cd $dir
input=$dir/$(echo $step.in)

datast1=$(echo $suffarea"_fcc"$orient"-"$temp".lmp")
cp  ../../systems/$datast1 $(echo "Fstep1."$zwf".data") 

cleavw=$(awk '{if($3 == "zlo"){print ($2-$1)/2.0}}' $(echo "Fstep1."$zwf".data") )

preparedir $dir


~/bin/inputs/s3zeroinpljcry.sh $(echo $suffarea"_fcc"$orient"-"$temp"-walls.lmp") $cleavw

mv step3.in s3.in

awk  -v clw=$cleavw   '{  if($2 == "clwall1"){print $1,$2,$3,clw;}
	else{print $0;} }' s3.in > step3.in

    rm -r ./dat 2> /dev/null
    mkdir ./dat

mpirun -np $PBS_NUM_PPN $lmp  < $input > $log

