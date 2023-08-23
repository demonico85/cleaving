#! /bin/bash

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

# -------- Load some modules for the simulations, set up env variables ---------
# cluster dependent




##PBS_NUM_PPN=20
##PBS_O_WORKDIR=$(pwd)


lmp="lmp_mpi_20Jun"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $PBS_O_WORKDIR
echo "Entering woring dir..."
echo $(pwd)
echo

lmpath=$(echo "./$lmp")
if [ ! -e  $lmpath ];
  then
    echo "ERROR: Lammps executable not found in $lmpath"
    exit
fi

temp="T0"
startstep="no"
orient="no"
temperature="no"
runsim=1
now=$(date +"%m_%d_%Y_%s")
log=$(echo $now.log)
cp $lmpath $PBS_O_WORKDIR

lmp=$(echo "../"$lmp)

suffarea=$(echo $PBS_O_WORKDIR | rev | cut -d/ -f1 | cut -d_ -f1 | rev )

# Some inputs that can be changed

zwf=0.62


while getopts s:o:t:w: option
do
 case "${option}"
 in
 s) startstep=${OPTARG}
    if [ "${#OPTARG}" -eq 0 ];
      then
        echo "ERROR: -s needs an argument:"
        echo "<step1> or <step2> or <step3> or <step4>"
        echo "exiting..."
        exit 2
    fi
    if [ "$startstep" != "step1" ] && [ "$startstep" != "step2" ] && [ "$startstep" != "step3" ] && [ "$startstep" != "step4" ];
      then
        echo "ERROR: the input should be stepX with X=1,2,3,4"
        echo
        exit 3
    else
        runsim=0
    fi
    ;;
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
 t) temperature=${OPTARG}
    if [ "${#OPTARG}" -eq 0 ];
      then
        echo "ERROR: -t needs an argument:"
        echo "Temperature"
        echo "exiting..."
        exit 6
    fi
    if  (( $(echo "$temperature < 0.0" | bc -l) ));
      then
        echo "ERROR: Temperature must be positive "
        echo "Exiting..."
        echo -1
    fi
    ;;
 w) zwf=${OPTARG}
    if [ "${#OPTARG}" -eq 0 ];
      then
        echo "ERROR: -w needs an argument:"
        echo "Temperature"
        echo "exiting..."
        exit 7
    fi
    if  (( $(echo "$zwf < 0.0" | bc -l) ));
      then
        echo "ERROR: zwf must be positive "
        echo "Exiting..."
        echo -1
    fi
    ;;
   \?)
      echo "ERROR: Invalid option: -$OPTARG" 
      echo "Exiting ..."
      exit 4
      ;;
 esac
done

if [ $temperature == "no" ];
  then
    echo "Temperature needs to be specified (option -t X)"
    echo "Exiting..."
    echo
    exit 5
fi



if [ $orient == "no" ];
  then
    echo "orientation needs to be specified (option -o X)"
    echo "Exiting..."
    echo
    exit 5
fi




if [ $runsim -eq 1 ];
  then
    for step in step1 step3 step4
      do
        rm -r $PBS_O_WORKDIR/$step 2> /dev/null
        cp -r ../wallscry/$step .
        if [ $type == "step1" ];
          then
            cp -r ../systems/$(echo $suffarea"_fcc"$orient"-"$temp".lmp") ./$step/.
        fi
        cp -r ../systems/$(echo $suffarea"_fcc"$orient"-"$temp"-walls.lmp") ./$step/.
    done
else
    nin=${startstep:4:1}
    for i in `seq $nin 4`;
      do
        step=$(echo "step"$i)
        rm -r $PBS_O_WORKDIR/$step 2> /dev/null
        cp -r ../wallscry/$step . 2> /dev/null
        if [ $step == "step1" ];
          then
	    cp -r ../systems/$(echo $suffarea"_fcc"$orient"-"$temp".lmp") ./$step/.
        fi
        cp -r ../systems/$(echo $suffarea"_fcc"$orient"-"$temp"-walls.lmp") ./$step/
    done

fi


#for step in step1 step3 step4
#      do
#        rm -r $PBS_O_WORKDIR/$step 2> /dev/null
#        cp -r ../wallscry/$step .
#        if [ $step == "step1" ];
#          then
#            cp -r ../systems/$(echo $suffarea"_fcc"$orient"-T1.lmp") ./$step/.
#	fi
#	cp -r ../systems/$(echo $suffarea"_fcc"$orient"-T1-walls.lmp") ./$step/.
#done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 1

if [ $startstep == "step1" ] || [ $runsim -eq 1 ];
  then
    runsim=1

step="step1"

echo $step
echo
cd $PBS_O_WORKDIR
dir=$PBS_O_WORKDIR/$step
cd $dir
input=$dir/$(echo $step.in)

preparedir $dir

datast1=$(echo $suffarea"_fcc"$orient"-"$temp".lmp")
cleavw=$(awk '{if($3 == "zlo"){print ($2-$1)/2.0}}' $datast1 )


mv step1.in step1.in-bck

awk -v namefile=$(echo $suffarea"_fcc"$orient"-"$temp".lmp") -v cl=$cleavw  -v t=$temperature -v wall=$(echo $suffarea"_fcc"$orient"-"$temp"-walls.lmp") '{
		if($1 == "read_data"){print $1, namefile }
        else if($2 == "clwall"){print $1,$2,$3,cl}
        else if($2 == "Tsyst"){print $1, $2, $3,t}
		else if($1 == "fix" && $2 == "f2"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall} 
		else{print $0;}}' step1.in-bck  > step1.in

cp in.loop/loop in.loop/loop_bck

awk -v wall=$(echo $suffarea"_fcc"$orient"-"$temp"-walls.lmp") '{
        if($1 == "fix" && $2 == "f2"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall}
        else {print $0;}
}'  in.loop/loop_bck > in.loop/loop



mpirun -np $PBS_NUM_PPN $lmp  < $input > $log

fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 3

if [ $startstep == "step3" ] || [ $runsim -eq 1 ] ;
  then
runsim=1
step="step3"

echo $step
echo
cd $PBS_O_WORKDIR
dir=$PBS_O_WORKDIR/$step
cd $dir
input=$dir/$(echo $step.in)

cp ../step1/data/$(echo "Fstep1."$zwf".data") .

preparedir $dir

cleavw=$(awk '{if($2 == "clwall"){print $4}}' ../step1/step1.in )

~/bin/s3inpljcry.sh $(echo $suffarea"_fcc"$orient"-"$temp"-walls.lmp") $cleavw $temperature $zwf


    rm -r ./dat 2> /dev/null
    mkdir ./dat

mpirun -np $PBS_NUM_PPN $lmp  < $input > $log

fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparation Step4

if [ $startstep == "step4" ] || [ $runsim -eq 1 ];
  then

runsim=1
step="step4"

cd $PBS_O_WORKDIR
cd $PBS_O_WORKDIR/step1
cleavw=$(awk '{if($2 == "clwall"){print $4}}' step1.in )
cd $PBS_O_WORKDIR

dir=$PBS_O_WORKDIR/$step
cd $dir

cntlam=$(wc -l < $(echo $PBS_O_WORKDIR"/step3/lambda.txt"))

file=$(echo $PBS_O_WORKDIR"/step3/data/step3."$cntlam".data")

if [ ! -e "$file" ];
  then
    echo "File  $file does not exist"
    echo "something went wrong in step3. Check it better"
    exit
fi

preparedir $dir

cp $file .
echo "STEP4 QUI"
echo $(echo $suffarea"_fcc"$orient"-$temp-walls.lmp") $cntlam $cleavw
~/bin/s4inpljcry.sh $(echo $suffarea"_fcc"$orient"-$temp-walls.lmp") $cntlam $cleavw $temperature


input=$dir/$(echo $step.in)

mpirun -np $PBS_NUM_PPN $lmp  < $input > $log

fi


