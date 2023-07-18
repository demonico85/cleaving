#! /bin/bash -l
#

# Parameters to be changed

dirscripts="/home/mmm1133/Scratch/LJ_sol_liq/utils"
lmp="lmp_mpi"

trap "exit 15" TERM
export TOP_PID=$$


# ------------ Init the systems ------------------------------------------------


exampledir="../examples/lj_SL"
TT="T1" #temperature



runsim=1
re='^[0-9]+([.][0-9]+)?$'
continue=0
d=$(date +%Y-%m-%d)

echo "Starting of the simulation" 
date +"%c"
echo

# Read Inputs

while getopts l:s:n:c option
do
 case "${option}"
 in
 l) lmp=${OPTARG}
    ;;
 n) nproc=${OPTARG}
    if ! [[ $nproc =~ $re ]] ; then
       echo "ERROR: Input for -n not an integer"
       echo "Exiting ..." 
       exit 1
    fi
    ;;
 \?)
    echo "ERROR: Invalid option: -$OPTARG" 
    echo "Exiting ..."
    exit 4
    ;;
 esac
done


 if ! [[ $nproc =~ $re ]] ; then
   echo "ERROR: Need to include number of processes (option -n)"
   echo "Exiting ..."
   exit 1
fi

echo "Number of MPI processes: $nproc"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*******************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check if all the directories are in the correct place
################################################################################



if [ ! -x $lmp ];
  then
    echo "ERROR: I cannot find executable $lmp" 
    echo "Exiting..."
    exit
fi

echo "I am working with executable $lmp"


if [ ! -d $exampledir ];
  then
    echo "I cannot find $exampledir and I cannot run the test"
    echo "Exiting..."
    echo -1
fi


now=$(date +"%m_%d_%Y_%s")
log=$(echo $now.log)

#################################################################################


testdir=$(echo "test_LJ_fcc111_T1_walls")
rm -r $testdir 2> /dev/null
mkdir $testdir

cd $testdir

exampledir=../$exampledir

echo "Entering woring dir..."
echo $(pwd)
echo


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparations Dirs


echo "Preparations Dirs"
echo

dum=$exampledir/cases/fcc111_T1_walls/input


if [ ! -f $dum ];
 then
    echo "File input not found in the dir $exampledir/cases/fcc111_T1_walls"
    echo "Check it better"
    echo
    exit
fi

cp $dum .

cleav=$(awk '{if($1 == "111" || $1 == "110" || $1 == "100"){print $2}}' input)
temp=$(awk '{if($1 == "Temperature:"){print $2}}' input)

pwd
echo $exampledir/walls
    for step in step1 step2 step3 step4
      do
        cp -r $exampledir/walls/$step .
    done


CWD=$(pwd)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 1

step="step1"

echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir


mkdir out
mkdir data

zwi=$(head -1 zwalls.dat)
zwf=$(tail -1 zwalls.dat)

echo $zwi $zwf

inpdata=$(echo "fcc111-T1.lmp")
inpwall=$(echo "fcc111-T1-walls.lmp")

cp  ../$exampledir/systems/$inpdata .
cp  ../$exampledir/systems/$inpwall .

inpscript=$(echo $step".in")

cat $inpscript > tmp

cat tmp | awk -v T=$temp -v S=$inpdata -v C=$cleav '{
        if($2 == "Tsyst"){print $1, $2, $3, T;}
        else if($2 == "clwall"){print $1,$2,$3,C;}
        else if($1 == "read_data"){print $1,S;}
        else if( $2 == "nts"){print $1, $2, $3, 10;}
        else if( $2 == "firsteqnts"){print $1, $2, $3, 10;}    
        else if( $2 == "eqnts"){print $1, $2, $3, 10;}  
        else if( $2 == "Nevery"){print $1, $2, $3, 2;}
        else if( $2 == "Nrepeat"){print $1, $2, $3, 5;}    
        else if( $2 == "Nfreq"){print $1, $2, $3, 10;}                      
        else{print $0;}
        	}'  > $inpscript

cat ./in.loop/loop > tmp
cat tmp |  awk -v W=$inpwall  '{
        if($1 == "fix" && $2 == "totW"){print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,W;}
        else{print $0;}
        	}'  > ./in.loop/loop



mpirun -np $nproc ../../$lmp  -in  $inpscript > $log



if grep -q 'ERROR' $log;
  then
    echo "Something went wrong in $step. Check it better"
    exit 
fi 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 2

step="step2"

echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir


mkdir out
mkdir data

inpdata=$(echo "fcc111-T1.lmp")
inpwall=$(echo "fcc111-T1-walls.lmp")

cp  ../$exampledir/systems/$inpdata .
cp  ../$exampledir/systems/$inpwall .

inpscript=$(echo $step".in")

cat $inpscript > tmp

cat tmp | awk -v T=$temp -v S=$inpdata -v C=$cleav '{
        if($2 == "Tsyst"){print $1, $2, $3, T;}
        else if($2 == "clwall"){print $1,$2,$3,C;}
        else if($1 == "read_data"){print $1,S;}
        else if( $2 == "nts"){print $1, $2, $3, 10;}
        else if( $2 == "firsteqnts"){print $1, $2, $3, 10;}    
        else if( $2 == "eqnts"){print $1, $2, $3, 10;}  
        else if( $2 == "Nevery"){print $1, $2, $3, 2;}
        else if( $2 == "Nrepeat"){print $1, $2, $3, 5;}    
        else if( $2 == "Nfreq"){print $1, $2, $3, 10;}                      
        else{print $0;}
        	}'  > $inpscript

cat ./in.loop/loop > tmp
cat tmp |  awk -v W=$inpwall  '{
        if($1 == "fix" && $2 == "totW"){print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,W;}
        else{print $0;}
        	}'  > ./in.loop/loop



mpirun -np $nproc ../../$lmp  -in  $inpscript > $log



if grep -q 'ERROR' $log;
  then
    echo "Something went wrong in $step. Check it better"
    exit 
fi 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparation Step3

step="step3"

inpdata="inputStep3.lmp"
cd $CWD


file1=$(ls  $CWD/step1/data/F* | sort -t. -n -k2 -r | head -1) # F* to pickup the forward file
namefiles1=$(echo $file | rev | cut -d/ -f1 | rev)

file2=$(ls  $CWD/step2/data/F* | sort -t. -n -k2 -r | head -1) # F* to pickup the forward file
namefiles2=$(echo $file | rev | cut -d/ -f1 | rev)


if [ ! -e "$file1" ];
  then
    echo "File  $namefile1 does not exist"
    echo "something went wrong in step1. Check it better"
    exit 
fi

if [ ! -e "$file2" ];
  then
    echo "File  $namefile2 does not exist"
    echo "something went wrong in step2. Check it better"
    exit 
fi


cp $file1 .
cp $file2 .

a=$(echo $file1 | rev | cut -d'/' -f 1 | rev)
b=$(echo $file2 | rev | cut -d'/' -f 1 | rev)


gfortran $exampledir/utils/step3IN.f90 -o s3inp
./s3inp $a $b $cleav >  tmp
cleav2=$(awk '{if ($3 == "Position2:"){print $4}}' tmp )







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 3



echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir

mkdir out
mkdir data
mkdir dat


inpwall=$(echo "fcc111-T1-walls.lmp")

cp  ../$inpdata .
cp  ../$exampledir/systems/$inpwall .

inpscript=$(echo $step".in")

cat $inpscript > tmp
cat tmp | awk -v T=$temp -v S=$inpdata -v C=$cleav  -v cc=$cleav2 -v zwf=$zwf '{
        if($2 == "Tsyst"){print $1, $2, $3, T;}
        else if($2 == "clwall1"){print $1,$2,$3,C;}
        else if($2 == "clwall2"){print $1,$2,$3,cc;}
        else if($1 == "read_data"){print $1,S;}
        else if( $2 == "nts"){print $1, $2, $3, 10;}
        else if( $2 == "firsteqnts"){print $1, $2, $3, 10;}    
        else if( $2 == "eqnts"){print $1, $2, $3, 10;}  
        else if( $2 == "Nevery"){print $1, $2, $3, 2;}
        else if( $2 == "Nrepeat"){print $1, $2, $3, 5;}    
        else if( $2 == "Nfreq"){print $1, $2, $3, 10;}                      
        else if( $2 == "zwf"){print $1, $2, $3, zwf;}
        else{print $0;}
        	}'  > $inpscript


mpirun -np $nproc ../../$lmp  -in  $inpscript > $log



if grep -q 'ERROR' $log;
  then
    echo "Something went wrong in $step. Check it better"
    exit 
fi 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparation Step4

step="step4"


inpdata="inputStep4.lmp"
cd $CWD



echo "Now working on " $step 

file=$(ls  $CWD/step3/data/F* | sort -t. -n -k2 -r | head -1) # F* to pickup the forward fil1e
namefiles3=$(echo $file | rev | cut -d/ -f1 | rev)


cp $file .


gfortran $exampledir/utils/step4IN.f90 -o s4inp

./s4inp $namefiles3 $inpdata


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 4

echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir

mkdir out
mkdir data


inpwall=$(echo "fcc111-T1-walls.lmp")

cp  ../$inpdata .
cp  ../$exampledir/systems/$inpwall .

inpscript=$(echo $step".in")

cat $inpscript > tmp
cat tmp | awk -v T=$temp -v S=$inpdata -v C=$cleav  -v cc=$cleav2 '{
        if($2 == "Tsyst"){print $1, $2, $3, T;}
        else if($2 == "clwall1"){print $1,$2,$3,C;}
        else if($2 == "clwall2"){print $1,$2,$3,cc;}
        else if($1 == "read_data"){print $1,S;}
        else if( $2 == "nts"){print $1, $2, $3, 10;}
        else if( $2 == "firsteqnts"){print $1, $2, $3, 10;}    
        else if( $2 == "eqnts"){print $1, $2, $3, 10;}  
        else if( $2 == "Nevery"){print $1, $2, $3, 2;}
        else if( $2 == "Nrepeat"){print $1, $2, $3, 5;}    
        else if( $2 == "Nfreq"){print $1, $2, $3, 10;}                      
        else{print $0;}
        	}'  > $inpscript


mpirun -np $nproc ../../$lmp  -in  $inpscript > $log



if grep -q 'ERROR' $log;
  then
    echo "Something went wrong in $step. Check it better"
    exit 
fi 


