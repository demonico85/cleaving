#! /bin/bash -l
#

function checkcommands () {

if ! command -v $1 &> /dev/null
then
    echo "ERROR: $1 could not be found"
    exit
fi

}


# Parameters to be changed

lmp="lmp_mpi"

trap "exit 15" TERM
export TOP_PID=$$

checkcommands "mpirun"
checkcommands "gfortran"

# ------------ Init the systems ------------------------------------------------



runsim=1
nproc=4
re='^[0-9]+([.][0-9]+)?$'
continue=0
d=$(date +%Y-%m-%d)

echo "Starting of the simulation" 
date +"%c"
echo

# Read Inputs

while getopts l:n:h option
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
 h) echo "cleav_regression.sh -l <path lammps executable> OPTIONAL -n <nproc, default=4>"
    ;;
 \?)
    echo "ERROR: Invalid option: -$OPTARG" 
    echo "Exiting ..."
    exit 4
    ;;
 esac
done

lmp=$(realpath $lmp)
checkcommands $lmp

if ! [[ $nproc =~ $re ]] ; then
   echo "ERROR: Need to include number of processes (option -n)"
   echo "Exiting ..."
   exit 1
fi

echo "Number of MPI processes: $nproc"


# Check if all the commands are available



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*******************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check if all the directories are in the correct place
################################################################################



echo "I am working with executable $lmp"


if [ ! -d $exampledir ]; then
    echo "ERROR: I cannot find $exampledir and I cannot run the test"
    echo "Exiting..."
    echo -1
fi


now=$(date +"%m_%d_%Y_%s")
log=$(echo $now.log)

#################################################################################

echo "Starting regression test 1"

if [ ! -d iRT1 ];
  then
   echo "I cannot find dir iRT1 to run the test"
   echo "Check your directory"
   echo "Exiting..."
   exit
fi

rm -rf test_iRT1
cp -r iRT1 test_iRT1

echo "Entering working dir..."
cd ./test_iRT1
echo $(pwd)
echo

CWD=$(pwd)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 1

step="step1"

echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir

inpscript=$(echo $step".in")



rm -r ./out/ 2> /dev/null
mkdir out


mpirun -np $nproc $lmp  -in  $inpscript > $log



if grep -q 'ERROR' $log; then
    echo "Something went wrong in $step"
    echo "Test cannot be completed"
    echo "Consult the log file $log"
    exit 
else
    echo "$step successfully run"
    echo 
fi 

cd ./out/

for i in `seq 1 500`
  do
  file=$(echo "ave.F."$i".out")
  if [ ! -f $file ];
     then
     break
   fi
   awk 'NR>3 {l=$5;a=$4} END {print "Wells interactions step1:", a}' $file > $CWD/results.dat
done



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 3


step="step3"
echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir
inpscript=$(echo $step".in")



rm -r ./dat 2> /dev/null
mkdir dat

mpirun -np $nproc $lmp  -in  $inpscript > $log



if grep -q 'ERROR' $log; then
    echo "Something went wrong in $step"
    echo "Test cannot be completed"
    echo "Consult the log file $log"
    exit 
else
    echo "$step successfully run"
    echo 
fi 


cd ./dat
for i in `seq 1 500`
  do
  file=$(echo "inters3."$i".dat")
  if [ ! -f $file ];
     then
     break
   fi
   awk 'NR>7{w+=$3+$4;} END {print "Switch-off interactions work step3:", w}' $file >> $CWD/results.dat
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 4


step="step4"
echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir
inpscript=$(echo $step".in")

rm -r ./out 2> /dev/null
mkdir out

mpirun -np $nproc $lmp  -in  $inpscript > $log

if grep -q 'ERROR' $log; then
    echo "Something went wrong in $step"
    echo "Test cannot be completed"
    echo "Consult the log file $log"
    exit 
else
    echo "$step successfully run"
    echo 
fi 

cd ./out/

for i in `seq 1 500`
  do
  file=$(echo "ave.F."$i".out")
  if [ ! -f $file ];
     then
     break
   fi
   awk 'NR>3{l=$5;a=$4;} END {print "Wells interactions step4:", a}' $file >> $CWD/results.dat
done

cd $CWD
cd ../


#################################################################################

if [ ! -d iRT2 ];
  then
   echo "I cannot find dir RT2 to run the test"
   echo "Check your directory"
   echo "Exiting..."
   exit
fi

rm -rf test_iRT2
cp -r iRT2 test_iRT2

echo "Entering working dir..."
cd ./test_iRT2
echo $(pwd)
echo

CWD=$(pwd)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 1

step="step1"

echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir

inpscript=$(echo $step".in")



rm -r ./out/ 2> /dev/null
mkdir out


mpirun -np $nproc $lmp  -in  $inpscript > $log

if grep -q 'ERROR' $log; then
    echo "Something went wrong in $step"
    echo "Test cannot be completed"
    echo "Consult the log file $log"
    exit 
else
    echo "$step successfully run"
    echo 
fi 

cd ./out/

for i in `seq 1 500`
  do
  file=$(echo "ave.F."$i".out")
  if [ ! -f $file ];
     then
     break
   fi
    awk 'NR>3{a=$5;} END {print "Walls interactions step1:", a}' $file > $CWD/results.dat
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 2

step="step2"

echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir

inpscript=$(echo $step".in")

rm -r ./out/ 2> /dev/null
mkdir out


mpirun -np $nproc $lmp  -in  $inpscript > $log



if grep -q 'ERROR' $log; then
    echo "Something went wrong in $step"
    echo "Test cannot be completed"
    echo "Consult the log file $log"
    exit 
else
    echo "$step successfully run"
    echo 
fi 

cd ./out/

for i in `seq 1 500`
  do
  file=$(echo "ave.F."$i".out")
  if [ ! -f $file ];
     then
     break
   fi
    awk 'NR>3{a=$5;} END {print "Walls interactions step2:", a}' $file >> $CWD/results.dat
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 3


step="step3"
echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir
inpscript=$(echo $step".in")

rm -r ./dat 2> /dev/null
mkdir dat

mpirun -np $nproc $lmp  -in  $inpscript > $log



if grep -q 'ERROR' $log; then
    echo "Something went wrong in $step"
    echo "Test cannot be completed"
    echo "Consult the log file $log"
    exit 
else
    echo "$step successfully run"
    echo 
fi 

cd ./dat
for i in `seq 1 500`
  do
  file=$(echo "inters3."$i".dat")
  if [ ! -f $file ];
     then
     break
   fi
    awk 'BEGIN{ntypes=8}{i+=1; if(i>13){
     for (i=1; i< ntypes+1; i++){
     	getline
     	  for(j=1; j< ntypes+1; j++){
     	      k=j+2
     	      M[i,j]=$k;
     		}
    	}
    }
  }  
    END{
     CC=M[1,5]+ M[2,7]+ M[3,8]+M[4,6]+ M[1,6]+ M[2,8] + M[4,5]+ M[3,7]
     DD = M[1,2]+ M[1,3]+ M[2,4]+M[3,4]
    print "Switch interaction step3:", CC+2.0*DD}' $file >> $CWD/results.dat
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 4


step="step4"
echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir
inpscript=$(echo $step".in")

rm -r ./out 2> /dev/null
mkdir out


mpirun -np $nproc $lmp  -in  $inpscript > $log


if grep -q 'ERROR' $log; then
    echo "Something went wrong in $step"
    echo "Test cannot be completed"
    echo "Consult the log file $log"
    exit 
else
    echo "$step successfully run"
    echo 
fi 

cd ./out/

for i in `seq 1 500`
  do
  file=$(echo "ave.F."$i".out")
  if [ ! -f $file ];
     then
     break
   fi
   awk 'NR>3{l=$5;a=$4;} END {print "Walls interactions step4:", a}' $file >> $CWD/results.dat
done

cd $CWD
cd ..

echo "#####################################"
echo "Results from iRT1 (solid-vacuum LJ):"

cat test_iRT1/results.dat

echo "#####################################"

echo

echo "#####################################"
echo "Results from iRT2 (solid-liquid LJ):"

cat test_iRT2/results.dat

echo "#####################################"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

