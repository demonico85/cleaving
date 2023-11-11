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

if [ ! -d Rtest1 ];
  then
   echo "I cannot find dir Rtest1 to run the test"
   echo "Check your directory"
   echo "Exiting..."
   exit
fi

cp -r Rtest1 RT1

echo "Entering working dir..."
cd ./RT1
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
    awk '{i+=1; if(i==3){print $5, $4}}' $file >> $(echo $step".dat")
done


cd $CWD
cp ./$step/out/$(echo $step".dat") .



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
    awk '{i+=1; if(i>4){w+=$3+$4}}END{print w}' $file >> s3.tmp
done

paste ../zdir.dat s3.tmp > $(echo $step".dat")

cd $CWD
cp ./$step/dat/$(echo $step".dat") .

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
    awk '{i+=1; if(i==3){print $5, $4}}' $file >> $(echo $step".dat")
done


cd $CWD
cp ./$step/out/$(echo $step".dat") .

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm work.dat 2> /dev/null

for st in step1 step3 step4
do
awk '{n+=1; x[n]=$1; y[n]=$2}
END{
      for (j=1;j<n;j++){
	    yterm = 0.5*(y[j+1]+y[j]); 
        xterm = (x[j+1]-x[j]); 
        z = yterm * xterm;
        intg +=z;
	}
	print intg;
}' $(echo $st".dat") >> work.dat

done

SFE=$(awk '{sfe+=$1}END{print sfe}' work.dat)

echo "Surface Free Energy is: " $SFE



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




cd ../

echo "Starting regression test 2"

#################################################################################

if [ ! -d Rtest2 ];
  then
   echo "I cannot find dir Rtest2 to run the test"
   echo "Check your directory"
   echo "Exiting..."
   exit
fi

rm -r RT2 2> /dev/null

cp -r Rtest2 RT2

echo "Entering working dir..."
cd ./RT2
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
    awk '{i+=1; if(i==3){print $4, $5}}' $file >> $(echo $step".dat")
done


cd $CWD
cp ./$step/out/$(echo $step".dat") .

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
    awk '{i+=1; if(i==3){print $4, $5}}' $file >> $(echo $step".dat")
done


cd $CWD
cp ./$step/out/$(echo $step".dat") .

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
    awk 'BEGIN{ntypes=8}{i+=1; if(i>4){
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
    print CC+2.0*DD}' $file >> s3.tmp
done

paste ../lambda.dat s3.tmp > $(echo $step".dat")

cd $CWD
cp ./$step/dat/$(echo $step".dat") .

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
    awk '{i+=1; if(i==3){print $5, $4}}' $file >> $(echo $step".dat")
done


cd $CWD
cp ./$step/out/$(echo $step".dat") .

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm work.dat 2> /dev/null

for st in step1 step2 step3 step4
do
awk '{n+=1; x[n]=$1; y[n]=$2}
END{
      for (j=1;j<n;j++){
	    yterm = 0.5*(y[j+1]+y[j]); 
        xterm = (x[j+1]-x[j]); 
        z = yterm * xterm;
        intg +=z;
	}
	print intg;
}' $(echo $st".dat") >> work.dat

done

SFE=$(awk '{sfe+=$1}END{print sfe}' work.dat)

echo "Surface Free Energy is: " $SFE



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

