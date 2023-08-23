#! /bin/bash -l
#

# Parameters to be changed

dirscripts="./"
lmp="lmp_mpi"


trap "exit 15" TERM
export TOP_PID=$$

# ---------------------- SOME FUNCTION DEFINITION ------------------------------

function preparedir {


    rm $1/ave.F.* 2> /dev/null
    rm $1/ave.B.* 2> /dev/null
    rm $1/force.dump 2> /dev/null
    rm -r $1/dump 2> /dev/null
    mkdir $1/dump

    rm -r $1/group 2> /dev/null
    mkdir $1/group

    rm -r $1/X.restart 2> /dev/null
    mkdir $1/X.restart


    rm -r $1/data 2> /dev/null
    mkdir $1/data

    rm -r $1/out 2> /dev/null
    mkdir $1/out

    rm $1/forces.log 2> /dev/null
}


function continuestep3 {

    step="step3"
    fulldir=$(echo $1/$step)
    delam=$(awk '{ if($1 == "variable" && $2 == "delam" ){print $4}}' $( echo $fulldir/$step".in") )
    for datafile in `ls $fulldir/data/ | sort -r -n -t . -k 2`
      do
        break
    done

#most recent file 

    recfile=$(ls -Art $fulldir/data/ | tail -n 1)



    if [ ${recfile:0:1} == "B" ];
      then
        echo "NON ancora implementato per partire dai file backwards Bstep..."
        echo "exit"
        kill -s TERM $TOP_PID
   fi


    cp $fulldir/data/$datafile .

# prendo il numero da data file invece che da out perche' potrebbe essere che out esista senza che quel 
# punto abbia finito. Invece data viene creato quando ha finito quel valore di lambda

    numb=$( echo $datafile  | cut -f2 -d.)
# leggo il valore attual di lambda dalla dir data
     lambda=$(awk '{i+=1; if(i==4){print $2}}' ./dat/$(echo "lambda."$numb".dat"))
#    lambda=$(echo "$delam * $numb" | bc )
#    let numb=numb-1


    if [ -z $lambda ];
      then
        echo "ERROR: problems in reading lambda at cnt=$numb"
        echo "Exiting..."
        kill -s TERM $TOP_PID
        exit
   fi

    
    cp $( echo $fulldir/$step".in") $( echo $fulldir/$step"-restart.in")
    
    awk -v lam=$lambda -v fil=$datafile -v nn=$numb 'BEGIN {lb=0;ct=0;}{
        if( $1 == "variable" && $2 == "lambda" && lb == 0){ print $1,$2,$3,lam; lb=1;}
        else if($1 == "read_data"){print $1,fil}
        else if($1 == "variable" && $2 == "cnt" && ct == 0){ print $1,$2,$3,nn; ct=1;}
        else {print $0;}
    }' $( echo $fulldir/$step"-restart.in") > $( echo $fulldir/$step".in")

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DEBUGGING

#CWD=$(pwd)
#PBS_NUM_PPN=4
#PBS_NUM_PPN=$NSLOTS


# ------------ Init the systems ------------------------------------------------



runsim=1
startstep="no"
pbcpressure=0
re='^[0-9]+$'
continue=0
d=$(date +%Y-%m-%d)

echo "Starting of the simulation" 
date +"%c"
echo

#PBS_NUM_PPN=$NSLOTS


# Read Inputs

while getopts l:s:n:c option
do
 case "${option}"
 in
 l) lmp=${OPTARG}
    ;;
 n) PBS_NUM_PPN=${OPTARG}
    if ! [[ $PBS_NUM_PPN =~ $re ]] ; then
       echo "ERROR: Input for -n not an integer"
       echo "Exiting ..." 
       exit 1
    fi
    ;;
# p) pbcpressure=1
#    ;;
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
 c) continue=1
    ;; 
 \?)
    echo "ERROR: Invalid option: -$OPTARG" 
    echo "Exiting ..."
    exit 4
    ;;
 esac
done

     if [ $continue -eq 1 ];
       then
        if [ "$startstep" != "step3" ];
          then
            echo "For now this option works only with step3."
            echo "Please select step3"
            echo
            exit 4
         fi
      fi


 if ! [[ $PBS_NUM_PPN =~ $re ]] ; then
   echo "ERROR: Need to include number of processes (option -n)"
   echo "Exiting ..."
   exit 1
fi

echo "Number of MPI processes: $PBS_NUM_PPN"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*******************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################################################################################



if [ ! -x $lmp ];
  then
    echo "ERROR: I cannot find executable $lmp" 
    echo "Exiting..."
    exit
fi

echo "I am working with executable $lmp"

#################################################################################

echo "Entering woring dir..."
echo $(pwd)
echo

CWD=$(pwd)


now=$(date +"%m_%d_%Y_%s")
log=$(echo $now.log)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparations Dirs


echo "Preparations Dirs"
echo

TT=$(echo $CWD | rev | cut -d/ -f1 | rev | cut -d_ -f2) # temperature T1, T2 or T3

if [ ! -e input ];
 then
    fff=$(pwd)
    echo "File input not found in the current directory"
    echo "Check it better"
    echo $fff
    echo
    exit
fi


re='^[0-9]+([.][0-9]+)?$'

type=$(awk '{if($1 == "walls"){print "wa";}else if($1 == "wells"){print "we";}}' input)
cryst=$(awk '{if($1 == "111" || $1 == "110" || $1 == "100"){print $1}}' input)
temp=$(awk '{if($1 == "Temperature:" ){print $2}}' input)
cleav=$(awk '{if($1 == "111" || $1 == "110" || $1 == "100"){print $2}}' input)

if [ $type == "wa" ];
  then 
    zwf=$(awk '{if($1 == "zwf:" ){print $2}}' input)
    zwi=$(awk '{if($1 == "zwi:" ){print $2}}' input)


    if [ ! -z "$zwf" ];
      then 
        if ! [[ "$zwf" =~ $re ]]; then
           echo "Error: zwf in input file Not a number"
           echo "zwf must be added when walls are used" >&2; exit 1    
        fi
        zwf=$(echo $zwf | rev | bc -l | rev) #eliminate the last zero if the number is in the form 0.40, because Lammps write only 0.4
        dew=$(awk '{if($1 == "dew:" ){print $2}}' input) # delta z for walls
    fi

fi

rhoc=$(awk '{if($1 == "rhoc:" ){print $2}}' input)
fact=$(awk '{if($1 == "fact:" ){print $2}}' input)
dw=$(awk '{if($1 == "dw:" ){print $2}}' input)   # delta r for wells
backloop=$(awk '{if($1 == "backloop:" ){print $2}}' input)
dew4=$(awk '{if($1 == "dew4:" ){print $2}}' input) #step4 can have a different spacing in zwall


a=$(awk '{if($1 == "cell:" ){print $2}}' input)
b=$(awk '{if($1 == "cell:" ){print $3}}' input)
c=$(awk '{if($1 == "cell:" ){print $4}}' input)


if ! [[ "$cleav" =~ $re ]];
  then
    echo "Error: cleav in input file is not a number or is not defined"
    echo "the position of the cleaving plane (cleav)  must be always added " >&2; exit 1 
fi


if [ $type == "wa" ];
  then
    zwi=$(awk '{print $1; exit;}' $(echo "zwall_"$cryst"_"$TT".dat")) 
    zcheckf=$(awk '{z=$1;} END{print z}' $(echo "zwall_"$cryst"_"$TT".dat"))

    zchecki=$(awk '{z=$1;} END{print z}' $(echo "zwall_back_"$cryst"_"$TT".dat"))
    zwf=$(awk '{print $1; exit;}' $(echo "zwall_back_"$cryst"_"$TT".dat"))

    if [ ! $zwi == $zchecki ] || [ ! $zwf == $zcheckf ];
      then
       echo "ERROR: boundaries in $(echo "zwall_"$cryst"_"$TT".dat")) and $(echo "zwall_back_"$cryst"_"$TT".dat")) do not coincide"
       echo "$zwi != $zchecki or $zwf != $zcheckf"
       echo "Check the files"
       echo "Exiting..."
       exit -10
    fi
fi




if [ $runsim -eq 1 ]; 
  then
    for step in step1 step2 step3 step4
      do
        rm -r $CWD/$step 2> /dev/null
        if [ $type == "we" ];
          then
            cp -r ../wells/$step .
        else
            cp -r ../walls/$step .
        fi
    done
else   
    nin=${startstep:4:1}
    if [ $continue -eq 1 ];
      then
        let nin=nin+1
    fi
    for i in `seq $nin 4`;
      do
        step=$(echo "step"$i)
        rm -r $CWD/$step 2> /dev/null
        if [ $type == "we" ];
          then
            cp -r ../wells/$step .
        else
            cp -r ../walls/$step .
        fi
    done

fi

#-z string True if the string is null (an empty string)


if [ -z "$rhoc" ];
 then
    rhoc=-100
fi
if [ -z "$zwi" ];
 then
    zwi=-100
fi
if [ -z "$zwf" ];
 then
    zwf=-100
fi
if [ -z "$dew" ];
 then
    dew=-100
fi
if [ -z "$dw" ];
 then
    dw=-100
fi
if [ -z "$fact" ];
 then
    fact=-100
fi
if [ -z "$a" ];
 then
    a=-100
fi
if [ -z "$b" ];
 then
    b=-100
fi
if [ -z "$c" ];
 then
    c=-100
fi
if [ -z "$dew4" ];
 then
    dew4=$dew
fi
if [ -z "$backloop" ];
 then
    backloop=0
fi

#echo "QUII" $cryst $type $temp $cleav $TT $rhoc $zwi $zwf $dew $fact $dw $a $b $c $backloop $dew4  #$pbcpressure
#QUII 111 we 0.617 T1 0.945 -100 -100 -100 1.0 6 -100 -100 -100 0 -100



$dirscripts/preparation.sh $cryst $type $temp $cleav $TT $rhoc $zwi $zwf $dew $fact $dw $a $b $c $backloop $dew4  #$pbcpressure


if [ $? -ne 0 ];
  then
   echo "preparation.sh failed. Check it out"
   echo "Exiting ..."
   exit -1 
fi 


#
if [ $startstep == "step1" ] || [ $runsim -eq 1 ];
  then
    runsim=1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 1

step="step1"

echo $step
echo
cd $CWD
dir=$CWD/$step

if [ $type == "wa" ];
  then
    cp $CWD/$(echo "zwall_"$cryst"_"$TT".dat") $dir/zwalls.dat
fi


cd $dir
preparedir $dir 
input=$dir/$(echo $step.in)



#mpirun -np $PBS_NUM_PPN $lmp  -in  $input > $log

gerun $lmp  -in  $input > $log


if grep -q 'ERROR' $log;
  then
    echo "Something went wrong in $step . Check it better"
    exit 
fi 

fi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 2

if [ $startstep == "step2" ] || [ $runsim -eq 1 ];
  then
    runsim=1


step="step2"

echo $step
echo
cd $CWD


dir=$CWD/$step

if [ $type == "wa" ];
  then
    cp $CWD/$(echo "zwall_"$cryst"_"$TT".dat") $dir/zwalls.dat
fi



cd $dir
preparedir $dir 
input=$dir/$(echo $step.in)



#mpirun -np $PBS_NUM_PPN $lmp  -in  $input > $log

gerun $lmp  -in  $input > $log

if grep -q 'ERROR' $log;
  then
    echo "Something went wrong in $step . Check it better"
    exit 
fi 

fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparation Step3

if [ $startstep == "step3" ] || [ $runsim -eq 1 ];
  then
    runsim=1

inputnames3="inputStep3.lmp"
step="step3"
cd $CWD


file1=$(ls  $CWD/step1/data/F* | sort -t. -n -k2 -r | head -1) # F* to pickup the forward file
namefiles1=$(echo $file | rev | cut -d/ -f1 | rev)

file2=$(ls  $CWD/step2/data/F* | sort -t. -n -k2 -r | head -1) # F* to pickup the forward file
namefiles2=$(echo $file | rev | cut -d/ -f1 | rev)


#if [ $type == "we" ];
#  then
#    file1=$(echo $CWD/step1/data/Fstep1.1.data)
#    file2=$(echo $CWD/step2/data/Fstep2.1.data)
#else
#    file1=$(echo $CWD/step1/data/Fstep1.$zwf.data)
#    file2=$(echo $CWD/step2/data/Fstep2.$zwf.data)
#fi


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

rm $CWD/step3/$inputnames3 2> /dev/null 

echo $file1
echo $file2

cp $file1 .
cp $file2 .

a=$(echo $file1 | rev | cut -d'/' -f 1 | rev)
b=$(echo $file2 | rev | cut -d'/' -f 1 | rev)

$dirscripts/compiled/s3inp $a $b $cleav >  tmp

cleav2=$(awk '{if ($3 == "Position2:"){print $4}}' tmp )
cp $CWD/$step/$(echo $step".in") $CWD/$step/$(echo $step".in-back")
awk -v c=$cleav2 -v nnm=$inputnames3 '{if($2 == "clwall2"){print $1,$2,$3,c}
            else if($1 == "read_data"){print $1,nnm}
            else{print $0}}' $CWD/$step/$(echo $step".in-back")  > $CWD/$step/$(echo $step".in")

mv $inputnames3 $CWD/step3




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 3



echo $step
echo
cd $CWD
dir=$CWD/$step
cd $dir

if [ $continue -eq 1 ];
  then
    continuestep3 $CWD
    cp info.log $(echo "info."$d".log")
    now=$(date +"%m_%d_%Y_%s")
    log=$(echo $now.log)
else
    preparedir $dir 
    rm -r $dir/dat 2> /dev/null
    mkdir $dir/dat
    rm info.log 2> /dev/null
fi



input=$dir/$(echo $step.in)

gerun $lmp  -in  $input > $log
#mpirun -np $PBS_NUM_PPN $lmp  -in  $input > $log

if tail -50 $log | grep -q 'ERROR';
  then
    echo "Something went wrong in $step . Check it better"
    exit -2 
fi 

fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparation Step4

if [ $startstep == "step4" ] || [ $runsim -eq 1 ];
  then
    runsim=1

    inputnames4="inputStep4.lmp"
    if [ $startstep == "step4" ];
      then
        if [ ! -e $CWD/step3/step3.in ];
          then
            echo "ERROR: Something went wrong"
            echo "I cannot find the step3.in file in step3 dir. Check it better"
            echo
            exit
        fi
        cleav2=$(awk '{if($2 == "clwall2"){print $4;}}' $CWD/step3/step3.in)
    fi

step="step4"
echo "Now working on " $step 

cd $CWD
file=$(ls  $CWD/step3/data/F* | sort -t. -n -k2 -r | head -1) # F* to pickup the forward fil1e
namefiles3=$(echo $file | rev | cut -d/ -f1 | rev)

echo "STEP3 $namefiles3"
echo "STEP3 $file"

if [ ! -e "$file" ];
  then
    echo "File  $file does not exist"
    echo "something went wrong in step3. Check it better"
    exit 
fi
rm $step/$inputnames4 2> /dev/null

cp $file .
$dirscripts/compiled/s4inp $namefiles3 $inputnames4
mv $inputnames4 $CWD/$step


cp $CWD/$step/$(echo $step".in") $CWD/$step/$(echo $step".in-back")
awk -v c=$cleav2  -v nnm=$inputnames4 '{if($2 == "clwall2"){print $1,$2,$3,c}
        else if($1 == "read_data"){print $1,nnm}
        else{print $0}}' $CWD/$step/$(echo $step".in-back")  > $CWD/$step/$(echo $step".in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STEP 4

dir=$CWD/$step

if [ $type == "wa" ];
  then
    cp $CWD/$(echo "zwall_back_"$cryst"_"$TT".dat") $dir/zwalls_back.dat
fi



cd $dir
preparedir $dir 
input=$dir/$(echo $step.in)
#mpirun -np $PBS_NUM_PPN $lmp  -in  $input > $log

gerun $lmp  -in  $input > $log

if grep -q 'ERROR' $log;
  then
    echo "Something went wrong in $step . Check it better"
    exit 
fi 

fi
