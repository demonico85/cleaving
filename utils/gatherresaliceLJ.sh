#! /bin/bash


currdir=$(pwd)

rm *dat 2> /dev/null
rm *log 2> /dev/null

log=$currdir/gath.log


star=0
end=50


while getopts s:e: option
do
 case "${option}"
 in
 s) star=${OPTARG}
    if [ "${#OPTARG}" -eq 0 ];
      then
        echo "ERROR: -s needs an argument:"
        echo "<step1> or <step2> or <step3> or <step4>"
        echo "exiting..."
        exit 2
    fi
    if [ "$star" -lt 0 ] || [ "$star" -gt 50 ]; 
      then
        echo "ERROR: the input should be between 0 and 50"
        echo
        exit 3
    fi
    ;;
 e) end=${OPTARG}
    if [ "${#OPTARG}" -eq 0 ];
      then
        echo "ERROR: -e needs an argument:"
        echo "<step1> or <step2> or <step3> or <step4>"
        echo "exiting..."
        exit 2
    fi
    if [ "$end" -lt 0 ] || [ "$end" -gt 50 ];
      then
        echo "ERROR: the input should be between 0 and 50"
        echo
        exit 3
    fi
    ;;
   \?)
      echo "ERROR: Invalid option: -$OPTARG" 
      echo "Exiting ..."
      exit 4
      ;;
 esac
done

if [ $star -gt $end ]; 
  then
    echo "ERROR: $star > $end "
    echo "Exiting..."
    exit 2
fi

direction=$(pwd | rev | cut -f1 -d/ | rev | cut -f1 -d_)
temp=$(pwd | rev | cut -f1 -d_ | rev)

echo "I am taking files from the following dir"
echo $currdir 
echo




for i in `seq $star $end`
  do
    area=$(echo "area_A"$i)
    echo "Area: " $area

    step="step1"
    echo "step: " $step
    cd $currdir/$area/$step/out
    echo $currdir/$area/$step/out
    if [ ! -d "$currdir/$area/$step/out" ];
      then
        echo "ERROR: No out dir in step1"
        echo "Check better"
        echo "Exiting..."
        exit
    fi

    work.sh F
    cp F-work.dat ../../../$(echo $area"_"$step"-work.dat") 2> $log

    step="step3"
    echo "step: " $step
    cd $currdir/$area/$step/dat
    c3cry
    work.sh F
    if [ ! -e F-work.dat ];
      then
	echo "ERROR: No .out files in step3 dir"
        echo "Check better"
	echo "Exiting..."
        exit
    fi
    cp F-work.dat ../../../$(echo $area"_"$step"-work.dat") 2> $log

    step="step4"
    echo "step: " $step
    cd $currdir/$area/$step/out
    work.sh F
    cp F-work.dat ../../../$(echo $area"_"$step"-work.dat") 2> $log
done

cd $currdir

rm -r results 2> /dev/null

mkdir results


mv  *work.dat ./results/.




