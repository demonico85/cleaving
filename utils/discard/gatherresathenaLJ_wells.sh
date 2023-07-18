#! /bin/bash


currdir=$(pwd)

rm *dat 2> /dev/null
rm *log 2> /dev/null

log=$currdir/gath.log

#athenadir=/gpfs/home/leics/ndp8/ljcrystal/cleavcry

direction=$(pwd | rev | cut -f1 -d/ | rev | cut -f1 -d_)
orient=$(pwd | rev | cut -f2 -d/ | rev)
temp=$(pwd | rev | cut -f1 -d_ | rev)
athenadir=$(echo "/gpfs/home/leics/ndp8/ljcrystal/fcc"$orient"/"$temp"/"$direction)


athenadir="/gpfs/home/leics/ndp8/ljcrystal/fcc111/wells/T01/exx"

echo "I am taking files from the following dir"
echo $athenadir
echo

##for area in A0 A1;

beg=0
end=50

#for area in A0 A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12;

for i in `seq $beg $end`
  do
    area=$(echo "A"$i)
    cd $currdir
    rm -r $area 2> /dev/null
    mkdir $area
    cd $currdir/$area

    echo "Area: " $area

    step="step1"
    echo "step: " $step
    scp -r ndp8@athena.hpc-midlands-plus.ac.uk:$athenadir/$(echo "area_"$area)/$step/out $currdir/$area/$step 2> $log
    cd $currdir/$area/$step
    work.sh F
    cp F-work.dat ../../$(echo $area"_"$step"-work.dat") 2> $log

    step="step3"
    echo "step: " $step
    scp -r ndp8@athena.hpc-midlands-plus.ac.uk:$athenadir/$(echo "area_"$area)/$step/dat $currdir/$area/$step 2> $log 
    cd $currdir/$area/$step
    c3crywells
    work.sh F
    if [ ! -e F-work.dat ];
      then
	echo "ERROR: No .out files in step3 dir"
        echo "Check better"
	echo "Exiting..."
        exit
    fi
    cp F-work.dat ../../$(echo $area"_"$step"-work.dat") 2> $log

    step="step4"
    echo "step: " $step
    scp -r ndp8@athena.hpc-midlands-plus.ac.uk:$athenadir/$(echo "area_"$area)/$step/out $currdir/$area/$step 2> $log
    cd $currdir/$area/$step
    work.sh F
    cp F-work.dat ../../$(echo $area"_"$step"-work.dat") 2> $log
done

