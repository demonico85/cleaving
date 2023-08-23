#! /bin/bash

or=$1
T=$2
ty=$3
i=$4
dir=$5

if [ -z $dir ];
  then
  dir="F"
fi

#or=111
#T="T1"
#ty="wells"
#i=3

step1=1
step3=1
step4=1

currdir=$(pwd)
inpdir="./utils"
nbins=5
replace=0

tags=()

outf=$( echo $dir".f_totW-work.dat")

tags=("f_totW" "v_totW" "f_totWA" "f_totWB" "v_totWc" "v_totWl" "f_LW_c" "f_RW_c" "f_LW_l" "f_RW_l")
#tags=("f_totW" "v_totW" "f_totWA" "f_totWB" "v_totWC" "v_totWL" "f_fc2" "f_fc3" "f_fl2" "f_fl3")


         fold=$(echo "fcc"$or"_"$T"_"$ty"_"$i)
         if [ ! -d $fold ];
           then
               echo "$fold not found" 
               echo "Skipping it"
               exit
         fi
         echo "Gathering from $fold"
         cd $fold
         if [ $step1 -eq 1 ];
          then 
         outf=$(echo $dir".f_totW-work.dat")
         for st in step1 step2  
           do
             cd $currdir/$fold
             cd $st/out
         echo "run work.sh "
             $inpdir/work_s3.sh $dir 
	     if [ ! -f $outf ];
	       then
		       echo "Results non produced in $fold for $st"
		       echo "Chek it better"
	       else
		       cp $outf $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st"_"$dir"_"$i".dat")
	       fi
	   done
           fi
           if [ $step4 -eq 1 ];
             then
           cd $currdir/$fold
           st=step4
           cd $st/out
           if [ -e $(echo "ave."$dir".liq.1.out") ];
             then
#echo "DDDDDDDDDDDDDDDDDDDDD"
               for syst in liq sol
                 do
                    $inpdir/work_s3.sh $dir $syst 
                    wfile=$(echo $syst".F")
                    for string in "${tags[@]}"
                      do
#                        echo "AAAAAAA" $(echo $wfile"."$string-work.dat )
                        if [ $string == "v_totW" ];
                          then
#echo "CCCCCCCCCCCCCC"
                            continue
                        fi
 #echo "BBBBBB" $(echo $wfile"."$string-work.dat )
                        outfs4=$(echo $wfile"."$string-work.dat )
                        if [ -e $outfs4 ];
                          then
                            cp $outfs4 $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st"_"$syst"_"${string:2}"_"$dir"_"$i".dat")
                          else
                            echo "$outfs4 non produced in $fold for $st"
                            echo "Chek it better"
                         fi
                   done

               done
           else
             $inpdir/work_s3.sh $dir 
             outf=$(echo $dir".v_totW-work.dat")
             if [ ! -f $outf ];
               then
                       echo "Results non produced in $fold for $st"
                       echo "Chek it better"
               else
                       cp $outf $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st"_"$dir"_"$i".dat")
               fi
           fi
           fi
          if [ $step3 -eq 1 ];
            then
           outf="F.f_totW-work.dat"
	   cd $currdir/$fold
	   st=step3
             cd $st/dat
	     $inpdir/compiled/s3work # _nodupl
	     echo $inpdir/compiled/s3work
	     if [ ! -x $inpdir/compiled/s3work ];
               then
                echo "ERROR: s3work not found or not executable"
		echo "Exiting..."
		exit -1
             fi
             $inpdir/work_s3.sh $dir 
             if [ ! -f $outf ];
               then
                       echo "Results non produced in $fold for $st"
                       echo "Chek it better"
               else
                       cp $outf $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st"_"$dir"_"$i".dat")
               fi
          fi
