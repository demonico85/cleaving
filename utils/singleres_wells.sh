#! /bin/bash




step1=1
step3=1
step4=1

insim=1
totsim=5
currdir=$(pwd)
inpdir="./utils"
nbins=5
replace=0

tags=()



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
         for st in step1 step2 
           do
             cd $currdir/$fold
             cd $st/out
         echo "run work.sh "
#             $inpdir/work.sh F 
$inpdir/work_temp.sh F 
	     if [ ! -f F-work.dat ];
	       then
		       echo "Results non produced in $fold for $st"
		       echo "Chek it better"
	       else
		       cp F-work.dat $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st"_"$i".dat")
	       fi
	   done
           fi
           if [ $step4 -eq 1 ];
             then
           cd $currdir/$fold
           st=step4
           cd $st/out
           if [ -e "ave.F.liq.1.out" ];
             then
           for syst in liq sol
           do
             $inpdir/work_s2.sh F $syst 
               tags[1]="totWc"
               tags[2]="totWl"
               tags[3]="LW_c"
               tags[4]="RW_c"
               tags[5]="LW_l"
               tags[6]="RW_l"
               wfile=$(echo $syst".F")
               for nmb in 1 2 3 4 5 6
                 do
                  outfs4=$(echo $wfile"."${tags[$nmb]}-work.dat )
                  if [ -e $outfs4 ];
                    then
                      cp $outfs4 $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st"_"$syst"_"${tags[$nmb]}"_"$i".dat")
                  else
                       echo "$outfs4 non produced in $fold for $st"
                       echo "Chek it better"
                  fi
              done

           done
           else
             $inpdir/work.sh F 
             if [ ! -f F-work.dat ];
               then
                       echo "Results non produced in $fold for $st"
                       echo "Chek it better"
               else
                       cp F-work.dat $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st"_"$i".dat")
               fi
           fi
           fi
          if [ $step3 -eq 1 ];
            then

	   cd $currdir/$fold
	   st=step3
             cd $st/dat
	     $inpdir/compiled/s3work_fin # _nodupl
	     echo $inpdir/compiled/s3work
	     if [ ! -x $inpdir/compiled/s3work ];
               then
                echo "ERROR: s3work not found or not executable"
		echo "Exiting..."
		exit -1
             fi
             $inpdir/work.sh F 
             if [ ! -f F-work.dat ];
               then
                       echo "Results non produced in $fold for $st"
                       echo "Chek it better"
               else
                       cp F-work.dat $currdir/results/$(echo "fcc"$or"_"$T"_"$ty"_"$st"_"$i".dat")
               fi
          fi
