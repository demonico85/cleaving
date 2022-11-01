#! /bin/bash


currdir=$(pwd)
inpdir=/mnt/iusers01/pp01/mjkssnd2/scratch/interface/inputs_8Jul2020
insim=2
totsim=2
replace=0


runs3=1

qstat > qstat.log

for or in 111  #100 110 
do
   for T in T1 #T2 T3
   do
     for ty in walls # wells
       do
         cd $currdir
         fold=$(echo "fcc"$or"_"$T"_"$ty)
         for i in `seq $insim $totsim `
           do
		cd $currdir
		newfold=$(echo $fold"_"$i)
                namerun=$( echo "fcc"$or"_"$T"_"$ty"_"$i)

               echo $newfold 
                    if grep -q $namerun $currdir/qstat.log
                      then
# per ora conta anche i casi con C completed
                        echo "Job < $namerun > is already running"
                      else
                        echo "No job called $namerun"
                        echo "Looking for the related dir"
                        if [ ! -d $newfold ];
                          then
                            echo "WARNING. I cannot find $cleavdir"
                            echo "Skipping it"
                        else
                            cd $newfold
                            cp $inpdir/run.sh .
                               lambda=75 # Per ora e' hardcoded
                            step=step1
                               zwi=$( awk '{if($1 == "variable" &&  $2 == "zwi"){print $4}}' ./$step/$(echo $step".in") )
                               zwf=$( awk '{if($1 == "variable" &&  $2 == "zwf"){print $4}}' ./$step/$(echo $step".in") )
                               dew=$( awk '{if($1 == "variable" &&  $2 == "dew"){print $4}}' ./$step/$(echo $step".in") )
                               nfiles=$(echo "t" | awk -v zi=$zwi -v zf=$zwf -v d=$dew '{print int((zi-zf)/d)+1}')
                               cd ./$step
                               if [ ! -e ./out/$(echo "ave.F."$nfiles".out") ];
                                 then
                                   cd $currdir/$newfold
                                   qsub -N $jobname -F $step ./run.sh 
                                   echo "Run $step" 
                               else
				  step=step2
                                  cd $currdir/$newfold/$step
                                  zwi=$(awk '{if($1 == "variable" &&  $2 == "zwi"){print $4}}' $(echo $step".in") )
                                  zwf=$(awk '{if($1 == "variable" &&  $2 == "zwf"){print $4}}' $(echo $step".in") )
                                  dew=$(awk '{if($1 == "variable" &&  $2 == "dew"){print $4}}' $(echo $step".in") )
                                  nfiles=$(echo "t" | awk -v zi=$zwi -v zf=$zwf -v d=$dew '{print int((zi-zf)/d)+1}')
				  if [ ! -e ./out/$(echo "ave.F."$nfiles".out") ];
                                    then
                                      cd $currdir/$newfold
                                   qsub -N $jobname -F $step ./run.sh 
                                      echo "Run $step" 
                                  else
                                    step=step3
                                    cd $currdir/$newfold/$step
                                    if [ ! -e ./dat/$(echo "lambda."$lambda".dat") ] || [ $runs3 -eq 1 ];
                                      then
                                        cd $currdir/$newfold
pwd
if [ -e run.sh ];
  then 
   echo "OK"
fi                                      
	                                qsub -N $jobname run.sh
                                        echo "Run $step" 
                                    else
                                    step=step4
                                    cd $currdir/$newfold/$step
                                      zwi=$(awk '{if($1 == "variable" &&  $2 == "zwi"){print $4}}' $(echo $step".in") )
                                      zwf=$(awk '{if($1 == "variable" &&  $2 == "zwf"){print $4}}' $(echo $step".in") )
                                      dew=$(awk '{if($1 == "variable" &&  $2 == "dew"){print $4}}' $(echo $step".in") )
                                      nfiles=$(echo "echo" | awk -v zi=$zwi -v zf=$zwf -v d=$dew '{print int((zi-zf)/d)+1}')
                                      if [ ! -e ./out/$(echo "ave.F."$nfiles".out") ];
                                        then
                                          cd $currdir/$newfold
                                          qsub -N $jobname -F $step ./run.sh 
                                          echo "Run $step" 
                                      else
					  echo "Case completed"
			              fi
                                   fi
                                 fi
                              fi
                          fi
                    fi
         done
      done
   done
done

