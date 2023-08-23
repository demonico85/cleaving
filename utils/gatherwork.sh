#! /bin/bash

# Function definition

function collectwork {

    absdir=$(pwd)

    cd ./$name

    for step in  step1 step2 step3 step4;
      do
        if [ -d $step ];
          then
            cd ./$step
            if [ $step != "step3" ];
              then
                dir=./out
                if [ -d $dir ];
                  then
                    cp -r $dir ../../$(echo "work_"$name)/$(echo $step)
                else
                    echo "Dir $dir not found moving to the next step"
                fi
            else
                dir=./dat
                if [ -d $dir ];
                  then
                    cp -r $dir ../../$(echo "work_"$name)/$(echo $step)
                else
                    echo "Dir $step/$dir not found moving to the next step"
                fi
            fi
            cd ../.
        else
            echo "Dir $step not found moving to the next"
        fi
    done 

    cd $absdir 

}



#*******************************************************************************

# Starting of the script

for name in `ls`
 do

    if [ ${#name} -gt 2 ];
      then
        root=${name:0:3}
    else
        root="wqptr" #lettere a caso
    fi

    if [ -d "$name" ] && [ $root == "fcc" ];
       then
        severedname=${name:0:8}
        systcr=${name:3:3}
        Tmp=$(echo $name | cut -f2 -d_)
        typ=$(echo $name | cut -f3 -d_)
        if [ $typ == "walls" ];
          then
            TT="WA"
        else
            TT="WE"
        fi
        jobname=$( echo "f"$systcr$Tmp$TT)

        running=$(squeue | grep ndp8 | awk -v  nm=$jobname 'BEGIN{a=0}{if($3 == nm){if($5 == "R" || $5 == "PD" ){a+=1;}}} END {if(a == 0) {print 0;} else if(a == 1){print 1;}else {print "ERROR more than two occurrence submission" >"sub.log"}}')

        if [ "$running" == "1" ];
          then
            echo $name "is running, collection work will be incomplete"
        else
            echo $name "is not running"
        fi

        if [ -d "$(echo "work_"$name)" ];
          then
            rm -r  $(echo "work_"$name)
        fi
        mkdir $(echo "work_"$name)
        collectwork $name

    fi
done
