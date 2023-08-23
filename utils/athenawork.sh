#! /bin/bash

# Starting of the script



workdir="tot_work"

if [ -d $workdir ];
  then
    rm -r $workdir
fi

    absdir=$(pwd)

mkdir $absdir/$workdir




for name in `ls`
 do
    root=$(echo $name | cut -f1 -d_)
    if [ -d "$name" ] && [ $root == "work" ];
      then
        cd $absdir/$name
        echo $name


        Tmp=$(echo $name | cut -f3 -d_) #Temperature
        typ=$(echo $name | cut -f4 -d_)
        severedname=$(echo $name | cut -f2 -d_  )
        systcr=${severedname:3:3}

        if [ $typ == "walls" ];
          then
             TT="WA"
          else
           TT="WE"
        fi
        namedir=$( echo "fcc"$systcr"_"$Tmp"_"$TT)

        if [ ! -d $absdir/$workdir/$namedir ];
          then
            mkdir $absdir/$workdir/$namedir
        fi

        for step in  step1 step2 step3 step4;
          do
            if [ -d $step ];
              then
                echo $step

                jobname=$( echo "f"$systcr$Tmp$TT)

                cd ./$step

                if [ $step == "step3" ];        
                  then
                    s3work
                fi 

                work.sh F
                cp F-work.dat $absdir/$workdir/$namedir/$(echo $step"_work.dat")
                work.sh B 
                cp B-work.dat $absdir/$workdir/$namedir/$(echo "B_"$step"_work.dat")
        
                cd ../.
            else
                echo "No Dir $name/$step"
            fi
        done
        cd ../.
    fi
done
