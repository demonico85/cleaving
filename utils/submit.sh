#! /bin/bash

function subm {

sbatch -J $1  ~/bin/inputs/runcleaving.sh 
}


abs=$(pwd)


for name in `cat "dirs"`
 do

    if [ ! -d "$name" ] && [  -z $1 ];
       then
            echo "You did not insert any name but I didn't find the dirs "
            echo "check your inputs, nothing will be done"
            exit
    fi

    if [ -d "$name" ] && [ ! -z $1 ]; 
      then
        echo  $(echo $name"_"$1)
        mv $name $(echo $name"_"$1)
        name=$(echo $name"_"$1)

    fi

   if [ ! -d "$name" ] && [ ! -z $1 ];
     then
        name=$(echo $name"_"$1)
        if [ ! -d "$name" ];
          then
            echo "The name of the dirs is already changed "
            echo "but the new name does not agree with the suffix you wrote."
            echo "check your inputs, nothing will be done"
            exit
        fi
   fi

    if [ -d $name ];
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



        if [ -s sub.log ]; then
            cat sub.log
        fi

        if [ "$running" == "1" ];
          then
            echo $name "is already running"
        elif [ "$running" == "0" ];
          then
           echo $name "job submitted"
           cd $abs/$name
           subm $jobname
          cd ../.
        else 
            echo "ERROR: something went wrong in the reading of the queue. Check it!"
        fi
        rm sub.log 2> /dev/null 
    fi
done

