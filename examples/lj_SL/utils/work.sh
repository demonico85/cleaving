#! /bin/bash

cnt=1

if [ ! "$1" == "B" ] && [ ! "$1" == "b" ] && [ ! "$1" == "F" ] && [ ! "$1" == "f" ];
  then
   echo "Usage: work.sh F (or f, B, b)" 
   exit 
fi 

side=${1^}
wfile=$(echo $side-work.dat)

echo $wfile

rm $wfile 2> /dev/null


nfiles=$(ls ave.* | wc -l)

while [ $cnt -le $nfiles ];
  do
    file=$( echo "ave."$side"."$cnt".out" ) 
    echo $file
        if [ -e $file ];
          then

            lam=$(awk '{i+=1; if(i == 5){print $5; exit}}' $file)
            awk -v lam=$lam '{i+=1; if(i > 3){
                    s+=$4; ssq+=$4*$4; n+=1;}
			}END{ssq=ssq/n; s=s/n; print  lam, s, sqrt(ssq - s*s); }' $file >> $wfile 
        fi
        let cnt=cnt+1    
done

