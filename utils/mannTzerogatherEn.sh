#! /bin/bash

int=50
dir=$(pwd)


rm $outfile 2> /dev/null

rm *work.dat 2> /dev/null

for tp in lj coul
  do
    for i in `seq 0 $int`;
      do
       if [ -f $(echo $tp"."$i".dat") ];
       	  then
	    awk '{i+=1; if(i>2){print $3,$2}}' $(echo $tp"."$i".dat") >> $(echo $tp"-work.dat")
	fi
     done
done
