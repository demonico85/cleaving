#! /bin/bash

int=50
dir=$(pwd)



for i in `seq 0 $int`;
  do
   area=$(echo "area_A"$i)
   echo $dir/$area
   cd $dir/$area
   if [ ! -d ./step3/dat ];
     then
       echo "Something wrong for $area"
       echo "Dir step3/dat not found"
       echo "Skipping to next area"
   else
      cd ./step3/dat
      rm -r g.dat 2> /dev/null
      echo $area
      for file in ` ls inters3.* | sort -n -t . -k 2 `
        do
#echo $file 
#cat $file 
        if [ $(echo $file | cut -f1 -d. ) == "inters3" ];
          then
            awk '{i+=1; if(i == 6){print $4}}' $file >> g.dat
         fi
       done

       paste ../lambda.txt g.dat | column -s $'\t' -t > gamma.dat

       mv gamma.dat $dir/$(echo $area"_gamma.dat")
   fi
done

cd $dir

direc=$(pwd | rev | cut -d/ -f1 | rev)
orient=$(pwd | rev | cut -d/ -f3 | rev)
temp=$(pwd | rev | cut -d/ -f2 | rev)

out=$(echo "gamma_"$direc"_"$orient"_"$temp)

if [ -d $out ];
  then
    rm -r $out
fi

mkdir $out

mv *dat ./$out/.


