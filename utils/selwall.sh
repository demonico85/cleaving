#! /bin/bash

#awk '{i+=1; if(i > 11 &&  i < 4764){print $1,1,$3,$4,$5}else{print $0;}}' InputCW.lmp > dataYang.dump

#awk '{i+=1; if(i > 2){work += $4; cnt += 1}}END{print work/cnt}' 

cnt=1

rm work.dat 2> /dev/null
rm work_a.dat 2> /dev/null
rm work_b.dat 2> /dev/null
rm work_tot.dat 2> /dev/null


while [ $cnt -le 1000 ]
  do


    file=$( echo "ave.F."$cnt".out" )
    if [ -e $file ];
      then
        awk '{i+=1; if(i > 2){work += $4; cnt += 1; zw=$5;}}END{print $5,"    ", work/cnt;}' $file >> work.dat
        awk '{i+=1; if(i > 2){work += $6; cnt += 1; zw=$5;}}END{print work/cnt;}' $file >> work_a.dat
        awk '{i+=1; if(i > 2){work += $7; cnt += 1; zw=$5;}}END{print work/cnt;}' $file >> work_b.dat

        paste work.dat work_a.dat work_b.dat > work_tot.dat

    else
        echo "File " $file " not found" 
        exit
    fi
    let cnt=cnt+1

done

