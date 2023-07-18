#! /bin/bash

#awk '{i+=1; if(i > 11 &&  i < 4764){print $1,1,$3,$4,$5}else{print $0;}}' InputCW.lmp > dataYang.dump

#awk '{i+=1; if(i > 2){work += $4; cnt += 1}}END{print work/cnt}' 


scriptdir="/mnt/iusers01/pp01/mjkssnd2/scratch/sandpit/utils"
minimumsize=5
maxblock=5
cnt=1
re='^[0-9]+$'



if [ ! "$1" == "B" ] && [ ! "$1" == "b" ] && [ ! "$1" == "F" ] && [ ! "$1" == "f" ];
  then
   echo "Usage work F (or f, B, b)" 
   exit 
fi 


pdf=$(echo $side"pdfs")

if [ -d $pdf ];
  then
    rm -r $pdf 2> /dev/null
fi
mkdir $pdf

nfiles=$(ls *out | wc -l | awk '{print $1/2}')

echo $(ls *out  | wc -l) $nfiles

 if ! [[ $nfiles =~ $re ]] ; then
   echo "ERROR: number of files in out is not an even number, check"
   echo "Exiting ..."
   exit 1
fi


for phase in sol liq
  do
    cnt=1
    side=${1^}
    wfile=$(echo $side"-work_"$phase".dat")

    echo $wfile

    rm $wfile 2> /dev/null



    while [ $cnt -le $nfiles ];
      do
        ex1=0
        ex2=0
        file=$( echo "ave."$side"."$phase"."$cnt".out" ) 
        forw="${file%*F*}"
        if [[ "$file" == *"$side"* ]];
          then
            if [ -e $file ];
              then
                if [ $phase == "sol" ];
                  then
                     awk '{i+=1; if(i > 3){print $6}}' $file > 1.tmp # Queste dipendono dall'output di Lammps
                     awk '{i+=1; if(i == 5){print $5}}' $file > 2.tmp
                     awk '{i+=1; if(i > 3){print $8}}' $file > a.tmp
                     awk '{i+=1; if(i > 3){print $9}}' $file > b.tmp
                  elif [ $phase == "liq" ];
                    then
                     awk '{i+=1; if(i > 3){print $7}}' $file > 1.tmp # Queste dipendono dall'output di Lammps
                     awk '{i+=1; if(i == 5){print $5}}' $file > 2.tmp
                     awk '{i+=1; if(i > 3){print $10}}' $file > a.tmp
                     awk '{i+=1; if(i > 3){print $11}}' $file > b.tmp

                   fi
                     fn="$(echo ${file%.*}.pdf)"
                     echo  $fn
                     python $scriptdir/analysis.py -f 1.tmp -n $maxblock -o $fn  |  awk '{if ($1 == "<x>"){av=$3; err=$5;}else if($1 == "95%"){conf=$7;}} END{print av, err, conf }' > 3.tmp 
            
                     grep . a.tmp > 4.tmp
                     grep . b.tmp > 5.tmp

                     actualsize=$(wc -c <"4.tmp")

                     if [ $actualsize -ge $minimumsize ];
                       then
                         fn="$(echo ${file%.*}.f2.pdf)"
                         python $scriptdir/analysis.py -f 4.tmp -n $maxblock -o $fn  | awk '{if ($1 == "<x>"){av=$3; err=$5;}else if($1 == "95%"){conf=$7;}} END{print av, err, conf }' > 6.tmp 
                         ex1=1
                     fi 

                     actualsize=$(wc -c <"5.tmp")

                     if [ $actualsize -ge $minimumsize ];
                       then
                         fn="$(echo ${file%.*}.f3.pdf)"
                         python $scriptdir/analysis.py -f 5.tmp -n $maxblock -o $fn  | awk '{if ($1 == "<x>"){av=$3; err=$5;}else if($1 == "95%"){conf=$7;}} END{print av, err, conf }' > 7.tmp
                         ex2=1
                     fi 
    
                     if [ $ex1 -ne $ex2 ];
                       then
                         echo "WARNING: f_2 or f_3 has zero size."
                         echo "something could be wrong, the calculation"
                         echo "won't be stopped but check"
                     fi 

                     if [ $actualsize -ge $minimumsize ];
                       then
                         paste  2.tmp 3.tmp 6.tmp 7.tmp | column -s $'\t' -t >> $wfile
                     else
                         paste  2.tmp 3.tmp  | column -s $'\t' -t >> $wfile
                     fi
         
                     rm *tmp 2> /dev/null 
                     mv $fn ./$pdf/.
                 fi
                 let cnt=cnt+1    
             fi
      done
done
