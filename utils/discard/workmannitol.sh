#! /bin/bash

#awk '{i+=1; if(i > 11 &&  i < 4764){print $1,1,$3,$4,$5}else{print $0;}}' InputCW.lmp > dataYang.dump

#awk '{i+=1; if(i > 2){work += $4; cnt += 1}}END{print work/cnt}' 

minimumsize=5
maxblock=2
#natoms=19968




pdf=$(echo $side"pdfs")

if [ -d $pdf ];
  then
    rm -r $pdf 2> /dev/null
fi
mkdir $pdf

nfiles=$(ls | wc -l)

for forcetype in coul lj
do

    wfile=$(echo $forcetype"-work.dat")

    echo $wfile
    rm *tmp 2> /dev/null 
    rm $wfile 2> /dev/null

    cnt=1

    while [ $cnt -le $nfiles ];
      do
        ex1=0
        ex2=0
        file=$( echo $forcetype"."$cnt".dat" ) 
        if [[ "$file" == *"$forcetype"* ]];
          then
            if [ -e $file ];
              then
                awk '{i+=1; if(i > 2){print $2}}' $file > 1.tmp # Queste dipendono dall'output di Lammps
                awk '{i+=1; if(i == 5){print $3}}' $file > 2.tmp
    
                fn="$(echo ${file%.*}.pdf)"
                echo  $fn
                python ~/bin/analysis.py -f 1.tmp -n $maxblock -o $fn  | awk '{if ($1 == "<x>"){av=$3; err=$5;}else if($1 == "95%"){conf=$7;}} END{print av, err, conf }' > 3.tmp
                
    
                paste  2.tmp 3.tmp  | column -s $'\t' -t >> 4.tmp 
       
    
    
                mv $fn ./$pdf/.
            fi
            let cnt=cnt+1    
        fi
    done

cp 4.tmp $wfile

#            awk -v NAT=$natoms 'BEGIN{Ang2m=1e-20; NAv=6.023e23; kcal2mJ=4.184e6; conversion=1.}{ print $1,$2*conversion ,$3*conversion}' 4.tmp > $wfile
#            rm *tmp 2> /dev/null 

done



rm *.tmp 2> /dev/null






#conversion=kcal2mJ/NAv/Ang2m/NAT























