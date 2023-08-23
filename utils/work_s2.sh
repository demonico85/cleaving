#! /bin/bash

function varcleav ()
{

 file=$2
 awk -v N=$1 '{i+=1; if(i == 5){print $N}}' $file  > var.tmp 
}


function createav ()
{
 
 file=$2
 awk -v N=$1 '{i+=1; if(i > 3){print $N}}' $file > dummy.tmp
 fn="$(echo ${file%.*}.pdf)"
 python $scriptdir/analysis.py -f dummy.tmp -n $maxblock -o $fn  |  awk '{if ($1 == "<x>"){av=$3; err=$5;}else if($1 == "95%"){conf=$7;}} END{print av, err, conf }' > $(echo $1".tmp")

 rm dummy.tmp 2> /dev/null
 rm $fn 2> /dev/null

}


#awk '{i+=1; if(i > 11 &&  i < 4764){print $1,1,$3,$4,$5}else{print $0;}}' InputCW.lmp > dataYang.dump

#awk '{i+=1; if(i > 2){work += $4; cnt += 1}}END{print work/cnt}' 


scriptdir="./utils"
minimumsize=5
maxblock=5
cnt=1

if [ ! "$1" == "B" ] && [ ! "$1" == "b" ] && [ ! "$1" == "F" ] && [ ! "$1" == "f" ];
  then
   echo "Usage work F (or f, B, b)" 
   exit 
fi 

syst=$2

side=${1^}

if [ -z  $syst ];
  then
    wfile=$(echo $side)
else
    wfile=$(echo $syst"."$side)
fi

tags=()


echo $wfile

rm $wfile 2> /dev/null

pdf=$(echo $side"pdfs")

if [ -d $pdf ];
  then
    rm -r $pdf 2> /dev/null
fi
mkdir $pdf


nfiles=$(ls | wc -l)

#nfiles=5


rm *.tmp 2> /dev/null
rm $(echo $syst"*-work.dat")  2> /dev/null

while [ $cnt -le $nfiles ];
  do
    ex1=0
    ex2=0
    if [ -z  $syst ];
      then
        file=$( echo "ave."$side"."$cnt".out" ) 
    else
        file=$( echo "ave."$side"."$syst"."$cnt".out" )
    fi
    echo $file
    forw="${file%*F*}"
    if [[ "$file" == *"$side"* ]];
      then
        if [ -e $file ];
          then
            varcleav 5 $file            
# These numbers depend on the column in the output of LAMMPS
            
            if [ ! -z  $syst ];
               then
               createav 6 $file  
               tags[6]="totWc"
               createav 7 $file
               tags[7]="totWl"
               createav 8 $file  
               tags[8]="LW_c"
               createav 9 $file
               tags[9]="RW_c"
               createav 10 $file 
               tags[10]="LW_l"
               createav 11 $file
               tags[11]="RW_l"
            else
               createav 4 $file
               tags[4]="totW"
            fi
           
            for nmb in 4 6 7 8 9 10 11
              do
                if [ -e $(echo $nmb".tmp") ];
                  then
                 #   echo $nmb $(echo $wfile"."${tags[$nmb]}-work.dat )
                    paste var.tmp $(echo $nmb".tmp") | column -s $'\t' -t >> $(echo $wfile"."${tags[$nmb]}-work.dat )
                fi
            done

        fi
        let cnt=cnt+1    
    fi
done

