#! /bin/bash

function varcleav ()
{

if [ $2 == "step3" ];
  then
   findw="v_lambda"
else
  if [ $1 == "walls" ];
    then
      findw="v_zwalls"
  else
      findw="v_lambda"
  fi
fi

#echo $findw
 file=$3
 N=$(awk -v F=$findw '{i+=1; if(i == 2){for(i=1; i<=NF; i++) {if($i==F){ print i-1; exit}}}}' $file) # -1 because the first field in the lammps file is the comment symbol #
 awk -v N=$N '{i+=1; if(i == 5){print $N}}' $file  > var.tmp 
#echo $N
# cat var.tmp
}


function createav ()
{
 
 file=$2
# echo $1
 N=$(awk -v F=$1 '{i+=1; if(i == 2){for(i=1; i<=NF; i++) {if($i==F){ print i-1; exit}}}}' $file) # -1 because the first field in the lammps file is the comment symbol #
 if [ -z $N ];
  then
   return $?
 fi
 awk -v N=$N '{i+=1; if(i > 3){print $N}}' $file > dummy.tmp
 fn="$(echo ${file%.*}.pdf)"
 python $scriptdir/analysis.py -f dummy.tmp -n $maxblock -o $fn  |  awk '{if ($1 == "<x>"){av=$3; err=$5;}else if($1 == "95%"){conf=$7;}} END{print av, err, conf }' > $(echo $1".tmp")

 rm dummy.tmp 2> /dev/null
 rm $fn 2> /dev/null

}

function fast_createav ()
{

 file=$2
# echo $1
 N=$(awk -v F=$1 '{i+=1; if(i == 2){for(i=1; i<=NF; i++) {if($i==F){ print i-1; exit}}}}' $file) # -1 because the first field in the lammps file is the comment symbol #
 if [ -z $N ];
  then
   return $?
 fi
 awk -v N=$N '{i+=1; if(i > 3){print $N}}' $file > dummy.tmp

 awk '{n+=1; sq+=$1*$1; s+=$1;}END{Av=s/n; stdv=sqrt(sq/n - Av*Av); print Av, stdv, stdv}' dummy.tmp >  $(echo $1".tmp") 

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

extp=$(pwd | rev | cut -d/ -f3 | cut -d_ -f2 | rev )
step=$(pwd | rev | cut -d/ -f2 | rev)


syst=$2

side=${1^}

if [ -z  $syst ];
  then
    wfile=$(echo $side)
else
    wfile=$(echo $syst"."$side)
fi



echo $wfile

rm $wfile 2> /dev/null


nfiles=$(ls | wc -l)

#nfiles=5


tags=("f_totW" "v_totW" "f_totWA" "f_totWB" "v_totWc" "v_totWl" "f_LW_c" "f_RW_c" "f_LW_l" "f_RW_l")
#tags=("f_totW" "v_totW" "f_totWA" "f_totWB" "v_totWC" "v_totWL" "f_fc2" "f_fc3" "f_fl2" "f_fl3")

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
            varcleav $extp $step $file            
# These numbers depend on the column in the output of LAMMPS
            
            for i in "${tags[@]}"
              do
             #   createav $i $file
                fast_createav  $i $file
                if [ -e $(echo $i".tmp") ];
                  then
                 #   echo $nmb $(echo $wfile"."$i-work.dat )
                    paste var.tmp $(echo $i".tmp") | column -s $'\t' -t >> $(echo $wfile"."$i-work.dat )
                fi
            done

        fi
        let cnt=cnt+1    
    fi
done

