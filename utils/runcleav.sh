
#
#
# NON TORNANO i timestep se ricominci. Non dovrebbe essere un problema
#
#
#

# INPUTS


step="step3"
dir=$(pwd)
cont=0


# Read Inputs

while getopts cs: option
do
 case "${option}"
 in
 c)  cont=1;;
 s) startstep=${OPTARG};;
 esac
done




if [ $cont -eq 0 ];then

rm $dir/ave.F.* 2> /dev/null
rm $dir/ave.B.* 2> /dev/null
rm $dir/force.dump 2> /dev/null

rm -r $dir/dump 2> /dev/null
mkdir $dir/dump

rm -r $dir/group 2> /dev/null
mkdir $dir/group

rm -r $dir/X.restart 2> /dev/null
mkdir $dir/X.restart


rm -r $dir/data 2> /dev/null
mkdir $dir/data

rm -r $dir/dat 2> /dev/null
mkdir $dir/dat

rm -r $dir/out 2> /dev/null
mkdir $dir/out

rm $dir/forces.log 2> /dev/null

elif [ $cont -eq 1 ]; then
    delam=$(awk '{ if($1 == "variable" && $2 == "delam" ){print $4}}' $( echo $step".in") )
    for datafile in `ls ./data/ | sort -r -n -t . -k 2`
      do
        break
    done

    cp ./data/$datafile .

# prendo il numero da data file invece che da out perche' potrebbe essere che out esista senza che quel 
# punto abbia finito. Invece data viene creato quando ha finito quel valore di lambda

    numb=$( echo $datafile  | cut -f2 -d.)
    let numb=numb-1
    lambda=$(echo "$delam * $numb" | bc )
    
    cp $( echo $step".in") $( echo $step"-restart.in")
    
    awk -v lam=$lambda -v fil=$datafile -v nn=$numb 'BEGIN {lb=0;ct=0;}{
        if( $1 == "variable" && $2 == "lambda" && lb == 0){ print $1,$2,$3,lam; lb=1;}
        else if($1 == "read_data"){print $1,fil}
        else if($1 == "variable" && $2 == "cnt" && ct == 0){ print $1,$2,$3,nn+1; ct=1;}
        else {print $0;}
    }' $( echo $step"-restart.in") > $( echo $step".in")

else
	echo "NO option selected"
	echo "Set the enviroment variable for cont"
	echo "sbatch -D `pwd` runcleav.sh --export=cont"

fi 


