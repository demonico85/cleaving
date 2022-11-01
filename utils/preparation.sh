#! /bin/bash

function changeinwe12 {

#    echo $1 $2 $3 $4 $5
    cp ./$1/$5 ./$1/$(echo $5"-back") 
#    echo ./$1/$(echo $5"-back") 
    awk -v T=$2 -v S=$3 -v W=$4 -v rho=$6 -v fact=$7 -v dw=$8  -v bloop=$9 '{
        if($2 == "Tsyst"){print $1, $2, $3, T;}
        else if($2 == "rhoc" ){ print $1, $2, $3, rho;}
        else if($2 == "fact" ){ print $1, $2, $3, fact;}
        else if($2 == "dw" ){ print $1, $2, $3, dw;}
        else if($1 == "read_data"){print $1,S;}
        else if($1 == "fix" && $2 == "f2"){print $1, $2,$3,$4,$5,$6,$7,$8,$9,W;}
        else if($1 == "variable" &&  $2 == "backloop"){print $1,$2,$3,bloop}
        else{print $0;}
}' ./$1/$(echo $5"-back")  > ./$1/$5
}

function changeinwa12 {

#    echo $1 $2 $3 $4 $5
    cp ./$1/$5 ./$1/$(echo $5"-back") 
#    echo ./$1/$(echo $5"-back") 
    awk -v T=$2 -v S=$3 -v W=$4 -v zwi=$6 -v zwf=$7 -v dew=$8 -v bloop=$9 '{
        if($2 == "Tsyst"){print $1, $2, $3, T;}
        else if($1 == "read_data"){print $1,S;}
        else if($1 == "fix" && $2 == "f2"){print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,W;}
        else if($1 == "variable" &&  $2 == "zwi"){print $1,$2,$3,zwi}
        else if($1 == "variable" &&  $2 == "zwf"){print $1,$2,$3,zwf}
        else if($1 == "variable" &&  $2 == "dew"){print $1,$2,$3,dew}
        else if($1 == "variable" &&  $2 == "backloop"){print $1,$2,$3,bloop}
        else{print $0;}
}' ./$1/$(echo $5"-back")  > ./$1/$5

}


function changeinwe34 {

#    echo $1 $2 $3 $4 
    cp ./$1/$4 ./$1/$(echo $4"-back") 
#    echo ./$1/$(echo $4"-back") 
    awk -v T=$2 -v W=$3 -v fix=$5 -v rho=$6 -v fact=$7 -v dw=$8  -v bloop=$9 '{
        if($2 == "Tsyst"){print $1, $2, $3, T;}
        else if($2 == "rhoc" ){ print $1, $2, $3, rho;}
        else if($2 == "fact" ){ print $1, $2, $3, fact;}
        else if($2 == "dw" ){ print $1, $2, $3, dw;}
        else if($1 == "fix" && $2 == fix){print $1, $2,$3,$4,$5,$6,$7,$8,$9,W;}
        else if($1 == "variable" &&  $2 == "backloop"){print $1,$2,$3,bloop}
        else{print $0;}
    }' ./$1/$(echo $4"-back")  > ./$1/$4
}

function changeinwa34 {

#    echo $1 $2 $3 $4
    cp ./$1/$4 ./$1/$(echo $4"-back") 
#    echo ./$1/$(echo $4"-back") 
    awk -v T=$2 -v W=$3 -v fix=$5 -v zwi=$6 -v zwf=$7 -v dew=$8  -v bloop=$9 '{
        if($2 == "Tsyst"){print $1, $2, $3, T;}
        else if($1 == "fix" && $2 == fix){print $1, $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,W;}
        else if($1 == "variable" &&  $2 == "zwi"){print $1,$2,$3,zwi}
        else if($1 == "variable" &&  $2 == "zwf"){print $1,$2,$3,zwf}
        else if($1 == "variable" &&  $2 == "dew"){print $1,$2,$3,dew}
        else if($1 == "variable" &&  $2 == "backloop"){print $1,$2,$3,bloop}
        else{print $0;}
    }' ./$1/$(echo $4"-back")  > ./$1/$4
}

################################################################################

dirscripts=/mnt/iusers01/pp01/mjkssnd2/scratch/inputs_8Jul2020

if [ ! "$1" == "111" ] && [ ! "$1" == "110" ] && [ ! "$1" == "100" ];
then
    echo "1usage: preparation.sh <system=111,100,110> <walls/wells> <temp> <cleaving pos>"
    exit 1
fi
 
syst=$1

z=${2^^}

if [ ! "$z" == "WE" ] && [ ! "$z" == "WA" ];
then
    echo "2usage: preparation.sh <system=111,100,110> <walls/wells> <temp> <cleaving pos>"
    exit 1
fi
 
type=$z

re='^[0-9]+([.][0-9]+)?$'
#echo $3
if ! [[ "$3" =~ $re ]] ; then
    echo "error: Not a number"
    echo "3usage: preparation.sh <system=111,100,110> <temp> <cleaving pos>" >&2; exit 1
fi

temp=$3

re='^[0-9]+([.][0-9]+)?$'
#echo $3
if ! [[ "$4" =~ $re ]] ; then
    echo "error: Not a number"
    echo "3usage: preparation.sh <system=111,100,110> <temp> <cleaving pos> " >&2; exit 1
fi

cleav=$4
TT=$5
rho=$6
zwi=$7
zwf=$8
dew=$9
fact=${10}
dw=${11}
a=${12}
b=${13}
c=${14}
backloop=${15}
dew4=${16}
#pbcpressure=${16}

dir="../systems"

if [ ! -d $dir ];
  then
    echo "error: I did not found the systems dir"
    exit 1
fi

# file in step 1/2

for step in step1 step2
  do
    if [ $step == "step1" ];
      then 
        names=$(echo "fcc"$syst"-"$TT".lmp")
    else
        names=$(echo "inputLIQ-"$syst"-"$TT".lmp")        
    fi
    cp $dir/$names ./$step/$names

    if [ "$type" == "WE" ];
      then
        namew=$(echo "fcc"$syst"-"$TT"-wells.lmp")
        cp $dir/$namew ./$step/.
        changeinwe12 "$step" $temp $names $namew "$step.in" $rho $fact $dw $backloop
        changeinwe12 "$step/in.loop" $temp $names $namew "loop" -100  -100 -100 -100 
        changeinwe12 "$step/in.loop" $temp $names $namew "loop.back" -100 -100 -100 -100
    else
        namew=$(echo "fcc"$syst"-"$TT"-walls.lmp")
        cp $dir/$namew ./$step/.
        changeinwa12 "$step" $temp $names $namew $(echo $step".in") $zwi $zwf $dew $backloop
        changeinwa12 "$step/in.loop" $temp $names $namew "loop"        -100 -100 -100
        changeinwa12 "$step/in.loop" $temp $names $namew "loop.back" -100 -100 -100
    fi 
done


# file in step 3/4

for step in step3 step4
  do
    if [ "$type" == "WE" ];
      then
        namew=$(echo "fcc"$syst"-"$TT"-hwells.lmp")
        cp $dir/$namew ./$step/.
        changeinwe34 "$step" $temp $namew "$step.in" "f3" $rho  $fact $dw $backloop
        changeinwe34 "$step/in.loop" $temp $namew "loop" "f3" -100 -100 -100 -100

        namew=$(echo "fcc"$syst"-"$TT"-wells.lmp")
        cp $dir/$namew ./$step/.
        changeinwe34 "$step" $temp $namew "$step.in" "f2" $rho  $fact $dw $backloop
        changeinwe34 "$step/in.loop" $temp $namew "loop" "f2" -100 -100 -100 -100
    else

        namew=$(echo "fcc"$syst"-"$TT"-walls.lmp")
        cp $dir/$namew ./$step/.
        changeinwa34 "$step" $temp $namew "$step.in" "f3"  $zwi $zwf $dew4 $backloop 
        changeinwa34 "$step/in.loop" $temp $namew "loop" "f3" -100 -100 -100 -100

        namew=$(echo "fcc"$syst"-"$TT"-walls.lmp")
        cp $dir/$namew ./$step/.
        changeinwa34 "$step" $temp $namew "$step.in" "f2" $zwi $zwf $dew4 $backloop
        changeinwa34 "$step/in.loop" $temp $namew "loop" "f2" -100 -100 -100 -100
    fi 

done


if [[ -z $a ]];
  then 
    $dirscripts/fxlyrs $syst $rho
else
    if [ $a -gt 0 ]; 
      then
        $dirscripts/fxlyrs $syst $rho -c $a $b $c
    else
        $dirscripts/fxlyrs $syst $rho
    fi
fi


#zhi1=$(awk '{if($1 == "zhi1"){print $2;}}' fxdlyr.tmp)


#if [ "$pbcpressure" == 1 ];
#  then
#   mv fxdlyr.tmp zeroT_fxdlyr.tmp
#   awk -v zl=$zhi1 'BEGIN{zl=zl+0.05*zl; }
#	{if($3 == "xlo"){print "xlo",$1;print "xhi",$2} 
#	  else if ($3 == "ylo"){print "ylo",$1;print "yhi",$2}
#	  else if ($3 == "zlo"){
#			zedge=$2-$1; 
#			print "zlo1",$1; print "zhi1",zl;
#			print "zlo2",$2-zl; print "zhi2",$2;
#		        print "zlo1a",$1+zedge; print "zhi1a",zl+zedge;	
#			print "zlo2a",$2+zedge-zl; print "zhi2a",$2+zedge
#		}
#	  }' ./step1/$(echo "fcc"$syst"-"$TT".lmp") > fxdlyr.tmp
#fi



xlo=$(awk '{if($1 == "xlo"){print $2;}}' fxdlyr.tmp)
xhi=$(awk '{if($1 == "xhi"){print $2;}}' fxdlyr.tmp)
ylo=$(awk '{if($1 == "ylo"){print $2;}}' fxdlyr.tmp)
yhi=$(awk '{if($1 == "yhi"){print $2;}}' fxdlyr.tmp)
zlo1=$(awk '{if($1 == "zlo1"){print $2;}}' fxdlyr.tmp)
zhi1=$(awk '{if($1 == "zhi1"){print $2;}}' fxdlyr.tmp)
zlo2=$(awk '{if($1 == "zlo2"){print $2;}}' fxdlyr.tmp)
zhi2=$(awk '{if($1 == "zhi2"){print $2;}}' fxdlyr.tmp)

zlo1a=$(awk '{if($1 == "zlo1a"){print $2;}}' fxdlyr.tmp)
zhi1a=$(awk '{if($1 == "zhi1a"){print $2;}}' fxdlyr.tmp)
zlo2a=$(awk '{if($1 == "zlo2a"){print $2;}}' fxdlyr.tmp)
zhi2a=$(awk '{if($1 == "zhi2a"){print $2;}}' fxdlyr.tmp)



for step in step1 step3 step4
  do
    inp=$(echo $step.in)
    cp ./$step/$inp ./$step/$(echo $inp"-back") 
    awk -v xlo=$xlo -v xhi=$xhi  -v ylo=$ylo -v yhi=$yhi -v zlo1=$zlo1 -v zhi1=$zhi1 -v zlo2=$zlo2 -v zhi2=$zhi2 -v zlo1a=$zlo1a -v zhi1a=$zhi1a -v zlo2a=$zlo2a -v zhi2a=$zhi2a '{
        if($1 == "region" && $2 == "Rfrz1"){print $1,$2,$3,xlo,xhi,ylo,yhi,zlo1,zhi1,0.0,0.0,0.0;}
        else if($1 == "region" && $2 == "Rfrz2"){print $1,$2,$3,xlo,xhi,ylo,yhi,zlo2,zhi2,0.0,0.0,0.0;}
        else if($1 == "region" && $2 == "Rfrz1a"){print $1,$2,$3,xlo,xhi,ylo,yhi,zlo1,zhi1,0.0,0.0,0.0;}
        else if($1 == "region" && $2 == "Rfrz2a"){print $1,$2,$3,xlo,xhi,ylo,yhi,zlo2,zhi2,0.0,0.0,0.0;}
        else if($1 == "region" && $2 == "Rfrz1b"){print $1,$2,$3,xlo,xhi,ylo,yhi,zlo1a,zhi1a,0.0,0.0,0.0;}
        else if($1 == "region" && $2 == "Rfrz2b"){print $1,$2,$3,xlo,xhi,ylo,yhi,zlo2a,zhi2a,0.0,0.0,0.0;}
else{print $0}}' ./$step/$(echo $inp"-back") > ./$step/$inp
done

# Per i walls la posizione del cleaving plane deve essere specificata, per i wells e' contenuta 
# nei file wells.lmp

if [ "$z" == "WA" ];
  then
    for step in step1 step2 step3 step4
      do
        inp=$(echo $step.in)
        cp ./$step/$inp ./$step/$(echo $inp"-back") 
        awk -v C=$cleav -v z=$zwf '{if($2 == "clwall"){print $1,$2,$3,C;}else if($2 == "clwall1"){print $1,$2,$3,C;} else{print $0}}' ./$step/$(echo $inp"-back") > ./$step/$inp
    done
fi


