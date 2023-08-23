#! /bin/bash

dumpparameter=9
percentraj=0.7


if [ ! -e "$1" ];
  then
    echo "File $1 not found"
    echo
    exit
fi

if [ ! -e "$2" ];
  then
    echo "File $2 not found"
    echo
    exit
fi

echo "First file must be dump"
echo "Second file must be data"
echo
echo

frames=$(awk 'BEGIN{n=0;}{if ($2 == "TIMESTEP"){n+=1;}}END{print n;}' $1)

newfr=$(echo "$frames*$percentraj" | bc )
newframes=${newfr%.*}

awk -v fr=$newframes -v dumpar=$dumpparameter 'BEGIN{n=0;} {
                if ($2 == "TIMESTEP"){n+=1;}
                if ( n == fr ){
                        for(i=1;i<dumpar;i++){getline; if(i==3){natoms=$1;}}
                        for(i=1;i<natoms+1;i++){getline; print $3,$4,$5;}
    }}' $1 > coord.tmp

norm=0
xedge=1.0
yedge=1.0
zedge=1.0

norm=$( awk '{if($2 == "ATOMS"){if($5=="xs"){print 1; exit}}}' $1)

if [ $norm == 1 ];
  then
    xedge=$( awk '{if($2 == "BOX"){getline; xlo=$1; xhi=$2; print xhi-xlo; exit }}' $1)
    yedge=$( awk '{if($2 == "BOX"){getline; getline; ylo=$1; yhi=$2; print yhi-ylo; exit }}' $1)
    zedge=$( awk '{if($2 == "BOX"){getline; getline; getline; zlo=$1; zhi=$2; print zhi-zlo; exit }}' $1)
fi

nline=$(awk 'BEGIN{n=0;}{n+=1; if ($1 == "Atoms"){cnt=n;}}END{print cnt;}' $2)


./dump2data coord.tmp $2 $nline $xedge $yedge $zedge
