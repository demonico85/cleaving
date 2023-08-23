#! /bin/bash



rm *log 2> /dev/null

# Inputs 

maxblock=500



# Look for data file to extract box dimensions
 
ndatafiles=$(ls *.data 2> /dev/null | wc -l)
re='^[0-9]+$'

if [ "$ndatafiles" -gt 1 ];
  then
    echo "There is more than one data file in the dir"
    echo "The script does not know which one to use"
    echo " leave just the data you need"
    echo "Exiting ... "
    echo
    exit
elif [ "$ndatafiles" -eq 0 ];
  then
    echo "There is no data file in the dir"
    echo "The script needs one data file to run"
    echo "Exiting ... "
    echo
    exit
elif ! [[ "$ndatafiles" =~ $re ]];
  then
    echo "There is no data file in the dir"
    echo "The script needs one data file to run"
    echo "Exiting ... "
    echo
    exit
fi


data=$(ls *.data)
includeortho=0
direction="a"
inputf="stress_profile.out"

# Read Inputs


while getopts on:d:f: option
do
 case "${option}"
 in
 n) maxblock=${OPTARG}
    if ! [[ $maxblock =~ $re ]] ; then
       echo "ERROR: Input for -n not an integer"
       echo "Exiting ..." 
       exit 1
    fi
    ;;
 o) includeortho=1
    ;;

 d) direction=${OPTARG} 
    if [ "$direction" != "x" ] && [ "$direction" != "y" ] && [ "$direction" != "z" ] ; 
      then
        echo "ERROR: the input should be x,y,z"
        echo
        exit 3
    fi
    ;;
 f) inputf=${OPTARG}
    here=0
    if [ ! -e $inputf ];
      then
       let here=here+1
    fi
    if  [ ! -e ../$inputf ];
      then
       let here=here+1
    fi
    if  [ $here -gt 1 ];
      then
        echo "There is no file stress_profile.oute either in this dir nor in parent dir"
        echo " check your inputs"
        echo "Exiting ... "
        echo
        exit 5
    elif [ $here -eq 0 ];
      then
      echo "ERROR: There are two  stress_profile.out files"
      echo "(one in this dir, the other in the parent dir)"
      echo "Choose one to keep"
      echo "Nothing will be done. Exiting..."
      exit 4 
    fi
    if  [ -e ../$inputf ];
      then
       inputf=../$inputf
    fi
    ;;
 \?)
      echo "ERROR: Invalid option: -$OPTARG" 
      echo "Exiting ..."
      exit -1
      ;;
 esac
done


awk 'BEGIN{ortho=1;}{if ($3 == "xlo"){xlo=$1; xhi=$2} else if($3 == "ylo"){ylo=$1; yhi=$2} else if($3 == "zlo"){zlo=$1; zhi=$2} else if ($4 == "xy"){ortho=0;}}END{print "lx",xhi-xlo; print "ly",yhi-ylo; print "zlo", zlo,zhi-zlo;print "ortho", ortho;}' $data > box.log

zlo=$(awk '{if($1 ==  "zlo"){print $2}}' box.log)
lx=$(awk '{if($1 ==  "lx"){print $2}}'   box.log)
ly=$(awk '{if($1 ==  "ly"){print $2}}'   box.log)
lz=$(awk '{if($1 ==  "zlo"){print $3}}' box.log)
ortho=$(awk '{if($1 ==  "ortho"){print $2}}' box.log)

if [ $includeortho == 1 ];
  then
   ortho=0
fi

# In non-orthogonal box lammps requires reduced units for chunk and stress calculation
# that means that I need to multiply the final results of the spatial coordinate by
# lz to recover the real length of the interval 


pdf=$(echo "profile_pdfs")

if [ -d $pdf ];
  then
    rm -r $pdf 2> /dev/null
fi
mkdir $pdf





# Calculation

deltaz=$( awk '{i+=1; if(i == 5){zl=$2}else if(i == 6){zh=$2; exit}}END{print zh-zl;}' $inputf)

if [ $ortho == 0 ];
  then
    deltaz=$(echo 't' | awk -v dzz=$deltaz -v lz=$lz '{print dzz*lz}')
fi

bvol=$( echo "t" | awk -v lx=$lx -v ly=$ly -v dz=$deltaz '{print lx*ly*dz}' )
nbin=$( awk '{i+=1; if(i == 4){print $2; exit}}' $inputf)

cnt=1

while [ $cnt -le $nbin ];
  do
    awk  -v delta=$deltaz -v lobound=$zlo -v nbins=$nbin -v currbin=$cnt '{ cnt+=1;
       if(cnt > 3){
        for(i=0;i<nbins;i++){
            getline;
            if( i == currbin-1) print $0
            }
            
    }
}' $inputf > $(echo $cnt"_bin.tmp")

if [ $direction == "x" ] || [ $direction == "a" ];
  then
    awk -v binvol=$bvol '{print -$4/binvol > "pxx.tmp" }' $(echo $cnt"_bin.tmp")
    nlines=$(wc -l < pxx.tmp)
    if [ "$nlines" -lt "$(echo "t" | awk -v n=$maxblock '{print n*4}' )" ];
      then
        echo "ERROR(x): analysis.py requires at least nobservation*4"
        echo "The files contains fewer of them"
        echo "reduce maxblock (option -n)"
        echo "Exiting..."
        exit
    fi
    mv pxx.tmp $(echo $cnt"_pxx.tmp")
fi

if [ $direction == "y" ] || [ $direction == "a" ];
  then
    awk  -v binvol=$bvol '{print -$5/binvol > "pyy.tmp" }' $(echo $cnt"_bin.tmp")
    nlines=$(wc -l < pyy.tmp)
    if [ "$nlines" -lt "$(echo "t" | awk -v n=$maxblock '{print n*4}' )" ];
      then
        echo "ERROR(y): analysis.py requires at least nobservation*4"
        echo "The files contains fewer of them"
        echo "reduce maxblock (option -n)"
        echo "Exiting..."
        exit
    fi

    mv pyy.tmp $(echo $cnt"_pyy.tmp")
fi

if [ $direction == "z" ] || [ $direction == "a" ];
  then
    awk  -v binvol=$bvol '{print -$6/binvol > "pzz.tmp" }' $(echo $cnt"_bin.tmp")
    nlines=$(wc -l < pzz.tmp)
    if [ "$nlines" -lt "$(echo "t" | awk -v n=$maxblock '{print n*4}' )" ];
      then
        echo "ERROR(z): analysis.py requires at least nobservation*4"
        echo "The files contains fewer of them"
        echo "reduce maxblock (option -n)"
        echo "Exiting..."
        exit
    fi

    mv pzz.tmp $(echo $cnt"_pzz.tmp")
fi  
    let cnt=cnt+1

done

if [ $direction == "x" ];
  then
    dircs=(pxx)
elif [ $direction == "y" ];
  then
    dircs=(pyy)
elif [ $direction == "z" ];
  then
    dircs=(pzz)
elif [ $direction == "a" ];
  then
    dircs=(pxx pyy pzz)
fi


for pdir in ${dircs[@]};
do
    cnt=1
    rm $(echo $pdir".dat") 2> /dev/null
    while [ $cnt -le $nbin ];
      do
        finp=$(echo $cnt"_"$pdir".tmp")
echo $cnt $finp
        fn=$(echo $cnt"_"$pdir".pdf")
        python ~/bin/analysis.py -f $finp -n $maxblock -o $fn  | awk -v n=$cnt -v z=$zlo -v delta=$deltaz -v ortho=$ortho -v lz=$lz  '{if ($1 == "<x>"){print z+(n-1/2)*delta,$3, $5}}' >>  $(echo $pdir".dat")
#        python ~/bin/analysis.py -f $finp -n $maxblock -o $fn  | awk -v n=$cnt -v z=$zlo -v delta=$deltaz -v ortho=$ortho -v lz=$lz  '{if ($1 == "<x>"){if(ortho == 1){print z+(n-1/2)*delta,$3, $5}else if(ortho == 0){print z+lz*(n-1/2)*delta,$3, $5}}}' >>  $(echo $pdir".dat")

        mv $fn ./$pdf/.

        let cnt=cnt+1
    done

done

# Controllare questa parte, non l'ho ancora usata

awk -v binvol=$bvol -v delta=$deltaz -v lobound=$zlo -v ortho=$ortho -v lz=$lz '{ cnt+=1;
       if(cnt > 3){
         nbins=$2;
        for(i=0;i<nbins;i++){
            getline;
            pxx[i]+=$4;
            pyy[i]+=$5;
            pzz[i]+=$6;
            }
        sample+=1;
    }
}END{   

        for(i=0;i<nbins;i++){
            if(ortho == 1){print lobound+(i*delta+(i+1)*delta)/2.0, pxx[i]/binvol/sample, pyy[i]/binvol/sample,pzz[i]/binvol/sample, (pxx[i]+pyy[i]+pzz[i])/3.0/binvol/sample;} 
            else if (ortho == 0){print lobound+lz*(i*delta+(i+1)*delta)/2.0, pxx[i]/binvol/sample, pyy[i]/binvol/sample,pzz[i]/binvol/sample, (pxx[i]+pyy[i]+pzz[i])/3.0/binvol/sample;}
            else {print "ERROR";}
            }
}' $inputf > stressZ.dat

rm *tmp 2> /dev/null


