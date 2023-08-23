#! /bin/bash

int=50
dir=$(pwd)

direc=$(pwd | rev | cut -d/ -f1 | rev)
orient=$(pwd | rev | cut -d/ -f3 | rev)
temp=$(pwd | rev | cut -d/ -f2 | rev)

outfile=$(echo "gamma_energy_"$direc"_"$orient"_"$temp".dat")

rm $outfile 2> /dev/null


for i in `seq 0 $int`;
  do
   area=$(echo "area_A"$i)
   echo $dir/$area
   cd $dir/$area/step3
   # Devo prendere il logfile ma ha un nome che cambia, quindi lo seleziono come 
   # il file piu' grande .log ma c'e' il rischio che non lo becchi se non e' il piu' grande
   #logfile=$(ls -S | grep "\.log") #


   natoms=$(awk '{if($2 == "atoms"){print $1}}' inputS3.lmp)
   lx=$(awk '{if($3 == "xlo"){print $2-$1}}' inputS3.lmp)
   ly=$(awk '{if($3 == "ylo"){print $2-$1}}' inputS3.lmp)
   
   #echo $natoms $lx $ly

   # Altrimenti uso il log.lammps
   logfile="log.lammps"
 
   # 2 interfaces
   gamma=$(awk -v nat=$natoms -v lx=$lx -v ly=$ly '{if($1 == "run" && $2 == 0){i+=1; getline;  getline; getline;  en[i]=$3;}}END{print (en[2]-en[1])*nat/lx/ly/2.0}' $logfile  ) 

   echo $gamma

   echo "t" | awk -v A=$i -v g=$gamma '{print A, g;}' >> $dir/$outfile

done

cd $dir 
mv $outfile ../../.
