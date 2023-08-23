#! /bin/bash 

wall=$1

if [ ! -e $wall ];
  then
   echo "ERROR: There is no wall file in the dir" $(pwd)
   echo "Exiting..."
   echo -1
fi

cleav=$2

if  (( $(echo "$cleav < 0.0" |bc -l) ));
  then
   echo "ERROR: Second input must be positive (cleav) "
   echo "Exiting..."
   echo -1
fi


echo "PROVA" $cleav

inputfile="Fstep1.0.62.data"
ntypes=2
s3input="inputS3.lmp"

if [ ! -e $inputfile ];
  then
	echo "ERROR: $inputfile not found"
	pwd
	echo "Exiting..."
	exit
fi

echo "PROVA" $cleav
pwd


awk -v cleavplane=$cleav -v ntyp=$ntypes '{
	regex="Atoms";
	where = match($0,regex);
	if(where > 0){
	print $0; getline;print " ";
        for(j=1;j<natoms+1;j++){
                getline; i+=1;
                indxA[i]=$1;
                indxT[i]=$2;
                x[i]=$3;
                y[i]=$4;
                z[i]=$5;
                a[i]=$6;
                b[i]=$7;
                c[i]=$8;
	}

    	for(j=1;j<natoms+1;j++){
	#i contiene il totale di campi letti
    		if(z[j] < cleavplane ){print indxA[j],indxM[j],1,x[j],y[j],z[j],a[j],b[j],c[j]; }
		else {print indxA[j],indxM[j],2,x[j],y[j],z[j],a[j],b[j],c[j];}
    }
}
else if($2 == "atoms") {natoms=$1; print $0;}
else if($3 == "types"){print ntyp,"atom types"}
else if($3 == "xlo"){print $0; print "xbox",$1,$2 > "box.log"}
else if($3 == "ylo"){print $0; print "ybox",$1,$2 > "box.log"}
else{print $0;}
}' $inputfile > $s3input

echo "PROVA3" $cleav

xlo=$(awk '{if($1 == "xbox" ){print $2}}' box.log)
xhi=$(awk '{if($1 == "xbox" ){print $3}}' box.log)
ylo=$(awk '{if($1 == "ybox" ){print $2}}' box.log)
yhi=$(awk '{if($1 == "ybox" ){print $3}}' box.log)

echo "PROVA4" $cleav

cp step3.in step3.in_bck 

awk -v wall=$wall -v s3inp=$s3input -v xlo=$xlo -v xhi=$xhi -v ylo=$ylo -v yhi=$yhi '{
	if($1 == "read_data"){print $1,s3inp}
	else if($2 == "Rfrz1a"){print $1,$2,$3,xlo,xhi,ylo,yhi,$8,$9,$10,$11,$12}
        else if($2 == "Rfrz2a"){print $1,$2,$3,xlo,xhi,ylo,yhi,$8,$9,$10,$11,$12}
	else if($2 == "f2"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall;}
        else {print  $0}
}' step3.in_bck > step3.in

echo "FINE"











