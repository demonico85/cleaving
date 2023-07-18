#! /bin/bash 

wall=$1

if [ ! -e "$wall" ];
  then
   echo "ERROR: There is no wall file in the dir" $(pwd)
   echo "Exiting..."
   echo -1
fi


cnt=$2

if ! [[ "$cnt" =~ $re ]] ; 
  then
   echo "ERROR: Second input must be an integer"
   echo "Exiting..."
   echo -1
fi

cleav=$3

if  (( $(echo "$cleav < 0.0" |bc -l) )); 
  then
   echo "ERROR: Third input must be positive "
   echo "Exiting..."
   echo -1
fi

temp=$4

if  (( $(echo "$temp < 0.0" |bc -l) )); 
  then
   echo "ERROR: Fourth input must be positive (temperature) "
   echo "Exiting..."
   echo -1
fi



inputfile="step3.50.data"
ntypes=3
zshift=11.5
vapshift=1.5
inputs4="inputS4.lmp"


awk -v cleavplane=$cleav -v ntyp=$ntypes -v shift=$zshift '{
        regex="Atoms";
        where = match($0,regex);
        if(where > 0){
        print $0; getline;print " ";
        NF=10;
        while(NF>0){
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
                if(indxT[j] == 1){if(z[j] < 0.47){typA=3; n+=1; print z[j] > "n.tmp";}else{typA=1;} print indxA[j],typA,x[j],y[j],z[j]-shift; }
		else{ if(z[j] > 33.1){typA=3;n+=1;  print indxA[j],z[j] > "n.tmp";}else{typA=2;} print indxA[j],typA,x[j],y[j],z[j]+shift; }
    }
    print cnt > "n.tmp"; 
    print "  ";
}
else if($1 == "Masses"){
        getline;
	NF=10;
	while(NF>0){
                getline;
	}
}
else if($2 == "atoms") {natoms=$1; print $0;}
else if($3 == "zlo"){print $1-shift,$2+shift,"zlo zhi";}
else if($3 == "types"){print ntyp,"atom types"}
else{print $0;}
} END{print "ciao", n > "n.tmp"}' $inputfile > $inputs4 

cp step4.in step4.in_bck

awk -v wall=$wall -v s4inp=$inputs4 -v zshift=$zshift -v cleav=$cleav -v vshift=$vapshift -v t=$temp '{
	if($1 == "read_data"){print $1,s4inp}
        else if($2 == "Tsyst"){print $1, $2, $3,t}
        else if($2 == "clwall1"){print $1,$2,$3, cleav-zshift }
        else if($2 == "clwall2"){print $1,$2,$3, cleav+zshift }
        else if($2 == "vapwall1"){print $1,$2,$3, cleav-vshift }
        else if($2 == "vapwall2"){print $1,$2,$3, cleav+vshift }
	else if($2 == "fv1"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall}
        else if($2 == "fv2"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall}
        else if($1 == "fix" && $2 == "fc1"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall}
        else if($1 == "fix" && $2 == "fc2"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall}
	else {print $0}
}' step4.in_bck > step4.in


cp in.loop/loop in.loop/loop_bck

awk -v wall=$wall '{
        if($1 == "fix" && $2 == "fc1"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall}
        else if($1 == "fix" && $2 == "fc2"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,wall}
	else {print $0;}
}'  in.loop/loop_bck > in.loop/loop 


















