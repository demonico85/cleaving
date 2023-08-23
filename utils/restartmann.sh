#! /bin/bash

file=$(ls -t ./data/mann.* | head -1)

file=$(echo $file | rev | cut -f1 -d/ | rev)
echo $file

count=$(echo $file | cut -f2 -d.)
delam=$( awk '{if($2 == "delam"){print $4}}' cleaving.in.init )
lambda=$(awk '{if($2 == "lambda"){print $4}}' cleaving.in.init )

#echo "AAA" $count $delam $lambda

newlam=$( echo $count $delam | awk '{print 1.0-$1*$2}')

count=$((count+1))

#echo $startlam $count


cp cleaving.in.init cleaving.in.init_bkp 

awk -v file=$file -v count=$count -v newlam=$newlam '{if ($1 == "read_data"){print $1,"\""file"\""}else if($2 == "cnt" ){print $1,$2,$3,count}else if($2 == "lambda"){print $1,$2,$3,newlam}else {print $0;}}' cleaving.in.init_bkp > cleaving.in.init

cp ./data/$file .


