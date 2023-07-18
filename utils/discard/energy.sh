#! /bin/bash

currdir=$(pwd)
temp=$(pwd | rev | cut -f2 -d/ | rev)
#athenadir=$(echo "/gpfs/home/leics/ndp8/mannitol/001/gromos/"$temp)

zwanzigdir=/home/cg/Desktop/nicodemo/mannitol/$temp

echo "I am taking files from the following dir"
echo $athenadir
echo

rm *dat 2> /dev/null

for syst in slab bulk
  do

#     scp -r  cg@delmat.ddns.warwick.ac.uk:$zwanzigdir/$syst/$(echo "ave."$syst"_pe.out") $currdir/.

#    scp -r ndp8@athena.hpc-midlands-plus.ac.uk:$athenadir/$syst/$(echo "ave."$syst"_pe.out") $currdir/.

    awk -v S=$syst '{i+=1; if(i>2){en+=$3; n+=1; }}END{print S,": ", en/n}' $(echo "ave."$syst"_pe.out")  >> a.tmp 



done

type=bulk
#scp -r  cg@delmat.ddns.warwick.ac.uk:$zwanzigdir/$syst/$(echo "nvt_"$type"_mannitol.data") $currdir/.
#scp -r ndp8@athena.hpc-midlands-plus.ac.uk:$athenadir/$syst/$(echo "nvt_"$type"_mannitol.data") $currdir/.

lx=$(awk '{if($3 == "xlo"){print $2-$1}}' $currdir/$(echo "nvt_"$type"_mannitol.data") )

ly=$(awk '{if($3 == "ylo"){print $2-$1}}' $currdir/$(echo "nvt_"$type"_mannitol.data") )

awk -v lx=$lx -v ly=$ly 'BEGIN{
	kcal2mJ=4.184e6;
	NA=6.023e23;
	Ang2m=1e-10;
}
{print $0; if($1 == "bulk"){b=$3}else if($1 == "slab"){s=$3}}END{SFE=0.5*(s-b)/lx/ly; print " "; print "SFE", ":",SFE*kcal2mJ/NA/Ang2m/Ang2m;}' a.tmp > $(echo "energy."$temp".dat")

rm *tmp 2> /dev/null
