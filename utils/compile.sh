#! /bin/bash

dir=compiled

rm -r $dir 2> /dev/null

if [ ! -d $dir ];
  then
    mkdir $dir
fi


gfortran fixedlayer.f90 -o ./$dir/fxlyrs
gfortran calcintegrals.f90 -o  ./$dir/calcint
gfortran calcworkstep3.f90 -o  ./$dir/s3work
gfortran -g step3IN.f90 -o ./$dir/s3inp
gfortran step4IN.f90 -o ./$dir/s4inp
gfortran dump2data.f90 -o ./$dir/dump2data
gfortran calcinterstep3.f90 -o ./$dir/s3int
gfortran -g c3cryst_wells.f90 -o ./$dir/s3CrySv


#cp *.sh ../.
#cp fxlyrs calcint s3work s3inp s4inp dump2data s3int ../.
