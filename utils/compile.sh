#! /bin/bash



gfortran fixedlayer.f90 -o fxlyrs
gfortran calcintegrals.f90 -o  calcint
gfortran calcworkstep3.f90 -o  s3work
gfortran -g step3IN.f90 -o s3inp
gfortran step4IN.f90 -o s4inp
gfortran dump2data.f90 -o dump2data
gfortran calcinterstep3.f90 -o s3int

#cp *.sh ../.
#cp fxlyrs calcint s3work s3inp s4inp dump2data s3int ../.
