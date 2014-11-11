#!bin/csh


ifort -c ODE.f90
ifort -c test.f90 
ifort test.o ODE.o -o a.out
./a.out

