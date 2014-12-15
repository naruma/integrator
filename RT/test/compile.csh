#!bin/csh

ifort -warn all -traceback -g time_evolution.f90 RT_test.f90 -mkl

