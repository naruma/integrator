FC = ifort
FCFLAGS= -O3 -trackback -warn all -g -mkl
OBJS=  gear_parameter.o ODE.o MD_EOM.o test_MD.o

a.out:$(OBJS)
	$(FC) $(FCFLAGS) -o $@ $(OBJS)

.SUFFIXES: .f90 .f03
.f90.o:
	$(FC) $(FCFLAGS) -c $<
.f03.o:
	$(FC) $(FCFLAGS) -c $<


clean:
	rm -f $(OBJS) a.out *.mod
