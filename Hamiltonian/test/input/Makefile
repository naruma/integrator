FC = ifort
FCFLAGS= -O3 -trackback -warn all -g
OBJS=  input.o test_input.o

a.out:$(OBJS)
	$(FC) $(FCFLAGS) -o $@ $(OBJS)

.SUFFIXES: .f90
.f90.o:
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f $(OBJS) a.out *.mod
