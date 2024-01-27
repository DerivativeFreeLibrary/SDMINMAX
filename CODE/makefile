FC = gfortran
RM = rm -f

FFLAGS = -O3

OBJS = problem.o main.o SDMINMAX.o

all:  $(OBJS) 
	$(FC) -o minmax $(OBJS)

.SUFFIXES : .f90 .o
	    .f   .o

.f90.o: $* ; $(FC) $(FFLAGS) -c $*.f90
.f.o:   $* ; $(FC) $(FFLAGS) -c $*.f

clean: 
	$(RM) *.o
	$(RM) minmax

