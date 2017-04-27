FC = gfortran -ffixed-line-length-none
optimize=-O
FFLAGS=$(optimize)
CFLAGS=$(optimize)
#if SAC library has been installed, uncomment the next two lines
CFLAGS=$(optimize) -DSAC_LIB
SACLIB=-L$(SACHOME)/lib -lsac -lsacio

all: trav2 


trav2: trav.o tau_p.o
	$(LINK.f) -o $@ trav.o tau_p.o

clean:
	rm -f *.o 
