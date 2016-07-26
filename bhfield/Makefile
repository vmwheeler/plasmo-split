FC=gfortran

LIBS = -lm libcuba.a

OBJ = suzuki.o

CUBALIBS = -lm -L/home/vmwheeler/Software/local/lib -lcuba
#CUBAINC = -I/home/vmwheeler/Software/local/include

# Fortran compiler flags. 
#FCFLAGS=
#FCLIBS=

DEFINES = -DNMAX_COEFFS=700 -DCHECK_SURFACE_MODE -DDATA_DIR="'/home/vmwheeler/Code/plasmo-split/bhfield/data/'"

# DBOPT = -DCHECK_UNDERFLOW -DCHECK_TANGENTIAL_CONTINUITY -DDEBUG_BESSEL -DOUTPUT_YBESSE -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -g
DBOPT = -DCHECK_UNDERFLOW -DCHECK_TANGENTIAL_CONTINUITY -DDEBUG_BESSEL -fbounds-check

MAIN = full_spectrum
#MAIN = spectral_shells

# clear up the default gnu make
.DEFAULT:
.SUFFIXES:
.SUFFIXES: .f .o .exe
# .PRECIOUS: %.o %-arp.o %-arp-bd.o
.PHONY: all clean allclean

%.exe: %.o

%.o: %.f
	$(FC) -cpp $(DEFINES) $(LIBLOCS) $(DBOPT) -c -O2 -Wall -o $@ $^ 

all: $(MAIN).exe


%.exe: %.f $(OBJ)
#	$(FC) -cpp $(DEFINES) $(LIBLOCS) $(LIBS) -g -static -O2 -Wall -o $@ $^ libcuba.a
#	$(FC) -cpp $(DEFINES) $(DBOPT) -g -O2 -Wall -static -o $@ $^ $(CUBALIBS)
	$(FC) -cpp $(DEFINES) $(DBOPT) -g -O2 -Wall -o $@ $^ $(CUBALIBS)
clean:
	rm -f *o *exe

allclean:
	rm -f *o *exe *dat *log
	

