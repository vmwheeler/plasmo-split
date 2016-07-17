FC=gfortran

#LIBLOCS = -L/usr/lib64
LIBS = -lm libcuba.a

OBJ = suzuki.o

#CUBASTUFF = -L/home/vmwheeler/Software/local/lib -lcuba -lm

# Fortran compiler flags. 
#FCFLAGS=
#FCLIBS=

DEFINES = -DNMAX_COEFFS=700 -DCHECK_SURFACE_MODE -DDATA_DIR="'/home/vmwheeler/Code/plasmo-split/bhfield/nkdata/'"

# DBOPT = -DCHECK_UNDERFLOW -DCHECK_TANGENTIAL_CONTINUITY -DDEBUG_BESSEL -DOUTPUT_YBESSE -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -g
DBOPT = -DCHECK_UNDERFLOW -DCHECK_TANGENTIAL_CONTINUITY -DDEBUG_BESSEL -fbounds-check

MAIN = test_main

# clear up the default gnu make
.DEFAULT:
.SUFFIXES:
.SUFFIXES: .f .o .exe
# .PRECIOUS: %.o %-arp.o %-arp-bd.o
.PHONY: all clean allclean

%.exe: %.o

%.o: %.f
	$(FC) -cpp $(DEFINES) $(LIBLOCS) $(DBOPT) -c -O2 -Wall -o $@ $+ 

# all: bhfield-arp.exe bhfield-std.exe bhfield-arp-db.exe bhfield-std-db.exe
all: $(MAIN).exe


%.exe: %.f $(OBJ)
	$(FC) -cpp $(DEFINES) $(LIBLOCS) $(LIBS) -g -static -O2 -Wall -o $@ $+ libcuba.a

clean:
	rm -f *o *exe

allclean:
	rm -f *o *exe *dat *log
	

