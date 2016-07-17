FC=gfortran

LIBLOCS = -L/usr/lib64 -L/home/vmwheeler/Software/local/lib
LIBS = -lm

CUBASTUFF = -L/home/vmwheeler/Software/local/lib -lcuba -lm

# Fortran compiler flags. 
#FCFLAGS=
#FCLIBS=

DEFINES = -DNMAX_COEFFS=700 -DCHECK_SURFACE_MODE -DDATA_DIR="'/home/vmwheeler/Code/plasmo-split/bhfield/nkdata/'"

# DBOPT = -DCHECK_UNDERFLOW -DCHECK_TANGENTIAL_CONTINUITY -DDEBUG_BESSEL -DOUTPUT_YBESSE -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -g
DBOPT = -DCHECK_UNDERFLOW -DCHECK_TANGENTIAL_CONTINUITY -DDEBUG_BESSEL -fbounds-check


# clear up the default gnu make
.DEFAULT:
.SUFFIXES:
.SUFFIXES: .f .o .exe
# .PRECIOUS: %.o %-arp.o %-arp-bd.o
.PHONY: all clean ftnchek allclean
%.o: %.f
%.exe: %.o

# all: bhfield-arp.exe bhfield-std.exe bhfield-arp-db.exe bhfield-std-db.exe
all: bhfield-std-db.exe

# standard version: one-step compilation
%-std.exe: %.f
	$(FC) -cpp $(DEFINES) $(GSLSTUFF)       -O2 -Wall -static -o $@ $+

%-std-db.exe: %.f
	$(FC) -cpp $(DEFINES) $(LIBLOCS) $(LIBS) $(DBOPT) -g -static -O2 -Wall -o $@ $+ libcuba.a
#	$(FC) -cpp $(DEFINES) $(CUBASTUFF) $(DBOPT) -static -O2 -Wall -o $@ $+

clean:
	rm -f *o *exe

allclean:
	rm -f *o *exe *dat *log
	

