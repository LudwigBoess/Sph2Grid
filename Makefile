## Sph2Grid Makefile ##
SHELL = /bin/bash		# This should always be present

# Target Computer #
ifndef $(SYSTYPE)
SYSTYPE	= $(shell hostname)
endif

#Std systype
CC	  = mpicc
OPTIMIZE  = -Wall -g -O3 -fopenmp
MPI_LIBS  = -lmpi -L/usr/local/lib 
MPI_INCL  = -I/usr/local/include
GSL_INCL  =
GSL_LIBS  = 
HDF5_INCL = 
HDF5_LIBS = -lhdf5
SILO_INCL = 
SILO_LIBS = 
FFTW_LIBS = -lfftw3_mpi -lfftw3 -lfftw3_omp
FFTW_INCL = 

ifeq ($(SYSTYPE),"INAF_IRA")
CC	  = mpicc
OPTIMIZE  = -std=c99 -strict-ansi -Wall -openmp -O3 -m64  -xhost -mkl -ipo4 -ansi-alias -ansi-alias-check -opt-matmul -g
MPI_LIBS  = -lmpich  -L/homes/donnert/Libs/lib
MPI_INCL  = -I/homes/donnert/Libs/include
GSL_INCL  = 
GSL_LIBS  = 
HDF5_INCL = 
HDF5_LIBS = 
SILO_INCL = -I/homes/donnert/Libs/include
SILO_LIBS = -lsilo -lhdf5  -lz -L/homes/donnert/Libs/lib
FFTW_LIBS = -lfftw3_mpi -lfftw3 -lfftw3_omp
FFTW_INCL = 
endif

ifeq ($(SYSTYPE),mach64.ira.inaf.it)
CC	  = mpicc
OPTIMIZE  = -O2 -std=c99 -Wall -g -m64 -fopenmp # -Ofast -mtune=native -march=bdver1  -mprefer-avx128 -fopenmp  -minline-all-stringops -fprefetch-loop-arrays --param prefetch-latency=300 -funroll-all-loops
MPI_LIBS  = -lmpich  -L/homes/donnert/Libs/lib
MPI_INCL  = -I/homes/donnert/Libs/include
GSL_INCL  = 
GSL_LIBS  = 
HDF5_INCL = 
HDF5_LIBS = 
SILO_INCL = -I/homes/donnert/Libs/include
SILO_LIBS =  -lhdf5 -lz -L/homes/donnert/Libs/lib
FFTW_LIBS = -lfftw3_mpi -lfftw3 -lfftw3_omp
FFTW_INCL = 
endif


ifeq ($(SYSTYPE),"CRUX")
CC       = /home/crux/abeck/progs/openmpi143/bin/mpicc   
OPTIMIZE = -O3 -Wall
GSL_INCL = -I/home/crux/abeck/progs/gsl115/include
GSL_LIBS = -L/home/crux/abeck/progs/gsl115/lib
FFTW_INCL= -I/home/crux/abeck/progs/fftw2/include
FFTW_LIBS= -L/home/crux/abeck/progs/fftw2/lib
MPI_INCL = -L/home/crux/abeck/progs/openmpi143/include
MPI_LIBS = -L/home/crux/abeck/progs/openmpi143/lib
HDF5_INCL = -I/home/crux/abeck/progs/hdf5/include
HDF5_LIBS  = -L/home/crux/abeck/progs/hdf5/lib
endif

ifeq ($(SYSTYPE),"DARWIN")
CC     	  =  mpicc
OPTIMIZE  = -O3 -Wall -m64 -g -fopenmp -mtune=native
MPI_LIBS  = -lmpich
MPI_INCL  = 
GSL_INCL  =  
GSL_LIBS  =  
HDF5_INCL = 
HDF5_LIBS =
SILO_INCL = 
SILO_LIBS = 
FFTW_LIBS = -lfftw3_mpi -lfftw3 -lfftw3_omp
FFTW_INCL = 
endif

ifeq ($(SYSTYPE),"MPA")
CC        =  mpicc
OPTIMIZE  = -O3 -Wall -g  -fipa-struct-reorg -combine
MPI_LIBS  = -lmpi 
MPI_INCL  = -lm
GSL_INCL  = 
GSL_LIBS  =
HDF5_INCL = -I/usr/common/pdsoft/appl/hdf5-1.8.4p1/include
HDF5_LIBS = -L/usr/common/pdsoft/appl/hdf5-1.8.4p1/lib
SILO_INCL = -I/usr/common/pdsoft/appl/silo-4.7.2/include
SILO_LIBS = -L/usr/common/pdsoft/appl/silo-4.7.2/lib -lsilo 
FFTW_LIBS = 
FFTW_INCL = 
endif


ifeq ($(SYSTYPE),"MPA64")	# fo 64bit support
CC        =  mpicc
OPTIMIZE  = -O3 -Wall -g -m64 -fipa-struct-reorg -combine
MPI_LIBS  = -lmpi 
MPI_INCL  = -lm
GSL_INCL  = 
GSL_LIBS  =
HDF5_INCL = -I/usr/common/pdsoft/appl/hdf5-1.8.6-64bit/include
HDF5_LIBS = -lz -L/usr/common/pdsoft/appl/hdf5-1.8.6-64bit/lib
SILO_INCL = -I/usr/common/pdsoft/appl/silo-4.7.2-64bit/include
SILO_LIBS = -L/usr/common/pdsoft/appl/silo-4.7.2-64bit/lib
FFTW_LIBS = 
FFTW_INCL = 
endif

ifeq ($(SYSTYPE),"RZG_OPA") 
CC       =  mpicc
OPTIMIZE = -cc=icc -O3  -Wall -g -m64 -openmp 
MPI_LIBS = -lm
MPI_INCL = 
GSL_INCL =  -I/afs/@cell/common/soft/gsl/1.14/@sys/include
GSL_LIBS =  -L/afs/@cell/common/soft/gsl/1.14/@sys/lib
HDF5_INCL = -I/afs/@cell/common/soft/hdf5/1.8.7/@sys/intel/11.1/impi/4.0.0/include
HDF5_LIBS = -L/afs/@cell/common/soft/hdf5/1.8.7/@sys/intel/11.1/impi/4.0.0/lib
SILO_INCL = 
SILO_LIBS = 
FFTW_LIBS = -lfftw3_mpi -lfftw3 -L/afs/ipp/home/j/jdonnert/Libs/lib
FFTW_INCL = -I/afs/ipp/home/j/jdonnert/Libs/include
endif

SRCDIR	= src/

# Sources and objects
SRCFILES = allvars.c main.c domain_decomp.c print_settings.c \
	  			timing.c setup.c fill_grid.c test_grid.c\
	  			input.c output.c ngb.c tree.c \
				unit.c cosmo.c	\
				sph.c ngp.c cic.c tsc.c dwavelet.c \
				fft.c powerspectrum.c
SOURCES = $(addprefix $(SRCDIR),$(SRCFILES))

OBJECTS	= $(SOURCES:.c=.o)

# Include files
INCLFILES= allvars.h  input.h ngb.h tree.h cosmo.h unit.h \
		   schemes.h proto.h config.h

INCLUDES= $(addprefix $(SRCDIR),$(INCLFILES)) Makefile

CFLAGS	= -std=c99 $(OPTIMIZE) $(OPT) $(SILO_INCL) $(HDF5_INCL) \
			  $(GSL_INCL) $(MPI_INCL) $(FFTW_INCL)  

LIBS	= $(MPI_LIBS)  $(HDF5_LIBS) $(GSL_LIBS) $(FFTW_LIBS) $(SILO_LIBS)\
		  -lhdf5 -lgsl -lgslcblas   


EXEC = Sph2Grid


## TARGETS  ##

default : D20 $(EXEC)	# why would you want something else ?

SPH : CFLAGS += -DSPH 	# we need those or the objects will fail
SPH : $(OBJECTS)
	$(call compile,SPH)

NGP : CFLAGS += -DNGP 
NGP : $(OBJECTS)
	$(call compile,NGP)

CIC : CFLAGS += -DCIC
CIC : $(OBJECTS)
	$(call compile,CIC)

TSC : CFLAGS += -DTSC 
TSC : $(OBJECTS)
	$(call compile,TSC)

D20 : CFLAGS += -DD20 
D20 : $(OBJECTS)
	$(call compile,D20)

clean	: 
	rm -f  $(OBJECTS) $(EXEC)_* $(SRCDIR)config.h 

install	: 
	cp -i $(EXEC)_* ~/bin

## RULES ##

$(OBJECTS) : $(INCLUDES)

$(EXEC) :
	ln -s Sph2Grid_D20 Sph2Grid

$(SRCDIR)config.h : Config		# This is essentially a BASH script
	sed '/^#/d; /^$$/d; s/^/#define /g' Config > $(SRCDIR)config.h	
	echo '#include "allvars.h"' >  $(SRCDIR)print_settings.c
	echo 'void print_compile_time_settings(){' >> $(SRCDIR)print_settings.c
	echo '	printf("Compiled with : \n " ' >> $(SRCDIR)print_settings.c
	sed '/^#/d; /^$$/d; s/^/"	/g; s/$$/ \\n"/g;' Config >>  $(SRCDIR)print_settings.c
	echo ');}' >> $(SRCDIR)print_settings.c

## FUNCTIONS ##

define compile
	$(CC) $(CFLAGS) $(OBJECTS)  $(LIBS) -o $(EXEC)_$1
	cd $(SRCDIR) && ctags *.[ch]
	rm $(SRCDIR)config.h	
endef

