################################################
##                                            ##
##   MAKEFILE FOR G4CF                        ##
##                                            ##
##  Four-body anisotropy polarisability       ##
##  correlation function.                     ##
##                                            ##
##  Author: Aleksander Dawid     (c) 2014     ##
##                                            ##
################################################



###############################
## USER VARIABLE INFORMATION ##
###############################
SEQEXEC = ../../Example/G4GPU
SEQOBJ = g4gpu.o
ALGFLAGS = --use_fast_math -arch sm_20 -DPSINGLE -DLX 
# sm_20 - CUDA ver 2.0 We need makefile to detect CUDA version instaled or not
POTFLAGS = 
PREFLAGS = 
CFILES = g4gpu.cu

########################################
## ARCHITECTURE DEPENDENT INFORMATION ##
########################################
NVCC = /usr/local/cuda/bin/nvcc
CC = g++ -fPIC
LINK = gcc

## LOCAL CUDA DIRECTORY TREE ## Not portable to all Linux based systems

INCDIR = /usr/local/cuda/sdk/C/common/inc/
LIBDIR = /usr/local/cuda/sdk/C/lib/

## SIMPLE FLAG SETTINGS ##
OPTIM    = -O2 #-fomit-frame-pointer
ARCHFLAG = #-march=i686
MATHFLAG = #-lm #-ffast-math#
LIBFLAG  = -lcudart /usr/local/cuda/sdk/C/lib/libcutil_x86_64.a
WARNFLAG = -Wall #-Wno-deprecated

LOCFLAG  = -L$(LIBDIR) $(LIBFLAG)
############################################
$(SEQOBJ): $(CFILES)
	$(NVCC) $(OPTIM) $(ALGFLAGS) -c $(CFILES) -I$(INCDIR)
	make $(SEQEXEC)
	
$(SEQEXEC): $(SEQOBJ)
	$(NVCC) $(LOCFLAG) $^ -o $(SEQEXEC)
	rm $(SEQOBJ)

## MAKE CLEAN ##
clean:
	rm $(SEQEXEC) $(SEQOBJ)

## MAKE ALL ##
all: $(SEQEXEC) 
