################################################
##                                            ##
##   MAKEFILE FOR G4GPU                       ##
##                                            ##
##  Four-body anisotropy polarisability       ##
##  correlation function.                     ##
##                                            ##
##  Author: Aleksander Dawid     (c) 2019     ##
##  CUDA 10.1                                 ##
################################################


###############################
## USER VARIABLE INFORMATION ##
###############################
SEQEXEC = G4GPU
SEQOBJ = g4gpu.o
CFILES = g4gpu.cu
ALGFLAGS = -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52  -gencode=arch=compute_52,code=compute_52 -DPSINGLE -DLX


########################################
## ARCHITECTURE DEPENDENT INFORMATION ##
########################################
NVCC = /usr/local/cuda/bin/nvcc
CC = g++ -fPIC
LINK = gcc

## LOCAL DIRECTORY TREE ## Not portable to all Linux based systems

INCDIR = /usr/local/cuda/samples/common/inc/

## SIMPLE FLAG SETTINGS ##
OPTIM    = -O2 
LIBFLAG  = -lcudart 
WARNFLAG = -Wall 

LOCFLAG  = -L$(LIBDIR) $(LIBFLAG)
############################################
$(SEQOBJ): $(CFILES)
	$(NVCC) $(OPTIM) $(ALGFLAGS) -c $(CFILES) -I$(INCDIR)
	make $(SEQEXEC)
	
$(SEQEXEC): $(SEQOBJ)
	$(NVCC) $^ -o $(SEQEXEC)
	rm $(SEQOBJ)

## MAKE CLEAN ##
clean:
	rm $(SEQEXEC) $(SEQOBJ)

## MAKE ALL ##
all: $(SEQEXEC) 
