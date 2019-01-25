################################################################
##                                                            ##
##   MAKEFILE                                                 ##
##                                                            ##
##  4-body anisotropy polarizability                          ##
##  correlation function                                      ##
##                                                            ##
##  Authors: A. Dawid, K. Górny, D. Wojcieszyk, Z. Gburski    ##     
##  University of Silesia 2014                                ##
##                                                            ##
################################################################

## DLX       -Linux headers

###############################
## USER VARIABLE INFORMATION ##
###############################
SEQEXEC = ../../Example/G4CPU
SEQOBJ = g4cpu.o
PREFLAGS = -DLX
CFILES = g4cpu.cc


########################################
## ARCHITECTURE DEPENDENT INFORMATION ##
########################################
CC = g++

## LOCAL DIRECTORY TREE ##
INCDIR = .
LIBDIR = lib

## FLAG SETTINGS ##
OPTIM    = -O2 
ARCHFLAG = -march=native
WARNFLAG = -Wall 


LOCFLAG  = -L$(LIBDIR) 

############################################
## ACTUAL MAKE CODE BELOW (DO NOT CHANGE) ##
############################################

## LOCAL SEQUENTIAL VERSION ##
$(SEQOBJ): $(CFILES)
	$(CC) $(OPTIM) $(WARNFLAG) $(PREFLAGS) $(ARCHFLAG) -c $(CFILES) -I$(INCDIR)
	make $(SEQEXEC)
	
$(SEQEXEC): $(SEQOBJ)
	$(CC) $(LOCFLAG) $^ -o $(SEQEXEC)
	rm $(SEQOBJ)
    

## MAKE CLEAN ##
clean:
	rm $(SEQEXEC) $(SEQOBJ)

## MAKE ALL ##
all: $(SEQEXEC) 
