/*==============================================================*/
/* Interaction induced scattering (IIS) calculation module for  */
/* RIGMD simulation application                                 */
/*                                                              */
/* Author: Aleksander Dawid                                     */
/* University of Silesia, Katowice, Poland                      */
/*==============================================================*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <signal.h>


float ver=1.00f;

#define real float

// Decorelation time
#define TCOR 200

#include "headers/structures.h"
#include "headers/globals.h"
#include "headers/ci_func.h"
#include "headers/inout.h"
//===========================================

int main(int argc, char *argv[])
{

	if(argc<2){

/*=======================================================================================================================*/
	
		printf("\nInteraction Induced Scattering calculation 4-body - ver %1.2f\nSilesia University 2014\nAuthors: A. Dawid, K. Górny, D. Wojcieszyk, Z. Gburski\n",ver);
        printf("=========================\nFunction\n=========================\n");
		printf("Usage: G4CPU Points.xyz [OPTIONS]\n");
		printf("[OPTIONS]\n");
		printf("-n <value> : # of points to calculate\n");
        printf("-s <value> : starting point\n");
        printf("-i <value> : interval between points\n");
		
		return -1;
	}

	if(argc>=2){
   if(argv[1][0]=='-'){
     return 0;
    }


//Starting point of 4-body CF; intended to use for many Graphics Cards
tmstart=0;
//Number of points to calculate
trwanie=1;
// The interval between steps 
steptime=1;
int i;
  for(i=2;i<argc;i++){

    if(argv[i][1]=='s'){
        char *nnp;
        nnp=&argv[i][2];   
        tmstart=atoi(nnp);    
     }

     if(argv[i][1]=='n'){
        char *nnp;
        nnp=&argv[i][2];   
        trwanie=atoi(nnp);    
     }

    if(argv[i][1]=='i'){
        char *nnp;
        nnp=&argv[i][2];   
        steptime=atoi(nnp);    
     }
       
   }

	signal(SIGINT,hand_sigint);
	signal(SIGKILL,hand_sigint);
	
	
	CalculationLoop(argv[1]);
  }
	
	
	return 0;
}
