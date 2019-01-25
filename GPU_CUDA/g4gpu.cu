/*==============================================================*/
/* Program: G4CF                                                */
/*                                                              *
/* Description: Program calculate four-body polarisability      */
/* anisotropy correlation function. The calculations of this    */
/* function are the most intensive in the view of interaction   */
/* induced light scattering phenomena.                          */
/*                                                              */
/* Authors: Aleksander Dawid, Krzysztof Górny,
/*          Daniel Wojcieszyk, Zygmunt Gburski                  */
/* University of Silesia, Katowice, Poland  (c) 2014            */
/*==============================================================*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <cutil.h>

float ver=1.0f;

#define TCOR 200

//===========================================
#include "headers/defines.h"
#include "headers/structures.h"
#include "headers/globals.h"
#include "headers/fgpu.h"
#include "headers/cuTab.h"
#include "headers/inout.h"
//===========================================

int main(int argc, char *argv[])
{


	if(argc<2){

/*=======================================================================================================================*/
	
		printf("\n4-Body Correlation Function - ver %1.2f(CUDA)\nSilesia University 2014\nAuthors: A. Dawid, K. Górny, D. Wojcieszyk, Z. Gburski\n",ver);
        printf("=========================\nFunction\n=========================\n");
		printf("G4GPU reads *.xyz files\n");
		printf("Usage: G4GPU Points.xyz <OPTIONS>\n");
        printf("<OPTIONS> \n");
        printf("-n<value> : # of points to calculate\n");
        printf("-s<value> : starting point\n");
        printf("-i<value> : interval between points\n");
		return -1;
	}

 if(argv[1][1]=='V')
            {
			  PrintDevicesInformation();
              return 0;			 
		    }
		    
 
	

/*===================================================*/
/*Main program                                       */
/*===================================================*/
  
  if(argc>=2){
   if(argv[1][0]=='-'){
     return 0;
    }

//----------------------------------------
//Set parameters
//----------------------------------------


/*=====================================================*/
/* Default params                                      */
/*=====================================================*/

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
	
/*==================================================*/	  
	  
	return 0;
}
