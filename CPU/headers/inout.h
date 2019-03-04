/*=============================================================*/
/* Input-Output functions                                      */
/*=============================================================*/

inline POS** InitTab(int SimTime, int NUMOL)
{	
	int p;
	POS **tTab;

	
          if ((tTab = new POS*[NUMOL+1]) == NULL)
		  {
            printf("Out of memory\n");
            exit(1);  
		  }
  
          for(p=0;p<NUMOL;p++){
            if((tTab[p] = new POS[SimTime+1]) == NULL)
			{
             printf("Out of memory\n");
             exit(1);  
			} 
		  }
    
          ile_mem+=(NUMOL+1)*(SimTime+1)*sizeof(double);
	  
        
		return tTab;
}

/*=============================================================*/
int inline CheckNumOfPoints(char *file)
{
	char AtomName[64],AtomNumber[32];
	int idx=0,n;
	float x,y,z;
	FILE *out;
	

	 out=fopen(file,"r");
  if(out==NULL){
	  printf("File not found: %s\n",file);
	  return -1;
  }
	
   
  do{

	  ff=fscanf(out,"%d\n",&N);
	  ff=fscanf(out,"%f %s %f %f %f %f\n",&dt,AtomName,&Lx,&Ly,&Lz,&sigmaAr);
          if(idx==0){psdt=dt;}
      

      for(n=0;n<N;n++)
	  {       
       ff=fscanf(out,"%s %f %f %f\n",AtomNumber,&x,&y,&z);	   
	  }
    
    idx++;
  }while(!feof(out)); 

	fclose(out);
	printf("Simulation parameters\n--------------------------------------\n Cell size: (%f,%f,%f)\n LJsigma=%f\n",Lx*2,Ly*2,Lz*2,sigmaAr);
return idx;
}
/*=============================================================*/
int inline GetPoints(char *file)
{
   char AtomName[64],AtomNumber[32];
	int idx=0,n;
	float x,y,z;
	FILE *out;

	 out=fopen(file,"r");
  if(out==NULL){
	  printf("File not found: %s\n",file);
	  return -1;
  }
	
   
  do{

	  ff=fscanf(out,"%d\n",&N);
	  ff=fscanf(out,"%f %s %f %f %f %f\n",&dt,AtomName,&Lx,&Ly,&Lz,&sigmaAr);
      

      for(n=0;n<N;n++)
	  {       
       ff=fscanf(out,"%s %f %f %f\n",AtomNumber,&x,&y,&z);	
	   
		  Traj[n+N*idx].x = x;
		  Traj[n+N*idx].y = y;
		  Traj[n+N*idx].z = z;
	  }

	  
    idx++;
  }while(!feof(out)); 

	fclose(out);  
	  
	return 0;
}

/*=============================================================*/
/* Main Calculation Loop                                       */
/*=============================================================*/
int inline CalculationLoop(char *file)
{
	
	int TotalPoints;
	
//------------------------------------------------------------


       TotalPoints=CheckNumOfPoints(file);  

	Lxx=Lx*2;
	Lyy=Ly*2;
	Lzz=Lz*2;
	vol=sigmaAr*sigmaAr*sigmaAr;

	SimTime=TotalPoints-TCOR;
    printf("# of particles: %d\n# of time points: %d\n",N,SimTime);
	
    ARGON=N;
//------------------- Init Trajectory Table -----------------------
       pa2Size=C2_idx(ARGON);
	   n4p=C4_pNum(ARGON);  
//==================================================================   
    printf("-----------------------\nNumber of pairs: %d\n-------------------------\n",pa2Size);


    
 int blockres=256;

   
    TPB=blockres;
    xsiz=blockres,ysiz=1;
 

  
    parSize=xsiz*ysiz*TPB;

    if(pa2Size==parSize){
		MulTot=1;
	}else{
        MulTot=pa2Size/parSize + 1;
	}

    sharedSize=(TPB+1)*sizeof(real);
     
   
   InitTraj(TotalPoints,N);

 
/*---------------------------------------------------------------------------*/
  
   TR2=InitTab(TotalPoints,N);

//-----------------------------------------------------------------

	TotalPoints=GetPoints(file);


  
         sdkCreateTimer(&timer);  
         cudaDeviceSynchronize();          
          sdkResetTimer(&timer);         
          sdkStartTimer(&timer);
           
           LoadBeta(ARGON,SimTime);
	
            cudaDeviceSynchronize();            
            sdkStopTimer(&timer);
            gpuTimeRun = sdkGetTimerValue(&timer);
            TimeInGPU=gpuTimeRun/1000;
            printf("Beta Table load time %lf [s] \n",TimeInGPU);



  
  C2_idx_pair(ARGON);
  CUDAStart();
  SaveCvv4(T,tmstart);
 
 
return 0;
}

/*-------------------------------------------------------------------------------*/
void hand_sigint(int t)
{

printf("\nUSER INTERUPTION OR KILL SIGNAL\nFree CUDA Memory\n");
bEnd=1;
}
