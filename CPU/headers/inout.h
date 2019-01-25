/*=============================================================*/
/* Input/outpu functions                                       */
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

	  if(fscanf(out,"%d\n",&N)>0){}
	  if(fscanf(out,"%f %s %f %f %f %f\n",&dt,AtomName,&Lx,&Ly,&Lz,&sigmaAr)>0){}
	  if(idx==0){ psdt=dt;}
      

      for(n=0;n<N;n++)
	  {       
       if(fscanf(out,"%s %f %f %f\n",AtomNumber,&x,&y,&z)>0){}	   
	  }
    
    idx++;
  }while(!feof(out)); 

	fclose(out);
	printf("------Simulation parameters------\nCell size: (%f,%f,%f)\nSigma: %f\n",Lx*2,Ly*2,Lz*2,sigmaAr);
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

	  if(fscanf(out,"%d\n",&N)>0){}
	  if(fscanf(out,"%f %s %f %f %f %f\n",&dt,AtomName,&Lx,&Ly,&Lz,&sigmaAr)>0){}
      

      for(n=0;n<N;n++)
	  {       
       if(fscanf(out,"%s %f %f %f\n",AtomNumber,&x,&y,&z)>0){}	
	   
		  TR2[n][idx].x = x;
		  TR2[n][idx].y = y;
		  TR2[n][idx].z = z;
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
    printf("# of particles: %d\n# of simulation points: %d\n",N,SimTime);
	
   
//------------------- Initialize trajectory table TR2 -----------------------  
    TR2=InitTab(TotalPoints,N);
//------------------- Read points from file to TR2 table --------------------
    TotalPoints=GetPoints(file);  

    GetScatterPoints(SimTime);

    CalculateScattering(SimTime);
  
return 0;
}

/*=====================================================*/
void hand_sigint(int t)
{

printf("\nUSER INTERUPTION OR KILL SIGNAL\n");
exit(1);
}

