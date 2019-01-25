/*====================================================================*/
/* Initialization of trajectory table                                 */
/*====================================================================*/

void InitTraj(int TotalPoints, int N)
{
 double TotMB;
 Traj_size = (TotalPoints+1) * N * sizeof(real4);
 Cvv_size = (SimTime+1) * sizeof(real);
 N_size = parSize * sizeof(real);
 calc_size = parSize * sizeof(ushort4);
 p2_size = pa2Size * sizeof(ushort4);
 BetaS = (SimTime+1) * pa2Size * sizeof(real);
	

 
TotMB=(Traj_size)/1000000;
//printf("TotalPoints=%d 	Traj_size=%.2lf MB\n",TotalPoints,TotMB);
 
 TotMB=(N_size + p2_size + calc_size + Cvv_size + BetaS)/1000000;

 printf("CUDA memory: %.2lf MB\n",TotMB);
	if(TotMB>3999){
     exit(1);
	}
 
	
 Traj = (real4*)malloc(Traj_size);
 Cvv = (real*)malloc(Cvv_size);
 ArTab = (real*)malloc(N_size);
 Beta = (real*)malloc(BetaS);

	
	

 PairTab = (ushort4*)malloc(calc_size);
 memset(PairTab,0,calc_size);

 P2T = (ushort4*)malloc(p2_size);
	
 memset(Traj,0,Traj_size);
 memset(Cvv,0,Cvv_size);
 memset(ArTab,0,N_size);

 memset(Beta,0,BetaS);

 memset(P2T,0,p2_size);

 
}


/*====================================================================*/
/* Number of pair indices                                             */
/*====================================================================*/
int C2_idx(int N)
{
	int Nm1=N-1;
	int sum=0;
        int i,j;

	for(i=0;i<Nm1;i++)
				  {
					  
				  for(j=i+1;j<N;j++)
					{
						sum++;
					}
					  
				  }
	
	return sum;
}
/*====================================================================*/
/* (delete)                                                           */
/*====================================================================*/
int inline C4_pNum(int N)
	{
	 unsigned long sum=0;
	 int Nm1;
	 register int k,i,j,l;

	 Nm1=N-1;

			 
				 i=4;j=5;
						 
					 for(k=0;k<Nm1;k++)
						{
						 for(l=k+1;l<N;l++)
							 {						 
						      if(i!=k && j!=k && i!=l && j!=l)
							  {
								  sum++;
								  //printf("sum=%d\n",sum);
                                
							  }//conditions
						 }
					 }						 			

return sum;
}


/*====================================================================*/
/* Number of pair indices used in GPU                                 */
/*====================================================================*/
void inline C2_idx_pair(int N)
{
	int Nm1=N-1;
	int idx=0;
        int i,j;

	for(i=0;i<Nm1;i++)
				  {
					  
				  for(j=i+1;j<N;j++)
					{
						P2T[idx].x=i;
						P2T[idx].y=j;
						P2T[idx].z=i;
						P2T[idx].w=j;
						
						idx++;
					}
					  
				  }
	
}

/*====================================================================*/
/* Polarizability anisotropy between two sites                        */
/*====================================================================*/
real inline CalcBeta(int N, int i,int j,int t){

	real4 Rij;
	real rv,IRij,IRij2;
	real Beta,VersRij;
	int it,jt;
	it=i+N*t;
	jt=j+N*t;

	
	
	Rij.x=Traj[it].x-Traj[jt].x;
	Rij.y=Traj[it].y-Traj[jt].y;
	Rij.z=Traj[it].z-Traj[jt].z;

//------------------------------------------------------------------------------
// Minimum Image Convention
//-----------------------------------------------------------------------------
              //X-axis
	                if (Rij.x>=Lx)         Rij.x-=Lxx;
	                else if (Rij.x<-Lx)    Rij.x+=Lxx;
              //Y-axis
	                if (Rij.y>=Ly)         Rij.y-=Lyy;
	                else if (Rij.y<-Ly)    Rij.y+=Lyy;
              //Z-axis
	                if (Rij.z>=Lz)        Rij.z-=Lzz;
	                else if (Rij.z<-Lz)   Rij.z+=Lzz;
//-----------------------------------------------------------------------------

	Rij.w=Rij.x*Rij.x+Rij.y*Rij.y+Rij.z*Rij.z;
    IRij2=1/Rij.w;
	
	rv=sqrt(Rij.w);
	IRij=1/rv;
	
     VersRij=vol*Rij.x*Rij.z*IRij;

	 Beta=VersRij*IRij2*IRij2;
	
	
	return Beta;
}

/*====================================================================*/
/* Load Beta_ij into 1dim Table                                       */
/*====================================================================*/
void inline LoadBeta(int N,int TimeSteps)
{
        int i,j;
	int Nm1=N-1;
	int idx=0;
	int t;

for(t=0;t<TimeSteps;t++)
	{
		idx=0;
	           for(i=0;i<Nm1;i++)
				  {	  
				  for(j=i+1;j<N;j++)
					{
						Beta[idx+pa2Size*t]=CalcBeta(N,i,j,t);
						idx++;
					}	  
				  }
	}

}

/*====================================================================*/
/* Allocation of GPU memory                                           */
/*====================================================================*/
void AllocateMemCUDA()
{
//------------------------------------------------------------
// Time measurement
//------------------------------------------------------------
cudaThreadSynchronize();
cutResetTimer(hTimer);
cutStartTimer(hTimer);

// Allocate positions in device memory


cudaMalloc((void**)&d_Cvv, Cvv_size);
	
cudaMalloc((void**)&d_Beta, BetaS);	
	
cudaMalloc((void**)&d_ArTab, N_size);

cudaMalloc((void**)&d_PairTab, calc_size);

cudaMalloc((void**)&d_P2T, p2_size);

cudaThreadSynchronize();
cutStopTimer(hTimer);
gpuTime = cutGetTimerValue(hTimer);
printf("GPU allocation time: %lf ms\n",gpuTime);
}

/*====================================================================*/
/* Free GPU memory                                                    */
/*====================================================================*/
void FreeCUDAMemory()
{
// Free device memory

cudaFree(d_Cvv);
cudaFree(d_Beta);
	
cudaFree(d_ArTab);
cudaFree(d_PairTab);
cudaFree(d_P2T);
}
//------------------------------------------------------------
//
//------------------------------------------------------------
void CopyHostToDevice()
{
cudaThreadSynchronize();
cutResetTimer(hTimer);
cutStartTimer(hTimer);
// Copy vectors from host memory to device memory


cudaMemcpy(d_Cvv, Cvv, Cvv_size, cudaMemcpyHostToDevice);	
	

cudaMemcpy(d_Beta, Beta, BetaS, cudaMemcpyHostToDevice);
	
cudaMemcpy(d_ArTab, ArTab, N_size, cudaMemcpyHostToDevice);
	

cudaMemcpy(d_PairTab, PairTab, calc_size, cudaMemcpyHostToDevice);


cudaMemcpy(d_P2T, P2T, p2_size, cudaMemcpyHostToDevice);

cudaThreadSynchronize();
cutStopTimer(hTimer);
gpuTime = cutGetTimerValue(hTimer);
printf("Transfer time to GPU: %lf ms\n",gpuTime);
//------------------------------------------------------------
free(Traj);
free(Beta);
}


/*===============================================================*/
/* OPTIMIZED 4-body PACF                                         */
/*===============================================================*/
void startC4(int tt,int rtime,int timestart)
{
FILE *sw;
double TTime;
int t,i,j;
char filename[256];
	
unsigned int start=0;



//printf("xsiz: %d\nysiz: %d\nTPB: %d\n",xsiz,ysiz,TPB);
dim3 grid(xsiz,ysiz,1);
dim3 grid2(ysiz,1,1);
	

TTime=0.0;
			   
	//printf("MulTot=%d\nn4p=%d\npa2Size=%d\nparSize=%d\n",MulTot,n4p,pa2Size,parSize);

	printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\nG4 function running from %d to %d every %d steps\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n",timestart,rtime,steptime);
    
	fflush(stdout);
	sprintf(filename,"G4CUDA_%d_%d_%d.swap",timestart,rtime,steptime);
	
   for(t=timestart;t<rtime;t+=steptime)
      {
		  cudaThreadSynchronize();
          cutResetTimer(hTimer);
          cutStartTimer(hTimer);
		  
       start=0;	
       
	   for(i=0;i<MulTot;i++)
		{	
		  for(j=0;j<n4p;j++)
			{	
				if(bEnd==1){ FreeCUDAMemory(); exit(1); }		
			// Create partial index			     
				
			  PairIndex<<<grid,TPB>>>(d_P2T, d_PairTab,start);
				
			  CalcFuncPair<<<grid,TPB>>>(d_Beta,d_ArTab,d_PairTab,t,start);								          	       							
                        }
                       // Sum anizotropy
			  Reducto<<<grid,TPB,sharedSize>>>(d_ArTab,d_ArTab);	   
	                  ReductoEnd<<<grid2,TPB,sharedSize>>>(d_ArTab,d_ArTab,d_Cvv,t);
                        
			//return; 		
			start += parSize;
			ZeroPairs<<<grid, TPB, sharedSize>>>(d_ArTab,d_PairTab);						
		}			 

		  
           cudaMemcpy(Cvv, d_Cvv, Cvv_size, cudaMemcpyDeviceToHost);
           
		    cudaThreadSynchronize();
            cutStopTimer(hTimer);
            gpuTimeRun = cutGetTimerValue(hTimer);
            TimeInGPU=gpuTimeRun/1000;
		  
			printf("Step %d G4=%5.10f Resolve time: %lf s\n",t,Cvv[t],TimeInGPU);
			TTime+=TimeInGPU;
			fflush(stdout);
			
			
			sw=fopen(filename,"a");
			fprintf(sw,"\n %lf %5.10f",(double)psdt*t,Cvv[t]);
			fclose(sw);
			
	  }
	
      
	printf("Time in GPU: %lf s\n",TTime);


}

/*====================================================*/
/* CUDA device information                            */
/*====================================================*/
void PrintDevicesInformation()
{
int deviceCount;
int nrdev;
int nMulProc;
cudaGetDeviceCount(&deviceCount);
printf("=========================================================\n");
printf("Number of CUDA devices: %d\n",deviceCount);
printf("---------------------------------------------------------\n");	

    for (nrdev = 0; nrdev < deviceCount; nrdev++) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, nrdev);
		printf("Device %d: %s\n", nrdev, deviceProp.name); 
		printf("maxThreadsPerBlock: %d\n", deviceProp.maxThreadsPerBlock); 
		printf("maxThreadsDim: {%d x %d x %d}\n", deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
		printf("maxGridSize: {%d x %d x %d}\n", deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
		
		nMulProc=deviceProp.multiProcessorCount;
		printf("multiProcessorCount: %d\n", nMulProc); 
		printf("computeMode: %d\n", deviceProp.computeMode); 
		printf("sharedMemPerBlock: %ld\n", deviceProp.sharedMemPerBlock); 		 			
	}
printf("=========================================================\n");	
}
/*------------------------------------------------------------------------------*/
void SetCurrentDevice(int nr){
	cudaSetDevice(nr);
	
int cdev;	
	
	cudaGetDevice(&cdev);
	cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, cdev);
		
	printf("Current device used: %d (%s)\n",cdev, deviceProp.name);
}
/*------------------------------------------------------------------------------*/
void CUDAInitParams(int Points)
{

	//Problem to be portable //Need to detect version of CUDA avaliable
    //sharedSize=4*4*4*256; //sm_20
	sharedSize=15000;   //sm_13 16384
	
printf("Block shared memory: %d Bytes\n",sharedSize);

	
cudaMemcpyToSymbol(gN, &pa2Size, sizeof(int));
cudaMemcpyToSymbol(gtt, &Points, sizeof(int));
cudaMemcpyToSymbol(gARGON, &ARGON, sizeof(int));
cudaMemcpyToSymbol(gVol, &vol, sizeof(real));

}


/*======================================================*/
/* Main CUDA procedure                                  */
/*======================================================*/
void CUDAStart()
{
  
  SetCurrentDevice(SelectDevice);
	
  AllocateMemCUDA();
  
  CopyHostToDevice();
 
  CUDAInitParams(SimTime);
  	
  startC4(SimTime,trwanie,tmstart);
	
  FreeCUDAMemory();
}


/*======================================================*/
/* Save PACF data to file                               */
/*======================================================*/
void Save_PACF(char *filename,real *Cvv,int trwanie,int tstart)
			 {
			 int nn;
			
			 FILE *dane_induk;


			 dane_induk = fopen(filename,"w");

				for(nn=tstart;nn<trwanie;nn++)
					 {
					  fprintf(dane_induk,"%f %5.10f\n",(real)psdt*nn,Cvv[nn]);
					 }

		fclose(dane_induk);

}

/*======================================================*/
/* Save 4-body PACF                                     */
/*======================================================*/
void inline SaveCvv4(double Ttot,int timestart)
{
 
 filename=new char[255];
 acf_file=new char[255];

 printf("\n=====================================\n");
 sprintf(filename,"fun_%d_%d_",timestart,trwanie);
	
	     strcpy(acf_file, filename);
	     strcat(acf_file,"g4cuda.dat");	
	     printf("Filename: %s\n---------------------------------\n",acf_file);
		 Save_PACF(acf_file,Cvv,trwanie,timestart);
	     	    
}
