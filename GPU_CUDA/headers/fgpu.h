/*=============================================================*/
/* GPU procedures                                              */
/* 4-body polarizability anisotropy correlation function       */
/*=============================================================*/




/*====================================================*/
/* Anizotropy for pairs                               */
/*====================================================*/
__global__ void CalcFuncPair(real* Beta,real* ArTab,ushort4* PairTab,int t,unsigned int start)
{
extern __shared__ real4 Rtab[]; 
	
unsigned int UniqueBlockIndex = blockIdx.y * gridDim.x + blockIdx.x;
unsigned int idx = UniqueBlockIndex * blockDim.x + threadIdx.x;
unsigned int idstart=idx+start;

//ArTab[idx]=0;
//__syncthreads();
	
//-----------------------------------------------------
if(idstart>=gN) return;
//----------------------------------------------------
unsigned int t_max,tau,t1,st1,st2;	
unsigned int Tx,Ty,N;
ushort4 P2r;
real ATA;
real Itmax;
real a,b;

	N=gN; //pa2Size number of ij pairs N*(N-1)/2

	t_max = gtt - t;
	Itmax =1/(real)t_max;
	
	P2r=PairTab[idx];

  
				     ATA = 0;
				     for(tau=1;tau<t_max;tau++)
					 { 	
						 st1=N*tau;
						 	 						 
						 Tx=idx + st1;
						 
						 t1=tau+t;
						 
                                                 st2=N*t1;
						 Ty=P2r.y-1 + st2;
						 
						  a=Beta[Tx];
						  b=Beta[Ty];
                        
						
						 ATA += a*b;						 
						
					 }

				    ArTab[idx] += ATA*Itmax;
			
}//G procedure



/*============================================================================*/
/* CUDA reduction procedure                                                   */
/*============================================================================*/
__global__ void Reducto(real* ArTab,real* RedTab)
{
	
extern __shared__ real lTab[];
	
unsigned int UniqueBlockIndex = blockIdx.y * gridDim.x + blockIdx.x;
unsigned int idx = UniqueBlockIndex * blockDim.x + threadIdx.x;

unsigned int tid = threadIdx.x;
	

lTab[tid] = ArTab[idx];
	
__syncthreads();

for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
if (tid < s) {
lTab[tid] += lTab[tid + s];
}
__syncthreads();
}	
	
	if(tid == 0) RedTab[UniqueBlockIndex] = lTab[0];
		
}//Reduction procedure

/*===========================================================================*/
__global__ void ReductoEnd(real* ArTab,real* RedTab,real* Cvv,int t)
{
	
extern __shared__ real lTab[];
	
unsigned int UniqueBlockIndex = blockIdx.y * gridDim.x + blockIdx.x;
unsigned int idx = UniqueBlockIndex * blockDim.x + threadIdx.x;

unsigned int tid = threadIdx.x;
	
lTab[tid] = ArTab[idx];
	
__syncthreads();

for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
if (tid < s) {
lTab[tid] += lTab[tid + s];
}
__syncthreads();
}	

	if(tid == 0){
		Cvv[t] += lTab[0];
	}
		
}//Reduction procedure




/*=========================================================================*/
/* Zero Pair in CUDA                                                       */
/*=========================================================================*/
__global__ void ZeroPairs(real* ArTab,ushort4* PairTab)
{
	
unsigned int UniqueBlockIndex = blockIdx.y * gridDim.x + blockIdx.x;
unsigned int idx = UniqueBlockIndex * blockDim.x + threadIdx.x;


	PairTab[idx].x=0;
    PairTab[idx].y=0;
	PairTab[idx].z=0;
	PairTab[idx].w=0;

    ArTab[idx]=0;
//__syncthreads();	
	__syncthreads();
}



/*=========================================================================*/
/* Beta Pair value calculation in CUDA                                     */
/*=========================================================================*/
__global__ void PairIndex(ushort4* P2, ushort4* PairTab,unsigned int start)
{
	
unsigned int UniqueBlockIndex = blockIdx.y * gridDim.x + blockIdx.x;
unsigned int idx = UniqueBlockIndex * blockDim.x + threadIdx.x;
unsigned int idstart=idx+start;
	
int Nm1,k,l,lk;
ushort Kstart,Lstart,cid;


if(idstart>=gN) return; // horizontal overflow check

Nm1=gARGON-1;
	
    ushort4 P2r;
	P2r=P2[idstart];

	cid=PairTab[idx].y;
	Kstart=PairTab[idx].z;
	Lstart=PairTab[idx].w;

	PairTab[idx].x=0;
    PairTab[idx].y=0;
	PairTab[idx].z=0;
	PairTab[idx].w=0;

	
	       for(k=Kstart;k<Nm1;k++)
						{
							if(k==Kstart)
							lk=Lstart+1;
							else
							lk=k+1;
							
						 for(l=lk;l<gARGON;l++)
							 {	
								 cid++; 
						        if(P2r.x!=k && P2r.y!=k && P2r.x!=l && P2r.y!=l)
							    {

								    PairTab[idx].x=idx;//P2r.x;
						            PairTab[idx].y=cid;//P2r.y;
						            PairTab[idx].z=k;
						            PairTab[idx].w=l;
									__syncthreads();
									return;							      							  							  								  
					            }
								
								
							 }	
							
						 }
					                       				 
					 
__syncthreads();				
}//G4 index procedure







