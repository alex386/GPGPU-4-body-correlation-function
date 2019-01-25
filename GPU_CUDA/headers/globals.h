/*===================================================*/
/* Global variables                                  */
/*===================================================*/
float sigmaAr;

real vol;
int RDO,RUP;
int Asize;
int N;
real dt,psdt;
real T;

int bEnd=0;
int tmstart=0;
int steptime=10;
int ff;

int VM_krok=0;
int trwanie=0;
int SimTime=0;
time_t Czas1,Czas2;
long ile_mem=0;

POS **TR2;
char *filename,*acf_file;

int SelectDevice=0;
float Lx,Ly,Lz,Lx2,Ly2,Lz2,Lxx,Lyy,Lzz;

int ARGON;
int CalcVer;
//============================================
// CUDA globals
//============================================
int n4p;
ushort4* PairTab;
ushort4* P2T;
ushort4* spair;
int pa2Size,pa4Size,parSize;

int sharedSize;
int threadsPerBlock, blocksPerGrid;
int MulTot;
unsigned int hTimer;
double gpuTime,gpuTimeRun,speedUP,TimeInCPU,TimeInGPU;
//-----------------------------------------------
int TPB,xsiz,ysiz;

size_t Traj_size,Cvv_size;
size_t N_size,p2_size;

size_t p4_size;
size_t calc_size;

size_t BetaS;
real* Cvv;
real* ArTab;
real4* Traj;
real* Beta;
// Device memory
real* d_Cvv;
real* d_ArTab;
real* d_Beta;
real4* d_Traj;
ushort4* d_PairTab;
ushort4* d_P2T;

