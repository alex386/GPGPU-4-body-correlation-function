/*====================================================================*/
/* Global variables                                                   */
/*====================================================================*/
POS* Traj;
float sigmaAr;
size_t Traj_size;
int steptime;
real vol;
real psdt;
real dt;
int N;
int tmstart=0;
int trwanie=0;
int SimTime=0;
long ile_mem=0;
int tt=0;


POS **TR2;
real *Cvv;
real ***TabWsk;

char *filename,*acf_file;

//============================================
float Lx,Ly,Lz,Lx2,Ly2,Lz2,Lxx,Lyy,Lzz;
