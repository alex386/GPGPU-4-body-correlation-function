
/*================================================================*/
/* Initialize table for trajectory                                */
/*================================================================*/

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
/*================================================================*/
/* Initialize table in memory                                     */
/*================================================================*/
int inline InitMemFunc(int SimTime, int NUM)
{
	long ile_mem;
	int p,q;


          if ((TabWsk = new real**[NUM]) == NULL)
		  {
            printf("Out of memory\n");
            exit(1);  /* terminate program if out of memory */
		  }

          for(p=0;p<NUM;p++){
            if((TabWsk[p] = new real*[NUM]) == NULL)
			{
             printf("Out of memory\n");
             exit(1);  /* terminate program if out of memory */
			} 
		  }

          for(p=0;p<NUM;p++){
              for(q=0;q<NUM;q++){
                 if((TabWsk[p][q] = new real[SimTime]) == NULL)
				 {
                  printf("Out of memory\n");
                  exit(1);  /* terminate program if out of memory */
				 }

                memset(TabWsk[p][q],0,(SimTime)*sizeof(real));
			  }
		  }

          ile_mem=(NUM+1)*(NUM+1)*(SimTime+1)*sizeof(real);
		  printf("Memory used: %ld bytes\n",ile_mem);

		  return 0;
}


/*===============================================================*/
/* Gather scattering information                                 */
/*===============================================================*/
void inline GetScatterPoints(int NrPoints)
{

int i,j,t;
POS Ra;
real Rijkwa,InversRijkwa,Rij,writeStan,VersRij;
timeval tv;
double cz1,cz2,cztot;
//------------------------------------------------	

//Initialize memory	
InitMemFunc(NrPoints, N);

// Get current time in miliseconds
gettimeofday (&tv, NULL);
cz1 = (tv.tv_sec) + 0.000001 * tv.tv_usec;

for(t=0;t<NrPoints;t++){	
	
     for(i=0;i<N;i++){
		  
	    for(j=0;j<N;j++){
			TabWsk[i][j][t]=0;
			
          if(i!=j){
          Ra.x = TR2[i][t].x-TR2[j][t].x;
          Ra.y = TR2[i][t].y-TR2[j][t].y;
          Ra.z = TR2[i][t].z-TR2[j][t].z;

//------------------------------------------------------------------------------
// Minimum Image Convention
//-----------------------------------------------------------------------------
              //X-axis
	                if (Ra.x>=Lx)         Ra.x-=Lxx;
	                else if (Ra.x<-Lx)    Ra.x+=Lxx;
              //Y-axis
	                if (Ra.y>=Ly)         Ra.y-=Lyy;
	                else if (Ra.y<-Ly)    Ra.y+=Lyy;
              //Z-axis
	                if (Ra.z>=Lz)         Ra.z-=Lzz;
	                else if (Ra.z<-Lz)   Ra.z+=Lzz;			  
//-----------------------------------------------------------------------------	
		

		  Rijkwa = Ra.x*Ra.x + Ra.y*Ra.y + Ra.z*Ra.z;
		  InversRijkwa=1.0/Rijkwa;
  
		  Rij=sqrt(Rijkwa);

		  VersRij=(vol*Ra.x*Ra.z)/Rij;
		  writeStan=VersRij*InversRijkwa*InversRijkwa;
			
		  TabWsk[i][j][t]=writeStan;
		  }
		
		}
		
	 }
	
	
 }			
gettimeofday (&tv, NULL);
cz2 = (tv.tv_sec) + 0.000001 * tv.tv_usec;
cztot = cz2-cz1;
printf("Anisotropy table creation time: %lf s\n",cztot);// TotalTest: %lf\n",TotalTest);
fflush(stdout);
	
}
/*========================================================================*/
void swap_data(char *filename,real Cvvt,int t)
			 {
			 		 
			 FILE *dane_swap;
			 dane_swap = fopen(filename,"a");				
			 fprintf(dane_swap,"%lf %5.14lf\n",(double)psdt*t,Cvvt);
					
		fclose(dane_swap);

}
/*========================================================================*/
/* 4-body anisotropy polarizability correlation function calculation      */
/*========================================================================*/
 void VACF_C4(int tt,int N,real *Cvv,int ttime,int timestart,int step)
	{
	 char filename[256];
	 int t_max,tau,t,Nm1;
	 real ATA,AVER;
	 register int k,i,j,l;
	 timeval tv;
	 double cz1,cz2,cztot;

	 Nm1=N-1;
	 sprintf(filename,"C4SWAP%dTO%d.swp",timestart,ttime);
	 printf("----------------------------------------\nCalcuating G4 function for %d atoms.\nFrom step %d to step %d every %d step\n----------------------------------------\n",N,timestart,ttime-1,steptime);
	 fflush(stdout);

	int id=0;	
		
    for(t=timestart;t<ttime;t+=step)
		 {
					 
             gettimeofday (&tv, NULL);
             cz1 = (tv.tv_sec) + 0.000001 * tv.tv_usec;
			 
			  AVER=0;
			  t_max = tt - t;
			  for(i=0;i<Nm1;i++)
				  {
				  for(j=i+1;j<N;j++)
					 {
					 for(k=0;k<Nm1;k++)
						{
						 for(l=k+1;l<N;l++)
							 {						 
						      if(i!=k && j!=k && i!=l && j!=l)
							  {
                                ATA=0;
							      for(tau=1;tau<t_max;tau++)
								  {
								   ATA += TabWsk[i][j][tau]*TabWsk[k][l][tau+t];
								  }

							    ATA/=t_max;
								
								id++;  
							    AVER+=ATA;
							  }//warunki
						 }
					 }
				  }
			  }
			  Cvv[t] = AVER;
			  

			 gettimeofday (&tv, NULL);
             cz2 = (tv.tv_sec) + 0.000001 * tv.tv_usec;
			 cztot = cz2-cz1;
			 
			  printf("step=%d t=%f [ps] Cvv=%5.10f time: %lf s\n",t,t*psdt,Cvv[t],cztot);
			  fflush(stdout);
			  swap_data(filename,Cvv[t],t);
        }
		
		}


/*==========================================================*/
/* Save correlation function to file                        */
/*==========================================================*/
void Save_acf(char *filename,real *Cvv,int trwanie,int timestart,int step)
			 {
			 int nn;
			 real czyn;
			
			 FILE *dane_induk;
			 
             czyn=(real)psdt;
                         
			 dane_induk = fopen(filename,"w");

				for(nn=timestart;nn<trwanie;nn++)
					 {
					  fprintf(dane_induk,"%f %5.10f\n",(real)czyn*nn,Cvv[nn]);
					 }

		fclose(dane_induk);

}
/*========================================================================*/
/* Main calculation procedure                                             */
/*========================================================================*/
void inline CalculateScattering(int Points)
{

 filename=new char[255];
 acf_file=new char[255];

 
 sprintf(filename,"funC4_%d_",trwanie);

         

          if ((Cvv = new real[Points]) == NULL)
		  {
 	       printf("\n Cvv-out of memory."); 
           exit(1);
		  }  	   

	 memset(Cvv,0,Points*sizeof(real));
	
	  		  
  /*====================================================*/
  /*   CALCULUS OF C4(t)                                */
  /*----------------------------------------------------*/
	memset(Cvv,0,Points*sizeof(real));
	
      // Calculate 4-body ACF    		
         VACF_C4(Points,N,Cvv,trwanie,tmstart,steptime);
	
	     strcpy(acf_file, filename);
	     strcat(acf_file,"c4.dat");
      // Save data to file
	     Save_acf(acf_file,Cvv,trwanie,tmstart,steptime);
	     printf("---------------------------------------\nCalculation completed\n---------------------------------------\n");
		
  // Free memory
		  delete[] filename;
		  delete[] acf_file;
		  
}

