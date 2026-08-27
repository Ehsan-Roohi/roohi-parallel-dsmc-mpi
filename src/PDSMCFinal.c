// DSMC.cpp : Defines the entry point for the console application.
//

#include "mpi.h"
#include <math.h>
#include <stdio.h>

#include <time.h>
     


//COMMON /CONST / 
const double PI=3.141592654, SPI=1.772454, BOLTZ=1.380622E-23;
//const int MNC=7501, MNSC=30001, MNSP=1, MNSG=1, MNSE=500;
#define MNM 9000000
#define MNC 7501
#define MNSC 30001
#define MNSP 1
#define MNSG 1
#define MNSE 500
struct MolPro{
	double X,Y;
	double PV[4];	//[1,2,3] : for velocity components of each particle, [0] is not used
	double PR;	
	int IPL,IR,IPS;
	int Part,IX,JY;	//i-index and j-index of the cell that this molecule belongs to
}Particle[MNM];
struct CellPro{
	double CG[7];	//1,2:x-comp., 4,5:y-comp., 3=dx; 6;dy;
	double CC;		//the area of the cell
}Cell[MNC];
struct SubDomain{
	double XStart,XEnd;
	int IStart,IEnd;
}Partition[128];

int NM;
int IC[3][MNC][2],ISC[MNSC],ISCG[3][MNSC][2],IG[2][2],NCX,NCY,IFCX,IFCY;
double CC[MNC],CG[7][MNC],CCG[3][MNC][2][2];
//COMMON /GAS   /
double SP[6][2],SPM[7][2][2];
int ISP[2];
//COMMON /GASR  / 
double SPR[4][2][2],CT[MNC];
int ISPR[4][2];
//COMMON /SAMP2 / 
double COL[2][2],SELT,MOVT,SEPT,CS[8][MNC][2],TIME,FND,FTMP,TIMI,FSP[2],VFX,VFY;
int NCOL,NPR,NSMP,ISPD;
//COMMON /SAMPR / 
double CSR[MNC][2];
//COMMON /SAMPS / 
double CSS[10][300][2];
//COMMON /COMP  / 
double FNUM,DTM;
int NIS,NSP,NPS,NPT;
//COMMON /GEOM2 / 
int NSCX,NSCY,IB[5],ISURF[3],LIMS[3][4],IIS,ISG;
double CB[5],TSURF[3],BMR[5][2],CW,FW,CH,FH,ALPI[3],ALPN[3],ALPT[3];
//COMMON /PDATA/ 
double PIN,POUT;
//COMMON /KN/ mfp,kn
double mfp,kn;
int IPLOT;
 
double GAMA,GASR;
//variables for parallel	
double XStP,XEnP, dxd;
int IStP,IEnP;
int partition[MNC],MNCP;
//partition[k] determines that the cell [k] belong to what pocessor
//MNMP determine maximum number of molecules in the partiotin [n]
//MNCP determine maximum number of cells in the partiotin [n]

//declaration
//of all the subprograms
void INIT2();
void startcode();
void DATA();
double GAM(double );
void mesh_set();
double RF(int);
void plot_grid();
void RVELC(double* ,double* ,double );
void SROT(double* ,double ,double );
void REMOVE(int );
void REFLECT2(int ,int ,int ,double ,double ,int );
void SAMPI2();
void SAMPLE2();
void MOVE2();
void INDEXM();
void COLLMR();
void ENTER2();
void PROPERTIES(int );
void OUT2();

int  namelen,myid, numprocs;

int main(int argc,char *argv[])
{	
//	FILE fp0;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double startwtime = 0.0, endwtime;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Get_processor_name(processor_name,&namelen);

    	fprintf(stdout,"Process %d of %d on %s\n",
	    myid, numprocs, processor_name);
	
	//array allocation for subdomains,
						// Partitioning \\;
	startcode();
	int n,ndebug;
	dxd=(CB[2]-CB[1])/numprocs;
	int i,j,k;
//	for(n=0;n<numprocs;n++){
	
	for(j=1;j<=NCY;j++){
		for(i=IStP;i<=IEnP;i++){
			k=(j-1)*NCX+i;
			partition[k]=myid;
		}
	}
	//if(myid==0){	
	if(myid==0){
		startwtime = MPI_Wtime();
	}
	int J,I;
	while (NPR<NPT){
		NPR=NPR+1;

		clock_t start, end;
     		double cpu_time_used;
		start = 0; //clock();

		if(NPR<=NPS) SAMPI2();     
		for( J=1;J<=NSP;J++){
			for( I=1;I<=NIS;I++){
				TIME=TIME+DTM;
				if(myid==0)				
					printf("Procs %d\tStep %d\t%d\t%d\t%d Mols %10d\tColls\t%g\n",myid,I,J,NPR,NM,NCOL,TIME);
				MOVE2();
				INDEXM();
				COLLMR();
			}
			SAMPLE2();
		}
		
		end = clock();
		cpu_time_used = ((double) (end - start)); ///CLOCKS_PER_SEC;		
		
		IPLOT=IPLOT+1;
		if((NPR-(NPR/10)*10)==0){
			OUT2(); 
		}
		if(myid==0) 	printf("\t\t\t< NPR=%d >\n",NPR);
		if(myid==0) 
		{
			//printf("\t\t\t< start=%d >\n",start/10000);
			//printf("\t\t\t< end=%d >\n",end/10000);
			printf("\t\t\t< cpu_time_used=%d >\n",(end - start)/10000);
		} 
		//if((NPR-(NPR/10000)*10000)==0){
			if(myid==0){
				//printf("id=%d,\tthe output file is created ",myid);
				endwtime=MPI_Wtime();
				printf("\t\t\t<CPU time=%f\n",endwtime-startwtime);
				//OUT2();		
				//scanf("%d",&ndebug);
			}
			
		//}
	}
	MPI_Finalize();
	return 0;
}


void INIT2()
{
    FND=0.0;
    FTMP=273.0;
    VFX=0.0;
    VFY=0.0;
    IFCX=0;
    IFCY=0;
    ALPI[1]=-1.0;
    ALPI[2]=-1.0;
	int N;	
	for( N=1;N<=4;N++){
		IB[N]=3;
        ISP[1]=1;
        FSP[1]=0.0;
//        BME[N][0]=0.0;
        BMR[N][1]=0.0;
	}
	DATA();

	int M;
    if(MNSP==1) ISPD=0;
	for( N=1;N<=MNSP;N++){
		for( M=1;M<=MNSP;M++){
			if((ISPR[3][N]==0)&&(M!=N)){
				SPR[1][N][M]=SPR[1][N][N];
				SPR[2][N][M]=SPR[2][N][N];
				SPR[3][N][M]=SPR[3][N][N];
			}
			if((ISPD==0)||(N==M)){
				SPM[1][N][M]=0.25*PI*pow((SP[1][N]+SP[1][M]),2.0);
				//*--the collision cross section is assumed to be given by eqn (1.35)
				SPM[2][N][M]=0.5*(SP[2][N]+SP[2][M]);
				SPM[3][N][M]=0.5*(SP[3][N]+SP[3][M]);
				SPM[4][N][M]=0.5*(SP[4][N]+SP[4][M]);
				//*--mean values are used for ISPD=0
			}else{
				SPM[1][N][M]=PI*pow(SPM[1][N][M],2.0);
				//*--the cross-collision diameter is converted to the cross-section
			}
			SPM[5][N][M]=(SP[5][N]/(SP[5][N]+SP[5][M]))*SP[5][M];
			//*--the reduced mass is defined in eqn (2.7)
			SPM[6][N][M]=GAM(2.5-SPM[3][N][M]);
		}
	}
    	TIME=0.0;
    	NM=0;
	NPR=0;
    	NCOL=0;
    	MOVT=0.0;
    	SELT=0.0;
    	SEPT=0.0;

	for( M=1;M<=MNSP;M++){
		for( N=1;N<=MNSP;N++){
			COL[M][N]=0.0;
		}
	}
//data set for partitions
/*	IStP=floor(myid*NCX/numprocs)+1;
	IEnP=floor((myid+1)*NCX/numprocs);
*/
	IStP=(myid)*NCX/numprocs+1;
	IEnP=(myid+1)*NCX/numprocs;

	mesh_set();


/*	XStP=CG[1][IStP];
	XEnP=CG[2][IEnP];
	int i;
	for(i=0;i<numprocs;i++){
		Partition[i].IStart=i*NCX/numprocs+1;
		Partition[i].IEnd=(i+1)*NCX/numprocs;
	}

	int prcnum;
	if(myid==0){
		printf("IStP=%d,\tIEnP=%d,\nXStP=%g,\tXEnP=%g\n",IStP,IEnP,XStP,XEnP);
		for(prcnum=0;prcnum<numprocs;prcnum++){
			printf("Processor %d\nIStart=%d,\tIEnd=%d,\nXStart=%g,\tXEnd=%g\n",
				prcnum,Partition[prcnum].IStart,Partition[prcnum].IEnd,Partition[prcnum].XStart,Partition[prcnum].XEnd);
		}
		printf("\nCB[0]=%g,\tCB[2]=%g\n",CB[1],CB[2]);
		scanf("%d",&prcnum);
	}*/

dxd=(CB[2]-CB[1])/(double)(numprocs);
	XStP=CB[1]+(double)(myid)*dxd;	//CG[1][IStP]; //
	//if(XStP!=CG[1][IStP])	printf("Error in XStP, Procs%d, XStP=%g, CG=%g\n",myid,XStP,CG[1][IStP]);
	XEnP=CB[1]+(double)(myid+1)*dxd;	//CG[2][IEnP];
	int prcnum;
	for(prcnum=0;prcnum<numprocs;prcnum++){
		Partition[prcnum].XStart=CB[1]+(double)(prcnum)*dxd;
		Partition[prcnum].XEnd=CB[1]+(double)(prcnum+1)*dxd;
		if(myid==0)printf("*^^^^^^^%g\t%g\t",Partition[prcnum].XStart,Partition[prcnum].XEnd);
		Partition[prcnum].IStart=floor(prcnum*NCX/numprocs)+1;
		Partition[prcnum].IEnd=floor((prcnum+1)*NCX/numprocs);
	}
	printf("\n");
	//if(XEnP!=CG[2][IEnP])	printf("Error in XEnP, Procs%d, XEnP=%g, CG=%g\n",myid,XEnP,CG[2][IEnP]);
	//printf("%d\t%d\t%g\t%g\n",IStP,IEnP,XStP,XEnP);
	//if((myid==(numprocs-1))&&(XEnP!=CB[2]))printf("Error in XEnP,%g\t%g\n",XEnP,CB[2]);
//data set for partitions, end

	
	double REM,VMP,A;
	int MM,NCOLM,NROW,L,K;
	int i,j;
	if((IIS>0)&&(ISG>0)){
	//*--if IIS=1 generate initial gas with temperature FTMP
		for( L=1;L<=MNSP;L++){
			REM=0;
			if(IIS==1) VMP=sqrt(2.0*BOLTZ*FTMP/SP[5][L]);
			//*--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
			//for( N=1;N<MNC;N++){
			for(j=1;j<=NCY;j++){
			    for(i=IStP;i<=IEnP;i++){
				N=(j-1)*NCX+i;
				A=FND*CC[N]*FSP[L]/FNUM+REM;
				//*--A is the number of simulated molecules of species L in cell N to
				//*--simulate the required concentrations at a total number density of FND
				if(N<MNC){
					MM=floor(A);
					REM=(A-MM);
					//if(myid==1)printf("%d\t%d\t%d\n",IStP,IEnP,MM);
					//*--the remainder REM is carried forward to the next cell
				}else{
					MM=floor(A);
				}
				if(MM>0){
					for( M=1;M<=MM;M++){
						if(NM<MNM){
							//*--round-off error could have taken NM to MNM+1
							NM=NM+1;
							Particle[NM].IPS=L;
							Particle[NM].X=CG[1][N]+RF(0)*(CG[2][N]-CG[1][N]);
				  			NCOLM=floor((Particle[NM].X-CG[1][N])*(NSCX-0.001)/CG[3][N])+1;

							Particle[NM].Y=CG[4][N]+RF(0)*(CG[5][N]-CG[4][N]);
				            		NROW=floor((Particle[NM].Y-CG[4][N])*(NSCY-0.001)/CG[6][N])+1;
                    
							Particle[NM].IPL=(N-1)*NSCX*NSCY+(NROW-1)*NSCX+NCOLM;
							//*--species, position, and sub-cell number have been set
							for( K=1;K<=3; K++){
								RVELC(&Particle[NM].PV[K],&A,VMP);
							}
							Particle[NM].PV[1]=Particle[NM].PV[1]+VFX;
							Particle[NM].PV[2]=Particle[NM].PV[2]+VFY;
							//*--velocity components have been set
							//*--set the rotational energy
							if(ISPR[1][L]>0) SROT(&Particle[NM].PR,FTMP,ISPR[1][L]);
						}
					} //end of for
				}
			    }//end of if(partition)
			}
		}
		printf("Number of initial molecules in Partiion %d := %d \n",myid,NM);
	}

}
void startcode()
{
    INIT2();
	SAMPI2();
}
void DATA()
{
    NCX=240;
    NCY=30;

    double AR=100.0;

    NSCX=2;
    NSCY=2;
    IIS=1;
    ISG=1;
    FTMP=300.0;
	PIN=2.0E5;
	POUT=PIN/2.0;
    FND=PIN/(BOLTZ*FTMP);
	//--FND is the number densty
    VFX=10.0;
    VFY=0.0;
    FSP[1]=1.0;
	//--FSP(N) is the number fraction of species N
    FNUM=0.175E10*20;
	//--FNUM  is the number of real molecules represented by a simulated mol.      
    SP[1][1]=2.33e-10;	//4.17E-10;
    SP[2][1]=273.0;
    SP[3][1]=0.66;	//0.74;
    SP[4][1]=1.0;
    SP[5][1]=6.65e-27;	//4.65E-26;

//      SP(1,1)=2.33E-10
//      SP(2,1)=273.
//      SP(3,1)=0.66
//      SP(4,1)=1.0
//      SP(5,1)=6.65E-27
	//--SP(1,N) is the molecular diameter of species N
	//--SP(2,N) is the reference temperature
	//--SP(3,N) is the viscosity-temperatire index
	//--SP(4,N) is the reciprocal of the VSS scattering parameter
	//--SP(5,N) is the molecular mass of species N
	mfp=1.0/(sqrt(2.0)*(PI*SP[1][1]*SP[1][1])*FND*
		pow((SP[2][1]/FTMP),(SP[3][1]-0.5)));
	//mfp=1.0/(SQRT(2.0)*(pi*sp(1,1)**2)*FND*(SP(2,1)/FTMP)
     	//	**(SP(3,1)-0.5))
	kn=0.02;
	CB[1]=0.0;
    CB[3]=0.0;
    CB[4]=mfp/(2.0*kn);	//	!0.45e-6*2.49d0 mfp/(2.*kn*HRATIO)
    CB[2]=AR*2.0*CB[4];	//	!/(1.-LRATIO)
	//*--the simulated region is from x=XB(1) to x=XB(2)	
      
	IB[1]=1; //!1:INFLOW,2:SYMMETRY
    IB[2]=1;
    IB[3]=2; 
    IB[4]=2;
	kn=SP[5][1]/(sqrt(2.0)*(PI*SP[1][1]*SP[1][1])*2.0*CB[4]);
    if(myid==0)
	printf("Kn=%g\n",kn);

    ISURF[1]=1;	//  !Whether there is a wall at sym. boundary, 1:positive y
	ISURF[2]=0;	//!3 FOR BUUFER !0:No wall

    LIMS[1][1]=1;
    LIMS[1][2]=1;
    LIMS[1][3]=NCX;
    LIMS[2][1]=0;
    LIMS[2][2]=0;
    LIMS[2][3]=0;

    TSURF[1]=300.0;
    TSURF[2]=-1.0;

	DTM=0.5E-13;
	//*--DTM is the time step
    ISPR[1][1]=2;
    SPR[1][1][1]=5.0;
    ISPR[2][1]=0;
    ISP[1]=1;
	//*--ISPR(1,N) is the number of degrees of freedom of species N
	//*--SPR(1,N,K) is the constant in the polynomial for the rotational
	//*--relaxation collision number of species N in collision with species K
	//*--ISPR(2,N) is 0,1 for constant, polynomial for collision number
    NIS=2;  
	//*--NIS is the number of time steps between samples
    NSP=10;
	//*--NSP is the number of samples between restart and output file updates
    NPS=30;
	//*--NPS is the number of updates to reach assumed steady flow
    NPT=800000;
	printf("%f\t%f\t%f\t%f\t%f\n",kn,PIN,POUT,CB[4],CB[2]);
	//*--NPT is the number of file updates to STOP
}
void RVELC(double *U,double *V,double VMP)
 {
	//*--generates two random velocity components U an V in an equilibrium
	//*--gas with most probable speed VMP  (based on eqns (C10) and (C12))

    double A,B;
	A=sqrt(-log(RF(0)));
    B=6.283185308*RF(0);
    *U=A*sin(B)*VMP;
    *V=A*cos(B)*VMP;
 }
void SROT(double *PR,double TEMP,double IDF)
{
	//*--selects a typical equuilibrium value of the rotational energy PR at
	//*----the temperature TEMP in a gas with IDF rotl. deg. of f.
    
	if(IDF==2){
		*PR=-log(RF(0))*BOLTZ*TEMP;
		//*--for 2 degrees of freedom, the sampling is directly from eqn (11.22)
	}else{
		//*--otherwise apply the acceptance-rejection method to eqn (11.23)
        double A, B, ERM;
		A=0.5*IDF-1.0;
Re2:    ERM=RF(0)*10.0;
		//*--the cut-off internal energy is 10 kT
        B=(pow((ERM/A),A))*exp(A-ERM);
        if (B<RF(0)) goto Re2;
        *PR=ERM*BOLTZ*TEMP;
	}
}
double GAM(double X)
{
	//*   GAM.FOR
	//*--calculates the Gamma function of X.
	double A,Y,Result;
	A=1.0;
    Y=X;
	if(Y<1.0){
		A=A/Y;
	}else{
		re:Y=Y-1;
		if(Y>=1.0){
			A=A*Y;
			goto re;
		}
	}
	Result=A*(1.0-0.5748646*Y+0.9512363*pow(Y,2.0)-0.6998588*pow(Y,3.0)+
		   0.4245549*pow(Y,4.0)-0.1010678*pow(Y,5.0));
    return Result;
}

void mesh_set()
{
	FW=CB[2]-CB[1];
    FH=CB[4]-CB[3];
    CG[1][1]=CB[1];
    CW=FW/NCX;
	//*--CW is the uniform cell width
    CG[4][1]=CB[3];
    CH=FH/NCY;
	//*--CH is the uniform cell height
	int N,M,MY,MX,L,K;
	for( MY=1;MY<=NCY;MY++){
		for( MX=1;MX<=NCX;MX++){	//IStP start index, IEnP end index
			M=(MY-1)*NCX+MX;
			//*--M is the cell number
			CT[M]=FTMP;
			//*--the macroscopic temperature is set to the freestream temperature
			//*--set the x coordinates
			if(MX==1) CG[1][M]=CG[1][1];
			if(MX >1) CG[1][M]=CG[2][M-1];
            		CG[2][M]=CG[1][M]+CW;
			CG[3][M]=CG[2][M]-CG[1][M];
			//*--set the y coordinates
			if(MY==1) CG[4][M]=CG[1][1];
			if((MY>1)&&(MX==1)) CG[4][M]=CG[5][M-1];
			if((MY>1)&&(MX >1)) CG[4][M]=CG[4][M-1];
		  	CG[5][M]=CG[4][M]+CH;
			CG[6][M]=CG[5][M]-CG[4][M];
			CC[M]=CG[3][M]*CG[6][M];
			for( L=1;L<=MNSG;L++){
				for( K=1;K<=MNSG;K++){
					CCG[2][M][L][K]=RF(0);
					CCG[1][M][L][K]=SPM[1][1][1]*300.0*sqrt(FTMP/300.0);
				}
			}
			//*--the maximum value of the (rel. speed)*(cross-section) is set to a
			//*--reasonable, but low, initial value and will be increased as necessary
		}
	}
	//int L;
	for( N=1;N<MNC;N++){
		for( M=1;M<=NSCY;M++){
			for( K=1;K<=NSCX;K++){
				L=(N-1)*NSCX*NSCY+(M-1)*NSCX+K;
				ISC[L]=N;
			}
		}
	}
	//printf("%d\n",ISC[18001]);
}
double RF(int IDUM)
{
	//*----IDUM will generally be 0, but negative values may be used to
	//*------re-initialize the seed
    static int MA[56],INEXT,INEXTP,IFF=0;
    const int MBIG=1000000000, MSEED=161803398, MZ=0;
	int MJ, MK, II;
    double RanF;

	const double FAC=1.0E-9;
	if((IDUM<0)||(IFF==0)){
		IFF=1;
        MJ=MSEED-fabs(IDUM);
        MJ=MJ%MBIG;
        MA[55]=MJ;
        MK=1;
	int I,K;
		for( I=1;I<=54;I++){
			II=(21*I)%55;
			MA[II]=MK;
			MK=MJ-MK;
			if(MK<MZ) MK=MK+MBIG;
			MJ=MA[II];
		}
		for( K=1;K<=4; K++){
			for( I=1;I<=55;I++){
				MA[I]=MA[I]-MA[1+((I+30)%55)];
				if(MA[I]<MZ) MA[I]=MA[I]+MBIG;
			}
		}
        INEXT=0;
        INEXTP=31;
	}
Re3:  INEXT=INEXT+1;
      if(INEXT==56) INEXT=1;
      INEXTP=INEXTP+1;
      if(INEXTP==56) INEXTP=1;
      MJ=MA[INEXT]-MA[INEXTP];
      if(MJ<MZ) MJ=MJ+MBIG;
      MA[INEXT]=MJ;
	  RanF=MJ*FAC;
      if((RanF>1.0E-8)&&(RanF<0.99999999)) return RanF;
      goto Re3;
	return RanF;
}
void plot_grid()
{
	FILE *fp;
	fp=fopen("Grid.plt","wt");
	fprintf(fp,"Variables=X,Y\nZone I=%d,\tJ=%d,\tF=Point\n",NCY,NCX);
	int M,i,j;
	double xm,ym;
	for( i=0;i<NCX; i++){
		for( j=0;j<NCY; j++){
			M=j*NCX+i;
			xm=0.5*(CG[0][M]+CG[1][M]);
			ym=0.5*(CG[3][M]+CG[4][M]);
			fprintf(fp,"%g\t%g\n",xm,ym);
		}
	}
	fclose(fp);
}
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\/\/\/\/\/\\/\/\\/\/\\/\/\\/\/\/\\/\/
void MOVE2()
{
	//*--the NM molecules are moved over the time interval DTM
	int IFT,N,MSC,MC,L1,L2,L3,MCS,K,MCX,MCY,MSCX,MSCY;
	double XI,X,XC,YI,Y,AT,YS,DX,DY,XSU,XSD;
	IFT=-1;  // Step w.o entering molecules
	//--a negative IFT indicates that molecules have not entered at this step
	N=0;
	int ndebug,nremoved,PartIndex,prc,reflected;
	int KS;
	nremoved=0;

Re4:    N=N+1;
	reflected=0;
	if (N<=NM){
	    if(IFT<0) AT=DTM;
	    if(IFT>0) AT=RF(0)*DTM;
	    //*--the time step is a random fraction of DTM for entering molecules
Re5:        MOVT=MOVT+1;
            MSC=Particle[N].IPL;
            MC=ISC[MSC];
	    //*--MC is the initial cell number
            XI=Particle[N].X;
	    if(((XI+0.00001*CG[3][1])<CB[1])||((XI-0.00001*CG[3][MNC-1])>CB[2])){
		printf(" MOL %d  X COORD OUTSIDE FLOW %g\n",N,XI);
		REMOVE(N);
		goto Re4;
	    }
	    /*if(reflected==0){
		if((XI+0.00001*CG[3][IStP]<XStP)||(XI-0.00001*CG[3][IEnP]>XEnP)){
			//printf(" MOL %d  X COORD OUTSIDE Partition %d,\t%g\t%g\t%g\n",N,myid,XI,XStP,XEnP);
			REMOVE(N);
			goto Re4;
		}
	    }*/
            YI=Particle[N].Y;
	    if(((YI+0.00001*CG[6][0])<CB[3])||((YI-0.00001*CG[6][MNC-1])>CB[4])){
		printf(" MOL %d  Y COORD OUTSIDE FLOW %g\n",N,YI);
		REMOVE(N);
		goto Re4;
	    }
	    DX=Particle[N].PV[1]*AT;
            DY=Particle[N].PV[2]*AT;
            X=XI+DX;
            Y=YI+DY;
    	    for( KS=1;KS<=2;KS++){
		//*--check the surfaces
		if(ISURF[KS]>0){
			if((ISURF[KS]==0)||(ISURF[KS]==1)){
				L1=LIMS[KS][1];
				if(L1<=NCY){
					YS=CG[4][(L1-1)*NCX+1];
				}else{
					YS=CB[4];
					L1=L1-1;
				}
				if(((ISURF[KS]==1)&&((YI>YS)&&(Y<YS)))||
				   ((ISURF[KS]==2)&&((YI>YS)&&(Y>YS)))){
					XC=XI+(YS-YI)*DX/DY;
					if((XC<=CB[1])&&(IB[1]==2)){	//!IB(1).EQ.2 SYMMETRY BOUNDARY
						XC=2.0*CB[1]-XC;
						Particle[N].PV[1]=-Particle[N].PV[1];
					}
					if((XC>=CB[2])&&(IB[2]==2)){	//!IB(1).EQ.2 SYMMETRY BOUNDARY
						XC=2.0*CB[2]-XC;
						Particle[N].PV[1]=-Particle[N].PV[1];
					}
					L2=LIMS[KS][2];
					L3=LIMS[KS][3];
					XSU=CG[1][L2];
					XSD=CG[2][L3];
					if((XC>XSU)&&(XC<XSD)){		//!WALL BOUNDARY
						//*--molecule collides with surface at XC
						MC=floor((XC-CB[1])/CW+0.99999);
						if(MC<1) MC=1;
						if(MC>NCX) MC=NCX;
						MCS=MC-(L2-1);
						if(ISURF[KS]==1) MC=MC+(L1-1)*NCX;
						if(ISURF[KS]==2){
							MC=MC+(L1-1)*NCX;			
							MCS=MC-L2;
						}
						//*--MC is the cell number for the reflected molecule
						if(KS==2) MCS=MCS+LIMS[1][3]-LIMS[1][2]+1;
						//*--MCS is the code number of the surface element
						AT=AT*(Y-YS)/DY;
						REFLECT2(N,KS,MCS,XC,YS,MC);
						reflected=1;
						goto Re5;
					}//end of if(XC)
				}//end of if(ISURF)
			}//end of if(ISURF)
		}//end of if ISURF
	    }//end of for(KS)
	    //if((myid==0)||(myid==(numprocs-1))){
	    	if((X<CB[1])||(X>CB[2])){
			if(X<CB[1]) K=1;
			if(X>CB[2]) K=2;
			//*--intersection with boundary K
			if(IB[K]==2){
				//*--specular reflection from the boundary (eqn (11.7)) !symmetry
				X=2.0*CB[K]-X;
				Particle[N].PV[1]=-Particle[N].PV[1];
			}else{
				//*--molecule leaves flow
				//printf("%d\t%g\t%g\t%g\t%d\n",myid,X,CB[1],CB[2],N);
				//scanf("%d",&ndebug);
				nremoved++;
				REMOVE(N);
				goto Re4;
			}
		}
	    //}
	    if((Y<CB[3])||(Y>CB[4])){
		if(Y<CB[3]) K=3;
		if(Y>CB[4]) K=4;
		//*--intersection with boundary K
		if(IB[K]==2){
			//*--specular reflection from the boundary (eqn (11.7))
			Y=2.0*CB[K]-Y;
			Particle[N].PV[2]=-Particle[N].PV[2];
		}else{
			//*--molecule leaves flow
			REMOVE(N);
			printf("Remove %d \n",myid);
			goto Re4;
		}
	    }
	    if((X<CG[1][MC])||(X>CG[2][MC])||(Y<CG[4][MC])||(Y>CG[5][MC])){
		//*--the molecule has moved from the initial cell
		MCX=floor((X-CB[1])/CW+0.99999);
		if(MCX<1) MCX=1;
		if(MCX>NCX) MCX=NCX;
		//*--MCX is the new cell column (note avoidance of round-off error)
        	MCY=floor((Y-CB[3])/CH+0.99999);
		if(MCY<1) MCY=1;
		if(MCY>NCY) MCY=NCY;
		//*--MCY is the new cell row (note avoidance of round-off error)
		MC=(MCY-1)*NCX+MCX;
	    }
            MSCX=floor(((X-CG[1][MC])/CG[3][MC])*(NSCX-0.001)+1);
	    MSCY=floor(((Y-CG[4][MC])/CG[6][MC])*(NSCY-0.001)+1);
            MSC=(MSCY-1)*NSCX+MSCX+NSCX*NSCY*(MC-1);
	    //*--MSC is the new sub-cell number
            if(MSC<1) MSC=1;
	    if(MSC>(MNSC-1)) MSC=MNSC-1;
            Particle[N].IPL=MSC;
	    Particle[N].X=X;
	    Particle[N].Y=Y; 
            goto Re4;
				//end of if (N<=NM)
	}else if((IFT<0)){	//clculating entering molecules
		if((myid==0)||(myid==(numprocs-1))){
		      	IFT=1;
			//*--new molecules enter
			if(NSMP!=0) ENTER2();
			N=N-1;
			//printf("N=%d,\tmyid=%d\n",N,myid);
			goto Re4;
		}
	}
	//scanf("%d",&ndebug);
}
void REFLECT2(int N,int KS,int K,double XC,double YC,int MC)
{
	//*--reflection of molecule N from surface KS, element K,
	//*----location XC,YC, cell MC
	int L;
	double VMP,VNI,UPI,WPI,ANG,VPI,ALPHAN,ALPHAT,
		ALPHAI,OM,CTH,X,A,R,TH,UM,VN,VP,WP;
	L=Particle[N].IPS;
	//*--sample the surface properies due to the incident molecules
    	CSS[1][K][L]=CSS[1][K][L]+1.0;
	if(ISURF[KS]==1){
		CSS[2][K][L]=CSS[2][K][L]-SP[5][L]*Particle[N].PV[2];
        	CSS[4][K][L]=CSS[4][K][L]+SP[5][L]*Particle[N].PV[1];
	}
	if(ISURF[KS]==2){
        	CSS[2][K][L]=CSS[2][K][L]+SP[5][L]*Particle[N].PV[2];
	        CSS[4][K][L]=CSS[4][K][L]+SP[5][L]*Particle[N].PV[1];
	}
	if(ISURF[KS]==3){
        	CSS[2][K][L]=CSS[2][K][L]-SP[5][L]*Particle[N].PV[1];
        	CSS[4][K][L]=CSS[4][K][L]+SP[5][L]*Particle[N].PV[2];
	}
	if(ISURF[KS]==4){
        	CSS[2][K][L]=CSS[2][K][L]+SP[5][L]*Particle[N].PV[1];
        	CSS[4][K][L]=CSS[4][K][L]+SP[5][L]*Particle[N].PV[2];
	}
	CSS[5][K][L]=CSS[5][K][L]+0.5*SP[5][L]*(pow(Particle[N].PV[1],2.0)+
			pow(Particle[N].PV[2],2.0)+pow(Particle[N].PV[3],2.0));
    	CSS[7][K][L]=CSS[7][K][L]+Particle[N].PR;

	if(TSURF[KS]<0.0){
		//*--specular reflection
        	if((ISURF[KS]==1)||(ISURF[KS]==2)) Particle[N].PV[2]=-Particle[N].PV[2];
        	if((ISURF[KS]==3)||(ISURF[KS]==4)) Particle[N].PV[1]=-Particle[N].PV[1];
	}else if(ALPI[KS]<0.0){
		//*--diffuse reflection
        	VMP=sqrt(2.0*BOLTZ*TSURF[KS]/SP[5][L]);
		//*--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
		if(ISURF[KS]==1){
			Particle[N].PV[2]=sqrt(-log(RF(0)))*VMP;
			RVELC(&Particle[N].PV[1],&Particle[N].PV[3],VMP);
		}
		if(ISURF[KS]==2){
			Particle[N].PV[2]=-sqrt(-log(RF(0)))*VMP;
			RVELC(&Particle[N].PV[1],&Particle[N].PV[3],VMP);
		}
		if(ISURF[KS]==3){
			Particle[N].PV[1]=sqrt(-log(RF(0)))*VMP;
			RVELC(&Particle[N].PV[2],&Particle[N].PV[3],VMP);
		}
		if(ISURF[KS]==4){
			Particle[N].PV[1]=-sqrt(-log(RF(0)))*VMP;
			RVELC(&Particle[N].PV[2],&Particle[N].PV[3],VMP);
		}
		//*--the normal velocity component has been generated (eqn(12.3))
		//*--a single call of RVELC generates the two tangential vel. components
        	if(ISPR[1][L]>0) SROT(&Particle[N].PR,TSURF[KS],ISPR[1][L]);
	}else if(ALPI[KS]>=0){
		//*--Cercignani-Lampis-Lord reflection model
        	VMP=sqrt(2.0*BOLTZ*TSURF[KS]/SP[5][L]);
		//*--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
		if((ISURF[KS]==1)||(ISURF[KS]==2)){
			if(ISURF[KS]==1) VNI=-Particle[N].PV[2]/VMP;
			if(ISURF[KS]==2) VNI=Particle[N].PV[2]/VMP;
			UPI=Particle[N].PV[1]/VMP;
		}
		if((ISURF[KS]==3)||(ISURF[KS]==4)){
			if(ISURF[KS]==3) VNI=-Particle[N].PV[1]/VMP;
			if(ISURF[KS]==4) VNI=Particle[N].PV[1]/VMP;
			UPI=Particle[N].PV[2]/VMP;
		}
        	WPI=Particle[N].PV[3]/VMP;
	        ANG=atan2(WPI,UPI);
        	VPI=sqrt(UPI*UPI+WPI*WPI);
		//*--VNI is the normalized incident normal vel. component (always +ve)
		//*--VPI is the normalized incident tangential vel. comp. in int. plane
		//*--ANG is the angle between the interaction plane and the x or y axis
		//*--first the normal component
        	ALPHAN=ALPN[KS];
        	R=sqrt(-ALPHAN*log(RF(0)));
        	TH=2.0*PI*RF(0);
        	UM=sqrt(1.0-ALPHAN)*VNI;
        	VN=sqrt(R*R+UM*UM+2.0*R*UM*cos(TH));
		//*--VN is the normalized magnitude of the reflected normal vel. comp.
		//*----from eqns (14.3)
		//*--then the tangential component
        	ALPHAT=ALPT[KS]*(2.0-ALPT[KS]);
        	R=sqrt(-ALPHAT*log(RF(0)));
        	TH=2.0*PI*RF(0);
        	UM=sqrt(1.0-ALPHAT)*VPI;
        	VP=UM+R*cos(TH);
        	WP=R*sin(TH);
		//*--VP,WP are the normalized reflected tangential vel. components in and
		//*----normal to the interaction plane, from eqns (14.4) and (14.5)
		if((ISURF[KS]==1)||(ISURF[KS]==2)){
			if(ISURF[KS]==1) Particle[N].PV[2]=VN*VMP;
			if(ISURF[KS]==2) Particle[N].PV[2]=-VN*VMP;
			Particle[N].PV[1]=(VP*cos(ANG)-WP*sin(ANG))*VMP;
		}
		if((ISURF[KS]==3)||(ISURF[KS]==4)){
			if(ISURF[KS]==3) Particle[N].PV[1]=VN*VMP;
			if(ISURF[KS]==4) Particle[N].PV[1]=-VN*VMP;
			Particle[N].PV[2]=(VP*cos(ANG)-WP*sin(ANG))*VMP;
		}
        	Particle[N].PV[3]=(VP*sin(ANG)+WP*cos(ANG))*VMP;
		if(ISPR[1][L]>0){
			//*--set CLL rotational energy by analogy with normal vel. component
			ALPHAI=ALPI[KS];
			OM=sqrt(Particle[N].PR*(1.0-ALPHAI)/(BOLTZ*TSURF[KS]));
			if(ISPR[1][L]==2){
				R=sqrt(-ALPHAI*log(RF(0)));
				CTH=cos(2.0*PI*RF(0));
			}else{
				//*--for polyatomic case, apply acceptance-rejection based on eqn (14.6)
Re10:			X=4.0*RF(0);
				A=2.7182818*X*X*exp(-X*X);
				if(A<RF(0)) goto Re10;
				R=sqrt(ALPHAI)*X;
				CTH=2.0*RF(0)-1.0;
			}
			Particle[N].PR=BOLTZ*TSURF[KS]*(R*R+OM*OM+2.0*R*OM*CTH);
		}
	}
	if(ISURF[KS]==1){
        	Particle[N].X=XC;
        	Particle[N].Y=YC+0.001*CG[6][MC];
	}
	if(ISURF[KS]==2){
        	Particle[N].X=XC;
        	Particle[N].Y=YC-0.001*CG[6][MC];
	}
	if(ISURF[KS]==3){
        	Particle[N].X=XC+0.001*CG[3][MC];
        	Particle[N].Y=YC;
	}
	if(ISURF[KS]==4){
        	Particle[N].X=XC-0.001*CG[3][MC];
        	Particle[N].Y=YC;
	}
    	Particle[N].IPL=(MC-1)*NSCX*NSCY+1;
	//*--sample the surface properties due to the reflected molecules
    	if(ISURF[KS]==1) CSS[3][K][L]=CSS[3][K][L]+SP[5][L]*Particle[N].PV[2];
    	if(ISURF[KS]==2) CSS[3][K][L]=CSS[3][K][L]-SP[5][L]*Particle[N].PV[2];
    	if(ISURF[KS]==3) CSS[3][K][L]=CSS[3][K][L]+SP[5][L]*Particle[N].PV[1];
    	if(ISURF[KS]==4) CSS[3][K][L]=CSS[3][K][L]-SP[5][L]*Particle[N].PV[1];
    	if(ISURF[KS]==1) CSS[9][K][L]=CSS[9][K][L]-SP[5][L]*Particle[N].PV[1];
    	if(ISURF[KS]==2) CSS[9][K][L]=CSS[9][K][L]-SP[5][L]*Particle[N].PV[1];
    	if(ISURF[KS]==3) CSS[9][K][L]=CSS[9][K][L]-SP[5][L]*Particle[N].PV[2];
    	if(ISURF[KS]==4) CSS[9][K][L]=CSS[9][K][L]-SP[5][L]*Particle[N].PV[2];
    	CSS[6][K][L]=CSS[6][K][L]-0.5*SP[5][L]*(pow(Particle[N].PV[1],2.0)+
		pow(Particle[N].PV[2],2.0)+pow(Particle[N].PV[3],2.0));
    	CSS[8][K][L]=CSS[8][K][L]-Particle[N].PR;
}
void REMOVE(int N)
{
	//*--remove molecule N and replace it by molecule NM
    Particle[N].X=Particle[NM].X;
    Particle[N].Y=Particle[NM].Y;
	int M;
	for( M=0;M<3;M++){
        Particle[N].PV[M]=Particle[NM].PV[M];
	}
    Particle[N].PR=Particle[NM].PR;
    Particle[N].IPL=Particle[NM].IPL;
    Particle[N].IPS=Particle[NM].IPS;
    NM=NM-1;
    N=N-1;
}


double ERF(double S)
{
	//*--calculates the error function of S
    double B,C,D,T;
	B=fabs(S);
	if(B>4.0){
        D=1.0;
	}else{
		C=exp(-B*B);
        T=1.0/(1.0+0.3275911*B);
        D=1.0-(0.254829592*T-0.284496736*T*T+1.421413741*T*T*T-
			1.453152027*T*T*T*T+1.061405429*T*T*T*T*T)*C;
	}
    if(S<0.0) D=-D;
	return D;
}
//COMMON /PLOTDATA/ 
double XC,YC,DEN,DENN,XM,P;
double VEL[4],TT,TROT,TEMP,DBOLTZ;
//data only for enter, properties and out2
void ENTER2()
{
	int L,NCS; 
	L=1;
	double VMP,SC,A;
	double BMEinlet[31]={0.0},BMEoutlet[31]={0.0};
	
	//BMEinlet=0.0;
	int n, ii,ndebug;	
	A=0;
	if(myid==0)	n=1;
	if(myid==(numprocs-1)) n=2;
	if(IIS>0){
		//for( n=1;n<=2; n++){
			L=1;
			for( ii=1;ii<=NCY;ii++){
				if(n==2) PROPERTIES(ii*NCX);
				if(n==1) PROPERTIES((ii-1)*NCX+1);

				VMP=sqrt(2.0*BOLTZ/SP[5][L]*TEMP);
				if(n==1) SC=VEL[1]/VMP;
				if(n==2) SC=-VEL[1]/VMP;
				//printf("%g\n",SC);
				if(fabs(SC)<10.1) A=(exp(-SC*SC)+SPI*SC*(1.0+ERF(SC)))/(2.0*SPI);
				//printf("%g\n",A);}
				//A=(exp(-SC*SC)+SPI*SC*(1.0+ERF(SC)))/(2.0*SPI);				
				if(SC>10.0) A=SC;
				/*if(A>0.0){
				printf("%g\t%g\t%g\t%g\n",SC,SPI,ERF(SC),exp(-SC*SC));//}*/
				if(SC<-10.0) A=0.0;
				//! for uniform cell
				if(n==1) BMEinlet[ii]=DENN*A*VMP*DTM*CG[6][(ii-1)*NCX+1]/FNUM;
				if(n==2) BMEoutlet[ii]=DENN*A*VMP*DTM*CG[6][(ii-1)*NCX+1]/FNUM; 
				/*if(A>0){
				printf("%g\t%g\t%g\t%g\t%g\n",DENN,A,VMP,DTM,BMEinlet[ii]);
				scanf("%d",&ndebug);}*/
			}
		//}
	}
	int N,NC,M,MC,MSCX,MSCY,MSC,K;
	double FS1,FS2,QA,U,UN,FTMP1;
//	for( N=1;N<=4;N++){
	if(myid==0)	N=1;
	if(myid==(numprocs-1)) N=2;
	//*--consider each boundary in turn
	if(IB[N]==1){
		//printf("%d\t%d\t%d\t%d\n",IB[1],IB[2],IB[3],IB[4]);			
		//scanf("%d",&NC);			
		if(N<3) NCS=NCY;
		if(N>2) NCS=NCX;

		for( NC=1;NC<=NCS;NC++){
			if(N==1) PROPERTIES((NC-1)*NCX+1);
			if(N==2) PROPERTIES(NC*NCX);
			for( L=1;L<=MNSP;L++){
				//*--consider each species in turn
				VMP=sqrt(2.0*BOLTZ/SP[5][L]*FTMP);
				if(N<3) VMP=sqrt(2.0*BOLTZ/SP[5][L]*TEMP);
				/*if(N==1){					
				printf("%g\t%g\n",VMP,TEMP);}*/
				if(N==1) A=BMEinlet[NC]+BMR[N][L];
				if(N==2) A=BMEoutlet[NC]+BMR[N][L];
				//if(N>2) A=BME[N][L]*CG[3][NC]/FW+BMR[N][L];
				M=floor(A);
				/*if((M>0)){					
				printf("%g\t%d\t\n",A,M);			
				scanf("%d",&M);}*/
				BMR[N][L]=A-(double)M;
				//*--M molecules enter, remainder has been reset
				if(M>0){
					if((N==1)||(N==2)){
						if(fabs(VEL[1])>1.0E-6){
							if(N==1) SC=VEL[1]/VMP;	 //-VFX/VMP
							if(N==2) SC=-VEL[1]/VMP; //!-VFX/VMP
						}
					}
					if((N==3)||(N==4)){
						if(fabs(VFY)>1.0E-6){
							if(N==3) SC=VFY/VMP;
							if(N==4) SC=-VFY/VMP;
						}
					}
					FS1=SC+sqrt(SC*SC+2.0);
					FS2=0.5*(1.0+SC*(2.0*SC-FS1));
					//* the above constants are required for the entering distn. of eqn (12.5)
					for( K=1;K<=M;K++){
						if(NM<MNM){
							NM=NM+1;
							//*--NM is now the number of the new molecule
							if(((N<3)&&(fabs(VEL[1])>1.0E-6))||((N>2)&&(fabs(VFY)>1.0E-6))){
								QA=3.0;
								if(SC<-3.0) QA=fabs(SC)+1.0;
Re2:								U=-QA+2.0*QA*RF(0);
								//*--U is a potential normalised thermal velocity component
								UN=U+SC;
								//*--UN is a potential inward velocity component
								if(UN<0.0) goto Re2;
								A=(2.0*UN/FS1)*exp(FS2-U*U);
								if(A<RF(0)) goto Re2;
								//*--the inward normalised vel. component has been selected (eqn (12.5))
								if(N==1) Particle[NM].PV[1]=UN*VMP;
								if(N==2) Particle[NM].PV[1]=-UN*VMP;
								if(N==3) Particle[NM].PV[2]=UN*VMP;
								if(N==4) Particle[NM].PV[2]=-UN*VMP;
							}else{
								if(N==1) Particle[NM].PV[1]=sqrt(-log(RF(0)))*VMP;
								if(N==2) Particle[NM].PV[1]=-sqrt(-log(RF(0)))*VMP;
								if(N==3) Particle[NM].PV[2]=sqrt(-log(RF(0)))*VMP;
								if(N==4) Particle[NM].PV[2]=-sqrt(-log(RF(0)))*VMP;
								//*--for a stationary external gas, use eqn (12.3)
							}
							if(N<3){
								RVELC(&Particle[NM].PV[2],&Particle[NM].PV[3],VMP);
								Particle[NM].PV[2]=Particle[NM].PV[2]+VEL[2]; //!VFY
							}
							if(N>2){
								RVELC(&Particle[NM].PV[1],&Particle[NM].PV[3],VMP);
								Particle[NM].PV[1]=Particle[NM].PV[1]+VEL[1]; //!VFX
							}
							//*--a single call of RVELC generates the two normal velocity components
							if(N==1) FTMP1=FTMP;
							if(N==2) FTMP1=TEMP;
							if(ISPR[1][L]>0) SROT(&Particle[NM].PR,FTMP1,ISPR[1][L]);
							if(N==1) Particle[NM].X=CB[1]+0.001*CG[3][1];
							if(N==2) Particle[NM].X=CB[2]-0.001*CG[3][MNC-1];
							if(N==3) Particle[NM].Y=CB[3]+0.001*CG[6][1];
							if(N==4) Particle[NM].Y=CB[4]-0.001*CG[6][MNC-1];
							//*--the molecule is moved just off the boundary
							Particle[NM].IPS=L;
							if(N<3){
								if(N==1) MC=(NC-1)*NCX+1;
								if(N==2) MC=NC*NCX;
								Particle[NM].Y=CG[4][MC]+RF(0)*CG[6][MC];
							}
							if(N>2){
								if(N==3) MC=NC;
								if(N==4) MC=(NCY-1)*NCX+NC;
								Particle[NM].X=CG[1][MC]+RF(0)*CG[3][MC];
							}
							MSCX=floor((((Particle[NM].X-CG[1][MC])/CG[3][MC])*(NSCX-0.001)+1));
							MSCY=floor((((Particle[NM].Y-CG[4][MC])/CG[6][MC])*(NSCY-0.001)+1));
							MSC=(MSCY-1)*NSCX+MSCX+NSCX*NSCY*(MC-1);
							//*--MSC is the new sub-cell number
							if(MSC<1) MSC=1;
							if(MSC>(MNSC-1)) MSC=MNSC-1;
							Particle[NM].IPL=MSC;
						}else{
							printf(" WARNING: EXCESS MOLECULE LIMIT RESTART WITH AN INCREASED FNUM\n");
						}
					}
				}
			}
		}
	//}
	}	
}

void PROPERTIES(int I)
{
	//*--Calculate cell properties
	int N;
	double COR_U,COR_U2,DEN_COR,T_COR;
	double A,SW,SUU,UU,SN,SM,SMU[4],SMCC,SRE,SRDF,SVEL[4],ASOUND;
	N=I;
	
	DBOLTZ=BOLTZ;
	
    A=FNUM/(CC[N]*NSMP);
    SN=0.0;
    SM=0.0;
	SW=0.0;	  
	int K,L;
	for( K=1;K<=3;K++){
		SMU[K]=0.0;
	}
	SUU=0.0;
	SMCC=0.0;
    SRE=0.0;
    SRDF=0.0;
	for( L=1;L<=MNSP;L++){
		SN=SN+CS[1][N][L];
		//*--SN is the number sum
	    SM=SM+SP[5][L]*CS[1][N][L];
		//*--SM is the sum of molecular masses
		for( K=1;K<=3;K++){
			SMU[K]=SMU[K]+SP[5][L]*CS[K+1][N][L];
			//*--SMU(1 to 3) are the sum of mu, mv, mw
		}
		SMCC=SMCC+(CS[5][N][L]+CS[6][N][L]+CS[7][N][L])*SP[5][L];
		//*--SMCC is the sum of m(u**2+v**2+w**2)

	    SRE=SRE+CSR[N][L];
		//*--SRE is the sum of rotational energy

		SRDF=SRDF+ISPR[1][L]*CS[1][N][L];
		//*--SRDF is the sum of the rotational degrees of freedom

		SUU=SUU+SP[5][L]*CS[5][N][L];
		//*--SUU is the sum of m*u*u

	}
	DENN=SN*A;
	//*--DENN is the number density, see eqn (1.34)
	mfp=1.0/(sqrt(2.0)*(PI*pow(SP[1][1],2.0))*DENN)  ;
    kn=mfp/CB[4];

	//* Cell Kundson No.

    DEN=DENN*SM/SN;
	//*--DEN is the density, see eqn (1.42)

	for( K=1;K<=3;K++){
		VEL[K]=SMU[K]/SM;
		SVEL[K]=VEL[K];
	}
	//*--VEL and SVEL are the stream velocity components, see eqn (1.43)

    UU=pow(VEL[1],2.0)+pow(VEL[2],2.0)+pow(VEL[3],2.0);
    TT=(SMCC-SM*UU)/(3.0*DBOLTZ*SN);
	//*--TT is the translational temperature, see eqn (1.51)

    if(SRDF>1.0E-6) TROT=(2.0/DBOLTZ)*SRE/SRDF;
	//*--TROT is the rotational temperature, see eqn (11.11)

    TEMP=(3.0*TT+(SRDF/SN)*TROT)/(3.0+SRDF/SN);
	P=DENN*BOLTZ*TEMP;
	//*--TEMP is the overall temperature, see eqn (11.12)

    CT[N]=TEMP;
    XC=0.5*(CG[1][N]+CG[2][N]);
    YC=0.5*(CG[4][N]+CG[5][N]);
	GAMA=(5.0+ISPR[1][1])/(3.0+ISPR[1][1]);
	GASR=BOLTZ/SP[5][1];
	ASOUND=sqrt(GAMA*GASR*TEMP);
	XM=sqrt(UU)/ASOUND;
	//*--XC,YC are the x,y coordinates of the midpoint of the cell
	//! Corrections FOR OUTPUT CELLS

	if((N-(N/NCX)*NCX)==0){
		//printf("%g\n",P);
		COR_U=(P-POUT)/(DEN*ASOUND);
		//printf("%g\n",COR_U);
		VEL[1]=VEL[1]+COR_U;
		DEN_COR=DEN+(POUT-P)/pow(ASOUND,2.0);
		T_COR=POUT/(DEN_COR*GASR);
		DEN=DEN_COR;
		TEMP=T_COR;
		DENN=DEN/SP[5][1];
		//P=DENN*BOLTZ*TEMP
	}

	//! Corrections FOR INPUT CELLS
	if((N-(N/NCX)*NCX)==1){
		//printf("%g\n",P);		
		COR_U2=(PIN-P)/(DEN*ASOUND);
		//printf("%g\n",COR_U2);
		VEL[1]=VEL[1]+COR_U2;
		TEMP=FTMP;
		//!P=PIN
		DEN=PIN/(GASR*FTMP);
		DENN=DEN/SP[5][1];
		//!P=DENN*BOLTZ*TEMP		 		 
	}
}
void OUT2()
{
     FILE *fp4;
	 int N,ndebug;

//--output a progressive set of results to file DSMC2.OUT.
 
	//include 'property.txt' 

//	COMMON /PLOTDATA/ XC,YC,DEN,DENN,TT,TROT,TEMP,VEL,XM,P

//      DOUBLE PRECISION VEL(3),SMU(3),SVEL(3,MNC),SN,SM,SMCC,SRDF,SRE,TT,
//     &                 TROT,DBOLTZ,SS(9)
//      DOUBLE PRECISION::ASOUND,GAMA,GASR,P,POUT,DEN_COR,T_COR,COR_U,TEMP
//	DOUBLE PRECISION::XC,YC,DEN,DENN,XM
	DBOLTZ=BOLTZ;
    	FILE *fp;
	if(myid==0){
		fp=fopen("SUB 0.xls","wt");
//		fprintf(fp,"Variables=X,Y,Density, TrTemp, RotTemp, T,U,V,W,Mach,Pressure \n Zone I=%d,\tJ=%d,\tF=Point\n",IEnP-IStP+1,NCY);
	}else if(myid==1){
		fp=fopen("SUB 1.xls","wt");
	}else if(myid==2){
		fp=fopen("SUB 2.xls","wt");
	}else if(myid==3){
		fp=fopen("SUB 3.xls","wt");
	}else if(myid==4){
		fp=fopen("SUB 4.xls","wt");
	}else if(myid==5){
		fp=fopen("SUB 5.xls","wt");
	}else if(myid==6){
		fp=fopen("SUB 6.xls","wt");
	}else if(myid==7){
		fp=fopen("SUB 7.xls","wt");
	}
	//OPEN (4,FILE='SUB.xls')
    //WRITE (4,*)'Variables=X,Y,Density, TrTemp, RotTemp, T,U,V,W,Mach,Pressure'
	//WRITE (4,*)'ZONE   I=',NCX,',  J=',NCY,',  F=POINT'

	//for (N=1;N<MNC;N++){
	fprintf(fp,"Variables=X,Y,Density,TrTemp,RotTemp,T,U,V,W,Mach,P/P<sub>in</sub>,Kn \n Zone I=%d,\tJ=%d,\tF=Point\n",IEnP-IStP+1,NCY);
	int i,j;
	for(j=1;j<=NCY;j++){
	for(i=IStP;i<=IEnP;i++){
		N=(j-1)*NCX+i;	
		PROPERTIES (N);
	    	fprintf(fp,"%g\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g\n",
		XC,YC,DEN,TT,TROT,TEMP,VEL[1],VEL[2],VEL[3],XM,P/1.0E6,kn);
	}}
      	fclose(fp);
	//scanf("%d",&ndebug);
}
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
struct MolPro SendBuffer[20000],RecvBuffer[20000];	//save the properties of molecules crossing partition bondaries
void Transmit(struct MolPro particle, int N, int Orig, int New, int RA)
{
    	MPI_Status status;
	if(RA==-1){
		MPI_Ssend(&particle.X ,1,MPI_DOUBLE,New,99,MPI_COMM_WORLD);
		MPI_Ssend(&particle.Y ,1,MPI_DOUBLE,New,99,MPI_COMM_WORLD);
		MPI_Ssend( particle.PV,4,MPI_DOUBLE,New,99,MPI_COMM_WORLD);
		MPI_Ssend(&particle.PR,1,MPI_DOUBLE,New,99,MPI_COMM_WORLD);
		MPI_Ssend(&particle.IPL,1,MPI_INT,New,99,MPI_COMM_WORLD);
		MPI_Ssend(&particle.IR ,1,MPI_INT,New,99,MPI_COMM_WORLD);
		MPI_Ssend(&particle.IPS,1,MPI_INT,New,99,MPI_COMM_WORLD);
		MPI_Ssend(&particle.IX,1,MPI_INT,New,99,MPI_COMM_WORLD);
		MPI_Ssend(&particle.JY,1,MPI_INT,New,99,MPI_COMM_WORLD);
//		MPI_Ssend(&particle.Part,1,MPI_INT,New,99,MPI_COMM_WORLD);
	}else{
		MPI_Recv(&particle.X ,1,MPI_DOUBLE, Orig, 99, MPI_COMM_WORLD,&status);
		if((particle.X+0.00001*CG[3][1]<XStP)||(particle.X-0.00001*CG[3][1]>XEnP))
			printf("****Error in move_Partition****\t%d\t%g\t%g\t%g\n",myid,particle.X,XStP,XEnP);
		MPI_Recv(&particle.Y ,1,MPI_DOUBLE, Orig, 99, MPI_COMM_WORLD,&status);
		MPI_Recv( particle.PV,4,MPI_DOUBLE, Orig, 99, MPI_COMM_WORLD,&status);
		MPI_Recv(&particle.PR,1,MPI_DOUBLE, Orig, 99, MPI_COMM_WORLD,&status);
		MPI_Recv(&particle.IPL,1,MPI_INT, Orig, 99, MPI_COMM_WORLD,&status);
		MPI_Recv(&particle.IR ,1,MPI_INT, Orig, 99, MPI_COMM_WORLD,&status);
		MPI_Recv(&particle.IPS,1,MPI_INT, Orig, 99, MPI_COMM_WORLD,&status);
		MPI_Recv(&particle.IX,1,MPI_INT, Orig, 99, MPI_COMM_WORLD,&status);
		MPI_Recv(&particle.JY,1,MPI_INT, Orig, 99, MPI_COMM_WORLD,&status);
//		MPI_Recv(&particle.Part,1,MPI_INT, Orig, 99, MPI_COMM_WORLD,&status);
		particle.Part=myid;
		RecvBuffer[N]=particle;
	}
}
void INDEXM()
{
	//checking each particle's partition
	MPI_Status status;
	int N,PFind,part;
	int i,j,ix,jy,recv,sent,ndebug;
	sent=0;	recv=0;
	//printf("%d\n",NM);
	int MovetoPartition[128];

	for(i=0;i<128;i++){
		MovetoPartition[i]=0;
	}
	//Syncronization//////////
/*	if(myid==(numprocs-1))
		MPI_Ssend(&i,1,MPI_INT,0,99,MPI_COMM_WORLD);
	if(myid==0)
		MPI_Recv(&j,1,MPI_INT,(numprocs-1),99,MPI_COMM_WORLD,&status);
	
	if(myid!=(numprocs-1)) MPI_Ssend(&i,1,MPI_INT,myid+1,99,MPI_COMM_WORLD);
	if(myid!=0)	 MPI_Recv(&j,1,MPI_INT,myid-1,99,MPI_COMM_WORLD,&status);
*/	//printf("Syncronization, id=%d\t%g\n",myid,CW);
	N=0;
RN:	N=N+1;	//N:NM loop
	if(N<=NM){
		Particle[N].IX=floor((Particle[N].X-CB[1])/CW+0.99999);
		if(Particle[N].IX<1) Particle[N].IX=1;
		if(Particle[N].IX>NCX) Particle[N].IX=NCX;
		ix=Particle[N].IX;
		if(((Particle[N].X-0.00001*CG[3][ix])>CG[2][ix])||((Particle[N].X+0.00001*CG[3][ix])<CG[1][ix])){
			printf("Error in INDEXM()\n");
		}
		//(XI+0.00001*CG[3][IStP]<XStP)||(XI-0.00001*CG[3][IEnP]>XEnP)
		if((ix<IStP)||(ix>IEnP)){	//the particle has moved from this partition (myid)
			//printf("%d\t%d\t%d\n",myid,N,ix);
			//scanf("%d",&ndebug);
			PFind=-1;
			sent++;
			for(part=0;part<numprocs;part++){
				if((ix>=Partition[part].IStart)&&(ix<=Partition[part].IEnd)){
					//printf("%d\t%d\t%d\n",myid,N,NM);
					//scanf("%d",&ndebug);
					PFind=part;
					MovetoPartition[PFind]++;
					Particle[N].Part=PFind;
					SendBuffer[sent]=Particle[N];
					REMOVE(N);
				}
			}
			if(PFind<0){
				printf("Error in index\t%d\t%d\t%d\t%d\n",myid,ix,Partition[myid].IStart,Partition[myid].IEnd);
				scanf("%d",&ndebug);
			}
		}
		goto RN;
	}

/*	//********************move the crossing molecules//*****************
	//nokteee dar ruze 14 Nov be an rasidam, yek processor nemitavanad 2 bar poshte sare ham be processore digari 		//		peigham befrestad dar hali ke hanuz peighame avvali daryaft nashode ast, 
	//					otherwise:deadlock
	//KARMA theory, vaghti to be kasi badi koni, ghtan shakhse sevvomi ham be to badi mikonad
*/	int mfap;

	recv=0;
 	for(part=0;part<numprocs;part++){
		if(part!=myid){
			if(myid<part){
				MPI_Ssend(&MovetoPartition[part],1,MPI_INT,part,99,MPI_COMM_WORLD);
				//printf("Procs%d\tsending %d to prc%d\n",myid,MovetoPartition[part],part);
			}
			if(myid>part){	
				MPI_Recv (&mfap,1,MPI_INT,part, 99, MPI_COMM_WORLD,&status);
				//printf("Procs%d\treceiving%d from prc%d\n",myid,mfap,part);
			}
			if((myid<part)&&(MovetoPartition[part]>0)){
				//printf("procs %d sending %d molecules to proces %d \n",myid,MovetoPartition[part],part);
				j=0;
				for(i=1;i<=sent;i++){
					if(part==SendBuffer[i].Part){
						Transmit(SendBuffer[i],-1,myid,SendBuffer[i].Part,-1);	
						j=j+1;
					}
				}
				if(j!=MovetoPartition[part])printf("Error in index 2\n");
			}
			if((myid>part)&&(mfap>0)){
				//printf("procs %d receiving %d molecules from proces %d \n",myid,mfap,part);
				for(i=1;i<=mfap;i++){
					recv++;
					Transmit(RecvBuffer[recv],recv,part,myid,+1);	
				}
			}
		}
	}
	for(part=0;part<numprocs;part++){
		if(part!=myid){
			if(myid>part){
				MPI_Ssend(&MovetoPartition[part],1,MPI_INT,part,99,MPI_COMM_WORLD);
				//printf("Procs%d\t**sending %d to prc%d\n",myid,MovetoPartition[part],part);
			}
			if(myid<part){	
				MPI_Recv (&mfap,1,MPI_INT,part, 99, MPI_COMM_WORLD,&status);
				//printf("Procs%d\t**receiving %d from prc%d\n",myid,mfap,part);
			}
			if((myid>part)&&(MovetoPartition[part]>0)){
				//printf("procs %d sending %d molecules to proces %d \n",myid,MtoPart[part],part);
				j=0;
				for(i=1;i<=sent;i++){
					if(part==SendBuffer[i].Part){
						Transmit(SendBuffer[i],-1,myid,SendBuffer[i].Part,-1);	
						j=j+1;
					}
				}
				if(j!=MovetoPartition[part])printf("Error in index 3\n");
			}
			if((myid<part)&&(mfap>0)){
				//printf("procs %d receiving %d molecules from proces %d \n\n",myid,mfap,part);
				for(i=1;i<=mfap;i++){
					recv++;
					Transmit(RecvBuffer[recv],recv,part,myid,+1);	
				}
			}
		}
	}
	//printf("total sent=%d,\ttatoal received=%d, myid=%d\n",sent,recv,myid);
	for(i=1;i<=recv;i++){
		NM++;
		Particle[NM]=RecvBuffer[i];
		if((Particle[NM].IX<IStP)||(Particle[NM].IX>IEnP))
			printf("Error in index 5\n");
	}


//--the NM molecule numbers are arranged in order of the molecule groups
//--and, within the groups, in order of the cells and, within the cells,in order of the sub-cells
	int MM, NN;
    	for( MM=1;MM<=MNSG;MM++){
		IG[2][MM]=0;      
		for( NN=1;NN<MNC;NN++)
         		 IC[2][NN][MM]=0;
		for( NN=1;NN<MNSC;NN++)
         		 ISCG[2][NN][MM]=0;
	}
   	int LS,MG,MSC,M,MC,K;

	int L;
	for( N=1;N<=NM;N++){
		LS=Particle[N].IPS;
        	MG=ISP[LS];
        	IG[2][MG]=IG[2][MG]+1;
	        MSC=Particle[N].IPL;
        	ISCG[2][MSC][MG]=ISCG[2][MSC][MG]+1;
        	MC=ISC[MSC];
        	IC[2][MC][MG]=IC[2][MC][MG]+1;
	}
	//--number in molecule groups in the cells and sub-cells have been counte
    	M=0;
    	for( L=1;L<=MNSG;L++){
        	IG[1][L]=M;
		//--the (start address -1) has been set for the groups
        	M=M+IG[2][L];
    	}
    	for( L=1;L<=MNSG;L++){
        	M=IG[1][L];
        	for( N=1;N<MNC;N++){
			IC[1][N][L]=M;
			M=M+IC[2][N][L];
	    	}
		//--the (start address -1) has been set for the cells
	        M=IG[1][L];
        	for( N=1;N<MNSC;N++){
			ISCG[1][N][L]=M;
			M=M+ISCG[2][N][L];
			ISCG[2][N][L]=0;
	        }
	}
	//--the (start address -1) has been set for the sub-cells
 	for( N=1;N<=NM;N++){
        	LS=Particle[N].IPS;
	        MG=ISP[LS];
        	MSC=Particle[N].IPL;
        	ISCG[2][MSC][MG]=ISCG[2][MSC][MG]+1;
        	K=ISCG[1][MSC][MG]+ISCG[2][MSC][MG];
        	Particle[K].IR=N;
		//--the molecule number N has been set in the cross-reference array
    	}
}    
void SAMPI2()
{
//--initialises all the sampling variables
 NSMP=0;
 TIMI=TIME;
	int L,M,N;
 for( L=1;L<=MNSP;L++){
	 for( N=1;N<MNC;N++){
          CS[1][N][L]=1.0E-6;
		  for( M=2;M<=7;M++){
            CS[M][N][L]=0.0;
		  }
          CSR[N][L]=0.0;
	 }

	 for( N=1;N<=MNSE;N++){
          CSS[1][N][L]=1.0E-6;
		  for( M=2;M<=9;M++){
            CSS[M][N][L]=0.0;
		}
	 }
 }

}
void SAMPLE2()
{
	int L,K,M,I;
//--sample the molecules in the flow.

    NSMP=NSMP+1;
	int NN,N,J,LL;	
	for( NN=1;NN<=MNSG;NN++){
		for( N=1;N<MNC;N++){
			L=IC[2][N][NN];
			if (L>=0) {
				for( J=1;J<=L;J++){
					K=IC[1][N][NN]+J;
					M=Particle[K].IR;
					I=Particle[M].IPS;
					CS[1][N][I]=CS[1][N][I]+1;
					for( LL=1;LL<=3;LL++){
						CS[LL+1][N][I]=CS[LL+1][N][I]+Particle[M].PV[LL];
						CS[LL+4][N][I]=CS[LL+4][N][I]+Particle[M].PV[LL]*Particle[M].PV[LL];
					}
					CSR[N][I]=CSR[N][I]+Particle[M].PR;
				}
			}
		}	  
	}
}

//*--VRC(3) are the pre-collision components of the relative velocity
double VRC[4],VRR,VR,CVR;
int L,M,MM,NN,NLS,MS,NMSC,NST,NSG,INC,N,MSC,LS;
//for collision 
void SELECT()
{
	//--selects a potential collision pair and calculates the product of the
	//--collision cross-section and relative speed
	//      include 'common.txt'
	int K,Ki;
	K=floor((RF(0)*(IC[2][N][NN]-0.001)))+IC[1][N][NN]+1;
    L=Particle[K].IR;
//--the first molecule L has been chosen at random from group NN in cell
Re11:   MSC=Particle[L].IPL;
    if(((NN==MM)&&(ISCG[2][MSC][MM]==1))||((NN!=MM)&&(ISCG[2][MSC][MM]==0))){
		//--if MSC has no type MM molecule find the nearest sub-cell with one
        NST=1;
        NSG=1;
Re5:    INC=NSG*NST;
        NSG=-NSG;
        NST=NST+1;
        MSC=MSC+INC;
        if((MSC<1)||(MSC>(MNSC-1))) goto Re5;
        if((ISC[MSC]!=N)||(ISCG[2][MSC][MM]<1)) goto Re5;
	}
	//--the second molecule M is now chosen at random from the group MM
	//--molecules that are in the sub-cell MSC
    K=floor(RF(0)*(ISCG[2][MSC][MM]-0.001))+ISCG[1][MSC][MM]+1;
    M=Particle[K].IR;
    if(L==M) goto Re11;
	//--choose a new second molecule if the first is again chosen
    for( Ki=1;Ki<=3;Ki++){
        VRC[Ki]=Particle[L].PV[Ki]-Particle[M].PV[Ki];
	}
	//--VRC(1 to 3) are the components of the relative velocity
    VRR=VRC[1]*VRC[1]+VRC[2]*VRC[2]+VRC[3]*VRC[3];
    VR=sqrt(VRR);
	//--VR is the relative speed
    LS=Particle[L].IPS;
    MS=Particle[M].IPS;
    CVR=VR*SPM[1][LS][MS]*(
		pow((2.0*BOLTZ*SPM[2][LS][MS]/(SPM[5][LS][MS]*VRR)),(SPM[3][LS][MS]-0.5)))/
		SPM[6][LS][MS];
	//--the collision cross-section is based on eqn (4.63)
}
void ELASTIC()
{
	//--generate the post-collision velocity components.
	//--VRCP(3) are the post-collision components of the relative velocity
	//--VCCM(3) are the components of the centre of mass velocity
	double RML,RMM,A,B,C,VRCP[4],VCCM[4];
	RML=SPM[5][LS][MS]/SP[5][MS];
    RMM=SPM[5][LS][MS]/SP[5][LS];
	int K;    
	for( K=1;K<=3;K++)	
        VCCM[K]=RML*Particle[L].PV[K]+RMM*Particle[M].PV[K];
	//--VCCM defines the components of the centre-of-mass velocity, eqn (2.1)
    if(fabs(SPM[4][LS][MS]-1.0)<1.0E-3){
		//--use the VHS logic
		B=2.0*RF(0)-1.0;
		//--B is the cosine of a random elevation angle
		A=sqrt(1.0-B*B);
        VRCP[1]=B*VR;
        C=2.0*PI*RF(0);
		//--C is a random azimuth angle
		VRCP[2]=A*cos(C)*VR;
        VRCP[3]=A*sin(C)*VR;
	}else{
		//--use the VSS logic
        B=2.0*pow( RF(0),SPM[5][LS][MS])-1.0;
		//--B is the cosine of the deflection angle for the VSS model, eqn (11.8)
        A=sqrt(1.0-B*B);
        double C=2.0*PI*RF(0);
        double OC=cos(C);
        double SC=sin(C);
        double D=sqrt(VRC[2]*VRC[2]+VRC[3]*VRC[3]);
        if(D>1.0E-6)
		{
			VRCP[1]=B*VRC[1]+A*SC*D;
			VRCP[2]=B*VRC[2]+A*(VR*VRC[3]*OC-VRC[1]*VRC[2]*SC)/D;
			VRCP[3]=B*VRC[3]-A*(VR*VRC[2]*OC+VRC[1]*VRC[3]*SC)/D;
		}else{
			VRCP[1]=B*VRC[1];
			VRCP[2]=A*OC*VRC[1];
			VRCP[3]=A*SC*VRC[1];
		}
		//--the post-collision rel. velocity components are based on eqn (2.22)
	}
	//--VRCP(1 to 3) are the components of the post-collision relative vel.
    for ( K=1;K<=3;K++){
        Particle[L].PV[K]=VCCM[K]+VRCP[K]*RMM;
        Particle[M].PV[K]=VCCM[K]-VRCP[K]*RML;
	}
}



void LBS(double XMA,double XMB,double ERM)
{
	//*--selects a Larsen-Borgnakke energy ratio using eqn (11.9)
Re100:   ERM=RF(0);
	if((XMA<1.0E-6)||(XMB<1.0E-6)){
		if((XMA<1.0E-6)&&(XMB<1.0E-6)) return;
        if(XMA<1.0E-6) P=pow((1.0-ERM),XMB);
        if(XMB<1.0E-6) P=pow((1.0-ERM),XMA);
	}else{
		P=pow(((XMA+XMB)*ERM/XMA),XMA)*(XMA+XMB)*pow(((1.0-ERM)/XMB),XMB);
	}
    if(P<RF(0)) goto Re100;
}


void INELR()
{
	//*--adjustment of rotational energy in a collision
    int IR[3];	//*--IR is the indicator for the rotational redistribution
	int IRT,K,KS,JS;
	double ETI,ECI,ECF,ECC,XIB,ATK,ERM,XIA,ETF,A;
	ETI=0.5*SPM[5][LS][MS]*VRR;	//*--ETI is the initial translational energy
    ECI=0.0;	//*--ECI is the initial energy in the active rotational modes
    ECF=0.0;	//*--ECF is the final energy in these modes
    ECC=ETI;	//*--ECC is the energy to be divided
    XIB=2.5-SPM[3][LS][MS];	//*--XIB is th number of modes in the redistribution
    IRT=0;		//*--IRT is 0,1 if no,any redistribution is made
	int NSP;
	for(NSP=1;NSP<=2;NSP++){
		//*--consider the molecules in turn
		if(NSP==1){
			K=L;
			KS=LS;
			JS=MS;
		}else{
			K=M;
			KS=MS;
			JS=LS;
		}
        IR[NSP]=0;
		if(ISPR[1][KS]>0){
			if(ISPR[2][KS]==0){
				ATK=1.0/SPR[1][KS][JS];
			}else{
				ATK=1.0/(SPR[1][KS][JS]+SPR[2][KS][JS]*CT[N]+SPR[3][KS][JS]*pow(CT[N],2.0));
			}
			//*--ATK is the probability that rotation is redistributed to molecule L
			if(ATK>RF(0)){
				IRT=1;
	            IR[NSP]=1;
		        ECC=ECC+Particle[K].PR;
			    ECI=ECI+Particle[K].PR;
				XIB=XIB+0.5*ISPR[1][KS];
			}
		}
	}
	//*--apply the general Larsen-Borgnakke distribution function
	if(IRT==1){
		for( NSP=1;NSP<=2;NSP++){
			if(IR[NSP]==1){
				if(NSP==1){
					K=L;
					KS=LS;
				}else{
					K=M;
					KS=MS;
				}
				XIB=XIB-0.5*ISPR[1][KS];
				//*--the current molecule is removed from the total modes
				if(ISPR[1][KS]==2){
					ERM=1.0-pow(RF(0),(1.0/XIB));
				}else{
					XIA=0.5*ISPR[0][KS];
					LBS(XIA-1.0,XIB-1.0,ERM);
				}
				Particle[K].PR=ERM*ECC;
				ECC=ECC-Particle[K].PR;
				//*--the available energy is reduced accordingly
				ECF=ECF+Particle[K].PR;
			}
		}
        ETF=ETI+ECI-ECF;
		//*--ETF  is the post-collision translational energy
		//*--adjust VR and, for the VSS model, VRC for the change in energy
        A=sqrt(2.0*ETF/SPM[5][LS][MS]);
		if(fabs(SPM[4][LS][MS]-1.0)<1.0E-3){
			VR=A;
		}else{
			for( K=1;K<=3;K++){
				VRC[K]=VRC[K]*A/VR;
			}
			VR=A;
		}
	}
}
void COLLMR()
{
	//*--calculates collisions appropriate to DTM in a gas mixture
	double AVN,ASEL,CVM,SN;
	int NSEL;
	int K,ISEL;
	//for(N=1;N<MNC;N++){
	int i,j;
	for(j=1;j<=NCY;j++){
	for(i=IStP;i<=IEnP;i++){
		N=(j-1)*NCX+i;
		//*--consider collisions in cell N
		for(NN=1;NN<=MNSG;NN++){
			for( MM=1;MM<=MNSG;MM++){
				SN=0.0;
				for( K=1;K<=MNSP;K++){
					if(ISP[K]==MM) SN=SN+CS[1][N][K];
				}
				if(SN>1.0){
					AVN=SN/(double)(NSMP);
				}else{
					AVN=IC[2][N][MM];
				}
				//*--AVN is the average number of group MM molecules in the cell
				ASEL=0.5*IC[2][N][NN]*AVN*FNUM*CCG[1][N][NN][MM]*DTM/CC[N]+CCG[2][N][NN][MM];
				//*--ASEL is the number of pairs to be selected, see eqn (11.5)
				NSEL=floor(ASEL);
				CCG[2][N][NN][MM]=ASEL-NSEL;
				if(NSEL>0){
					if(((NN!=MM)&&((IC[2][N][NN]<1)||(IC[2][N][MM]<1)))
						||((NN==MM)&&(IC[2][N][NN]<2))){
						CCG[2][N][NN][MM]=CCG[2][N][NN][MM]+NSEL;
						//*--if there are insufficient molecules to calculate collisions,
						//*--the number NSEL is added to the remainer CCG(2,N,NN,MM)
					}else{
						CVM=CCG[1][N][NN][MM];
						SELT=SELT+NSEL;
						for( ISEL=1;ISEL<=NSEL;ISEL++){
							SELECT();
							if(CVR>CVM) CVM=CVR;
							//*--if necessary, the maximum product in CVM is upgraded
							if(RF(0)<(CVR/CCG[1][N][NN][MM])){
								//*--the collision is accepted with the probability of eqn (11.6)
								NCOL=NCOL+1;
								SEPT=SEPT+
									sqrt(pow((Particle[L].X-Particle[M].X),2.0)+pow((Particle[L].Y-Particle[M].Y),2.0));
								COL[LS][MS]=COL[LS][MS]+1.00;
								COL[MS][LS]=COL[MS][LS]+1.00;
								if((ISPR[1][LS]>0.0)||(ISPR[1][MS]>0)) INELR();
								//*--bypass rotational redistribution if both molecules are monatomic
								ELASTIC();
							}
						}
						CCG[1][N][NN][MM]=CVM;
					}
				}
			}
		}
	}}
}
