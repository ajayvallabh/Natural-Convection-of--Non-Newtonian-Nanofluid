/*Nomenclature

Kf            - thermal conductivity of base fluid
Kp            - thermal conductivity of Nano-particle
Kef           - Effective thermal conductivity of Nanofluid
phi           - Volume fraction
beta_f        - Volumetric thermal expansion of base fluid
beta_s        - Volumetric thermal expansion of solid
beta_nf       - Volumetric thermal expansion of  Nanofluid
Muf           - Viscosity of non-Newtonian fluid
Mu_nf         - Effective Viscosity of Non-Newtonian nannofluid
rho_f         - Density of base fluid
rho_s         - Density of solid nano-particle
rho_nf        - Density of nanofluid
Cpf           - Heat capacity of base fluid
Cps           - Heat capacity of solid nano-particle
rhoCp_nf      - Volumetric Heat capacity of nanofluid
alpha_nf      - Thermal Diffusivity of nanofluid
N             - Consistency index
nu_nf         - Kinematic viscosity of nanofluid
Pr_nf         - Prandtl number of nanofluid
PrStar        - prandtl number of non-Newtonian nanofluid
Ra_nf         - Rayleigh number of non-Newtonian nanofluid
Nu_nf         - Nusselt number
n             - power law index
L             - Length of Square Geometry
Th            - Temperature of hot wall
Tc            - Temperature of cold wall
g             - gravity
u             - velocity of nanofluid in x direction
v             - velocity of nanofluid in y direction
tau           - stress
psi           - stream function
zeta          - Vorticity
theta         - Temperature

*/
/*
Code for Non-Newtonian Nanofluid

Solution of Natural convection steady state nanofluid flow inside enclosure

Ajay vallabh (16205402)

Task performed
1. Computation of stream function.
2. Computation of temperature field
3. Velocity profiles
4. Local Nusselt Number on hot wall

Files:
input.txt:
          first line  :  Kf, Kp, phi
          Second line : beta_f, beta_s
          Third line  : consistency index(N) , power law index (n)
          fourth line : Density of base fluid(rho_f) , Density of solid nano-particle (rho_s)
          fifth line  : Heat capacity of base fluid(Cpf),Heat capacity of solid nano-particle(Cps)
          sixth line  : Length(L),Height(H),Temperature of hot wall(Th),Temperature of cold wall(Tc),gravity(g)
          seventh line: delta_x, delta_y ,delta_t (time step for unsteady solution)
          eight line  : Precision, Relaxation_parameter, number of time steps to skip while printing the status on screen during execution

For e.g. :

    0.56 200 0.04
    0.87e-4 42e-6
    2.1 1.0
    1000 19300
    4185.5 125.6
    1.0 1.0 45 25 9.81
    0.025 0.025  0.005
    1e-4 1.5 50

Note: Commas are not included

Output files:

     1. output1 - psi
     2. output2 -zeta
     3. output3 -theta
     4. output4 - u
     5  output5 - v

Format:
	x
	y
	u(x,y)
	theta(x,y)
	psi(x,y)
	v(x,y)
	zeta(x,y)
  Algorithm:
     Gauss Seidel with successive over-relaxation, and transient/pseudo-transient method to reach steady state solution

*/
#include"nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>
#include "NanofluidProperties.h"
int main()
{
    int k, l, skip_itr ;
	// k is the number of nodes along x-axis and l is the same along y-axis
    // skip_iter - number of iterations to skip before priting running status on screen

double **psi, **psi_old,**psi_guess,**tau_xx , **tau_xy, **tau_yy, **zeta, **zeta_old, **zeta_guess ,**theta, **theta_old,**theta_guess,**u,**v,**Nu_nf ;
/*     zeta     - vorticity
	   theta 	- temperature
	   u 		- u_x, x-component of velocity
	   v 		- u_y y-component of velocity
	   psi 		- stream function

	   delta_x = delta_y = grid size
	   relaxation - relaxation factor in Gauss Seidel algorithm
	   Precision  - precision required in solution i.e. criterion for convergence of solution

	   Ra_nf 		- Rayleigh_number
	   Pr_nf 		- Prandtl_number
	   PrStar       - Prandtl number for non-Newtonian fluid
	   Nu_nf        - Nusselt number
	   tau_xx       - normal mean stress in x direction
	   tau_yy       - normal mean stress in y direction
	   tau_xy       - mean shear stress in y direction with x plane
*/
     int i, j  ;  // loop counters
     double delta_x, delta_y, delta_t ;
     /*fluid properties*/
     double Kf, Kp, Kef, phi, beta_f, beta_s, beta_nf ;
     double N, n, rho_f, rho_s, rho_nf, Cpf, Cps,rhoCp_nf ;
     double Muf, Mu_nf , alpha_nf, nu_nf;
	 double relaxation, Precision ;
     double L,H,Th, Tc, gravity ;
     double  PrS, Ra_nf,NuAvg ;

int    f,itr_time, itr_space,sume,sumo;
char outputFileName[120] = "output_u_theta_psi_" ;
char values[30] ; // used for form filename with input variables Ra, Pr
double error1,error2,error3,error4,error5,error6,r1,r2,r3,r4,r5,r6,errors,errort,time;
double ae,aw,as,an,ap;
double ae1,as1,aw1,an1,ap1;
double ae2,as2,aw2,an2,ap2,RHS2,STRESS,STRESS1,STRESS2,STRESS3,STRESS4;
double iue,iuw,ivs,ivn;
double q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13;

FILE *input,*output,*output1,*output2,*output3,*output4,*output5,*output6,*output7,*output8;
output1=fopen("psi.txt","w");
output2=fopen("zeta.txt","w");
output3=fopen("theta.txt","w");
output4=fopen("Uvelocity.txt","w");
output5=fopen("Vvelocity.txt","w");
output6=fopen("nusselt.txt","w");
output7=fopen("AvgNusselt.txt","w");
output8=fopen("Allerrors.txt","w");
     input=fopen( "input.txt", "r" ) ;
     fscanf(input, "%lf %lf %lf", &Kf, &Kp, &phi ) ;
     fscanf(input, "%lf %lf", &beta_f,&beta_s ) ;
     fscanf(input, "%lf %lf", &N, &n ) ;
     fscanf(input, "%lf %lf", &rho_f, &rho_s ) ;
     fscanf(input, "%lf %lf", &Cpf, &Cps ) ;
     fscanf(input, "%lf %lf %lf %lf %lf", &L, &H,&Th, &Tc, &gravity ) ;
	 fscanf(input, "%lf %lf %lf", &delta_x, &delta_y, &delta_t ) ;
	 fscanf(input, "%lf %lf %d", &Precision, &relaxation, &skip_itr ) ;

k=ceil(1.0/delta_x+1);
l=ceil(1.0/delta_y+1);
f=(k-1)*(k-1);

     Kef = Effective_Thermal_Conductivity(Kf,Kp,phi);
     rho_nf = Nanofluid_Density(rho_f,rho_s,phi);
     rhoCp_nf = Heat_Capacity(rho_f, Cpf, rho_s, Cps, phi);
     beta_nf = Volumetric_Thr_Exp( rho_f, beta_f, rho_s, beta_s, rho_nf, phi);
     alpha_nf = Thermal_Diffusivity( Kef, rhoCp_nf);
     nu_nf = Kinematic_Viscosity( N, rho_nf );
    // Pr_nf = Prandtl_Number( nu_nf, alpha_nf, N  );
     PrS = PrandtlStar_Number( N, n, L, alpha_nf, rho_nf);
     Ra_nf = Rayleigh_Number( alpha_nf, nu_nf, beta_nf , Th, Tc, gravity, L,n);
     printf("Rayleigh=%lf\n",Ra_nf);
     printf("PrandtlS=%lf\n",PrS);

    sprintf( values, "Ra%.1f_Pr%.1f.txt", Ra_nf, PrS ) ;
	strcat(outputFileName, values ) ;
/*define matrix*/
psi        =dmatrix(1,k,1,l);
psi_old    =dmatrix(1,k,1,l);
psi_guess  =dmatrix(1,k,1,l);
zeta       =dmatrix(1,k,1,l);
zeta_old   =dmatrix(1,k,1,l);
zeta_guess =dmatrix(1,k,1,l);
theta      =dmatrix(1,k,1,l);
theta_old  =dmatrix(1,k,1,l);
theta_guess=dmatrix(1,k,1,l);
tau_xx     =dmatrix(1,k,1,l);
tau_xy     =dmatrix(1,k,1,l);
tau_yy     =dmatrix(1,k,1,l);
u          =dmatrix(1,k,1,l);
v          =dmatrix(1,k,1,l);
Nu_nf      =dmatrix(1,k,1,l);

    output 	= fopen(outputFileName, "w") ;
	WriteXY_inTopRows( delta_x, delta_x, k, l, output) ; // writes X,Y coordinates in top two rows of the output file, see the header file Quiz.h

/*initialize 2 D matrix*/
for(i=1;i<=k;i++)
{
    for(j=1;j<=l;j++)
    {
    psi[i][j]        = 0;
    psi_old[i][j]    = 0;
    psi_guess[i][j]  = 0;
    zeta[i][j]       = 0;
    zeta_old[i][j]   = 0;
    zeta_guess[i][j] = 0;
    theta[i][j]      = 0;
    theta_old[i][j]  = 0;
    theta_guess[i][j]= 0;
    tau_xx[i][j]     = 0;
    tau_xy[i][j]     = 0;
    tau_yy[i][j]     = 0;
    u[i][j]          = 0;
    v[i][j]          = 0;

    }
}
/*main loop start here*/

/*----------------------------*/

/*time    loop  begins */
error1=1.0; // error for theta time loop
error2=1.0; // error for vorticity time loop
error3=1.0; // error for stream function time loop
itr_time=0;
do
{
error1=0;
error2=0;
error3=0;
itr_time++;
time=(itr_time-1)*delta_t;
/*solution of energy,vorticity,stream function*/
error4=1.0; // error for theta space convergence
error5=1.0; // error for zeta space convergence
error6=1.0; // error for psi space convergence
itr_space=0;
/*space loop iteration begin*/
do
{
error4=0;
error5=0;
error6=0;
itr_space++;
/*left wall & right wall*/
for(j=2;j<=l-1;j++)
{
 theta[1][j] = 1.0;
 theta[k][j] = 0;
}
/*upper & bottom*/
for(i=1;i<=k;i++)
 {
     theta[i][1]=(2.0*theta[i][4]-9.0*theta[i][3]+18.0*theta[i][2])/11.0;       //third order forward difference
     theta[i][l]=(2.0*theta[i][l-3]-9.0*theta[i][l-2]+18.0*theta[i][l-1])/11.0;
 }
 /*interior points*/
 for(i=2;i<=k-1;i++)
{
  for(j=2;j<=l-1;j++)
  {
   if(u[i][j]>0)
   {
       iuw=1.0;
       iue=0;
   }
   else
   {
       iuw=0;
       iue=1.0;
   }
   if(v[i][j]>0)
   {
       ivs=1.0;
       ivn=0;
   }
   else
   {
       ivs=0;
       ivn=1.0;
   }
    ae1 =-(delta_t*f+iue*fabs(u[i][j]*delta_t*(k-1)));
    aw1 =-(delta_t*f+iuw*fabs(u[i][j]*delta_t*(k-1)));
    as1 =-(delta_t*f+ivs*fabs(v[i][j]*delta_t*(l-1)));
    an1 =-(delta_t*f+ivn*fabs(v[i][j]*delta_t*(l-1)));
    ap1 = (1.0-(ae1+aw1+as1+an1));

    theta[i][j]=(theta_old[i][j]-(ae1*theta[i+1][j]+aw1*theta[i-1][j]+as1*theta[i][j-1]+an1*theta[i][j+1]))/ap1;
  }
 }
 /*calculation of shear stress*/
 /*interior points*/
  for(i=2;i<=k-1;i++)
 {
     for(j=2;j<=l-1;j++)
     {
      q1=(u[i+1][j]-u[i-1][j])*(k-1);
      q2=(u[i][j+1]-u[i][j-1])*(l-1)/2.0;
      q3=(v[i+1][j]-v[i-1][j])*(k-1)/2.0;
      q4=(v[i][j+1]-v[i][j-1])*(l-1);
      q5=q2+q3;
     if(fabs(q5)<0.001)
        tau_xy[i][j]=0;
     else
        tau_xy[i][j]=pow(fabs(q5),(n-1))*q5-q5;

    if(fabs(q1)<0.001)
        tau_xx[i][j]=0;
     else
        tau_xx[i][j]=pow(fabs(q1),(n-1))*q1-q1;

    if(fabs(q4)<0.001)
        tau_yy[i][j]=0;
     else
        tau_yy[i][j]=pow(fabs(q4),(n-1))*q4-q4;
     }
 }
 /*boundary condition for shear stress*/
 /*left wall*/
 for(j=2;j<=l-1;j++)
 {
     q6=(-3.0*u[1][j]+4.0*u[2][j]-u[3][j])*(k-1);
     q7=(-3.0*v[1][j]+4.0*v[2][j]-v[3][j])*(k-1)/2.0;

    if(fabs(q6)<0.001)
        tau_xx[1][j]=0;
     else
        tau_xx[1][j]=pow(fabs(q6),(n-1))*q6-q6;

    if(fabs(q7)<0.001)
        tau_xy[1][j]=0;
     else
        tau_xy[1][j]=pow(fabs(q7),(n-1))*q7-q7;
 }
  /*right wall*/
 for(j=2;j<=l-1;j++)
 {
     q8=(u[k-2][j]-4.0*u[k-1][j]+3.0*u[k][j])*(k-1);
     q9=(v[k-2][j]-4.0*v[k-1][j]+3.0*v[k][j])*(k-1)/2.0;

    if(fabs(q8)<0.001)
        tau_xx[k][j]=0;
     else
        tau_xx[k][j]=pow(fabs(q8),(n-1))*q8-q8;

    if(fabs(q9)<0.001)
        tau_xy[k][j]=0;
     else
        tau_xy[k][j]=pow(fabs(q9),(n-1))*q9-q9;
 }
 /*lower wall*/
 for(i=1;i<=k;i++)
 {
  q10=(-3.0*v[i][1]+4.0*v[i][2]-v[i][3])*(l-1);
  q11=(-3.0*u[i][1]+4.0*u[i][2]-u[i][3])*(l-1)/2.0;

     if(fabs(q10)<0.001)
        tau_yy[i][1]=0;
     else
        tau_yy[i][1]=pow(fabs(q10),(n-1))*q10-q10;

    if(fabs(q11)<0.001)
        tau_xy[i][1]=0;
     else
        tau_xy[i][1]=pow(fabs(q11),(n-1))*q11-q11;
 }
 /*upper wall*/
 for(i=1;i<=k;i++)
 {
  q12=(v[i][l-2]-4.0*v[i][l-1]+3.0*v[i][l])*(l-1);
  q13=(u[i][l-2]-4.0*u[i][l-1]+3.0*u[i][l])*(l-1)/2.0;

     if(fabs(q12)<0.001)
        tau_yy[i][k]=0;
     else
        tau_yy[i][k]=pow(fabs(q12),(n-1))*q12-q12;

    if(fabs(q13)<0.001)
        tau_xy[i][k]=0;
     else
        tau_xy[i][k]=pow(fabs(q13),(n-1))*q13-q13;
 }
 /*vorticity equation start*/

 /*vorticity boundary condition */

 /*left & right wall*/
 for(j=2;j<=l-1;j++)
{
 zeta[1][j]  = -2.0*( psi[2][j] )*f;
 zeta[k][j]  = -2.0*( psi[k-1][j] )*f;
}
 /*lower & upper wall*/
 for(i=1;i<=k;i++)
 {
 zeta[i][l]=-2.0*(psi[i][l-1])*f;
 zeta[i][1]=-2.0*(psi[i][2])*f;
 }
/*vorticity calculation at inner place*/

/*vorticity equation*/
for(i=2;i<=k-1;i++)
{
  for(j=2;j<=l-1;j++)
  {
   if(u[i][j]>0)
   {
       iuw=1.0;
       iue=0;
   }
   else
   {
       iuw=0;
       iue=1.0;
   }
   if(v[i][j]>0)
   {
       ivs=1.0;
       ivn=0;
   }
   else
   {
       ivs=0;
       ivn=1.0;
   }

    ae2 =-((PrS)*delta_t*f+iue*fabs(u[i][j]*delta_t*(k-1)));
    aw2 =-((PrS)*delta_t*f+iuw*fabs(u[i][j]*delta_t*(k-1)));
    as2 =-((PrS)*delta_t*f+ivs*fabs(v[i][j]*delta_t*(l-1)));
    an2 =-((PrS)*delta_t*f+ivn*fabs(v[i][j]*delta_t*(l-1)));
    ap2 = (1.0-(ae2+aw2+as2+an2));

    STRESS1 = 0.25 * (tau_yy[i+1][j+1] - tau_yy[i-1][j+1] - tau_yy[i+1][j-1] + tau_yy[i-1][j-1]);
    STRESS2 = 0.25 * (tau_xx[i+1][j+1] - tau_xx[i-1][j+1] - tau_xx[i+1][j-1] + tau_xx[i-1][j-1]);
    STRESS3 = tau_xy[i+1][j] - 2.0 * tau_xy[i][j] + tau_xy[i-1][j];
    STRESS4 = tau_xy[i][j+1] - 2.0 * tau_xy[i][j] + tau_xy[i][j-1];
    STRESS  = STRESS1-STRESS2+STRESS3-STRESS4 ;
    RHS2=(Ra_nf*PrS)*(theta[i+1][j]-theta[i-1][j])*(k-1)/2.0 ;
    zeta[i][j] = (RHS2*delta_t+zeta_old[i][j] - (ae2 * zeta[i+1][j] + aw2*zeta[i-1][j] + as2 * zeta[i][j-1] + an2 * zeta[i][j+1])+(delta_t*PrS*f)*(STRESS))/ap2;

  }
}
/*stream function boundary condition */
for(j=2;j<=l-1;j++)
{
 psi[1][j]   = 0;
 psi[k][j]   = 0;
}
for(i=1;i<=k;i++)
 {
 psi[i][l]=0;
 psi[i][1]=0;
 }

 /*inner iteration begin */
  for(i=2;i<=k-1;i++)
 {
     for(j=2;j<=l-1;j++)
     {
      ae=-delta_t*f ;
      aw=-delta_t*f ;
      as=-delta_t*f ;
      an=-delta_t*f ;
      ap=(1-(ae+aw+as+an));

      psi[i][j]=(delta_t*zeta[i][j] + psi_old[i][j] - (ae * psi[i+1][j] + aw * psi[i-1][j] + as * psi[i][j-1] + an * psi[i][j+1] ))/ap;
     }
 }

 /*calculate velocity profile*/
 for(i=2;i<=k-1;i++)
 {
     for(j=2;j<=l-1;j++)
     {
         u[i][j]=(psi[i][j+1] - psi[i][j-1])*(l-1)/2.0 ;
         v[i][j]=(psi[i-1][j] - psi[i+1][j])*(k-1)/2.0 ;
     }
 }
 /*space convergence for theta,vorticity & stream function */

 for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
  r4=fabs(theta[i][j]-theta_guess[i][j]);
  r5=fabs(zeta[i][j]-zeta_guess[i][j]);
  r6=fabs(psi[i][j]-psi_guess[i][j]);
  if(r4>error4)
    error4=r4;
  if(r5>error5)
    error5=r5;
  if(r6>error6)
    error6=r6;
  /*find maximum error in space criteria*/
  if(error4>error5)
    {
    if(error4>error6)
    errors=error4;
    }
  else
    errors=error6;
  if(error5>error4)
  {
      if(error5>error6)
        errors=error5;
  }
  else
    errors=error6;
 }
 }
/*copy array*/
 for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
      if(error4>Precision)
            theta_guess[i][j]=theta[i][j];
      if(error5>Precision)
        zeta_guess[i][j]=zeta[i][j];
      if(error6>Precision)
        psi_guess[i][j]=psi[i][j];
     }
 }
 if(itr_space%skip_itr==0)
 printf("error space=%lf \n",errors);
if(itr_space>10000)
    break;
}
while(errors>Precision);
/*time convergence criteria*/
for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
  r1=fabs(theta[i][j]-theta_old[i][j]);
  r2=fabs(zeta[i][j]-zeta_old[i][j]);
  r3=fabs(psi[i][j]-psi_old[i][j]);
  if(r1>error1)
    error1=r1;
  if(r2>error2)
    error2=r2;
  if(r3>error3)
    error3=r3;
    if(itr_time==100)
        fprintf(output8,"errorTheta=%lf \t errorZeta=%lf \t errorPsi=%lf",error1,error2,error3);
  /*find maximum error in space criteria*/
  if(error1>error2)
    {
    if(error1>error3)
    errort=error1;
    }
  else
    errort=error3;
  if(error2>error1)
  {
      if(error2>error3)
        errort=error2;
  }
  else
    errort=error3;
 }
 }
 for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
      theta_old[i][j]=relaxation*theta[i][j]+(1-relaxation)*theta_old[i][j];
      zeta_old[i][j] =relaxation*zeta[i][j] + (1-relaxation)*zeta_old[i][j];
      psi_old [i][j] =relaxation*psi[i][j] + (1-relaxation) * psi_old[i][j];
     }
 }
 //if(itr_time % skip_itr == 0)
    printf("time: %lf \t iteration space:%d \t error:%lf \n ",time,itr_space,errort);
    if( itr_time % 20 == 0 )
		{
			Write2Darray2File( u, k, l, output) ;	//  writes a mxn 2D matrix along a row of length mxn in an output file
			Write2Darray2File( theta, k, l, output) ;
			Write2Darray2File( psi, k, l, output) ;
			Write2Darray2File( v, k, l, output) ;
			Write2Darray2File( zeta, k, l, output) ;
		}
    if(itr_time>10000)
        break;
}
while(errort > Precision);
/*calculate nusselt number on hot wall*/
 for(j=1; j<=l; j++)
 {
     Nu_nf[1][j]=(k-1)*(2.0*theta[4][j]-9.0*theta[3][j]+18.0*theta[2][j]-11.0*theta[1][j])/6.0;
 }

/*calculate avg nusselt number*/
for(j=2;j<=l-1;j=j+2)
{
    sume=sume+Nu_nf[1][j];
}
for(j=3;j<=l-1;j=j+2)
{
    sumo=sumo+Nu_nf[1][j];
}
NuAvg=1.0/((k-1)*3.0)*(Nu_nf[1][1]+Nu_nf[1][l]+4.0*sumo+2.0*sume);

/*print values */
 for(i=1;i<=k;i++)
 {
     for(j=1;j<=l;j++)
     {
         fprintf(output1,"%lf,",psi[i][j]);
         fprintf(output2,"%lf,",zeta[i][j]);
         fprintf(output3,"%lf,",theta[i][j]);
         fprintf(output4,"%lf,",u[i][j]);
         fprintf(output5,"%lf,",v[i][j]);

     }
     fprintf(output1,"\n");
     fprintf(output2,"\n");
     fprintf(output3,"\n");
     fprintf(output4,"\n");
     fprintf(output5,"\n");
 }
 for(j=1;j<=l;j++)
 fprintf(output6,"%lf \n",Nu_nf[1][j]);
 fprintf(output7,"avgNusselt=%lf\n",NuAvg);
 fclose(output);
 fclose(output1);
 fclose(output2);
 fclose(output3);
 fclose(output4);
 fclose(output5);
 fclose(output6);
 fclose(output7);
 fclose(output8);
    return 0;
}

