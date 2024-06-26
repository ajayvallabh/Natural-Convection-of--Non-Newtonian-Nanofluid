//correationfor power law index  of Y-Al2O3/CMC-water nanofluid
double powerLawInd_Al2O3_CMC_Water(double Th,double Tc,double phi)
{
 double Tavg,n; //where  n is power law index
 Tavg = (Th+Tc)/2;
    if(phi==0)
    n = (2.71e-3)*Tavg + 0.47545;
    if(phi==0.001)
    n=(2.5e-5)*Tavg*Tavg + (5.4e-4)*Tavg+0.493175;
    if(phi==0.005)
    n=(5.0e-5)*Tavg*Tavg -(8.0e-4)*Tavg+0.4841;
    if(phi==0.015)
    n=(8.5e-5)*Tavg*Tavg -(2.2e-3)*Tavg+0.498375;
    return n;
}
//Correlation for power-law index of CuO/CMC- Water nanofluid
double powerLawInd_CuO_CMC_Water(double Th,double Tc,double phi)
{
 double Tavg,n; //where  n is power law index
 Tavg = (Th+Tc)/2;
    if(phi==0)
    n = (2.71e-3)*Tavg + 0.47545;
    if(phi==0.001)
    n=(6.5e-5)*Tavg*Tavg - (1.78e-3)*Tavg + 0.539275;
    if(phi==0.005)
    n=(4.0e-5)*Tavg*Tavg + (2.2e-4) * Tavg + 0.5069;
    if(phi==0.015)
    n=-(3.25e-5)*Tavg*Tavg + (5.58e-3)*Tavg+0.4386625;
    if(phi==0.03)
    n =-(3.42625e-5)*Tavg*Tavg + (5.13e-3)*Tavg+0.5400292;
    if(phi==0.04)
    n=2.87791e-3*Tavg + 0.5068315;
    return n;
}
//Correlation for power-law index of TiO2/CMC- Water nanofluid
double powerLawInd_TiO2_CMC_Water(double Th,double Tc,double phi)
{
 double Tavg,n; //where  n is power law index
 Tavg = (Th+Tc)/2;
    if(phi==0)
    n = (2.71e-3)*Tavg + 0.47545;
    if(phi==0.001)
    n = (2.5e-5)*Tavg*Tavg + (5.0e-4)*Tavg + 0.524875;
    if(phi==0.005)
    n = 2.36e-3*Tavg + 0.4812 ;
    if(phi==0.015)
    n = 1.48e-3*Tavg + 0.4846 ;
    if(phi==0.03)
    n=(2.51334e-3)*Tavg + 0.4405998;
    if(phi==0.040)
    n = -(2.8335e-5)*Tavg*Tavg + 2.6801e-3*Tavg + 0.4358069 ;
    return n;
}
//Calculation for Consistency index of Y-Al2O3/CMC- Water nanofluid
double ConsistencyInd_Al203_CMC_Water(double Th,double Tc,double phi)
{
    double N,Tavg;
    Tavg = (Th+Tc)/2;
    if(phi==0)
    N=(7.125e-5)*Tavg*Tavg - (7.98e-3)*Tavg + 0.2964938;
    if(phi==0.001)
    N=(4.75e-5)*Tavg*Tavg - (6.28e-3)*Tavg + 0.2924625;
    if(phi==0.005)
    N=-(4.18e-3)*Tavg + 0.3149;
    if(phi==0.015)
    N=-(5.98e-3)*Tavg + 0.3759 ;
    return N ;
}

 //Calculation for Consistency index of CuO/CMC- Water nanofluid
double ConsistencyInd_CuO_CMC_Water(double Th,double Tc,double phi)
{
    double N,Tavg;
    Tavg = (Th+Tc)/2;
    if(phi==0)
    N =(7.125e-5)*Tavg*Tavg - (7.98e-3)*Tavg + 0.2964938;
    if(phi==0.001)
    N =-(3.191e-3)*Tavg + 0.225805;
    if(phi==0.005)
    N = 3.225e-5*Tavg*Tavg-(5.708e-3)*Tavg + 0.2704088;
    if(phi==0.015)
    N =(1.1525e-4)*Tavg*Tavg - 0.010856*Tavg + 0.3759;
    if(phi==0.030)
    N = 3.4824e-4*Tavg*Tavg + 4.975e-3*Tavg + 0.2128093;
    if(phi==0.040)
    N= 8.333475e-4*Tavg*Tavg + 0.008816012*Tavg + 0.2901319;
    return N ;
}

//Calculation for Consistency index of TiO2/CMC- Water nanofluid
double ConsistencyInd_TiO2_CMC_Water(double Th,double Tc,double phi)
{
    double N,Tavg;
    Tavg = (Th+Tc)/2;
    if(phi==0)
    N=(7.125e-5)*Tavg*Tavg - (7.98e-3)*Tavg + 0.2964938;
    if(phi==0.001)
    N=(3.675e-5)*Tavg*Tavg - (5.124e-3)*Tavg + 0.2332263;
    if(phi==0.005)
    N = 9.4e-5*Tavg*Tavg-(9.472e-3)*Tavg + 0.33171;
    if(phi==0.015)
    N =(6.5e-5)*Tavg*Tavg - 8.02e-3*Tavg + 0.363975;
    if(phi==0.030)
    N = 1.121425e-4*Tavg*Tavg - 0.01270156*Tavg + 0.4868425;
    if(phi==0.040)
    N= 2.2242825e-4*Tavg*Tavg - 0.01895778*Tavg + 0.7000581;
    return N ;
}
//calculate thermal conductivity of base fluid which depend on temperature
//corelation is given by Carezzato
double Basefluid_Conductivity(double Th,double Tc)
{
    double Kf,Tavg;
    Tavg = (Th+Tc)/2;
    Kf=-7.74e-6*Tavg*Tavg+0.00186*Tavg+0.574;
    return Kf ;
}
//calculate viscosity of base fluid by  correlation
long double Basefluid_Viscosity(double Th,double Tc)
{
    double Muf,Tavg ;
    Tavg = (Th+Tc)/2;
    Muf = 3.12531e-11*pow((Tavg+273),4) - 4.28147e-8*pow((Tavg+273),3) + 2.20399e-5*pow((Tavg+273),2)-5.0587e-3*(Tavg+273)+0.437721 ;
   return Muf;
}
//calculate thermal conductivity of Nanofluid by Maxwell model
double Effective_Thermal_Conductivity( double Kf, double Kp, double phi)
{
 double Kef;
  Kef = Kf*( ( Kp + 2 * Kf + 2 *(Kp - Kf) * phi)/(Kp + 2*Kf - 2*(Kp - Kf) * phi));
 return Kef ;

}

// calculate effective viscosity of nanofluid by Brinkman Model
double Effective_Viscosity( double Muf, double phi)
{
    double Mu_nf;
    Mu_nf= (Muf / pow((1-phi) , 2.5 ));
    return Mu_nf ;
}

// Calculate nanofluid density by correlation
double Nanofluid_Density(double rho_f, double rho_s, double phi)
{
    double rho_nf ;
    rho_nf = ((1-phi) * rho_f + phi * rho_s) ;
    return rho_nf ;
}

// Calculate heat capacity of Nanofluid  by Correlation
double Heat_Capacity(double rho_f, double Cpf,double rho_s, double Cps, double phi )
{
    double rhoCp_nf ;
    rhoCp_nf = (1-phi) * (rho_f * Cpf) + phi * (rho_s * Cps) ;
    return rhoCp_nf ;
}
// Calculate Volumetric thermal expansion Coefficient of nanofluid  by Correlation
double Volumetric_Thr_Exp( double rho_f, double beta_f, double rho_s, double beta_s,double rho_nf, double phi)
{
    double rho_beta_nf, beta_nf ;
    rho_beta_nf = (1-phi) * ( rho_f * beta_f) + phi * (rho_s * beta_s) ;
    beta_nf = rho_beta_nf/rho_nf ;
    return beta_nf ;
}
// Calculate thermal diffusivity
double Thermal_Diffusivity(double Kef,double rhoCp_nf)
{
    double alpha_nf ;
    alpha_nf=( Kef / rhoCp_nf );
    return alpha_nf ;
}
// Calculate kinematic viscosity of Newtonian fluid
double Kinematic_Viscosity( double N,double rho_nf )
{
    double nu_nf;
    nu_nf= N / rho_nf ;
    return nu_nf ;
}

// calculation for Dimensionless number Prandtl number for nanofluid
double Prandtl_Number(double nu_nf,double alpha_nf,double N  )
{
    double Pr_nf ;
    Pr_nf = nu_nf / alpha_nf ;
    return Pr_nf ;
}

// calculation for Dimensionless Prandtl number for Non-Newtonian nanofluid
double PrandtlStar_Number(double N, double n, double L,double alpha_nf,double rho_nf)
{
    double  PrStar ;
    PrStar = ( pow(alpha_nf,(n-2)) * N) / ( pow(L,(2*n-2) ) * rho_nf) ;
    return PrStar ;
}

// calculation for Dimensionless Rayleigh Number
double Rayleigh_Number(double alpha_nf, double nu_nf,double beta_nf , double Th, double Tc, double gravity, double L,double n)
{
    double Ra_nf  ;
    Ra_nf = ((Th-Tc) * pow(L,(2.0*n+1)) * gravity * beta_nf ) / (pow(alpha_nf,n)* nu_nf);
    return Ra_nf ;
}
// writes 2D kxl array to a file in a single row, row-wise
void Write2Darray2File( double **x, int k, int l, FILE *output)
{
	int i, j ;
	for(j = 1 ; j <= l ; j++)
	{
		for( i = 1; i <= k ; i++ )
		{
 			fprintf(output, "%f \t" , x[i][j] ) ;
		}
	}
	fprintf( output, "\n") ;
}
// writes 1D 1 x k array to a file in a row
void Write1Darray2File( double *x, int k, FILE *output)
{
	int i ;
	for( i = 1; i <= k ; i++ )
	{
		fprintf(output, "%f \t" , x[i] ) ;
	}
	fprintf( output, "\n") ;
}

//  write x coordinates in a row
void WriteX_inTopRow( double delta_x, int k, FILE *output)
{
	int i ;
	double x;

	for( i = 1; i <= k ; i++ )
	{
		x = ( i-1 ) * delta_x ;
		fprintf(output, "%f \t" , x ) ;
	}
	fprintf( output, "\n") ;
}

//  write x and y coordinates in two rows respectively, x, the first
void WriteXY_inTopRows( double delta_x, double delta_y, int k, int l, FILE *output)
{
	int i, j ;
	double x, y ;

	for(j = 1 ; j <= k ; j++)
	{
		for( i = 1; i <= l ; i++ )
		{
			x = ( i-1 ) * delta_x ;
 			fprintf(output, "%f \t" , x ) ;
		}
	}
	fprintf( output, "\n") ;


	for(j = 1 ; j <= l ; j++)
	{
		y = ( j-1 ) * delta_y ;  // delta_y = delta_x
		for( i = 1; i <= k ; i++ )
		{
 			fprintf(output, "%f \t" , y ) ;
		}
	}
	fprintf( output, "\n") ;
}
