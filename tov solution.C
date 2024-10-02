// Call this routine as ./tov_solver_polytrope_mvsr gam nstep

// gam is tha adiabatic index for the equation of state, i.e., P = kappa rho^(gam)
//       the code will work out the proper value of kappa

// nstep determines the stepsize, which is set to 10^{-nstep), so 10^nstep steps per unit distance


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <math.h>

// Code tolerance level on compactness
double TOL=1.0e-6;

// Physical constants
double msun = 1.98847e33;
double kilometer=1.0e5;
double ggrav=6.67408e-8;
double speedlight=2.9979e10;


using namespace std;

void rk4(double x, double* y, double step, double k, double gam);
void derivs(double x, double* v, double* y, double k, double gam);

int main(int argc, char** argv) {

  // Rather than aim for a specific mass and radius, we range over central densities/enthalpies
  // and find compactnesses; searching for specific models happens in a different routine

  
  double mns,rns,compactns,rns_units,rrk4,mbrk4,mrk4;
  int nchoice;
  //  cout<<"Mass of NS in solar masses?"<<endl;
  //  cin >>mns;
  //  mns=1.4;
  mns=atof(argv[1]);

  //  cout<<"1 to set the compactness, 2 to set the radius"<<endl;
  //  cin >> nchoice;
  nchoice=atoi(argv[2]);
  
  // Use the radius to calculate the compactness

  if(nchoice==1) {
   //Use the compactness to calculate the radius
    
    //    cout<<"Compactness GM/rc^2?"<<endl;
    //    cin>>compactns;
    compactns=atof(argv[3]);

    rns=ggrav*mns*msun/compactns/speedlight/speedlight;
    rns_units = mns/compactns;

  } else if (nchoice==2) {

    //    cout<<" Radius of NS in km?"<<endl;
    //    cin>> rns;
    //    rns=14;
    rns=atof(argv[3]);
    
    rns *= kilometer;
    compactns = ggrav*mns*msun / rns / speedlight/speedlight;
    rns_units = mns/compactns;
  }
    
  cout<<"Mass: "<<mns<<" Msun;  Radius: "<<rns<<" in cm, "<<rns/kilometer<<" km, code units:"<<rns_units<<" ; compactness="<<compactns<<endl;
  
  double gam = atof(argv[4]);
  cout<<"Gamma value:"<<gam<<endl;
									      
  double step=pow(10,-1*atof(argv[5]));
  // The next line allows the code to go out to R=100
  int np = (int) (100.0/step);

  double* rho = new double [np];
  double* massg = new double [np];
  double* massb = new double [np];
  
  int i,j;

  //initial guess
  double k=1;
  double kold=-1.0;
  double compact,compactold;

  int nloop;
  int nloopmax=50;

  double error=1.0;

  double rhoc=1.0;
  
  for (nloop=0; nloop<nloopmax && error>TOL; nloop++) {
    
    for(i=0; i<np; i++) {
      rho[i] = 0.0;
      massg[i]=0.;
      massb[i]=0.;
    }
    
    rho[0]=rhoc;
    massg[0]=0.0;
    massb[0]=0.0;

    double entc=log(1.0+gam/(gam-1.0)*k*pow(rhoc,gam-1));
    
    double y[3];
    
    int done=0;
    double xsurf;

    for (i=1; i<np && done==0; i++) {
      double xpt=(i-1)*step;
      y[0]=rho[i-1];
      y[1]=massg[i-1];
      y[2]=massb[i-1];
	
      rk4(xpt,y,step,k,gam);
      
      if(y[0]<=0) {
	done=1;
	rho[i] = 0.0;
	//We need a better extrapolation approach than linear extrapolation
	double olderrhopow = pow(rho[i-3],gam-1.0);
	double oldrhopow = pow(rho[i-2],gam-1.0);
	double newrhopow = pow(rho[i-1],gam-1.0);

	double fracstep = newrhopow/(oldrhopow-newrhopow);

	//	cout<<"Fracstep:"<<fracstep<<" "<<newrhopow<<" "<<oldrhopow<<" "<<olderrhopow<<" rhovals:"<<rho[i-1]<<" "<<rho[i-2]<<" "<<rho[i-3]<<endl;
	xsurf=xpt+step*fracstep;
	
	double compactcheck=massg[i-1]/xpt;
	double diffcheck=(gam-1.0)/k/gam*massg[i-1]/pow(xpt,2)/(1.0-2.0*compactcheck);
	
	//cout<<"Check:"<<diffcheck<<" "<<diffcheck*step<<" "<<oldrhopow-newrhopow<<" "<<olderrhopow-oldrhopow<<endl;
        double massfrac=pow(fracstep,gam/(gam-1.0))/(pow(1.0+fracstep,gam/(gam-1.0))-pow(fracstep,gam/(gam-1.0)));
        //cout<<"massfrac:"<<massfrac<<endl;

        massg[i] = massg[i-1]+massfrac*(massg[i-1]-massg[i-2]);
	massb[i] = massb[i-1]+massfrac*(massb[i-1]-massb[i-2]);
	
	cout<<"r:"<<xsurf<<" massg:"<<massg[i]<<" "<<" massb:"<<massb[i]<<" compact:"<<massg[i]/xsurf<<" rhoc:"<<rhoc<<" entc:"<<entc<<endl;

	//	cout<<"checks: rratio"<<" "<<rprop[i]/xsurf<<" mratio:"<<massb[i]/massg[i]<<endl;

	rrk4=xsurf;
	mrk4=massg[i];
	mbrk4=massb[i];
	
      } else {
	rho[i]=y[0];
	massg[i]=y[1];
	massb[i]=y[2];
	xsurf = xpt+step;

      }
      
      compact = mrk4/xsurf;

    }
    
    
    error=fabs(compact-compactns);
    cout<<setprecision(12)<<"Current run: kappa:"<<k<<" rhoc:"<<rhoc<<" compact:"<<compact<<" error:"<<error<<endl;

    //    cout<<"Look here: log k:"<<log10(k)<<" "<<log10(compact)<<endl;

    cout<<"masses: target mns:"<<mns<<" actual:"<<mrk4<<"    radii -- target rns:"<<rns_units<<" actual:"<<rrk4<<endl;
    
    k=k*pow(mns/mrk4,2.0-gam)/pow(rns_units/xsurf,4.0-3.0*gam);
    rhoc=rhoc*mns/mrk4/pow(rns_units/rrk4,3.0);

    cout<<endl;
    cout<<endl;
  }
  /* This part of the code is now obsolete, we hope
    // Now check, by rescaling:
  cout<<"Rescaling with current M,R:"<<mrk4<<" "<<rrk4<<endl;
  double factor = mns/mrk4;
  double rhoc_new = 1.0/factor/factor;
  k *= pow(factor, 2*gam-2.0);
  cout<<"Factor:"<<factor<<endl;  */

  compact = mrk4/rrk4;

  double compval=(1.0-compact+sqrt(1.0-2.0*compact))/2.0;
  double rareal=rrk4*ggrav*msun/speedlight/speedlight/kilometer;
  double rcoord=compval*rareal;
  cout<<"In normal units -- areal radius in km:"<<rareal<<" coord radius in km:"<<rcoord<<endl;
  
  //3720.8 is the ratio between a solar mass per [solar mass in geometrized distance units]^3
  // and Lorene's density reference, nuclear density of 1.66e14 g/cm^3
  double klorene=k*pow(3720.8,1.0-gam);
  double denslorene=rhoc*3720.8;

  // enthalpy is defined in Lorene as
  //log(enthalpy) = log (1 + gam*epsilom) = log (1 + gam / (gam-1) P_c / rho_c)
  double entlorene=log(1.0+ gam / (gam-1.0)*klorene*pow(denslorene,gam-1.0));
  cout<<"Use the following with Lorene: k="<<klorene<<";   enthalpy log(gam/gam-1)*k*rho^(gam-1)):"<<entlorene<<endl;
  cout<<endl;
  cout<<endl;


  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //-------------------------------------------------------------

  int iloop;
  int iloopmax=10;
  double densfactor=10.0;
  
  for(iloop=0; iloop<=iloopmax; iloop++) {
    double densloop = rhoc*pow(densfactor,1.0-2.0*(iloopmax-iloop)/(1.0*iloopmax));

    for(i=0; i<np; i++) {
      rho[i] = 0.0;
      massg[i]=0.;
      massb[i]=0.;
    }
    
    rho[0]=densloop;
    massg[0]=0.0;
    massb[0]=0.0;

    double entc=log(1.0+gam/(gam-1.0)*k*pow(densloop,gam-1));
    
    double y[3];
    
    int done=0;
    double xsurf;

    for (i=1; i<np && done==0; i++) {
      double xpt=(i-1)*step;
      y[0]=rho[i-1];
      y[1]=massg[i-1];
      y[2]=massb[i-1];
	
      rk4(xpt,y,step,k,gam);
      
      if(y[0]<=0) {
	done=1;
	rho[i] = 0.0;
	//We need a better extrapolation approach than linear extrapolation
	double olderrhopow = pow(rho[i-3],gam-1.0);
	double oldrhopow = pow(rho[i-2],gam-1.0);
	double newrhopow = pow(rho[i-1],gam-1.0);

	double fracstep = newrhopow/(oldrhopow-newrhopow);

	//	cout<<"Fracstep:"<<fracstep<<" "<<newrhopow<<" "<<oldrhopow<<" "<<olderrhopow<<" rhovals:"<<rho[i-1]<<" "<<rho[i-2]<<" "<<rho[i-3]<<endl;
	xsurf=xpt+step*fracstep;
	
	double compactcheck=massg[i-1]/xpt;
	double diffcheck=(gam-1.0)/k/gam*massg[i-1]/pow(xpt,2)/(1.0-2.0*compactcheck);
	
	//cout<<"Check:"<<diffcheck<<" "<<diffcheck*step<<" "<<oldrhopow-newrhopow<<" "<<olderrhopow-oldrhopow<<endl;
        double massfrac=pow(fracstep,gam/(gam-1.0))/(pow(1.0+fracstep,gam/(gam-1.0))-pow(fracstep,gam/(gam-1.0)));
        //cout<<"massfrac:"<<massfrac<<endl;

        massg[i] = massg[i-1]+massfrac*(massg[i-1]-massg[i-2]);
	massb[i] = massb[i-1]+massfrac*(massb[i-1]-massb[i-2]);

	double rinkm = xsurf*ggrav*msun/speedlight/speedlight/kilometer;
	cout<<"r:"<<xsurf<<" r in km:"<<rinkm<<" massg:"<<massg[i]<<" "<<" massb:"<<massb[i]<<" compact:"<<massg[i]/xsurf<<" rhoc:"<<densloop<<" entc:"<<entc<<" kappa:"<<k<<endl;

	//	cout<<"checks: rratio"<<" "<<rprop[i]/xsurf<<" mratio:"<<massb[i]/massg[i]<<endl;

	rrk4=xsurf;
	mrk4=massg[i];
	mbrk4=massb[i];
	
      } else {
	rho[i]=y[0];
	massg[i]=y[1];
	massb[i]=y[2];
	xsurf = xpt+step;

      }
      
      compact = mrk4/xsurf;

    }
    
    
    //    cout<<setprecision(12)<<"compact:"<<compact<<endl;

  }
    

  
  
  return 0;

  }

void rk4(double x, double* y, double step, double k, double gam) {

  double h=step/2.0;    
  double t1[3], t2[3], t3[3], k1[3], k2[3], k3[3], k4[3], v[3], dv[3];
  
  int i;

  int done = 0;

  double xtry=x;

  derivs(xtry,v,y,k,gam);
  for (i=0; i<3; i++) {
    k1[i]=step*v[i];
    t1[i] = y[i]+0.5*k1[i];
  }
  if(t1[0]<0) {
    done=1;
    y[0]=0.;
    y[1]=0.;
  }
  xtry+=h;
  derivs(xtry,v,t1,k,gam);
  for (i=0; i<3; i++) {
    k2[i]=step*v[i];
    t2[i] = y[i]+0.5*k2[i];
  }
  if(t2[0]<0) {
    done=1;
    y[0]=0.;
    y[1]=0.;
  }
  derivs(xtry,v,t2,k,gam);	
  for (i=0; i<3; i++) {
    k3[i]=step*v[i];
    t3[i] = y[i]+k3[i];
  }
  if(t3[0]<0) {
    done=1;
    y[0]=0.;
    y[1]=0.;
  }
  xtry+=h;
  derivs(xtry,v,t3,k,gam);
  for (i=0; i<4; i++) {
    k4[i] = step*v[i];
  }

  for (i=0; i<3; i++) y[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;

  if(done==1)y[0]=0.;
}

void derivs(double x, double* v, double* y, double k, double gam) {

  double press,dpdrho,eps;
  double rho = y[0];
  double mass = y[1];
  if(y[0]>0) {
    press = k*pow(rho,gam);
    // This is rho*epsilon, not the specific version
    eps = press/(gam-1.0);
    dpdrho = k*gam*pow(rho,gam-1.0);
  } else {
    press = 0;
    eps=0;
    dpdrho = 0;
  }

  if(x>0 && rho>0) {
      v[0] = -1.0/x/x*(rho+eps+press)*(mass+4.0*M_PI*x*x*x*press)/(1.0-2.0*mass/x)/dpdrho;
      v[1] = 4.0*M_PI*x*x*(rho+eps);
      v[2] = 4.0*M_PI*x*x*rho/sqrt(1.0-2.0*mass/x);
  } else {
      v[0]=0.;
      v[1]=0.;
      v[2]=0.;
  }
  
}
