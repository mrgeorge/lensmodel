/*  Contra: Adiabatic Contraction of Dark Matter Halos

    Contra calculates contraction of a dark matter halo in response to
    the condensation of baryons in its center, based on the modified
    contraction model of Gnedin et al. 2004, ApJ, 616, 16.

    The following assumptions are made: The mass distribution of a
    dark matter halo is spherically symmetric and the velocity
    distribution is isotropic. The final baryon distribution does not
    need to be spherical. The code also calculates a line-of-sight
    velocity dispersion for a tracer population with a given density
    profile and velocity anisotropy.

    Units:  M_vir=1 (arbitrary in numerical profile), R_vir=1, G=1

    Written by Oleg Gnedin

    Last modified: January 20, 2008
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // for memset

#define sqr(x)  ((x)*(x))

#define nmax 1001
#define pi   3.14159265358979323846

int BAR, MAC, n=81, n3=101;   //n=81 n3=101
double c, fb, rb, nu=1., A=0.85, w=0.8, r0=0.03;
double ri[nmax], rf[nmax], mhi[nmax], mhf[nmax],
  rhohi[nmax], rhohf[nmax], mhi_av[nmax], g[nmax],
  logri[nmax], logrf[nmax],
  logmhi[nmax], logrhohi[nmax], logrhohf[nmax],
  Ms[nmax], Msmhi[nmax];


/* BEGIN CODE FROM SPLINE.C */
/* CUBIC SPLINE */

#define Nsi 999       /* required by Spline_Install */


void Spline_Install(int n,double x[],double y[],double Ms[] ) 
{
  int i;
  double c1[Nsi], c2[Nsi], c3[Nsi], d[Nsi], p[Nsi], q[Nsi];

  if( n > Nsi )
    { printf("Spline_Install: n=%d > Nsi=%d\n", n, Nsi); exit(1); }

  for( i=1; i<n-1; i++ )
  {
    c1[i]=(x[i]-x[i-1])/6.0;
    c2[i]=(x[i+1]-x[i-1])/3.0;
    c3[i]=(x[i+1]-x[i])/6.0;
    if( x[i] == x[i-1])
      d[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
    else
    { if( x[i] == x[i+1])
        d[i]=(y[i+2]-y[i+1])/(x[i+2]-x[i+1])-(y[i]-y[i-1])/(x[i]-x[i-1]);
      else
        d[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
    }
  };
  c1[1]=0.0;   c3[n-2]=0.0;
  p[1]=-c3[1]/c2[1];   q[1]=d[1]/c2[1];
  for( i=2; i<n-2; i++ )
  {
    p[i]=-c3[i]/(c2[i]+c1[i]*p[i-1]);
    q[i]=(d[i]-c1[i]*q[i-1])/(c2[i]+c1[i]*p[i-1]);
  };

  Ms[0] = 0.0;  Ms[n-1] = 0.0;
  Ms[n-2] = (d[n-2]-c1[n-2]*q[n-3]) / (c2[n-2]+c1[n-2]*p[n-3]);
  for( i=n-3; i>0; i-- )
    Ms[i] = p[i]*Ms[i+1] + q[i];
}


double DSpline(double z,int n,double x[],double y[],double Ms[],double *Derivative )
{
  double x1, x2, y1, y2, h, dx1, dx2, s, Ms1=0., Ms2=0.;
  int i, i1=0, i2=0, l;
  enum {outleft, in, outright} mode;

  if(x[0]<x[n-1])
  {
    l=1; i=0;
    if(z<=x[0]) mode=outleft;
    else { if(z>=x[n-1]) mode=outright; else mode=in; }
  }
  else
  {
    l=-1; i=n-1;
    if(z>=x[0]) mode=outleft;
    else { if(z<=x[n-1]) mode=outright; else mode=in; }
  }

  switch (mode)
  {
    case outleft : i1=0; i2=1; Ms1=0.0; Ms2=0.0; break;
    case outright: i1=n-2; i2=n-1; Ms1=0.0; Ms2=0.0; break;
    case in      : while(z >= x[i]) i+=l;
			    i1=i-l; i2=i; Ms1=Ms[i1]; Ms2=Ms[i2]; break;
  }

  x1=x[i1]; x2=x[i2]; y1=y[i1]; y2=y[i2];
  h=x2-x1; dx1=z-x1; dx2=x2-z;

  s=Ms1*dx2*dx2*dx2/6.0/h + Ms2*dx1*dx1*dx1/6.0/h
    + (y1-Ms1*h*h/6.0)*dx2/h + (y2-Ms2*h*h/6.0)*dx1/h;

  *Derivative = Ms2*dx1*dx1/2.0/h - Ms1*dx2*dx2/2.0/h
                + (y2-y1)/h - (Ms2-Ms1)*h/6.0;
  return( s );
}


/* FIND ROOT OF AN EQUATION */

#define MAXIT 100

double rtsafe(void (*funcd)(double, double *, double *, double, double), double x1, double x2,
	double xacc, double par1, double par2)
{
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;

	(*funcd)(x1,&fl,&df,par1,par2);
	(*funcd)(x2,&fh,&df,par1,par2);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
	  fprintf(stderr,"Root must be bracketed in rtsafe"); exit(1);}
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	(*funcd)(rts,&f,&df,par1,par2);
	for (j=1;j<=MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		(*funcd)(rts,&f,&df,par1,par2);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	fprintf(stderr,"Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}
#undef MAXIT

/* END CODE FROM SPLINE.C */


double mdm( double x, double c )    /* initial dark matter mass distribution */
{
  double f=0.0;
  
  /* Navarro, Frenk & White (1997) model */
  f = (1.-fb)*(log(1.+c*x) - c*x/(1.+c*x))/(log(1.+c) - c/(1.+c));
  return( f );
}


double y( double x, double *d )
{
  double f;
  if(MAC==1) {                   /* Modified model of AC - Gnedin */
    f = r0*A*pow(x/r0,w);
    *d = A*w*pow(x/r0,w-1.);
  } else if((MAC==0) || (MAC==2)){                    /* Standard model of AC - Blumenthal or Dutton's model*/ 
    f = x;
    *d = 1.0;
  }
  return( f );
}


double mb( double x, double *d )         /* final baryon mass distribution - Mass shells */
{
  double f=0.0, p, a, arg, argc;

  if(BAR==1) {                           /* Exponential disk */
    f = fb*(1.-(1.+x/rb)*exp(-x/rb))/(1.-2.*exp(-1./rb));
    *d = fb*x/sqr(rb)*exp(-x/rb)/(1.-2.*exp(-1./rb));
  }
  if(BAR==2) {                           /* Hernquist model */
    f = fb*sqr(x/(x+rb)*(1.+rb));
    *d = fb*sqr(1.+rb)*2.*x*rb/pow(x+rb,3.);
  }
  return( f );
}


void funcd( double r, double *f, double *df, double mhi, double g )
{
  double x, mbx, dy, dmb;
  x = y(r,&dy);
  mbx = mb(x,&dmb);
  *f = r*pow(mhi + mbx,nu) - g;
  *df = mhi + mbx + r*dmb*dy;
}


int pymain(int MACin, int BARin, double cin, double fbin, double rbin, double nuin, double Ain, double win, int nrad, double rfout[nmax], double rhofout[nmax])
{
  int i;
  double x, d, logmhf, dlogmhf, dlogmhi;

  // set output arrays to zero
  memset(rfout, 0, nmax*sizeof(double));
  memset(rhofout, 0, nmax*sizeof(double));

  // set global vars from inputs
  MAC=MACin;
  BAR=BARin;
  c=cin;
  fb=fbin;
  rb=rbin;
  nu=nuin;
  A=Ain;
  w=win;
  n=nrad;

  if(BAR < 1 || BAR > 2) { fprintf(stderr, ": baryon distribution must be 1 or 2\n"); return 1; }
  if((MAC!=2) && (nu!=1.)) { nu=1.; } // silently set nu=1 if we're not in Dutton mode

  /* Set up initial radial grid */
  for(i=0; i<n3; i++) {
    ri[i] = pow(10., 5.*(double)(i-n+1)/(double)(n-1));
    logri[i] = log(ri[i]);
    mhi[i] = mdm(ri[i],c);
    logmhi[i] = log(mhi[i]);
  }
  Spline_Install(n, logri, logmhi, Msmhi);

  /* Remap initial mass distributions on the average orbital radii */
  for(i=0; i<n; i++) {
    x = y(ri[i],&d);
    mhi_av[i] = mdm(x,c);
    g[i] = pow(mhi_av[i]/(1.-fb),nu)*ri[i]; // Mtot,i * ri
  }

  /* Solve for new radii of spherical dark matter shells */
  for(i=0; i<n; i++) {
    if(MAC == -1) { // No contraction
      rf[i] = ri[i];
      logrf[i] = logri[i];
    } else {
      rf[i] = rtsafe(funcd, 1.e-3*ri[i], 1.e4*ri[i], 1.e-10, mhi_av[i], g[i]);
      logrf[i] = log(rf[i]);
    }
    rfout[i]=rf[i];
  }
  Spline_Install(n, logrf, logmhi, Ms);

  /* Initial and final dark matter densities */
  for(i=0; i<n; i++) {
    logmhf = DSpline(logri[i], n, logrf, logmhi, Ms, &dlogmhf);
    mhf[i] = exp(logmhf);
    rhohf[i] = dlogmhf*mhf[i]/(4.*pi*pow(ri[i],3));
    rhofout[i]=rhohf[i];
    logrhohf[i] = log(rhohf[i]);
    DSpline(logri[i], n, logri, logmhi, Msmhi, &dlogmhi);
    rhohi[i] = dlogmhi*mhi[i]/(4.*pi*pow(ri[i],3));
    logrhohi[i] = log(rhohi[i]);
  }
  return 0;
}
