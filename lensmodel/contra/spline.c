/* CUBIC SPLINE */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Nsi 999       /* required by Spline_Install */


void Spline_Install( n, x, y, Ms ) 
int n;
double x[], y[], Ms[];
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


double DSpline( z, n, x, y, Ms, Derivative )
double z, x[], y[], Ms[], *Derivative;
int n;
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
