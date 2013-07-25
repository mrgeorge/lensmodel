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

#define sqr(x)  ((x)*(x))

#define nmax 1001
#define pi   3.14159265358979323846

int DM, BAR, MAC, NUM, ANIS, n=81, n3=101;   //n=81 n3=101
//double c, fb, rb, ra, ser_dm, ser_b, TRACE, Rmax = 3.0; // replaced by MRG following contra_aw.c
double c, fb, rb, ra, ser_dm, ser_b, TRACE, Rmax=3., A=0.85, w=0.8, r0=0.03;
double ri[nmax], rf[nmax], mhi[nmax], mhf[nmax],
  rhohi[nmax], rhohf[nmax], mhi_av[nmax], g[nmax],
  rhot[nmax], mtot[nmax], sigma_los[nmax],
  logri[nmax], logrf[nmax], logrbf[nmax], logmbi[nmax], logmbf[nmax],
  logmhi[nmax], logrhohi[nmax], logrhohf[nmax], dlogrhohi[nmax], dlogrhohf[nmax],
  Ms[nmax], Msrhoi[nmax], Msrhof[nmax], Msmhi[nmax], Msmbi[nmax], Msmbf[nmax];

//extern double rtsafe(); // Removed by MRG since spline.c is copied into this file
//extern void Spline_Install(); // Removed by MRG since spline.c is copied into this file
//extern double DSpline(); // Removed by MRG since spline.c is copied into this file
//extern void losdispersion_(); // Removed by MRG to avoid linking to los.f
//extern double gammp_(); // Removed by MRG to avoid linking to los.f

char *dmprofile[] = { "", "Navarro, Frenk & White (1997)", 
  "Spherical Sersic", "Kazantzidis et al. (2004)" };
char *barprofile[] = { "", "Exponential", "Hernquist (1990)", 
  "Jaffe (1983)", "Sersic (1968)" };
char *trprofile[] = { "power law", "Exponential", "Hernquist (1990)", 
  "Jaffe (1983)", "Plummer (1911)", "NFW (1997)", "Sersic (1968)", 
  "Spherical Sersic" };


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
  double f=0.0, a, arg, argc, logmh, dlogmh;
  
  if(DM==1) {                       /* Navarro, Frenk & White (1997) model */
    //    f = (log(1.+c*x) - c*x/(1.+c*x))/(log(1.+c) - c/(1.+c)); // replaced by MRG following contra_aw.c
    f = (1.-fb)*(log(1.+c*x) - c*x/(1.+c*x))/(log(1.+c) - c/(1.+c));
  }
  if(DM==2) {                       /* Navarro et al. (2004): */
    printf("ERROR: DM=2 option disabled"); // Added by MRG to avoid linking to los.f
    exit(1);
    //    a = 3.*ser_dm;                  /* spherical Sersic model */
    //    arg = 2.*ser_dm*pow(x*c, 1./ser_dm);
    //    argc = 2.*ser_dm*pow(c, 1./ser_dm);
    //    //    f = gammp_(&a, &arg)/gammp_(&a, &argc); // replaced by MRG following contra_aw.c
    //    f = (1.-fb)*gammp_(&a, &arg)/gammp_(&a, &argc);
  }
  if(DM==3) {                       /* Kazantzidis et al. (2004) profile */
    printf("ERROR: DM=3 option disabled"); // Added by MRG to avoid linking to los.f
    exit(1);
    //    a = 2.0;
    //    arg = x*c;
    //    argc = c;
    //    //    f = gammp_(&a, &arg)/gammp_(&a, &argc); // replaced by MRG following contra_aw.c
    //    f = (1.-fb)*gammp_(&a, &arg)/gammp_(&a, &argc);
  }
  if(NUM) {
    logmh = DSpline(log(x), n, logri, logmhi, Msmhi, &dlogmh);
    f = exp(logmh);
  }
  return( f );
}


double y( double x, double *d )
{
  double f;
  if(MAC) {                   /* Modified model of AC */
    //    f = 0.85*pow(x,0.8); // replaced by MRG following contra_aw.c
    //    *d = 0.85*0.8*pow(x,-0.2); // replaced by MRG following contra_aw.c
    f = r0*A*pow(x/r0,w);
    *d = A*w*pow(x/r0,w-1.);
  } else {                    /* Standard model of AC */ 
    f = x;
    *d = 1.0;
  }
  return( f );
}


double mb( double x, double *d )         /* final baryon mass distribution */
{
  double f=0.0, p, a, arg, argc, logmb, dlogmb;

  if(BAR==1) {                           /* Exponential disk */
    f = fb*(1.-(1.+x/rb)*exp(-x/rb))/(1.-2.*exp(-1./rb));
    *d = fb*x/sqr(rb)*exp(-x/rb)/(1.-2.*exp(-1./rb));
  }
  if(BAR==2) {                           /* Hernquist model */
    f = fb*sqr(x/(x+rb)*(1.+rb));
    *d = fb*sqr(1.+rb)*2.*x*rb/pow(x+rb,3.);
  }
  if(BAR==3) {                           /* Jaffe model */
    f = fb*x/(x+rb)*(1.+rb);
    *d = fb*(1.+rb)*rb/sqr(x+rb);
  }
  if(BAR==4) {                           /* Sersic model */
    printf("ERROR: BAR=4 option disabled"); // Added by MRG to avoid linking to los.f
    exit(1);
    //    p = 1.0 - 0.6097/ser_b + 0.05563/sqr(ser_b);
    //    a = (3.-p)*ser_b;
    //    arg = pow(x/rb, 1./ser_b);
    //    argc = pow(1./rb, 1./ser_b);
    //    f = fb*gammp_(&a, &arg)/gammp_(&a, &argc);
    //    *d = fb/gammp_(&a, &argc)*exp(-arg)*pow(arg,a)/(ser_b*x);
  }
  if(NUM) {
    logmb = DSpline(log(x), n, logrbf, logmbf, Msmbf, &dlogmb);
    f = exp(logmb);
    *d = dlogmb*f/x;
  }
  return( f );
}


void funcd( double r, double *f, double *df, double mhi, double g )
{
  double x, mbx, dy, dmb;
  x = y(r,&dy);
  mbx = mb(x,&dmb);
  *f = r*(mhi + mbx) - g;
  *df = mhi + mbx + r*dmb*dy;
}


double tracer_density( double x )              /* tracer population density */
{
  double f=0.0, p;
  if(TRACE == 1) f = exp(-x/rb);               /* Exponential disk */
  if(TRACE == 2) f = 1./x/pow(x+rb,3.);        /* Hernquist model */
  if(TRACE == 3) f = 1./sqr(x)/sqr(x+rb);      /* Jaffe model */
  if(TRACE == 4) f = 1./pow(1.+sqr(x/rb),2.5); /* Plummer model */
  if(TRACE == 5) f = 1./x/pow(x+rb,2.);        /* NFW model */
  if(TRACE == 6) {                             /* Sersic model */
    p = 1.0 - 0.6097/ser_b + 0.05563/sqr(ser_b);
    f = pow(x/rb, -p)*exp(-pow(x/rb, 1./ser_b));
  }
  if(TRACE == 7)                               /* spherical Sersic model */
    f = exp(-2.*ser_dm*pow(x/rb, 1./ser_dm));
  if(TRACE <= 0) f = pow(x, (double)TRACE);    /* power law */
  return( f );
}


int pymain(int MACin, int DMin, int BARin, double TRACEin, int ANISin, double cin, double ser_dmin, double fbin, double rbin, double ser_bin, double rain, double Ain, double win)
{
  int i, i0, tr;
  double x, d, mi_av, mbi, mbi_av, mbf, rbf,
    logmhf, dlogmhf, logmhi_av, logmbi_av, dlogmhi;
  char infile[80];
  FILE *in;

  int NUM=0;

  // set global vars from inputs
  MAC=MACin;
  DM=DMin;
  BAR=BARin;
  TRACE=TRACEin;
  ANIS=ANISin;
  c=cin;
  ser_dm=ser_dmin;
  fb=fbin;
  rb=rbin;
  ser_b=ser_bin;
  ra=rain;
  A=Ain;
  w=win;

  if((DM < 1 || DM > 3) && !NUM) { fprintf(stderr, ": DM distribution must be 1, 2 or 3\n"); return 1; }
  if((BAR < 1 || BAR > 4) && !NUM) { fprintf(stderr, ": baryon distribution must be 1, 2, 3 or 4\n"); return 1; }
  if(ANIS < 0 || ANIS > 3) { fprintf(stderr, ": anisotropy parameter must be 0, 1, 2 or 3\n"); return 1; }
  if(TRACE > 7) { fprintf(stderr, ": tracer population must be 1, 2, 3, 4, 5, 6, 7 or negative as a power law slope\n"); return 1; }

  //  printf("#"); for(i=0; i<argc; i++) printf(" %s", argv[i]); printf("\n");
  if(MAC >= 1) printf("# Modified model of Adiabatic Contraction\n");
  if(MAC == 0) printf("# Standard model of Adiabatic Contraction\n");
  if(MAC == -1) printf("# No adiabatic contraction\n");
  if(NUM) {
    printf("# reading numerical profile from <%s>\n", infile);
  } else {
    if(DM==1 || DM==3) printf("# Dark matter profile: %s model, c=%g\n", dmprofile[DM], c);
    if(DM==2) printf("# Dark matter profile: %s model, c=%g n=%g\n", dmprofile[DM], c, ser_dm);
    if(BAR==1 || BAR==2 || BAR==3) 
      printf("# Baryon distribution: %s model, fb=%g rb=%g\n", barprofile[BAR], fb, rb);
    if(BAR==4) 
      printf("# Baryon distribution: %s model, fb=%g rb=%g n=%g\n", barprofile[BAR], fb, rb, ser_b);
    if(TRACE>0) tr=(int)(TRACE+0.5); else tr=0;
    if(ANIS==0) printf("# Tracer population: %s model with isotropic velocity distribution\n", trprofile[tr]);
    if(ANIS==1) printf("# Tracer population: %s model with constant velocity anisotropy, beta=%g\n", trprofile[tr], ra);
    if(ANIS==2) printf("# Tracer population: %s model with Mamon-Lokas velocity anisotropy, ra=%g\n", trprofile[tr], ra);
    if(ANIS==3) printf("# Tracer population: %s model with Osipkov-Merritt velocity anisotropy, ra=%g\n", trprofile[tr], ra);
    if(TRACE<6 && BAR<4 && ser_b>0.0) Rmax = ser_b;
    printf("# los velocity dispersion integrated out to %g R_vir\n", Rmax);
  }

  if(NUM) {
    /* Read numerical profiles */
    if( (in = fopen( infile, "r" )) == NULL )
      { printf("Can't open input file <%s>\n", infile); exit(1); }
    i=0;
    while(!feof(in)) {
      fscanf(in,"%le %le %le %le %le %le",ri+i, &mbi, mhi+i, &rbf, &mbf, rhot+i);
      logri[i] = log(ri[i]);
      logmhi[i] = log(mhi[i]);
      logmbi[i] = log(mbi);
      logrbf[i] = log(rbf);
      logmbf[i] = log(mbf);
      i++;
      if(i>nmax) { fprintf(stderr, ": not enough memory for input arrays, increase nmax (currently nmax=%d)\n", nmax); return(1); }
    }
    fclose(in);
    n = n3 = i-1;
    Spline_Install(n, logri, logmbi, Msmbi);
    Spline_Install(n, logrbf, logmbf, Msmbf);   
  } else {
    /* Set up initial radial grid */
    for(i=0; i<n3; i++) {
      ri[i] = pow(10., 4.*(double)(i-n+1)/(double)(n-1));
      logri[i] = log(ri[i]);
      mhi[i] = (1.-fb)*mdm(ri[i],c);
      logmhi[i] = log(mhi[i]);
      rhot[i] = tracer_density(ri[i]);
    }
  }
  Spline_Install(n, logri, logmhi, Msmhi);

  /* Remap initial mass distributions on the average orbital radii */
  for(i=0; i<n; i++) {
    x = y(ri[i],&d);
    if(NUM) {
      logmhi_av = DSpline(log(x), n, logri, logmhi, Msmhi, &d);
      logmbi_av = DSpline(log(x), n, logri, logmbi, Msmbi, &d);
      mhi_av[i] = exp(logmhi_av);
      mbi_av = exp(logmbi_av);
      g[i] = (mhi_av[i]+mbi_av)*ri[i];
    } else {
      mi_av = mdm(x,c);
      mhi_av[i] = (1.-fb)*mi_av;
      g[i] = mi_av*ri[i];
    }
  }

  /* Solve for new radii of spherical dark matter shells */
  for(i=0; i<n; i++) {
    if(MAC == -1) {
      rf[i] = ri[i];
      logrf[i] = logri[i];
    } else {
      rf[i] = rtsafe(funcd, 1.e-3*ri[i], 2.*ri[i], 1.e-10, mhi_av[i], g[i]);
      logrf[i] = log(rf[i]);
    }
  }
  Spline_Install(n, logrf, logmhi, Ms);

  /* Initial and final dark matter densities */
  for(i=0; i<n; i++) {
    logmhf = DSpline(logri[i], n, logrf, logmhi, Ms, &dlogmhf);
    mhf[i] = exp(logmhf);
    rhohf[i] = dlogmhf*mhf[i]/(4.*pi*pow(ri[i],3));
    logrhohf[i] = log(rhohf[i]);
    DSpline(logri[i], n, logri, logmhi, Msmhi, &dlogmhi);
    rhohi[i] = dlogmhi*mhi[i]/(4.*pi*pow(ri[i],3));
    logrhohi[i] = log(rhohi[i]);
  }
  Spline_Install(n, logri, logrhohi, Msrhoi);
  Spline_Install(n, logri, logrhohf, Msrhof);

  /* Logarithmic slopes of density distributions */
  for(i=0; i<n; i++) {
    DSpline(logri[i], n, logri, logrhohi, Msrhoi, dlogrhohi+i);
    DSpline(logri[i], n, logri, logrhohf, Msrhof, dlogrhohf+i);
  }

  /* Total mass profile (baryons + dark matter) */
  for(i=0; i<n; i++)
    mtot[i] = mhf[i] + mb(ri[i],&d);
  for(i=n; i<n3; i++)
    mtot[i] = mdm(ri[i],c);

  /* Line-of-sight velocity dispersion of a tracer population
     integrated from 0 to Rmax, following Mamon & Lokas 2005 */

  if(Rmax > ri[n3-1]) Rmax = ri[n3-1];

  if(TRACE)
    { printf("ERROR: TRACE disabled"); exit(1); } // Added by MRG to avoid linking to los.f
    //    losdispersion_(ri, mtot, rhot, &ANIS, &ra, &n3, &Rmax, sigma_los);

  /* suppress output of unphysical solutions */
  //i0 = n-1;
  //while(-dlogrhohf[i0] >= 0.0 && i0 >= 0) i0--; i0++;
  i0 = 0;

  printf("#\n# r_i       r_f       m_i       m_f       rho_i     rho_f   gam_i gam_f   v_circ  sigma_los\n");
  for(i=i0; i<n; i++)
    printf("%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %5.3f %5.3f %9.3e %9.3e\n",
	   ri[i], rf[i], mhi[i], mhf[i], rhohi[i], rhohf[i], 
	   -dlogrhohi[i], -dlogrhohf[i], sqrt(mtot[i]/ri[i]), sigma_los[i]);
  return 0;
}

int main( int argc, char *argv[] )
{
  int i, i0, tr;
  double x, d, mi_av, mbi, mbi_av, mbf, rbf,
    logmhf, dlogmhf, logmhi_av, logmbi_av, dlogmhi;
  char infile[80];
  FILE *in;

  if(argc < 12) {
    fprintf(stderr, "Adiabatic Contraction of Dark Matter\n\n");
    fprintf(stderr, "Syntax: contra MAC DM BAR TRACE ANIS c ser_dm fb rb ser_b ra [in_file]\n\n");
    fprintf(stderr, "  MAC (model of adiabatic contraction) :\n");
    fprintf(stderr, "    0 = standard, 1 = modified (Gnedin et al. 2004), -1 = no contraction\n");
    fprintf(stderr, "  DM (initial dark matter profile) :\n");
    fprintf(stderr, "    1 = NFW, 2 = spherical Sersic (Navarro et al. 2004), 3 = satellite (Kazantzidis et al. 2004)\n");
    fprintf(stderr, "  BAR (final baryon profile) :\n");
    fprintf(stderr, "    1 = exponential disk, 2 = Hernquist, 3 = Jaffe, 4 = Sersic\n");
    fprintf(stderr, "  TRACE (tracer population) :\n");
    fprintf(stderr, "    1 = exponential, 2 = Hernquist, 3 = Jaffe, 4 = Plummer, 5 = NFW, 6 = Sersic, 7 = spherical Sersic, <0 = power law slope\n");
    fprintf(stderr, "  ANIS (type of anisotropy distribution) :\n");
    fprintf(stderr, "    0 = isotropic, 1 = constant beta, 2 = Mamon-Lokas, 3 = Osipkov-Merritt\n");
    fprintf(stderr, "  c = initial halo concentration parameter, R_vir/r_s\n");
    fprintf(stderr, "  ser_dm = index of a spherical Sersic profile for dark matter\n");
    fprintf(stderr, "  fb = baryon fraction, M_b/M_vir\n");
    fprintf(stderr, "  rb = final baryon scalelength, r_b/R_vir\n");
    fprintf(stderr, "  ser_b = index of a Sersic profile for baryons (if TRACE=6 or 7), otherwise tracer truncation radius\n");
    fprintf(stderr, "  ra = radius of velocity anisotropy for tracer population, r_a/R_vir (or beta, for ANIS=1)\n");
    fprintf(stderr, "  in_file = file with numerical profiles for initial dm and final baryons\n\n");
    fprintf(stderr, "more info at http://www.astro.lsa.umich.edu/~ognedin/contra/\n");
    return 1;
  }
  MAC = atoi(argv[1]);     /* modified or standard model of AC */
  DM = atoi(argv[2]);      /* initial dark matter profile */
  BAR = atoi(argv[3]);     /* final baryon profile */
  TRACE = atof(argv[4]);   /* type of tracer population */
  ANIS = atoi(argv[5]);    /* type of anisotropy distribution */
  c = atof(argv[6]);       /* halo concentration parameter, R_vir/r_s */
  ser_dm = atof(argv[7]);  /* index of a Sersic profile for dark matter */
  fb = atof(argv[8]);      /* baryon fraction */
  rb = atof(argv[9]);      /* baryon scale radius, R_b/R_vir */
  ser_b = atof(argv[10]);  /* index of a Sersic profile for baryons */
  ra = atof(argv[11]);     /* anisotropy radius for tracer population */
  NUM = 0;                 /* numerical profile */
  if(c*fb*rb == 0.) {
    NUM = 1;
    if(argc > 12)
      sprintf(infile, "%s", argv[12]);
    else {
      fprintf(stderr, " filename for a numerical profile must be provided\n");
      return 1;
    }
  }

  if((DM < 1 || DM > 3) && !NUM) { fprintf(stderr, ": DM distribution must be 1, 2 or 3\n"); return 1; }
  if((BAR < 1 || BAR > 4) && !NUM) { fprintf(stderr, ": baryon distribution must be 1, 2, 3 or 4\n"); return 1; }
  if(ANIS < 0 || ANIS > 3) { fprintf(stderr, ": anisotropy parameter must be 0, 1, 2 or 3\n"); return 1; }
  if(TRACE > 7) { fprintf(stderr, ": tracer population must be 1, 2, 3, 4, 5, 6, 7 or negative as a power law slope\n"); return 1; }

  printf("#"); for(i=0; i<argc; i++) printf(" %s", argv[i]); printf("\n");
  if(MAC >= 1) printf("# Modified model of Adiabatic Contraction\n");
  if(MAC == 0) printf("# Standard model of Adiabatic Contraction\n");
  if(MAC == -1) printf("# No adiabatic contraction\n");
  if(NUM) {
    printf("# reading numerical profile from <%s>\n", infile);
  } else {
    if(DM==1 || DM==3) printf("# Dark matter profile: %s model, c=%g\n", dmprofile[DM], c);
    if(DM==2) printf("# Dark matter profile: %s model, c=%g n=%g\n", dmprofile[DM], c, ser_dm);
    if(BAR==1 || BAR==2 || BAR==3) 
      printf("# Baryon distribution: %s model, fb=%g rb=%g\n", barprofile[BAR], fb, rb);
    if(BAR==4) 
      printf("# Baryon distribution: %s model, fb=%g rb=%g n=%g\n", barprofile[BAR], fb, rb, ser_b);
    if(TRACE>0) tr=(int)(TRACE+0.5); else tr=0;
    if(ANIS==0) printf("# Tracer population: %s model with isotropic velocity distribution\n", trprofile[tr]);
    if(ANIS==1) printf("# Tracer population: %s model with constant velocity anisotropy, beta=%g\n", trprofile[tr], ra);
    if(ANIS==2) printf("# Tracer population: %s model with Mamon-Lokas velocity anisotropy, ra=%g\n", trprofile[tr], ra);
    if(ANIS==3) printf("# Tracer population: %s model with Osipkov-Merritt velocity anisotropy, ra=%g\n", trprofile[tr], ra);
    if(TRACE<6 && BAR<4 && ser_b>0.0) Rmax = ser_b;
    printf("# los velocity dispersion integrated out to %g R_vir\n", Rmax);
  }

  if(NUM) {
    /* Read numerical profiles */
    if( (in = fopen( infile, "r" )) == NULL )
      { printf("Can't open input file <%s>\n", infile); exit(1); }
    i=0;
    while(!feof(in)) {
      fscanf(in,"%le %le %le %le %le %le",ri+i, &mbi, mhi+i, &rbf, &mbf, rhot+i);
      logri[i] = log(ri[i]);
      logmhi[i] = log(mhi[i]);
      logmbi[i] = log(mbi);
      logrbf[i] = log(rbf);
      logmbf[i] = log(mbf);
      i++;
      if(i>nmax) { fprintf(stderr, ": not enough memory for input arrays, increase nmax (currently nmax=%d)\n", nmax); return(1); }
    }
    fclose(in);
    n = n3 = i-1;
    Spline_Install(n, logri, logmbi, Msmbi);
    Spline_Install(n, logrbf, logmbf, Msmbf);   
  } else {
    /* Set up initial radial grid */
    for(i=0; i<n3; i++) {
      ri[i] = pow(10., 4.*(double)(i-n+1)/(double)(n-1));
      logri[i] = log(ri[i]);
      mhi[i] = (1.-fb)*mdm(ri[i],c);
      logmhi[i] = log(mhi[i]);
      rhot[i] = tracer_density(ri[i]);
    }
  }
  Spline_Install(n, logri, logmhi, Msmhi);

  /* Remap initial mass distributions on the average orbital radii */
  for(i=0; i<n; i++) {
    x = y(ri[i],&d);
    if(NUM) {
      logmhi_av = DSpline(log(x), n, logri, logmhi, Msmhi, &d);
      logmbi_av = DSpline(log(x), n, logri, logmbi, Msmbi, &d);
      mhi_av[i] = exp(logmhi_av);
      mbi_av = exp(logmbi_av);
      g[i] = (mhi_av[i]+mbi_av)*ri[i];
    } else {
      mi_av = mdm(x,c);
      mhi_av[i] = (1.-fb)*mi_av;
      g[i] = mi_av*ri[i];
    }
  }

  /* Solve for new radii of spherical dark matter shells */
  for(i=0; i<n; i++) {
    if(MAC == -1) {
      rf[i] = ri[i];
      logrf[i] = logri[i];
    } else {
      rf[i] = rtsafe(funcd, 1.e-3*ri[i], 2.*ri[i], 1.e-10, mhi_av[i], g[i]);
      logrf[i] = log(rf[i]);
    }
  }
  Spline_Install(n, logrf, logmhi, Ms);

  /* Initial and final dark matter densities */
  for(i=0; i<n; i++) {
    logmhf = DSpline(logri[i], n, logrf, logmhi, Ms, &dlogmhf);
    mhf[i] = exp(logmhf);
    rhohf[i] = dlogmhf*mhf[i]/(4.*pi*pow(ri[i],3));
    logrhohf[i] = log(rhohf[i]);
    DSpline(logri[i], n, logri, logmhi, Msmhi, &dlogmhi);
    rhohi[i] = dlogmhi*mhi[i]/(4.*pi*pow(ri[i],3));
    logrhohi[i] = log(rhohi[i]);
  }
  Spline_Install(n, logri, logrhohi, Msrhoi);
  Spline_Install(n, logri, logrhohf, Msrhof);

  /* Logarithmic slopes of density distributions */
  for(i=0; i<n; i++) {
    DSpline(logri[i], n, logri, logrhohi, Msrhoi, dlogrhohi+i);
    DSpline(logri[i], n, logri, logrhohf, Msrhof, dlogrhohf+i);
  }

  /* Total mass profile (baryons + dark matter) */
  for(i=0; i<n; i++)
    mtot[i] = mhf[i] + mb(ri[i],&d);
  for(i=n; i<n3; i++)
    mtot[i] = mdm(ri[i],c);

  /* Line-of-sight velocity dispersion of a tracer population
     integrated from 0 to Rmax, following Mamon & Lokas 2005 */

  if(Rmax > ri[n3-1]) Rmax = ri[n3-1];

  if(TRACE)
    { printf("ERROR: TRACE disabled"); exit(1); } // Added by MRG to avoid linking to los.f
    //    losdispersion_(ri, mtot, rhot, &ANIS, &ra, &n3, &Rmax, sigma_los);

  /* suppress output of unphysical solutions */
  //i0 = n-1;
  //while(-dlogrhohf[i0] >= 0.0 && i0 >= 0) i0--; i0++;
  i0 = 0;

  printf("#\n# r_i       r_f       m_i       m_f       rho_i     rho_f   gam_i gam_f   v_circ  sigma_los\n");
  for(i=i0; i<n; i++)
    printf("%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %5.3f %5.3f %9.3e %9.3e\n",
	   ri[i], rf[i], mhi[i], mhf[i], rhohi[i], rhohf[i], 
	   -dlogrhohi[i], -dlogrhohf[i], sqrt(mtot[i]/ri[i]), sigma_los[i]);
  return 0;
}
