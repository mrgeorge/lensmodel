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

    This version of Contra fits the best values of (A,w).

    Written by Oleg Gnedin

    Last modified: September 3, 2010
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sqr(x)  ((x)*(x))

#define nmax 1001
#define pi   3.14159265358979323846
int DM, BAR, MAC, NUM, ANIS, n=81, n3=101;   //n=81 n3=101
double c, fb, rb, ra, ser_dm, ser_b, TRACE, Rmax=3., A=0.85, w=0.8, Ndof, Nb, sigma2,
  r0=0.03, mrms;
double ri[nmax], rf[nmax], mhi[nmax], mhf[nmax], mhf0[nmax],
  rhohi[nmax], rhohf[nmax], mhi_av[nmax], g[nmax],
  rhot[nmax], mtot[nmax], sigma_los[nmax],
  logri[nmax], logrf[nmax], logrbf[nmax], logmbi[nmax], logmbf[nmax],
  logmhi[nmax], logrhohi[nmax], logrhohf[nmax], dlogrhohi[nmax], dlogrhohf[nmax],
  Ms[nmax], Msrhoi[nmax], Msrhof[nmax], Msmhi[nmax], Msmbi[nmax], Msmbf[nmax];

extern double rtsafe();
extern void Spline_Install();
extern double DSpline();
extern void losdispersion_();
extern double gammp_();

char *dmprofile[] = { "", "Navarro, Frenk & White (1997)", 
  "Spherical Sersic", "Kazantzidis et al. (2004)" };
char *barprofile[] = { "", "Exponential", "Hernquist (1990)", 
  "Jaffe (1983)", "Sersic (1968)" };
char *trprofile[] = { "power law", "Exponential", "Hernquist (1990)", 
  "Jaffe (1983)", "Plummer (1911)", "NFW (1997)", "Sersic (1968)", 
  "Spherical Sersic" };


double mdm( double x, double c )    // initial dark matter mass distribution
{
  double f=0.0, a, arg, argc, logmh, dlogmh;

  if(DM==1) {                       /* Navarro, Frenk & White (1997) model */
    f = (1.-fb)*(log(1.+c*x) - c*x/(1.+c*x))/(log(1.+c) - c/(1.+c));
  }
  if(DM==2) {                       /* Navarro et al. (2004): */
    a = 3.*ser_dm;                  /* spherical Sersic model */
    arg = 2.*ser_dm*pow(x*c, 1./ser_dm);
    argc = 2.*ser_dm*pow(c, 1./ser_dm);
    f = (1.-fb)*gammp_(&a, &arg)/gammp_(&a, &argc);
  }
  if(DM==3) {                       /* Kazantzidis et al. (2004) profile */
    a = 2.0;
    arg = x*c;
    argc = c;
    f = (1.-fb)*gammp_(&a, &arg)/gammp_(&a, &argc);
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
    p = 1.0 - 0.6097/ser_b + 0.05563/sqr(ser_b);
    a = (3.-p)*ser_b;
    arg = pow(x/rb, 1./ser_b);
    argc = pow(1./rb, 1./ser_b);
    f = fb*gammp_(&a, &arg)/gammp_(&a, &argc);
    *d = fb/gammp_(&a, &argc)*exp(-arg)*pow(arg,a)/(ser_b*x);
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


double CalculateChi2( double *rob )
{
  int i;
  double x, d, mbi_av, logmhf, dlogmhf, logmhi_av, logmbi_av, dlogmhi,
    Mmodel, Mdata, np, f, frob, sigmai2;

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
      mhi_av[i] = mdm(x,c);
      g[i] = mhi_av[i]/(1.-fb)*ri[i];
    }
  }

  /* Solve for new radii of spherical dark matter shells */
  for(i=0; i<n; i++) {
    if(MAC == -1) {
      rf[i] = ri[i];
      logrf[i] = logri[i];
    } else {
      rf[i] = rtsafe(funcd, 1.e-19*ri[i], 2.*ri[i], 1.e-10, mhi_av[i], g[i]);
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

  /* Calculate global error estimate */
  DM = 0.;
  if(*rob < 0) {
    sigma2=0.; i=0;
    while( ri[i]<=0.1 && i<n ) {
      if(i) {
	Mmodel = mhf[i] - mhf[i-1]*DM;  if(Mmodel<1.e-30) Mmodel=mhf[i];
	Mdata = mhf0[i] - mhf0[i-1]*DM; if(Mdata<1.e-30) Mdata=mhf0[i];
      } else {
	Mmodel = mhf[0];
	Mdata = mhf0[0];
      }
      sigma2 += sqr(Mmodel/Mdata - 1.);
      //sigma2 += sqr(Mmodel-Mdata)/(Mmodel*Mdata);
      i++;
    }
    Nb = (double)i;
    if(MAC==2) { Ndof = Nb-1.; } else { Ndof = Nb-2.; }
    sigma2 = sigma2/Ndof/100.;
  }
  
  /* Calculate chi^2 and a robust estimator at r < 0.1 Rvir */
  f=frob=mrms=0.; i=0;
  while( ri[i]<=0.1 && i<n ) {
    if(i) {
      Mmodel = mhf[i] - mhf[i-1]*DM;  if(Mmodel<1.e-30) Mmodel=mhf[i];
      Mdata = mhf0[i] - mhf0[i-1]*DM; if(Mdata<1.e-30) Mdata=mhf0[i];
    } else {
      Mmodel = mhf[0];
      Mdata = mhf0[0];
    }
    np = rhot[i];                   // number of particles in a bin
    sigmai2 = 1./np + sigma2;       // shot noise + global error
    //printf("r=%5.3f sigma=%6.4f sigmai=%6.4f dMd=%9.3e dMm=%9.3e z=%5.3f\n", ri[i], sigma2, sigmai2, dMdata, dMmodel, sqr(dMmodel-dMdata)/(dMmodel*dMdata));

    mrms += fabs(mhf[i]/mhf0[i] - 1.);
    d = sqr(Mmodel/Mdata - 1.)/sigmai2;
    //d = sqr(dMmodel-dMdata)/(dMmodel*dMdata)/sigmai2;
    f += d;
    // robust estimator
    frob += sqrt(d);
    i++;
  }
  *rob = frob;
  mrms /= Nb;

  return f;
}


int main( int argc, char *argv[] )
{
  int i, j, i0, tr;
  double d, mbi, mbf, rbf, chi2, chi2min=1.e30, rob, robmin=1.e30, x, f01, m01, logmhf, bulge,
    Amin=1., wmin=1., Arob=1., wrob=1.;
  char infile[80];
  FILE *in, *out;

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
      { printf("Can't open input file <%s>\n", infile); return 1; }
    i=0;
    while(!feof(in)) {
      fscanf(in,"%le %le %le %le %le %le %le", 
	     ri+i, &mbi, mhi+i, &rbf, &mbf, mhf0+i, rhot+i);
      logri[i] = log(ri[i]);
      logmhi[i] = log(mhi[i]);
      logmbi[i] = log(mbi);
      logrbf[i] = log(rbf);
      logmbf[i] = log(mbf);
      i++;
      if(i>nmax) { fprintf(stderr, ": not enough memory for input arrays, increase nmax (currently nmax=%d)\n", nmax); return 1; }
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
      mhi[i] = mdm(ri[i],c);
      logmhi[i] = log(mhi[i]);
      rhot[i] = tracer_density(ri[i]);
    }
  }
  Spline_Install(n, logri, logmhi, Msmhi);

  // calculate standard model and estimate global model error
  A=1.; w=1.; rob=-1.;
  CalculateChi2(&rob);

  if(MAC) {
    // Calculate chi^2 for a grid of (A,w): 0 < A < 3, 0 < w < 2
    if( (out = fopen( "chi2.dat", "w" )) == NULL )
      { printf("Can't open ouput file\n"); return 1; }

    int nA=291, nw=190, imin=1, jmin=1;
    double A_[nA+1], w_[nw+1], chi2_[nA+1][nw+1], d, dA1=0, dA2=0, dw1, dw2;
    
    for(i=0; i<=nA; i++) {
      if(MAC==2) i = 150;
      for(j=0; j<=nw; j++) {
	A_[i] = A = 0.1 + (double)i*0.01;
	w_[j] = w = 0.1 + (double)j*0.01;
	chi2_[i][j] = chi2 = CalculateChi2(&rob);
	if(chi2 < chi2min) { chi2min=chi2; Amin=A; wmin=w; imin=i; jmin=j; }
	if(rob < robmin) { robmin=rob; Arob=A; wrob=w; }
	fprintf(out, "%4.2f %5.2f %9.4f %9.4f\n", A, w, chi2, rob);
      }
      if(MAC==2) break;
    }
    fclose(out);

    d = 2.71;    // 90% confidence limits on the parameters
    if(MAC==1) {
    i=imin+1; while(chi2_[i][jmin]-chi2min < d && i<nA) i++; dA2 = A_[i]-A_[imin];
    i=imin-1; while(chi2_[i][jmin]-chi2min < d && i>0) i--;  dA1 = A_[imin]-A_[i];
    }
    j=jmin+1; while(chi2_[imin][j]-chi2min < d && j<nw) j++; dw2 = w_[j]-w_[jmin];
    j=jmin-1; while(chi2_[imin][j]-chi2min < d && j>0) j--;  dw1 = w_[jmin]-w_[j];

    // baryon domination factor at 0.01 Rvir
    x = 0.01;
    f01 = mb(x,&d)/mdm(x,c);
    // halo mass enhancement factor at 0.01 Rvir
    logmhf = DSpline(log(x), n, logrf, logmhi, Ms, &d);
    m01 = exp(logmhf)/mb(x,&d);
    // baryon bulge fraction at 0.004 Rvir
    bulge = mb(x,&d)/mb(1.,&d);

    /* output best-fitting model */
    A=Amin; w=wmin; //A=Arob; w=wrob;
    CalculateChi2(&rob);
    fprintf(stderr, "%7.2f  %4.2f %4.2f %4.2f  %4.2f %4.2f %4.2f %5.1f %4.2f %5.3f %6.3f %d\n",
	    chi2min, Amin, dA2, dA1, wmin, dw2, dw1, f01, m01, bulge, mrms, (int)Ndof);
    printf("# chi2=%g A=%g+%g-%g w=%g+%g-%g fb01=%3.1f m01=%4.2f bulge=%5.3f mrms=%5.3f Ndof=%g\n",
	   chi2min, Amin, dA2, dA1, wmin, dw2, dw1, f01, m01, bulge, mrms, Ndof);
  }

  // proceed as in main contra
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
    mtot[i] = mdm(ri[i],c)/(1.-fb);

  /* Line-of-sight velocity dispersion of a tracer population
     integrated from 0 to Rmax, following Mamon & Lokas 2005 */

  if(Rmax > ri[n3-1]) Rmax = ri[n3-1];

  if(TRACE)
    losdispersion_(ri, mtot, rhot, &ANIS, &ra, &n3, &Rmax, sigma_los);

  /* suppress output of unphysical solutions */
  //i0 = n-1; while(-dlogrhohf[i0] >= 0.0 && i0 >= 0) i0--; i0++;
  i0 = 0;

  printf("#\n# r_i       r_f       m_i       m_f       rho_i     rho_f   gam_i gam_f   v_circ  sigma_los\n");
  for(i=i0; i<n; i++)
    printf("%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %5.3f %5.3f %9.3e %9.3e\n",
	   ri[i], rf[i], mhi[i], mhf[i], rhohi[i], rhohf[i], 
	   -dlogrhohi[i], -dlogrhohf[i], sqrt(mtot[i]/ri[i]), sigma_los[i]);
  return 0;
}
