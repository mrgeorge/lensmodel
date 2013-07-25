%module contra
%{
// Headers from contra.c
double mdm( double x, double c );
double y( double x, double *d );
double mb( double x, double *d );
void funcd( double r, double *f, double *df, double mhi, double g );
double tracer_density( double x );
int pymain(int MAC, int DM, int BAR, double TRACE, int ANIS, double c, double ser_dm, double fb, double rb, double ser_b, double ra, double A, double w);
int main( int argc, char *argv[] );


// Headers from spline.c
#define Nsi 999       /* required by Spline_Install */

void Spline_Install(int n,double x[],double y[],double Ms[] ) ;
double DSpline(double z,int n,double x[],double y[],double Ms[],double *Derivative );

#define MAXIT 100
double rtsafe(void (*funcd)(double, double *, double *, double, double), double x1, double x2,
	      double xacc, double par1, double par2);

%}

// Headers from contra.c
double mdm( double x, double c );
double y( double x, double *d );
double mb( double x, double *d );
void funcd( double r, double *f, double *df, double mhi, double g );
double tracer_density( double x );
int pymain(int MAC, int DM, int BAR, double TRACE, int ANIS, double c, double ser_dm, double fb, double rb, double ser_b, double ra, double A, double w);
int main( int argc, char *argv[] );


// Headers from spline.c
#define Nsi 999       /* required by Spline_Install */

void Spline_Install(int n,double x[],double y[],double Ms[] ) ;
double DSpline(double z,int n,double x[],double y[],double Ms[],double *Derivative );

#define MAXIT 100
double rtsafe(void (*funcd)(double, double *, double *, double, double), double x1, double x2,
	      double xacc, double par1, double par2);

