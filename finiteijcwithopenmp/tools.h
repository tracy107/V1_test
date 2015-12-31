/* This is the header file for the tool box "tools.c"

  First Version: May 23, 2008
  Second Version: July 4, 2008
  Updated
*/



time_t gettimeofday_sec();
void nrerror(char error_text[]);

/*========================================================================
**  Utilities
========================================================================*/

//static float sqrarg;
//#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
//static double dsqrarg;
//#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
//static double dmaxarg1,dmaxarg2;
//#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?(dmaxarg1) : (dmaxarg2))
//static double dminarg1,dminarg2;
//#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?(dminarg1) : (dminarg2))
//static float maxarg1,maxarg2;
//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?(maxarg1) : (maxarg2))
//static float minarg1,minarg2;
//#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?(minarg1) : (minarg2))
//static long lmaxarg1,lmaxarg2;
//#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?(lmaxarg1) : (lmaxarg2))
//static long lminarg1,lminarg2;
//#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?(lminarg1) : (lminarg2))
//static int imaxarg1,imaxarg2;
//#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?(imaxarg1) : (imaxarg2))
//static int iminarg1,iminarg2;
//#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?(iminarg1) : (iminarg2))
//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/*========================================================================
**  Random Number Generators and probablity functions
========================================================================*/
// Uniform random number generator
double ran1(long *idum);
// Standard normal random number generator
double gasdev(long *idum);
// Standard Normal given uniform draw input (not internally generated)
double ltqnorm(double p);
// Return the value ln Gamma(xx) for xx>0
double gammln(double xx);
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double gammp(double a, double x);
double gammq(double a, double x);
double erffc(double x);
double gascdf(double x);
double dlanor ( double *x );
double eval_pol ( double a[], int *n, double *x );
double alnrel ( double *a );
//========================================================================
// Linear Algebra
//========================================================================
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void mxinv(double **source, double **target, int dimension);
void mx_inv_det(double **matrix, double **inv_matrix, double *det, int dimension);
void matmul(double **A, double **B, double **AB, long m, long n, long l);
void transpose(double **A, double **A_t, int m, int n);

/*========================================================================
**  Memory Allocation for Arrays
========================================================================*/
double *dvector(long nl, long nh);
int *ivector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double ***darray3d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h);
int ***iarray3d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h);
double ****darray4d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
                		long n4l, long n4h);
int ****iarray4d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		long n4l, long n4h);
double *****darray5d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
            		     long n4l, long n4h, long n5l, long n5h);
int *****iarray5d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		  long n4l, long n4h, long n5l, long n5h);
void free_dvector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_darray3d(double ***m, long n1l, long n1h, long n2l, long n2h,
             		   long n3l, long n3h);
void free_iarray3d(int ***m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h);
void free_darray4d(double ****m, long n1l, long n1h, long n2l, long n2h,
            		   long n3l, long n3h, long n4l, long n4h);
void free_iarray4d(int ****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h);
void free_darray5d(double *****m, long n1l, long n1h, long n2l, long n2h,
            		   long n3l, long n3h, long n4l, long n4h, long n5l, long n5h);
void free_iarray5d(int *****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h, long n5l, long n5h);
