/* This tool box contains useful numerical algorithms not 
  defined in ANSI C standard libraries. Some of the algorithms 
  are taken from Numerical Recipes in C.

  First Version: May 23, 2008
  Second Version: July 4, 2008
  Updated
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <sys/time.h>
#include <errno.h>
#include "tools.h"

/*========================================================================
  Get time
========================================================================*/
time_t gettimeofday_sec()
{
	time_t now = 0;
	now = time(NULL);
	return now;
   /* struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;*/
}


/*========================================================================
  nrerror
========================================================================*/
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

/*========================================================================
  Random Number Generators
========================================================================*/

// Uniform

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
	int j;
	long k;
	 static long iy=0;
	 static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


// Standard Normal

double gasdev(long *idum)
{
  /*  double ran1(long *idum); */
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if (*idum < 0) iset=0;
  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}	
	


// Standard Normal given uniform draw input (not internally generated)

/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 */

/* Coefficients in rational approximations. */
static const double a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00
};

static const double d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

#define LOW 0.02425
#define HIGH 0.97575

double ltqnorm(double p)
{
	double q, r;

	errno = 0;

	if (p < 0 || p > 1){
		errno = EDOM;
		return 0.0;
	}else if (p == 0){
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	}else if (p == 1){
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	}else if (p < LOW){
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}else if (p > HIGH){
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}else{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}


double gammln(double xx)
{
	double x,y,tmp,ser;
	const static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
  /*	double gammln(double xx);
	void nrerror(char error_text[]);  */
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS


#define ITMAX 200
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
  /*	double gammln(double xx);
	void nrerror(char error_text[]);  */
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) {
	  printf("ITMAX = %d.\n", ITMAX);
	  printf("a = %f.\n", a);
	  printf("x = %f.\n", x);
	  nrerror("a too large, ITMAX too small in dgcf");
	}
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN



double gammp(double a, double x)
{
  /*	void gcf(double *gammcf, double a, double x, double *gln);
 	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);  */
  double gamser,gammcf,gln;

  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf;
  }
}

double gammq(double a, double x)
{
  /*	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);  */
  double gamser,gammcf,gln;
  
  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}

double erffc(double x)
{
  return x < 0.0 ? 1.0+gammp(0.5,x*x) : gammq(0.5,x*x);
}

double gascdf(double x)
{
  return erffc(-x/sqrt(2))/2;
}

//***************************
// dlanor, eval_pol, alnrel are used to compute the log (normal cdf)
//***************************
/*
double dlanor ( double *x )

//****************************************************************************80
//
//  Purpose:
// 
//    DLANOR evaluates the logarithm of the asymptotic Normal CDF.
//
//  Discussion:
//
//    This routine computes the logarithm of the cumulative normal distribution
//    from abs ( x ) to infinity for  5 <= abs ( X ).
//
//    The relative error at X = 5 is about 0.5D-5.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.2.12.
//
//  Parameters:
//
//    Input, double *X, the value at which the Normal CDF is to be
//    evaluated.  It is assumed that 5 <= abs ( X ).
//
//    Output, double DLANOR, the logarithm of the asymptotic
//    Normal CDF.
//
{
# define dlsqpi 0.91893853320467274177e0

  static const double coef[12] = {
    -1.0e0,3.0e0,-15.0e0,105.0e0,-945.0e0,10395.0e0,-135135.0e0,2027025.0e0,
    -34459425.0e0,654729075.0e0,-13749310575.e0,316234143225.0e0
  };
  static int K1 = 12;
  static double dlanor,approx,correc,xx,xx2,T2;

  xx = fabs(*x);
  if ( xx < 5.0e0 ) nrerror(" Argument too small in DLANOR");
  approx = -dlsqpi-0.5e0*xx*xx-log(xx);
  xx2 = xx*xx;
  T2 = 1.0e0/xx2;
  correc = eval_pol ( coef, &K1, &T2 ) / xx2;
  correc = alnrel ( &correc );
  dlanor = approx+correc;
  return dlanor;
# undef dlsqpi
}


double eval_pol ( double a[], int *n, double *x )

//****************************************************************************80
// 
//  Purpose:
// 
//    EVAL_POL evaluates a polynomial at X.
//
//  Discussion:
//
//    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
//
//  Modified:
//
//    15 December 1999
//
//  Parameters:
//
//    Input, double precision A(0:N), coefficients of the polynomial.
//
//    Input, int *N, length of A.
//
//    Input, double *X, the point at which the polynomial 
//    is to be evaluated.
//
//    Output, double EVAL_POL, the value of the polynomial at X.
//
{
  static double devlpl,term;
  static int i;

  term = a[*n-1];
  for ( i = *n-1-1; i >= 0; i-- )
  {
    term = a[i]+term**x;
  }

  devlpl = term;
  return devlpl;
}

double alnrel ( double *a )

//****************************************************************************80
//
//  Purpose:
// 
//    ALNREL evaluates the function ln ( 1 + A ).
//
//  Modified:
//
//    17 November 2006
//
//  Reference:
//
//    Armido DiDinato, Alfred Morris,
//    Algorithm 708: 
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double *A, the argument.
//
//    Output, double ALNREL, the value of ln ( 1 + A ).
//
{
  double alnrel;
  static double p1 = -0.129418923021993e+01;
  static double p2 =  0.405303492862024e+00;
  static double p3 = -0.178874546012214e-01;
  static double q1 = -0.162752256355323e+01;
  static double q2 =  0.747811014037616e+00;
  static double q3 = -0.845104217945565e-01;
  double t;
  double t2;
  double w;
  double x;

  if ( fabs ( *a ) <= 0.375e0 )
  {
    t = *a / ( *a + 2.0e0 );
    t2 = t * t;
    w = (((p3*t2+p2)*t2+p1)*t2+1.0e0)
      / (((q3*t2+q2)*t2+q1)*t2+1.0e0);
    alnrel = 2.0e0 * t * w;
  }
  else
  {
    x = 1.0e0 + *a;
    alnrel = log ( x );
  }
  return alnrel;
}
*/


//========================================================================
//	Linear Algebra
//========================================================================

void matmul(double **A, double **B, double **AB, long m, long n, long l)
{
  int i,j,k;
  for (i=1; i<=m; i++){
    for (j=1; j<=l; j++){
      AB[i][j] = 0.0;
      for (k=1; k<=n; k++){
        AB[i][j] += A[i][k]*B[k][j];
      }
    }
  }
}

void transpose(double **A, double **A_t, int m, int n)
{
  int i,j;
  double temp_copy;
  for (i=1; i<=m; i++){
    for (j=1; j<=n; j++){
      temp_copy = A[i][j];
      A_t[j][i] = temp_copy;
    }
  }
}


#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv = (double *)dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}
#undef TINY


void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}

}



void mxinv(double **source, double **target, int dimension)
/* Calculates the inverse of the matrix source, and puts it in target */
{
  double d, *col;
  int i, j, *indx;

  col = (double *)dvector(1, dimension);
  indx = (int *)ivector(1, dimension);

  ludcmp(source, dimension, indx, &d);
  for(j=1; j<=dimension; j++) {
    for(i=1; i<=dimension; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(source, dimension, indx, col);
    for(i=1; i<=dimension; i++) target[i][j]=col[i];
  }

  free_dvector(col, 1, dimension); free_ivector(indx, 1, dimension);

}


/* This routine takes a matrix and dimension, and produce the inverse matrix
and its determinant. The initial matrix must be non-singluar. Otherwise, we
will have error in the subroutine "ludcmp". */

void mx_inv_det(double **matrix, double **inv_matrix, double *det,
		int dimension)
{
  double *col, **source;
  int i, j, *indx;

  /* printf("\nInside mx_inv_det(1).\n"); */
  /* printf("dimension = %d.\n", dimension); */

  col = dvector(1, dimension);
  indx = ivector(1, dimension);

  /* printf("\nInside mx_inv_det(2).\n"); */

  source = dmatrix(1, dimension, 1, dimension);

  /* printf("\nInside mx_inv_det(3).\n"); */

  for (i=1; i<=dimension; i++)
    for (j=1; j<=dimension; j++)
      source[i][j] = matrix[i][j];

  /* source = matrix; */

  /* for (i=1; i<=dimension; i++)
     for (j=1; j<=dimension; j++)
     printf("source[%d][%d] = %f.\n", i, j, source[i][j]); */

  ludcmp(source, dimension, indx, det);
  for(j=1; j<=dimension; j++) {
    for(i=1; i<=dimension; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(source, dimension, indx, col);
    for(i=1; i<=dimension; i++) inv_matrix[i][j]=col[i];
  }

  for(j=1; j<=dimension; j++)
    (*det) *= source[j][j];

  free_dvector(col, 1, dimension);
  free_ivector(indx, 1, dimension);
  free_dmatrix(source, 1, dimension, 1, dimension);

  /* printf("Leaving mx_inv_det\n"); */

}






/*========================================================================
  Memory Allocation for Arrays
========================================================================*/

#define NR_END 1
#define FREE_ARG char*

// vector ----------------------------------------------------------------

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

// matrix ----------------------------------------------------------------

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

// array3d ---------------------------------------------------------------

double ***darray3d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h)
/* allocate a double 3-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h] */
{
  long i, j, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1;
  double ***m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(double ***) malloc(((n1+NR_END)*sizeof(double**)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(double **)malloc(((n1*n2+NR_END)*sizeof(double *)));
  if (!m[n1l]) nrerror("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(double *)malloc(((n1*n2*n3+NR_END)*sizeof(double)));
  if (!m[n1l][n2l]) nrerror("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }

  /* return pointer to array of pointers to array of pointers to a row of
     memories */
  return m;

}

int ***iarray3d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h)
/* allocate a double 3-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h] */
{
  long i, j, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1;
  int ***m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(int ***) malloc(((n1+NR_END)*sizeof(int**)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(int **)malloc(((n1*n2+NR_END)*sizeof(int *)));
  if (!m[n1l]) nrerror("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(int *)malloc(((n1*n2*n3+NR_END)*sizeof(int)));
  if (!m[n1l][n2l]) nrerror("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }

  /* return pointer to array of pointers to rows */
  return m;

}

// array4d ---------------------------------------------------------------

double ****darray4d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		long n4l, long n4h)
/* allocate a double 3-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h][n4l..n4h][n5l..n5h] */
{
  long i, j, k, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1, n4=n4h-n4l+1;
  double ****m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(double ****) malloc(((n1+NR_END)*sizeof(double***)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(double ***)malloc(((n1*n2+NR_END)*sizeof(double **)));
  if (!m[n1l]) nrerror("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(double **)malloc(((n1*n2*n3+NR_END)*sizeof(double *)));
  if (!m[n1l][n2l]) nrerror("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }


  m[n1l][n2l][n3l]=(double *)malloc(((n1*n2*n3*n4+NR_END)*sizeof(double)));
  if (!m[n1l][n2l][n3l]) nrerror("allocation failure 3 in matrix()");
  m[n1l][n2l][n3l] += NR_END;
  m[n1l][n2l][n3l] -= n4l;

  for(k=n3l+1; k<=n3h; k++) m[n1l][n2l][k] = m[n1l][n2l][k-1] + n4;

  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++) {

      if (k==n3l)
	m[n1l][j][k] = m[n1l][j-1][n3h] + n4;
      else
	m[n1l][j][k] = m[n1l][j][k-1] + n4;

    }


  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++) {

	if ((j==n2l) & (k==n3l))
	  m[i][n2l][n3l] = m[i-1][n2h][n3h] + n4;
	else {
	  if (k==n3l)
	    m[i][j][k] = m[i][j-1][n3h] + n4;
	  else
	    m[i][j][k] = m[i][j][k-1] + n4;

	}

      }

  /* return pointer to array of pointers to rows */
  return m;

}


int ****iarray4d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		long n4l, long n4h)
/* allocate a double 3-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h][n4l..n4h][n5l..n5h] */
{
  long i, j, k, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1, n4=n4h-n4l+1;
  int ****m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(int ****) malloc(((n1+NR_END)*sizeof(int***)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(int ***)malloc(((n1*n2+NR_END)*sizeof(int **)));
  if (!m[n1l]) nrerror("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(int **)malloc(((n1*n2*n3+NR_END)*sizeof(int *)));
  if (!m[n1l][n2l]) nrerror("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }


  m[n1l][n2l][n3l]=(int *)malloc(((n1*n2*n3*n4+NR_END)*sizeof(int)));
  if (!m[n1l][n2l][n3l]) nrerror("allocation failure 3 in matrix()");
  m[n1l][n2l][n3l] += NR_END;
  m[n1l][n2l][n3l] -= n4l;

  for(k=n3l+1; k<=n3h; k++) m[n1l][n2l][k] = m[n1l][n2l][k-1] + n4;

  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++) {

      if (k==n3l)
	m[n1l][j][k] = m[n1l][j-1][n3h] + n4;
      else
	m[n1l][j][k] = m[n1l][j][k-1] + n4;

    }


  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++) {

	if ((j==n2l) & (k==n3l))
	  m[i][n2l][n3l] = m[i-1][n2h][n3h] + n4;
	else {
	  if (k==n3l)
	    m[i][j][k] = m[i][j-1][n3h] + n4;
	  else
	    m[i][j][k] = m[i][j][k-1] + n4;

	}

      }

  /* return pointer to array of pointers to rows */
  return m;

}


// array5d ---------------------------------------------------------------

double *****darray5d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		     long n4l, long n4h, long n5l, long n5h)
/* allocate a double 5-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h][n4l..n4h][n5l..n5h] */
{
  long i, j, k, l, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1, n4=n4h-n4l+1,
    n5=n5h-n5l+1;
  double *****m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(double *****) malloc(((n1+NR_END)*sizeof(double ****)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(double ****)malloc(((n1*n2+NR_END)*sizeof(double ***)));
  if (!m[n1l]) nrerror("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(double ***)malloc(((n1*n2*n3+NR_END)*sizeof(double **)));
  if (!m[n1l][n2l]) nrerror("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }


  m[n1l][n2l][n3l]=(double **)malloc(((n1*n2*n3*n4+NR_END)*sizeof(double *)));
  if (!m[n1l][n2l][n3l]) nrerror("allocation failure 4 in matrix()");
  m[n1l][n2l][n3l] += NR_END;
  m[n1l][n2l][n3l] -= n4l;

  for(k=n3l+1; k<=n3h; k++) m[n1l][n2l][k] = m[n1l][n2l][k-1] + n4;

  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++) {

      if (k==n3l)
	m[n1l][j][k] = m[n1l][j-1][n3h] + n4;
      else
	m[n1l][j][k] = m[n1l][j][k-1] + n4;

    }


  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++) {

	if ((j==n2l) & (k==n3l))
	  m[i][j][k] = m[i-1][n2h][n3h] + n4;
	else {
	  if (k==n3l)
	    m[i][j][k] = m[i][j-1][n3h] + n4;
	  else
	    m[i][j][k] = m[i][j][k-1] + n4;

	}

      }


  m[n1l][n2l][n3l][n4l]=(double *)malloc(((n1*n2*n3*n4*n5+NR_END)*sizeof(double)));
  if (!m[n1l][n2l][n3l][n4l]) nrerror("allocation failure 5 in matrix()");
  m[n1l][n2l][n3l][n4l] += NR_END;
  m[n1l][n2l][n3l][n4l] -= n5l;


  for(l=n4l+1; l<=n4h; l++) m[n1l][n2l][n3l][l] = m[n1l][n2l][n3l][l-1] + n5;

  for(k=n3l+1; k<=n3h; k++)
    for(l=n4l; l<=n4h; l++) {

      if (l==n4l)
	m[n1l][n2l][k][l] = m[n1l][n2l][k-1][n4h] + n5;
      else
	m[n1l][n2l][k][l] = m[n1l][n2l][k][l-1] + n5;

    }


  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++)
      for(l=n4l; l<=n4h; l++) {

	if ((k==n3l) & (l==n4l))
	  m[n1l][j][k][l] = m[n1l][j-1][n3h][n4h] + n5;
	else {
	  if (l==n4l)
	    m[n1l][j][k][l] = m[n1l][j][k-1][n4h] + n5;
	  else
	    m[n1l][j][k][l] = m[n1l][j][k][l-1] + n5;

	}

      }



  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++)
	for(l=n4l; l<=n4h; l++) {

	  if (((j==n2l) & (k==n3l)) & (l==n4l))
	    m[i][j][k][l] = m[i-1][n2h][n3h][n4h] + n5;
	  else if ((k==n3l) & (l==n4l))
	    m[i][j][k][l] = m[i][j-1][n3h][n4h] + n5;
	  else if (l==n4l)
	    m[i][j][k][l] = m[i][j][k-1][n4h] + n5;
	  else
	    m[i][j][k][l] = m[i][j][k][l-1] + n5;

	}

  /* return pointer to array of pointers to rows */
  return m;

}


int *****iarray5d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		  long n4l, long n4h, long n5l, long n5h)
/* allocate a double 5-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h][n4l..n4h][n5l..n5h] */
{
  long i, j, k, l, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1, n4=n4h-n4l+1,
    n5=n5h-n5l+1;
  int *****m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(int *****) malloc(((n1+NR_END)*sizeof(int****)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(int ****)malloc(((n1*n2+NR_END)*sizeof(int ***)));
  if (!m[n1l]) nrerror("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(int ***)malloc(((n1*n2*n3+NR_END)*sizeof(int **)));
  if (!m[n1l][n2l]) nrerror("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }


  m[n1l][n2l][n3l]=(int **)malloc(((n1*n2*n3*n4+NR_END)*sizeof(int *)));
  if (!m[n1l][n2l][n3l]) nrerror("allocation failure 4 in matrix()");
  m[n1l][n2l][n3l] += NR_END;
  m[n1l][n2l][n3l] -= n4l;

  for(k=n3l+1; k<=n3h; k++) m[n1l][n2l][k] = m[n1l][n2l][k-1] + n4;

  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++) {

      if (k==n3l)
	m[n1l][j][k] = m[n1l][j-1][n3h] + n4;
      else
	m[n1l][j][k] = m[n1l][j][k-1] + n4;

    }


  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++) {

	if ((j==n2l) & (k==n3l))
	  m[i][j][k] = m[i-1][n2h][n3h] + n4;
	else {
	  if (k==n3l)
	    m[i][j][k] = m[i][j-1][n3h] + n4;
	  else
	    m[i][j][k] = m[i][j][k-1] + n4;

	}

      }


  m[n1l][n2l][n3l][n4l]=(int *)malloc(((n1*n2*n3*n4*n5+NR_END)*sizeof(int)));
  if (!m[n1l][n2l][n3l][n4l]) nrerror("allocation failure 5 in matrix()");
  m[n1l][n2l][n3l][n4l] += NR_END;
  m[n1l][n2l][n3l][n4l] -= n5l;


  for(l=n4l+1; l<=n4h; l++) m[n1l][n2l][n3l][l] = m[n1l][n2l][n3l][l-1] + n5;

  for(k=n3l+1; k<=n3h; k++)
    for(l=n4l; l<=n4h; l++) {

      if (l==n4l)
	m[n1l][n2l][k][l] = m[n1l][n2l][k-1][n4h] + n5;
      else
	m[n1l][n2l][k][l] = m[n1l][n2l][k][l-1] + n5;

    }


  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++)
      for(l=n4l; l<=n4h; l++) {

	if ((k==n3l) & (l==n4l))
	  m[n1l][j][k][l] = m[n1l][j-1][n3h][n4h] + n5;
	else {
	  if (l==n4l)
	    m[n1l][j][k][l] = m[n1l][j][k-1][n4h] + n5;
	  else
	    m[n1l][j][k][l] = m[n1l][j][k][l-1] + n5;

	}

      }



  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++)
	for(l=n4l; l<=n4h; l++) {

	  if (((j==n2l) & (k==n3l)) & (l==n4l))
	    m[i][j][k][l] = m[i-1][n2h][n3h][n4h] + n5;
	  else if ((k==n3l) & (l==n4l))
	    m[i][j][k][l] = m[i][j-1][n3h][n4h] + n5;
	  else if (l==n4l)
	    m[i][j][k][l] = m[i][j][k-1][n4h] + n5;
	  else
	    m[i][j][k][l] = m[i][j][k][l-1] + n5;

	}


  /* return pointer to array of pointers to rows */
  return m;

}



// free vector ----------------------------------------------------

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

// free matrix ----------------------------------------------------

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

// free array3d ----------------------------------------------------------

void free_darray3d(double ***m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h)
/* free a double 3-dimensional array allocated by darray3d() */
{
  free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
  free((FREE_ARG) (m[n1l]+n2l-NR_END));
  free((FREE_ARG) (m+n1l-NR_END));
}

void free_iarray3d(int ***m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h)
/* free a int 3-dimensional array allocated by darray3d() */
{
        free((m[n1l][n2l]+n3l-NR_END));
	free((m[n1l]+n2l-NR_END));
	free((m+n1l-NR_END));
}


// free array4d ----------------------------------------------------------

void free_darray4d(double ****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h)
/* free a double 4-dimensional array allocated by darray3d() */
{
  free((FREE_ARG) (m[n1l][n2l][n3l]+n4l-NR_END));
  free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
  free((FREE_ARG) (m[n1l]+n2l-NR_END));
  free((FREE_ARG) (m+n1l-NR_END));
}

void free_iarray4d(int ****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h)
/* free a double 4-dimensional array allocated by darray3d() */
{
        free((FREE_ARG) (m[n1l][n2l][n3l]+n4l-NR_END));
	free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
	free((FREE_ARG) (m[n1l]+n2l-NR_END));
	free((FREE_ARG) (m+n1l-NR_END));
}


// free array5d ----------------------------------------------------------

void free_darray5d(double *****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h, long n5l, long n5h)
/* free a double 5-dimensional array allocated by darray3d() */
{
  free((FREE_ARG) (m[n1l][n2l][n3l][n4l]+n5l-NR_END));
  free((FREE_ARG) (m[n1l][n2l][n3l]+n4l-NR_END));
  free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
  free((FREE_ARG) (m[n1l]+n2l-NR_END));
  free((FREE_ARG) (m+n1l-NR_END));
}

void free_iarray5d(int *****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h, long n5l, long n5h)
/* free a double 5-dimensional array allocated by darray3d() */
{

        free((FREE_ARG) (m[n1l][n2l][n3l][n4l]+n5l-NR_END));
        free((FREE_ARG) (m[n1l][n2l][n3l]+n4l-NR_END));
	free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
	free((FREE_ARG) (m[n1l]+n2l-NR_END));
	free((FREE_ARG) (m+n1l-NR_END));
}

