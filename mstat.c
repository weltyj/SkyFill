/*
	FILE mstat.c
	PURPOSE - multiple linear regression
	CONTENTS -
*/

#include <stdio.h>
#include <stdlib.h>

#include "mstat.h"
#define MAX_P MSTAT_MAX_P
#define ROW_P MSTAT_ROW_P
#define ROW_P2 (MSTAT_ROW_P*2)
#define TINY (1.e-12)

// in case Windows MSVC, or 
#define _USE_MATH_DEFINES
#include <math.h>

#include "skyfill_tif.h"

// sample standard deviation
double stddev_calc(int n, double sum, double sumsq)
{
    if(n < 2)
	return 0. ;

    double nd=n ;
    double sq = (sumsq - sum*sum/nd)/(nd-1.) ;

    if(sq < 0.0)
	return 0. ;

    return sqrt(sq) ;
}

double cv_calc(int n, double sum, double sumsq)
{
    if(n < 2)
	return 0. ;

    double mean = sum/(double)n ;
    return stddev_calc(n,sum,sumsq) / mean ;
}

float  means_t_test(int n1, float sum1,  float sumsq1, int n2, float sum2, float sumsq2, int verbose)
{
    if(n1 < 2 || n1 < 1)
	return 0. ;

    float sd1 = stddev_calc(n1,(double)sum1,(double)sumsq1) ;
    float m1=sum1/(float)n1 ;

    if(n2 == 1) {
	// sample point t-test
	float SE_Xbar = sd1 / sqrt((float)n1) ;
	return (sum2 - m1) / SE_Xbar ;
    }

    float sd2 = stddev_calc(n2,(double)sum2,(double)sumsq2) ;
    float m2=sum2/(float)n2 ;

    if(verbose)
	fprintf(stderr, "MEANS T TEST: mean %f %f, sd %f %f\n", m1/MAX16f, m2/MAX16f, sd1/MAX16f, sd2/MAX16f) ;

    // t-test with pooled std dev
    return (m1-m2) / sqrtf( sd1*sd1/(float)n1+ sd2*sd2/(float)n2) ;
}

void gauss( double [][ROW_P2], short) ;
void lu(double [][ROW_P], short) ;
void lubksb( double [][ROW_P], short, short *) ;
void ludcmp(double [][ROW_P], short, short *, double *) ;
void matrix_solve(double [][ROW_P], short) ;

double y, f ;
int row, col, irow, i, j, k ;



#define abs(x) ((x) > 0 ? (x) : -(x))
#define tiny 1e-200
#define LU


void dv_exit(char *msg)
{
}

void cholesky(double mat[][ROW_P], short N)
{
    double lower[ROW_P+1][ROW_P+1] ;
    double c[ROW_P] ;
    double x[ROW_P] ;

#ifdef OTHER_DECOMP
    for(int row=0 ; row < N ; row++) {
	c[row]=0. ; // initialize for following step
	x[row]=0. ; // initialize for following step

	for(int col=0 ; col <= row ; col++) {
	    double sum=0. ;

	    if(col == row) {
		// diagonal
		for(int k=0 ; k < col ; k++) {
		    sum += lower[col][k] * lower[col][k] ;
		}

		lower[col][col] = sqrt(mat[col][col]-sum) ;

	    } else {
		for(int k=0 ; k < col ; k++) {
		    sum += lower[row][k] * lower[col][k] ;
		}

		lower[row][col] = (mat[row][col]-sum)/lower[col][col] ;
		lower[col][row] = 0. ;
	    }
	}
    }

    // Cholesky-Banachiewicz
    for(int i=0 ; i < N ; i++) {
	c[i]=0. ; // initialize for following step
	x[i]=0. ; // initialize for following step

	for(int j=0 ; j <= i ; j++) {
	    double sum=0. ;

	    for(int k=0 ; k < j ; k++)
		sum += lower[i][k] * lower[i][k] ;

	    if(col == row) {
		// diagonal
		lower[i][j] = sqrt(mat[i][i]-sum) ;

	    } else {
		lower[i][j] = (mat[i][j]-sum)/lower[j][j] ;
/*  		lower[j][i] = 0. ;  */
	    }
	}
    }
#endif

    // Cholesky-Crout
    for(int j=0 ; j < N ; j++) {
	c[j]=0. ; // initialize for following step
	x[j]=0. ; // initialize for following step

	double sum=0. ;

	for(int k=0 ; k < j ; k++)
	    sum += lower[j][k] * lower[j][k] ;

	lower[j][j] = sqrt(mat[j][j]-sum) ;

	for(int i=j+1 ; i < N ; i++) {
	    sum=0. ;

	    for(int k=0 ; k < j ; k++)
		sum += lower[i][k] * lower[i][k] ;

	    lower[i][j] = (mat[i][j]-sum)/lower[j][j] ;
/*  	    lower[j][i] = 0. ;  */
	}
    }

/*      fprintf(stderr, "************ cholesky lower ************ \n") ;  */
/*      for(int row = 0 ; row < N ; row++) {  */
/*  	for(int col = 0 ; col < N ; col++) {  */
/*  	    fprintf(stderr, " %10g", lower[row][col]) ;  */
/*  	}  */
/*  	fprintf(stderr, "\n") ;  */
/*      }  */

    // efficiently compute c as lower*c=b, where b is at column N of input matrix
    // only use lower triangular elements of lower

    for(int row=0 ; row < N ; row++) {
	double sum = 0. ;
/*  	fprintf(stderr, "c[%d] = (%5g", row, mat[row][N]) ;  */
	for(int col=0 ; col < row ; col++) {
/*  	    fprintf(stderr, " - %5g",lower[row][col]*c[col]) ;  */
	    sum += lower[row][col]*c[col] ;
	}
/*  	fprintf(stderr, ") /  %5g\n",lower[row][row]) ;  */
	c[row] = (mat[row][N]-sum) / lower[row][row] ;
    }

/*      fprintf(stderr, "************ cholesky c ************ \n") ;  */
/*      for(int row = 0 ; row < N ; row++) {  */
/*  	fprintf(stderr, " %10g", c[row]) ;  */
/*      }  */
/*      fprintf(stderr, "\n") ;  */

    // efficiently compute x as c/transpose of lower
    // only use upper triangular elements of transpose of lower

    for(int row=N-1 ; row >= 0 ; row--) {
/*  	fprintf(stderr, "x[%d] = (%5g", row, c[row]) ;  */
	double sum = 0. ;
	for(int col=N-1 ; col > row ; col--) {
/*  	    fprintf(stderr, " - %5g*%5g",lower[col][row],x[col]) ;  */
	    sum += lower[col][row]*x[col] ;
	}
/*  	fprintf(stderr, ") /  %5g\n",lower[row][row]) ;  */
	x[row] = (c[row] - sum) / lower[row][row] ;
	// store result in last column of mat, where the calling routine expects to find it.
	mat[row][N] = x[row] ;
    }

/*      fprintf(stderr, "************ cholesky x ************ \n") ;  */
/*      for(int row = 0 ; row < N ; row++) {  */
/*  	fprintf(stderr, " %10g", x[row]) ;  */
/*      }  */
/*      fprintf(stderr, "\n") ;  */
}





void matrix_solve(double coef[][ROW_P], short N)
/*  double coef[MAX_P+1][MAX_P+2] ;  */
/*  short N ;  */
{

#ifdef GAUSS
    gauss(coef, N) ;
#else
/*      cholesky(coef, N) ;  */
    lu(coef, N) ;
#endif
}

void lu(double coef[][ROW_P], short N)
/*  double coef[MAX_P+1][MAX_P+2] ;  */
/*  short N ;  */
{
    double d ;
    short indx[MAX_P] ;
    ludcmp(coef, N, indx, &d) ;
    lubksb(coef, N, indx) ;
}

void ludcmp(double a[][ROW_P], short N, short *indx, double *d)
{
    double vv[MAX_P+2], dum, sum, aamax ;
    short imax = 0 ;

/** print_mat(a, N) ; **/

    *d = 1. ;

    for(i = 0 ; i < N ; i++) {
	aamax = 0. ;
	for(j = 0 ; j < N ; j++)
	    if(abs(a[i][j]) > aamax)
		aamax = abs(a[i][j]) ;
	if(aamax <= 1.e-200)
	    dv_exit("singular matrix") ;
	vv[i] = 1./aamax ;
    }

    for(j = 0 ; j < N ; j++) {
	aamax = 0. ;

	if(j > 0) {
	    for(i = 0 ; i < j ; i++) {
		sum = a[i][j] ;
		if(i > 0) {
		    for(k = 0 ; k < i ; k++)
			sum -= a[i][k]*a[k][j] ;
		    a[i][j] = sum ;
		}
	    }
	}

	for(i = j ; i < N ; i++) {
	    sum = a[i][j] ;
	    if(j > 0) {
		for(k = 0 ; k < j ; k++)
		    sum -= a[i][k]*a[k][j] ;
		a[i][j] = sum ;
	    }

	    dum = vv[i]*abs(sum) ;

	    if(dum >= aamax) {
		imax = i ;
		aamax = dum ;
	    }
	}

/** imax = j ; **/
/** printf("Pivot %d\n", imax) ; **/

	if(j != imax) {
	    for(k = 0 ; k < N ; k++) {
		dum = a[imax][k] ;
		a[imax][k] = a[j][k] ;
		a[j][k] = dum ;
	    }
	    *d = -*d ;
	    vv[imax] = vv[j] ;
	}

	indx[j] = imax ;

	if( j != N-1 ) {
	    if(abs(a[j][j]) < tiny)
		a[j][j] = tiny ;
	    dum = 1./a[j][j] ;
	    for(i = j+1 ; i < N ; i++)
		a[i][j] *= dum ;
	}
    }

    if(abs(a[N-1][N-1]) < tiny)
	a[N-1][N-1] = tiny ;
}

#define b(x) a[x][N]

void lubksb( double a[][ROW_P], short N, short *indx )
{
    short ii = -1 ;

    for(i = 0 ; i < N ; i++) {
	short ll = indx[i] ;
	double sum = a[ll][N] ;

	a[ll][N] = a[i][N] ;

	if(ii != -1) {
	    for(j = ii ; j < i ; j++) {
		sum -= a[i][j]*a[j][N] ;
	    }
	} else if(sum != 0.) {
	    ii = i ;
	}

	a[i][N] = sum ;
    }

    for(i = N-1 ; i >= 0 ; i--) {
	double sum = a[i][N] ;

	if(i < N-1) {
	    for(j = i+1 ; j < N ; j++)
		sum -= a[i][j]*a[j][N] ;
	}
	a[i][N] = sum/a[i][i] ;
    }
}

void gauss_inv( double coef[][ROW_P2], short N )
{
    /** gauss elimination **/
    for(row = 0 ; row < N ; row++) {
	for(irow = row+1 ; irow < N ; irow++) {
	    f = coef[row][row] / coef[irow][row] ;
	    for(col = row ; col < N*2 ; col++) {
		coef[irow][col] *= f ;
		coef[irow][col] -= coef[row][col] ;
	    }
	}
    }

    /** work back up leaving 1's on the diagonal **/
    for(row = N-1 ; row > -1 ; row--) {
	for(irow = row-1 ; irow > -1 ; irow--) {
	    f = coef[row][row] / coef[irow][row] ;
	    for(col = 0 ; col < N*2 ; col++) {
		coef[irow][col] *= f ;
		coef[irow][col] -= coef[row][col] ;
	    }
	}
	f = coef[row][row] ;
	for(col = 0 ; col < N*2 ; col++)
	    coef[row][col] /= f ;
    }
}

void gauss( double coef[][ROW_P2], short N )
{
    /** gauss elimination **/
    for(row = 0 ; row < N+1 ; row++) {
	for(irow = row+1 ; irow < N+2 ; irow++) {
	    f = coef[row][row] / coef[irow][row] ;
	    for(col = row ; col < N+2 ; col++) {
		coef[irow][col] *= f ;
		coef[irow][col] -= coef[row][col] ;
	    }
	}
    }

    /** work back up leaving 1's on the diagonal **/
    for(row = N ; row > -1 ; row--) {
	for(irow = row-1 ; irow > -1 ; irow--) {
	    f = coef[row][row] / coef[irow][row] ;
	    for(col = 0 ; col < N+2 ; col++) {
		coef[irow][col] *= f ;
		coef[irow][col] -= coef[row][col] ;
	    }
	}
	f = coef[row][row] ;
	for(col = 0 ; col < N+2 ; col++)
	    coef[row][col] /= f ;
    }
}

struct mstat_75 init_reg_75(int P)
{
    struct mstat_75 m ;
    m.P = P ;

    if(m.P > MAX_P) {
	fprintf(stderr, "reg: more than %d independent variables\n", MAX_P) ;
	exit(1) ;
    }

    m.nobs=0 ;
    m.sse=0. ;
    m.syy=0 ;
    m.sy=0 ;
    m.wghtn=0. ;
    m.rbarsize=(P+1)*k/2+1 ;
    m.rsq=0. ;
    m.adjrsq=0. ;
    m.sigmasq=0. ;

    for(int i=0 ; i < P ; i++) {
	m.d[i]=0. ;
	m.thetabar[i]=0. ;
	m.rbar[i]=0. ;
	m.theta[i]=0. ;
	m.standard_error[i]=0. ;
    }

    return m ;

}

int sum_reg_75(struct mstat_75 *m, double data[], double yelement, double weight)
{
    double xcopy[MAX_P] ;
    for(int i=1 ; i <= m->P ; i++) {
	xcopy[i] = data[i-1] ;
    }

    m->syy += weight*yelement*yelement ;
    m->sy += weight*yelement ;

    if(weight >= 0.0)
	m->nobs++ ;
    else
	m->nobs-- ;

    m->wghtn += weight ;

    for(int i=1 ; i <= m->P ; i++) {
	if(weight < 1.e-30) return m->nobs ;
	if(fabsf(xcopy[i]) > TINY) {
	    double xi = xcopy[i] ;

	    double di = m->d[i] ;
	    double dprimei = di+weight*(xi*xi) ;
	    double cbar = di/dprimei ;
	    double sbar=weight*xi/dprimei ;
	    weight *= cbar ;
	    m->d[i] = dprimei ;
	    int nextr = (int)( (float)((i-1)*(2.*m->P-i))/2.+1.) ;
	    if( !(nextr<=m->rbarsize)) {
		fprintf(stderr, "FATAL nextr (%d) <= rbarsize (%d)\n", nextr, m->rbarsize) ;
		exit(1) ;
	    }
	    double xP ;
	    for(int Pc=i+1 ; Pc <= m->P ; Pc++) {
		xP=xcopy[Pc] ; xcopy[Pc] = xP-xi*m->rbar[nextr] ;
		m->rbar[nextr] = cbar * m->rbar[nextr] + sbar*xP ;
		++nextr ;
	    }
	    xP=yelement ; yelement -= xi*m->thetabar[i] ;
	    m->thetabar[i]=cbar*m->thetabar[i]+sbar*xP ;
	}
    }
    m->sse += weight*yelement*yelement ;

    return m->nobs ;
}

int rbrindex(int a, int b, int K)
{
    return (a == b) ? -9 :
	((((float)a-1.)* (2.*(float)K-(float)a)/2.+1.) + b - 1 - a);
}

double rbr(int a,int b, double *rbar, int K)
{
    return (a == b) ? 1 : ( (a>b) ? 0 : (rbar[rbrindex(a,b,K)]));
}

void estimate_reg_75(struct mstat_75 *m, double *B)
{

    // compute coefficients
    for(int i=m->P ; i >= 1 ; i--) {
	m->theta[i] = m->thetabar[i] ;
	int nextr = (int)( (double)((i-1)*(2.*m->P-i))/2.+1.) ;
	for(int Pc=i+1 ; Pc <= m->P ; Pc++) {
	    m->theta[i] -= m->rbar[nextr]*m->theta[Pc] ;
	    nextr++ ;
	}
    }

    for(int i=m->P ; i >= 1 ; i--) {
	B[i-1] = m->theta[i] ;
    }

    m->ybar = m->sy / m->wghtn ;
    m->sst = m->syy - m->wghtn*m->ybar*m->ybar ;
    m->rsq = 1. - m->sse/m->sst ;
    m->sigmasq = m->nobs <= m->P ? 1.e300 : m->sse/((double)(m->nobs-m->P)) ;

    // compute standard errors
    // note, in orginal algorithm, K is used instead of P for number of independent variables, so the code below will keep with that convention
    int K = m->P ;

    double U[MAX_P*MAX_P+1] ;

#define ui(a,b) ( (K*(a-1))+(b-1) )
#define u(a,b) ( U[ui(a,b)])
#define setu(a,b,c) ( U[ui(a,b)] = c )
#define add2u(a,b,c) ( U[ui(a,b)] += c )
#define mult2u(a,b,c) ( U[ui(a,b)] *= c )

    for(int i=0 ; i <= K*K ; i++) {
	U[i] = 0. ;
    }

    for (int j=K; j>=1; --j) {
	setu(j,j, 1.0/(rbr(j,j,m->rbar,K)));
	for (int k=j-1; k>=1; --k) {
	    setu(k,j,0);
	    for (i=k+1; i<=j; ++i) { add2u(k,j, rbr(k,i,m->rbar,K)*u(i,j)); }
	    mult2u(k,j, (-1.0)/rbr(k,k,m->rbar,K));
	}
    }

#ifdef DEBUG75
    fprintf(stderr, "inverse of R is:\n") ;
    for(int i=1 ; i <= K ; i++) {
	fprintf(stderr, "[%d]:\t", i) ;
	for(int j=1 ; j <= K ; j++) {
	    fprintf(stderr, "%f\t", U[ui(i,j)]) ;
	}
	fprintf(stderr, "\n") ;
    }
#endif

    for (int i=1;i<=K;++i) { 
	for (int j=1;j<=K;++j) {
	    if (abs(m->d[j])<TINY) {
		mult2u(i,j, sqrt(1.0/TINY));
		if (abs(m->d[j])==0.0) {
		    fprintf(stderr,"FATAL cannot compute theta covariance matrix for variable %d\n", j) ;
		    exit(1) ;
		}
	    } else {
		mult2u(i,j, sqrt(1.0/m->d[j]));
	    }
	}
    }

    m->sigmasq= (m->nobs<=K) ? 1.e300 : (m->sse/(double)(m->nobs - K));

    double xpxinv[MAX_P+1] ;

    for (int i=1;i<=K; ++i) {
	for (int j=i;j<=K;++j) {
	    int indexij= ui(i,j);
	    xpxinv[indexij]= 0.0;
	    for(int k=1;k<=K;++k) {
		xpxinv[indexij] += U[ui(i,k)]*U[ui(j,k)];
	    }
	    xpxinv[ui(j,i)]= xpxinv[indexij]; // this is symmetric
	}
    }

#ifdef DEBUG75
    fprintf(stderr, "Full inverse of X'X is:\n") ;
    for(int i=1 ; i <= K ; i++) {
	fprintf(stderr, "[%d]:\t", i) ;
	for(int j=1 ; j <= K ; j++) {
	    fprintf(stderr, "%f\t", xpxinv[ui(i,j)]) ;
	}
	fprintf(stderr, "\n") ;
    }

    printf("sigma^2 is %f\n", m->sigmasq) ;

    double secoefs[MAX_P+1] ;
    for (int i=1; i<=K; ++i) {
	int j=i ;
	int indexij= ui(i,j);
	printf("X'X at %d,%d (%d) is %f\n", i, j, indexij, xpxinv[indexij]) ;
	m->standard_error[i] = sqrt(xpxinv[indexij] * m->sigmasq);
    }





    printf("%-15s\t%12s\t%12s\t%7s\n", "Name", "Theta", "StdErr", "T-stat") ;

    for(int i=1; i<= m->k; ++i) {
	char name[20] ;
	sprintf(name, "B[%d]", i) ;
	printf("%-15s\t", name);
	printf("%12.4f\t", m->theta[i]);
	printf("%12.4f\t", m->standarderrors[i]);
	printf("%7.2f", (m->theta[i]/m->standard_error[i]) );
	printf("\n");
    }

    printf("\nR^2= %.3f N=%d K=%d", m->rsq, m->nobs, m->k) ;
    printf("****************************************************************\n");
#endif


}

struct mstat init_reg(int n_indep)
{
    struct mstat m ;
    m.P = n_indep ;
    m.sum_y=0. ;
    m.sum_y_sq=0. ;
    m.nobs=0 ;
    m.sum_wgts=0. ;

    if(m.P > MAX_P) {
	fprintf(stderr, "reg: more than %d independent variables\n", MAX_P) ;
	exit(1) ;
    }

    /** zero the coef array which will hold the sums **/
    for(row = 0 ; row <= m.P  ; row++)
	for(col = 0 ; col <= m.P  ; col++)
	    m.coef[row][col] = 0 ;

    return m ;

}

void print_reg_matrix(struct mstat *m, char *name)
{
    fprintf(stderr, " mstat matrix for %s\n", name) ;

    for(row = 0 ; row <= m->P ; row++) {
	for(col = 0 ; col <= m->P ; col++) {
	    fprintf(stderr, " %10g", m->coef[row][col]) ;
	}
	fprintf(stderr, "\n") ;
    }

    fprintf(stderr, "\n") ;

}

void sum_reg_int(struct mstat *m, double data[], double y)
{
    double row_data[MAX_P+1] ;

    for(i = 1 ; i <= m->P  ; i++)
        row_data[i] = data[i-1] ;

    row_data[m->P+1] = y ;
    row_data[0] = 1. ;

    for(row = 0 ; row <= m->P ; row++)
        for(col = 0 ; col <= m->P+1 ; col++)
            m->coef[row][col] += row_data[row] * row_data[col] ;
}


void sum_reg(struct mstat *m, double data[], double y, double wgt)
{
    double row_data[MAX_P+1] ;

    for(i = 0 ; i < m->P  ; i++) {
	row_data[i] = data[i] ;
    }

    row_data[m->P] = y ; 

    if(wgt > 0.) 
	m->nobs++ ;
    else
	m->nobs-- ;

    m->sum_wgts += wgt ;
    m->sum_y += y*wgt ;
    m->sum_y_sq += y*y*wgt ;

    for(row = 0 ; row <= m->P ; row++) {
	for(col = 0 ; col <= m->P ; col++) {
	    m->coef[row][col] += wgt * row_data[row] * row_data[col] ;
	}
    }

}

void estimate_reg_cholesky(struct mstat *m, double *B)
{
    double XprimeY[MAX_P+1] ;

    for(col = 0 ; col <= m->P ; col++)
	XprimeY[col] = m->coef[m->P][col] ;

    cholesky(m->coef, m->P) ;

    for(row = 0 ; row < m->P ; row++) {
	B[row] = m->coef[row][m->P] ;
    }

    double B_XprimeY=0 ;
    for(col = 0 ; col < m->P ; col++)
	B_XprimeY += B[col] * XprimeY[col] ;

/*      SSE = Y'Y - B' * X'Y ;  */
    m->sse = m->sum_y_sq - B_XprimeY ;
}

void estimate_reg(struct mstat *m, double *B)
{
    double XprimeY[MAX_P+1] ;

    for(col = 0 ; col <= m->P ; col++) {
	XprimeY[col] = m->coef[m->P][col] ;
    }

    double mat[MSTAT_MAX_P+1][MSTAT_ROW_P] ;

    for(row = 0 ; row <= m->P ; row++) {
	for(col = 0 ; col <= m->P ; col++) {
	    mat[row][col] = m->coef[row][col] ;
	}
    }

    //fprintf(stderr, "\n") ;
    matrix_solve(mat, m->P) ;


    for(row = 0 ; row < m->P ; row++) {
	B[row] = mat[row][m->P] ;
	//fprintf(stderr, "est_reg B%d=%f\n", row,b[row]) ;
    }

    double B_XprimeY=0. ;

    for(col = 0 ; col < m->P ; col++) {
	B_XprimeY += B[col] * XprimeY[col] ;
    }

/*      SSE = Y'Y - B' * X'Y ;  */
    m->sse = m->sum_y_sq - B_XprimeY ;
}

void estimate_reg_full(struct mstat *m, double *B)
{
    double XprimeY[MAX_P+1] ;

    for(col = 0 ; col <= m->P ; col++)
	XprimeY[col] = m->coef[m->P][col] ;

    //fprintf(stderr, "\n") ;
    matrix_solve(m->coef, m->P) ;


    for(row = 0 ; row < m->P ; row++) {
	B[row] = m->coef[row][m->P] ;
	m->B[row] = m->coef[row][m->P] ;
	//fprintf(stderr, "est_reg B%d=%f\n", row,b[row]) ;
    }

    double B_XprimeY=0 ;
    for(col = 0 ; col < m->P ; col++)
	B_XprimeY += B[col] * XprimeY[col] ;

/*      SSE = Y'Y - B' * X'Y ;  */
    m->sse = m->sum_y_sq - B_XprimeY ;

    m->sum_y_sq = m->coef[m->P][m->P] ;
    m->mean_y = m->sum_y/m->sum_wgts ;
    m->ssto = m->sum_y_sq - m->sum_y*m->sum_y/m->sum_wgts ;

    double XprimeX_plusI[ROW_P][ROW_P2] ;
    double XprimeX[ROW_P][ROW_P] ;
    double XprimeX_inv[ROW_P][ROW_P] ;

    fprintf(stderr, "X\'X\n") ;
    for(row = 0 ; row < m->P ; row++) {
	for(col = 0 ; col < m->P ; col++) {
	    XprimeX[row][col] = m->coef[row][col] ;
	    XprimeX_plusI[row][col] = m->coef[row][col] ;
	    fprintf(stderr, "  %5g,", XprimeX[row][col] ) ;
	}
	for(col = m->P ; col < m->P*2 ; col++) {
	    if((col-m->P) == row) {
		XprimeX_plusI[row][col] = 1. ;
	    } else {
		XprimeX_plusI[row][col] = 0. ;
	    }
	    fprintf(stderr, "  %5g,", XprimeX_plusI[row][col] ) ;
	}
	fprintf(stderr, "\n") ;
    }

    gauss_inv(XprimeX_plusI, m->P) ;

    fprintf(stderr, "X\'X inverse\n") ;
    for(row = 0 ; row < m->P ; row++) {
	for(col = m->P ; col < m->P*2 ; col++) {
	    XprimeX_inv[row][col-m->P] = XprimeX_plusI[row][col] ;
	    fprintf(stderr, "  %5g,", XprimeX_plusI[row][col] ) ;
	}
	fprintf(stderr, "\n") ;
    }







    m->ssr = m->ssto - m->sse ;
    m->rsq = 1.-m->sse/m->ssto ;
    m->sigmasq = m->sse/((double)(m->sum_wgts-m->P)) ;
    printf("sum_wgts is %f\n", m->sum_wgts) ;
    printf("sigma^2 is %f\n", m->sigmasq) ;

    for(row = 0 ; row < m->P ; row++)
	m->standard_error[row] = sqrt(XprimeX_inv[row][row] * m->sigmasq);

    printf("%-15s\t%12s\t%12s\t%7s\n", "Name", "Beta", "StdErr", "T-stat") ;

    for(int i=0; i< m->P; ++i) {
	char name[20] ;
	sprintf(name, "B[%d]", i) ;
	printf("%-15s\t", name);
	printf("%12.4f\t", m->B[i]);
	printf("%12.4f\t", m->standard_error[i]);
	printf("%7.2f", (m->B[i]/m->standard_error[i]) );
	printf("\n");
    }

    printf("\nR^2= %.3f N=%d P=%d", m->rsq, m->nobs, m->P) ;
    printf("****************************************************************\n");
}


#include <math.h>                    // required for powl(), fabsl(),
                                     // expl() and logl().
#include <float.h>                   // required for LDBL_EPSILON.
#include <limits.h>    // required for LONG_MAX

static const long double ln_LDBL_MAX =  1.13565234062941435e+4L;

static long double const e =  2.71828182845904523536028747L;
static long double const pi = 3.14159265358979323846264338L;
static long double const g =  9.65657815377331589457187L;
static long double const exp_g_o_sqrt_2pi = +6.23316569877722552586386e+3L;
static double max_double_arg = 171.0;
static long double max_long_double_arg = 1755.5L;

static long double const a[] = {
                                 +1.14400529453851095667309e+4L,
                                 -3.23988020152318335053598e+4L,
                                 +3.50514523505571666566083e+4L,
                                 -1.81641309541260702610647e+4L,
                                 +4.63232990536666818409138e+3L,
                                 -5.36976777703356780555748e+2L,
                                 +2.28754473395181007645155e+1L,
                                 -2.17925748738865115560082e-1L,
                                 +1.08314836272589368860689e-4L
                              };

static const long double log_sqrt_2pi = 9.18938533204672741780329736e-1L;
// Bernoulli numbers B(2),B(4),B(6),...,B(20).  Only B(2),...,B(6) currently //
// used.                                                                     //

static const long double Bernoulli_numbers[] = {   1.0L / (long double)(6 * 2 * 1),
                                  -1.0L / (long double)(30 * 4 * 3),
                                   1.0L / (long double)(42 * 6 * 5),
                                  -1.0L / (long double)(30 * 8 * 7),
                                   5.0L / (long double)(66 * 10 * 9),
                                -691.0L / (long double)(2730 * 12 * 11),
                                   7.0L / (long double)(6 * 14 * 13),
                               -3617.0L / (long double)(510 * 16 * 15),
                               43867.0L / (long double)(796 * 18 * 17),
                             -174611.0L / (long double)(330 * 20 * 19) 
                           };

/*  static const int nBernoulli = sizeof(Bernoulli_numbers) / sizeof(long double);  */

static long double Beta_Continued_Fraction( long double x, long double a, long double b);
static long double xBeta_Distribution(double x, double a, double b);
static long double xLnGamma_Asymptotic_Expansion( long double x ) ;
static long double Duplication_Formula( long double two_x ) ;
static long double xGamma(long double x) ;

//                         Internally Defined Routines                        //
double Beta_Distribution(double x, double a, double b);
double Ln_Beta_Function(double a, double b);
long double xBeta_Function(long double a, long double b) ;
double Gamma_Function(double x) ;
long double xLn_Beta_Function(long double a, long double b) ;
long double xGamma_Function(long double x) ;
double Gamma_Function_Max_Arg( void ) ;
long double xLn_Gamma_Function(long double x) ;
long double xGamma_Function_Max_Arg( void ) ;

double Beta_Distribution(double x, double a, double b)
{

   if ( x <= 0.0 ) return 0.0;
   if ( x >= 1.0 ) return 1.0;

   return (double) xBeta_Distribution( x, a, b);
}

static long double xBeta_Distribution(double xx, double aa, double bb) 
{
   long double x = (long double) xx;
   long double a = (long double) aa;
   long double b = (long double) bb;

           /* Both shape parameters are strictly greater than 1. */

   if ( aa > 1.0 && bb > 1.0 ) {
      if ( x <= (a - 1.0L) / ( a + b - 2.0L ) )
         return Beta_Continued_Fraction(x, a, b);
      else
         return 1.0L - Beta_Continued_Fraction( 1.0L - x, b, a );
    }
  
             /* Both shape parameters are strictly less than 1. */

   if ( aa < 1.0 && bb < 1.0 )  
      return (a * xBeta_Distribution(xx, aa + 1.0, bb) 
                      + b * xBeta_Distribution(xx, aa, bb + 1.0) ) / (a + b); 
   
              /* One of the shape parameters exactly equals 1. */

   if ( aa == 1.0 )
      return 1.0L - powl(1.0L - x, b) / ( b * xBeta_Function(a,b) );

   if ( bb == 1.0 ) return powl(x, a) / ( a * xBeta_Function(a,b) );

      /* Exactly one of the shape parameters is strictly less than 1. */

   if ( aa < 1.0 )  
      return xBeta_Distribution(xx, aa + 1.0, bb)
            + powl(x, a) * powl(1.0L - x, b) / ( a * xBeta_Function(a,b) );
 
                   /* The remaining condition is b < 1.0 */

   return xBeta_Distribution(xx, aa, bb + 1.0)
            - powl(x, a) * powl(1.0L - x, b) / ( b * xBeta_Function(a,b) );
}

static long double Beta_Continued_Fraction( long double x, long double a,
                                                                 long double b)
{
   long double Am1 = 1.0L;
   long double A0 = 0.0L;
   long double Bm1 = 0.0L;
   long double B0 = 1.0L;
   long double e = 1.0L;
   long double Ap1 = A0 + e * Am1;
   long double Bp1 = B0 + e * Bm1;
   long double f_less = Ap1 / Bp1;
   long double f_greater = 0.0L;
   long double aj = a;
   long double am = a;
   static long double eps = 10.0L * LDBL_EPSILON;
   int j = 0;
   int m = 0;
   int k = 1;

   if ( x == 0.0L ) return 0.0L;
   
   while ( (2.0L * fabsl(f_greater - f_less) > eps * fabsl(f_greater + f_less)) ) {
      Am1 = A0;
      A0 = Ap1;
      Bm1 = B0;
      B0 = Bp1;
      am = a + m;
      e = - am * (am + b) * x / ( (aj + 1.0L) * aj );
      Ap1 = A0 + e * Am1;
      Bp1 = B0 + e * Bm1;
      k = (k + 1) & 3;
      if (k == 1) f_less = Ap1/Bp1;
      else if (k == 3) f_greater = Ap1/Bp1;
      if ( fabsl(Bp1) > 1.0L) {
         Am1 = A0 / Bp1;
         A0 = Ap1 / Bp1;
         Bm1 = B0 / Bp1;
         B0 = 1.0;
      } else {
         Am1 = A0;
         A0 = Ap1;
         Bm1 = B0;
         B0 = Bp1;
      }
      m++;
      j += 2;
      aj = a + j;
      e = m * ( b - m ) * x / ( ( aj - 1.0L) * aj  );
      Ap1 = A0 + e * Am1;
      Bp1 = B0 + e * Bm1;
      k = (k + 1) & 3;
      if (k == 1) f_less = Ap1/Bp1;
      else if (k == 3) f_greater = Ap1/Bp1;
   }
   return expl( a * logl(x) + b * logl(1.0L - x) + logl(Ap1 / Bp1) ) /
                                                ( a * xBeta_Function(a,b) );
}

double Beta_Function(double a, double b)
{
   long double beta = xBeta_Function( (long double) a, (long double) b);
   return (beta < DBL_MAX) ? (double) beta : DBL_MAX;
}

long double xBeta_Function(long double a, long double b)
{
   long double lnbeta;

     // If (a + b) <= Gamma_Function_Max_Arg() then simply return //
     //  gamma(a)*gamma(b) / gamma(a+b).                          //

   if ( (a + b) <= Gamma_Function_Max_Arg() )
      return xGamma_Function(a) / (xGamma_Function(a + b) / xGamma_Function(b));

     // If (a + b) > Gamma_Function_Max_Arg() then simply return //
     //  exp(lngamma(a) + lngamma(b) - lngamma(a+b) ).           //

   lnbeta = xLn_Gamma_Function(a) + xLn_Gamma_Function(b)
                                                 - xLn_Gamma_Function(a + b);
   return (lnbeta > ln_LDBL_MAX) ? (long double) LDBL_MAX : expl(lnbeta);
}

double Gamma_Function(double x)
{
   long double g;

   if ( x > max_double_arg ) return DBL_MAX;
   g = xGamma_Function( (long double) x);
   if (fabsl(g) < DBL_MAX) return (double) g;
   return (g < 0.0L) ? -DBL_MAX : DBL_MAX;

}

long double xGamma_Function(long double x)
{
   long double sin_x;
   long double rg;
   long int ix;

             // For a positive argument (x > 0)                 //
             //    if x <= max_long_double_arg return Gamma(x)  //
             //    otherwise      return LDBL_MAX.              //

   if ( x > 0.0L ) {
      if (x <= max_long_double_arg)
	  return xGamma(x);
      else
	  return LDBL_MAX;
   }

                   // For a nonpositive argument (x <= 0) //
                   //    if x is a pole return LDBL_MAX   //

   if ( x > -(long double)LONG_MAX) {
      ix = (long int) x;
      if ( x == (long double) ix) return LDBL_MAX;
   }
   sin_x = sinl(pi * x);
   if ( sin_x == 0.0L ) return LDBL_MAX;

          // if x is not a pole and x < -(max_long_double_arg - 1) //
          //                                     then return 0.0L  //

   if ( x < -(max_long_double_arg - 1.0L) ) return 0.0L;

            // if x is not a pole and x >= -(max_long_double - 1) //
            //                               then return Gamma(x) //

   rg = xGamma(1.0L - x) * sin_x / pi;
   if ( rg != 0.0L ) return (1.0L / rg);
   return LDBL_MAX;
}

static long double xGamma(long double x)
{

   long double xx = (x < 1.0L) ? x + 1.0L : x;
   long double temp;
   int const n = sizeof(a) / sizeof(long double);
   int i;

   if (x > 1755.5L) return LDBL_MAX;

   if (x > 900.0L) return Duplication_Formula(x);

   temp = 0.0L;
   for (i = n-1; i >= 0; i--) {
      temp += ( a[i] / (xx + (long double) i) );
   }
   temp += 1.0L;
   temp *= ( powl((g + xx - 0.5L) / e, xx - 0.5L) / exp_g_o_sqrt_2pi );
   return (x < 1.0L) ?  temp / x : temp;
}

static long double Duplication_Formula( long double two_x )
{
   long double x = 0.5L * two_x;
   long double g;
   int n = (int) two_x - 1;

   g = powl(2.0L, two_x - 1.0L - (long double) n);
   g = ldexpl(g,n);
   g /= sqrt(pi);
   g *= xGamma_Function(x);
   g *= xGamma_Function(x + 0.5L);

   return g;
}

double Gamma_Function_Max_Arg( void ) { return max_double_arg; }

long double xGamma_Function_Max_Arg( void ) { return max_long_double_arg; }

double Ln_Beta_Function(double a, double b)
{
   return (double) xLn_Beta_Function( (long double) a, (long double) b );
}

long double xLn_Beta_Function(long double a, long double b)
{

     // If (a + b) <= Gamma_Function_Max_Arg() then simply return //
     //  log(gamma(a)*gamma(b) / gamma(a+b)).                     //

   if ( (a + b) <= (long double) Gamma_Function_Max_Arg() ) {
      if ( a == 1.0L && b == 1.0L )
	  return 0.0L;
      else
	  return logl( xGamma_Function(a) / ( xGamma_Function(a + b) / xGamma_Function(b) ));
   }

     // If (a + b) > Gamma_Function_Max_Arg() then simply return //
     //  lngamma(a) + lngamma(b) - lngamma(a+b).                 //

   return xLn_Gamma_Function(a) + xLn_Gamma_Function(b) - xLn_Gamma_Function(a+b);
}

double Ln_Gamma_Function(double x)
{

       // For a positive argument, 0 < x <= Gamma_Function_Max_Arg() //
       // then  return log Gamma(x).                                 //

   if (x <= Gamma_Function_Max_Arg()) return log(Gamma_Function(x));

    // otherwise return result from asymptotic expansion of ln Gamma(x). //

   return (double) xLnGamma_Asymptotic_Expansion( (long double) x );
}

long double xLn_Gamma_Function(long double x)
{

       // For a positive argument, 0 < x <= Gamma_Function_Max_Arg() //
       // then  return log Gamma(x).                                 //

   if (x <= Gamma_Function_Max_Arg()) return logl(xGamma_Function(x));

    // otherwise return result from asymptotic expansion of ln Gamma(x). //

   return xLnGamma_Asymptotic_Expansion( x );
}


static long double xLnGamma_Asymptotic_Expansion( long double x ) {
   const int  m = 3;
   long double term[3];
   long double sum = 0.0L;
   long double xx = x * x;
   long double xj = x;
   long double lngamma = log_sqrt_2pi - xj + (xj - 0.5L) * logl(xj);
   int i;

   for (i = 0; i < m; i++) { term[i] = Bernoulli_numbers[i] / xj; xj *= xx; }
   for (i = m - 1; i >= 0; i--) sum += term[i]; 
   return lngamma + sum;
}


double Student_t_Density( double x, int n )
{
   double ln_density;

   ln_density = -(double)(n+1)/2.0 * log(1.0 + x * x /(double)n)
                - 0.5*log((double)n)
                - Ln_Beta_Function(0.5 * (double)n, 0.5);

   return exp(ln_density);
}

//                         Externally Defined Routines                        //
extern double Beta_Distribution(double x, double a, double b);

double Student_t_Distribution(double x, int n)
{
   double a = (double) n / 2.0;
   double beta = Beta_Distribution( 1.0 / (1.0 + x * x / n), a, 0.5);

   if ( x > 0.0 ) return 1.0 - 0.5 * beta;
   else if ( x < 0.0) return 0.5 * beta;
   return 0.5;
}

double Student_t_probability(double t, int df)
{
    return Student_t_Distribution(t, df) ;
}
