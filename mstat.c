/*
	FILE mstat.c
	PURPOSE - multiple linear regression
	CONTENTS -
*/

#include <stdio.h>
#include <stdlib.h>

#define MAX_N MSTAT_MAX_N
#define ROW_N MSTAT_ROW_N

#include "mstat.h"

void gauss( double [][ROW_N], short) ;
void lu(double [][ROW_N], short) ;
void lubksb( double [][ROW_N], short, short *) ;
void ludcmp(double [][ROW_N], short, short *, double *) ;
void matrix_solve(double [][ROW_N], short) ;

double y, f ;
int row, col, irow, i, j, k ;



#define abs(x) ((x) > 0 ? (x) : -(x))
#define tiny 1e-200


void dv_exit(char *msg)
{
}

void matrix_solve(double coef[][ROW_N], short N)
/*  double coef[MAX_N+1][MAX_N+2] ;  */
/*  short N ;  */
{
/*      gauss(coef, N) ;  */
    lu(coef, N) ;
}

void lu(double coef[][ROW_N], short N)
/*  double coef[MAX_N+1][MAX_N+2] ;  */
/*  short N ;  */
{
    double d ;
    short indx[MAX_N] ;
    ludcmp(coef, N, indx, &d) ;
    lubksb(coef, N, indx) ;
}

void ludcmp(double a[][ROW_N], short N, short *indx, double *d)
{
    double vv[MAX_N+2], dum, sum, aamax ;
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

void lubksb( double a[][ROW_N], short N, short *indx )
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

#ifdef GAUSS
void gauss( double coef[][ROW_N], short N )
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
#endif

struct mstat init_reg(int n)
{
    struct mstat m ;
    m.N = n ;

    if(m.N > MAX_N) {
	fprintf(stderr, "reg: more than %d independent variables\n", MAX_N) ;
	exit(1) ;
    }

    /** zero the coef array which will hold the sums **/
    for(row = 0 ; row <= m.N  ; row++)
	for(col = 0 ; col <= m.N  ; col++)
	    m.coef[row][col] = 0 ;

    return m ;

}

void print_reg_matrix(struct mstat *m, char *name)
{
    fprintf(stderr, " mstat matrix for %s\n", name) ;

    for(row = 0 ; row <= m->N ; row++) {
	for(col = 0 ; col <= m->N ; col++) {
	    fprintf(stderr, " %10g", m->coef[row][col]) ;
	}
	fprintf(stderr, "\n") ;
    }

    fprintf(stderr, "\n") ;

}

void sum_reg(struct mstat *m, double data[], double y, double wgt)
{
    double row_data[MAX_N+1] ;

    for(i = 0 ; i < m->N  ; i++) {
	row_data[i] = data[i] ;
    }

    row_data[m->N] = y ; 

    for(row = 0 ; row <= m->N ; row++)
	for(col = 0 ; col <= m->N ; col++)
	    m->coef[row][col] += wgt * row_data[row] * row_data[col] ;

}

void estimate_reg(struct mstat *m, double *b)
{
    //fprintf(stderr, "est_reg\n") ;
    for(row = 0 ; row <= m->N ; row++) {
	for(col = 0 ; col <= m->N+1 ; col++) {
	    //fprintf(stderr, "  %f,", m->coef[row][col] ) ;
	}
	//fprintf(stderr, "\n") ;
    }

    //fprintf(stderr, "\n") ;
    matrix_solve(m->coef, m->N) ;

    for(row = 0 ; row < m->N ; row++) {
	b[row] = m->coef[row][m->N] ;
	//fprintf(stderr, "est_reg B%d=%f\n", row,b[row]) ;
    }
}
