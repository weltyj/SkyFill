#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amoeba_06.h"

#define NMAX 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0
#define GRID_ITER_LIMIT 20

/*** !!! IMPORTANT NOTE:  Indexing starts at 1, not zero, based on old FORTRAN code!!! ***/
/*** ... and, this is some pretty awful code which is based on the original Numerical Recipies in C */

#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++)\
						sum += p[i][j]; psum[j]=sum;}

static float constraints_lo[100], constraints_hi[100] ;
static int n_iter=0 ;

float EvaluateVertex(float (*function)(), float *p, int ndim)
{
    int j ;

    /* check constraints */
    for (j=1;j<=ndim;j++) {
	if(p[j] < constraints_lo[j]) p[j] = constraints_lo[j] ;
	if(p[j] > constraints_hi[j]) p[j] = constraints_hi[j] ;
    }

    n_iter++ ;

    return (*function)(p);
}


float amotry(float **p,float *y, float *psum,int ndim,float (*function)(),int ihi,float fac)
{
	int j;
	float fac1,fac2,ytry,*ptry ;

	ptry=(float *)calloc(ndim+1, sizeof(float));

	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;

	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;

	ytry=EvaluateVertex(function,ptry,ndim);

	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}

	free(ptry);

	return ytry;
}




int amoeba(float **p, float y[], int ndim,float ftol, int max_iter, float (*function)(), int verbose)
{
    int i,j,ilo,ihi,inhi,mpts=ndim+1;
    float ytry,ysave,sum,rtol,*psum ;
    int iter=0 ;
    n_iter=0 ;

    psum=(float *)calloc(ndim+1, sizeof(float));

    GET_PSUM

	while(1) {

	    // find lowest, highest and next highest vertices ;
	    ilo=1;

	    if(y[1] > y[2]) {
		ihi=1 ;
		inhi=2 ;
	    } else {
		ihi=2 ;
		inhi=1 ;
	    }

	    for (i=1;i<=mpts;i++) {
		if (y[i] < y[ilo]) ilo=i;
		if (y[i] > y[ihi]) {
		    inhi=ihi;
		    ihi=i;
		} else if (y[i] > y[inhi])
		    if (i != ihi) inhi=i;
	    }

#define TINY 1.e-16
	    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]+TINY));

	    if(verbose) {

		printf("iter=%d, rtol=%g, ftol=%g\n", iter, rtol, ftol) ;
		for(i = 1 ; i <= ndim+1 ; i++) {
		    printf("   ") ;

		    for(j = 1 ; j <= ndim ; j++)
			printf(" %10f", p[i][j]) ;

		    printf(" = %10g\n", y[i]) ;
		}
	    }

	    iter++ ;

	    if (rtol < ftol) break;

	    if (n_iter >= max_iter*2) {
		fprintf(stderr, "Iterations limit reached in AMOEBA\n");
		break ;
	    }

	    float a = ALPHA ;

	    // don't contract too fast if early in search
	    if(! (iter % ndim) && iter < ndim*GRID_ITER_LIMIT) {
		a *= 10.f ;
	    }

	    // Reflection of simplex
	    ytry=amotry(p,y,psum,ndim,function,ihi,-a);

	    if (ytry <= y[ilo]) {

		float g = GAMMA ;

		// don't contract too fast if early in search
		if(! (iter % ndim) && iter < ndim*GRID_ITER_LIMIT) {
		    g *= 10.f ;
		}

		// Extrapolation of simplex
		ytry=amotry(p,y,psum,ndim,function,ihi,g);

	    } else if (ytry >= y[inhi]) {
		// reflected point worse than 2nd highest

		ysave=y[ihi];
		ytry=amotry(p,y,psum,ndim,function,ihi,BETA);

		if (ytry >= ysave) {
		    for (i=1;i<=mpts;i++) {
			if (i != ilo) {
			    for (j=1;j<=ndim;j++) {
				psum[j]=0.5*(p[i][j]+p[ilo][j]);
				p[i][j]=psum[j];
			    }
			    y[i]=EvaluateVertex(function,psum,ndim);
			}
		    }
		    GET_PSUM
		}
	    }
	}

    free(psum) ;

    return ilo ;
}

void print_parms(float *p, float y, int n)
{
    int i ;
    for(i = 1 ; i <= n ; i++)
	printf(" %10f", p[i]) ;

    printf(" = %10g\n", y) ;
}

float AmoebaFit(float (*function)(), float tol, int itmax, int ndim, float guesses[] ,float guess_delta[], float c_lo[], float c_hi[], int verbose)
{
    float **p ;  // holds all the vertices of the simplex
    float *y ;  // holds function result of parameters at a simplex point
    int nvert = ndim+1 ;
    int i, j ;

    if(verbose) printf("Initial guesses\n") ;

    /* put constraints in global array */
    for (j=1;j<=ndim;j++) {
	constraints_lo[j] = c_lo[j] ;
	constraints_hi[j] = c_hi[j] ;
	if(verbose) printf("%2d:%f + %f (range %f to %f)\n", j, guesses[j], guess_delta[j], c_lo[j], c_hi[j]) ;
    }

// the +1 is for the indexing starting at 1 ...
    p = (float **)calloc(nvert+1, sizeof(float *)) ;
    y = (float *)calloc(nvert+1, sizeof(float)) ;

    for (j=1;j<=ndim+1;j++) {
	p[j] = (float *)calloc(ndim+1, sizeof(float)) ;
    }

    // first vertex is just the guesses
    for (j=1;j<=ndim;j++) {
	p[1][j] = guesses[j] ;
    }


    if(verbose) printf("=======  INITIAL SIMPLEX ==============\n") ;
    y[1]=(*function)(p[1]);
    if(verbose) print_parms(p[1],y[1],ndim) ;

    // now we don't need the +1 on "i" in the loop
    for (i=1;i<=ndim;i++) {

	for (j=1;j<=ndim;j++) {
	    p[i+1][j] = guesses[j] ;
	}

	p[i+1][i] = guesses[i]+guess_delta[i] ;
	y[i+1]=(*function)(p[i+1]);

	if(verbose) print_parms(p[i+1],y[i+1],ndim) ;
    }
    if(verbose) printf("=======================================\n") ;

    int ilo = amoeba(p,y,ndim,tol,itmax,function,verbose) ;

    for (j=1;j<=ndim;j++) {
	guesses[j] = p[ilo][j] ;
    }


    for (j=1;j<=ndim+1;j++) {
	free(p[j]) ;
    }

    free(p) ;
    float best = y[ilo] ;

    free(y) ;

    return best ;
}

#define MAIN_NOT
#ifdef MAIN

float test_func(float *p)
{
    float d1 = p[1]-7 ;
    float d2 = p[2]+3 ;
    float feval = (d1*d1 + d2*d2) ;

    //printf("Test func %f, %f = %f\n", p[1], p[2], feval) ;
    if(n_iter > 200) exit(1) ;
    return(d1*d1 + d2*d2) ;
}
int main(int argc, char **argv)
{
    float guess[3] = {0,1.,1.} ;
    float guess_delta[3] = {0,1.,1.} ;
#define FMAX 1.e30f
    float c_lo[3] = {-FMAX,-FMAX,-FMAX} ;
    float c_hi[3] = { FMAX, FMAX, FMAX} ;

    float feval = AmoebaFit(test_func,1.e-7,100,2,guess,guess_delta,c_lo,c_hi) ;

    printf("%f, %f = %f\n", guess[1], guess[2], feval) ;
}
#endif

#undef ALPHA
#undef BETA
#undef GAMMA
#undef NMAX
