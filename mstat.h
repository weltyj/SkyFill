#ifndef MSTAT_HEADER
#define MSTAT_HEADER

#define MSTAT_MAX_P 10
#define MSTAT_ROW_P (MSTAT_MAX_P+2)

struct mstat {
    double coef[MSTAT_MAX_P+1][MSTAT_ROW_P] ;
    double B[MSTAT_MAX_P] ;
    double standard_error[MSTAT_MAX_P] ;
    double sum_wgts ;
    double sum_y ;
    double sum_y_sq ;
    double mean_y ;
    double ssto ;
    double sse ;
    double ssr ;
    double mse ;
    double msr ;
    double sigmasq ;
    double rsq ;
    int P ; // number if independent variables, including intercept
    int nobs ;
} ;

struct mstat_75 {
    // for Gentleman's AS 75 algorithm
    int P ;
    int nobs ;
    double sse ;
    double syy ;
    double sy ;
    double wghtn ;
    double d[MSTAT_MAX_P+1] ;
    double thetabar[MSTAT_MAX_P+1] ;
    int rbarsize ;
    double rbar[MSTAT_MAX_P+1] ;
    double theta[MSTAT_MAX_P+1] ;
    double sigmasq ;
    double standard_error[MSTAT_MAX_P+1] ;
    double rsq ;
    double adjrsq ;
    double ybar  ;
    double sst ;
} ;

void estimate_reg(struct mstat *, double *) ;
void estimate_reg_cholesky(struct mstat *, double *) ;
struct mstat init_reg(int) ;
void sum_reg(struct mstat *, double x[], double y, double weight) ;

struct mstat_75 init_reg_75(int) ;
int sum_reg_75(struct mstat_75 *, double x[], double y, double weight) ;
void estimate_reg_75(struct mstat_75 *, double *) ;

void sum_reg_int(struct mstat *, double x[], double y) ;
void print_reg_matrix(struct mstat *m, char *name) ;
double stddev_calc(int n, double sum, double sumsq) ;
double cv_calc(int n, double sum, double sumsq) ;
float  means_t_test(int n1, float sum1,  float sumsq1, int n2, float sum2, float sumsq2, int verbose) ;
double Student_t_probability(double t, int df) ;
void cholesky(double mat[][MSTAT_ROW_P], short N) ;
#endif
