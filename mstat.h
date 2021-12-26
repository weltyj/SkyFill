
#define MSTAT_MAX_N 7
#define MSTAT_ROW_N (MSTAT_MAX_N+2)

struct mstat {
    double coef[MSTAT_MAX_N+1][MSTAT_ROW_N] ;
    int N ;
} ;

void estimate_reg(struct mstat *, double *) ;
struct mstat init_reg(int) ;
void sum_reg(struct mstat *, double x[], double y, double weight) ;
void print_reg_matrix(struct mstat *m, char *name) ;
