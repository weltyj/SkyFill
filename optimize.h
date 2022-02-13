#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#define OPTIMIZE_V_from_HSV 0x01
#define OPTIMIZE_FULL_RGB 0x02

void set_xy_constraints(struct OPT_PARM opt_parms[]) ;
void set_cie_parms(int CIE_index) ;
float optimize_grid_function(int CIE_index) ;
float optimize_horizon_function(int CIE_index, float B) ;
float find_horizon_py(int n_samples, SKYFILL_DATA_t *pData,struct OPT_PARM opt_parms[]) ;
float optimize_grid(int n_samples, int verbose, int allowed_sky_type, int CIE_sky_index,struct OPT_PARM opt_parms[]) ;

float smart_optimize_function(float *pParams) ;
void smart_optimize(int n_samples, SKYFILL_DATA_t *pData,struct OPT_PARM opt_parms[]) ;

void optimize_sky_model(int force_grid_optimization, int n_custom_optimizations, int optimization_var[10][MAX_OPT_PARMS+1], int n_samples, SKYFILL_DATA_t *pData,struct OPT_PARM opt_parms[]) ;
#endif
