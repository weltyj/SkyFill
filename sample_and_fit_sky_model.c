#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "tiffio.h"

#include "skyfill_tif.h"
#include "mstat.h"
#include "colorspace_conversions.h"
#include "pixel_tests.h"
#include "sample_and_fit_sky_model.h"
#include "amoeba_06.h"
#include "mpfit.h"

#define USE_AMOEBA 0

void output_model_results(int n_samples, SKYFILL_DATA_t *pData) ;
void initialize_sample(int n_samples, uint16_t x, uint16_t y, float px, float py, uint16_t *rgb, float *hsv) ;
void find_sky_outliers(int n_samples) ;
void output_normalized_sample_data(int n_normalized_samples, char filename[]) ;

struct sample_point *samples = NULL ;
struct sample_point *normalized_samples = NULL ; // has value from HSV normalized by px location
int n_samples_to_optimize = 0 ;
int n_normalized_samples = 0 ;
float *val_base_knots ;
float *sat_base_knots ;

int g_full_hsv_model_level ;
int g_location_index ;
int g_calc_sample_error=0 ;
int g_sat_model_full=1 ;
int g_val_model_full=1 ;

void reset_sample_errors(int n_samples)
{
    for(int i = 0 ; i < n_samples ; i++) {
	for(int l = 0 ; l < 2 ; l++) {
	    samples[i].abs_error[l] = 1. ;
	    samples[i].abs_h_error[l] = 1. ;
	    samples[i].abs_s_error[l] = 1. ;
	    samples[i].abs_v_error[l] = 1. ;
	}
    }
}


// coefs to model sky HSV from input file sky area
// model is H|S|V = B0 + B1*px + B2*py + B3*px*py ;
struct HSV_model_coefs raw_hsv_coefs[2] ;

float sat_optimize_function_reduced(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;
    float coefs_all[7] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = 0. ;
    coefs_all[4] = 0. ;
    coefs_all[5] = 0. ;

    float sse=0. ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;

    wgt_func = px_wgt_single ;

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float s_hat = local_s_from_pxpy(samples[i].px, py_hc, coefs_all) ;

	float err = samples[i].s - s_hat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
/*  	if(i < 5) printf("NMERR WGT %d %f %f %f\n", i, samples[i].s, err, wgt) ;  */

	if(samples[i].is_outlier == 0) {
	    sse += err*err*wgt ;
	}

	if(g_calc_sample_error) {
	    samples[i].abs_s_error[g_location_index] = fabsf(err) ;
/*  	    if(fabsf(err) > 0.1) {  */
/*  		samples[i].is_outlier |= SAT_OUTLIER ;  */
/*  	    }  */
	}
    }

    return sse ;
}

float lm_sat_func(float *coefs, double *dy, void *vars, int wgt_flag_single)
{
    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;

    if(g_full_hsv_model_level == 1 || wgt_flag_single == 1) {
	wgt_func = px_wgt_single ;
    } else if(g_location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }

    float sse=0. ;

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float s_hat = local_s_from_pxpy(samples[i].px, py_hc, coefs) ;

	float err = samples[i].s - s_hat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
/*  	if(i < 5) printf("LMERR WGT %d %f %f %f\n", i, samples[i].s, err, wgt) ;  */

	if(samples[i].is_outlier == 0) {
	    sse += err*err*wgt ;
	    dy[i] = -err*wgt ; // for levenberg marquardt divide the error by the wgt
	} else {
	    wgt *= .001 ;
	    dy[i] = -err*wgt ; // for levenberg marquardt divide the error by the wgt
	}


	if(g_calc_sample_error) {
	    samples[i].abs_s_error[g_location_index] = fabsf(err) ;
	    if(fabsf(err) > 0.1) {
#ifdef SATURDAY
		samples[i].is_outlier |= SAT_OUTLIER ;
#endif
	    }
	}
    }

    return sse ;
}

int lm_sat_func_full(int m, int n, double *coefs, double *dy, double **dvec, void *vars)
{
    int i;

#ifdef QUADRATIC_SAT_VERSION
// this is actually a misnomer, it has a coefficient for a quadratic term of px, but it is set to 0.
    float coefs_all[6] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = coefs[3] ;
    coefs_all[4] = 0. ;
    coefs_all[5] = coefs[4] ;
#else
    float coefs_all[5] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = coefs[3] ;
    coefs_all[4] = coefs[4] ;
#endif

    //lm_sat_func(coefs_all, dy, vars, 0) ;
    float sse = lm_sat_func(coefs_all, dy, vars, 0) ;
/*      printf("LMSFF: %f %f %f %f %f = %f\n", coefs[0], coefs[1], coefs[2], coefs[3], coefs[4], sse) ;  */

    if(isnan(coefs[0])) exit(1) ;
    return 0 ;
}

int lm_sat_func_reduced(int m, int n, double *coefs, double *dy, double **dvec, void *vars)
{
    int i;

    float coefs_all[6] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = 0. ;
    coefs_all[4] = 0. ;
    coefs_all[5] = 0. ;

    float coefs_one_based[7] ;
    coefs_one_based[1] = coefs[0] ;
    coefs_one_based[2] = coefs[1] ;
    coefs_one_based[3] = coefs[2] ;
    coefs_one_based[4] = 0. ;
    coefs_one_based[5] = 0. ;
    coefs_one_based[6] = 0. ;

    float sse_old = sat_optimize_function_reduced(coefs_one_based) ;
    float sse = lm_sat_func(coefs_all, dy, vars, 1) ;
/*      printf("LMSFR: %f %f %f = %f old=%f\n", coefs[0], coefs[1], coefs[2], sse, sse_old) ;  */

    if(isnan(coefs[0])) exit(1) ;

    return 0 ;
}

float sat_optimize_function(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;

#ifdef QUADRATIC_SAT_VERSION
// this is actually a misnomer, it has a coefficient for a quadratic term of px, but it is set to 0.
    float coefs_all[6] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = coefs[3] ;
    coefs_all[4] = 0. ;
    coefs_all[5] = coefs[4] ;
#else
    float coefs_all[5] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = coefs[3] ;
    coefs_all[4] = coefs[4] ;
#endif


    float sse=0. ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;

    if(g_full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(g_location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float s_hat = local_s_from_pxpy(samples[i].px, py_hc, coefs_all) ;

	float err = samples[i].s - s_hat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;

	if(samples[i].is_outlier == 0) {
	    sse += err*err*wgt ;
	}

	if(g_calc_sample_error) {
	    samples[i].abs_s_error[g_location_index] = fabsf(err) ;
	    if(fabsf(err) > 0.1) {
#ifdef SATURDAY
		samples[i].is_outlier |= SAT_OUTLIER ;
#endif
	    }
	}
    }

    return sse ;
}

void lm_get_saturation_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift, float best_horizon, float best_damp, float best_rate)
{

/*      py - hpy - c[0] > 0  */
/*      -c[0] > hpy - py  */
/*       c[0] < lowest_sky_py - hpy  */


     float upper_0 = pData->lowest_sky_py_found - pData->horizon_py ;
     upper_0 = 1. ;
     if(best_horizon > upper_0)
	best_horizon = upper_0 - .1 ;
    fprintf(stderr, "============ Leveberg Marquardt SAT START: lowest sky py:%f, horizon py:%f upper_0:%f\n", pData->lowest_sky_py_found, pData->horizon_py, upper_0) ;

    fprintf(stderr, "============ Leveberg Marquardt SAT START: location %d SAT horizon:%f, damp:%f rate:%f\n", location_index, best_horizon, best_damp, best_rate) ;
#define NSP 5
    // this little struct just makes it easier to see the parameters and their constraints
    struct valparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
    } parms[NSP] = {
	{best_horizon, .05, 0., upper_0}, // "horizon"
	{best_damp, .02, .001, 10.},  // 
	{best_rate, .2, .001, 10.},  // rate
	{.05, .01, -.3, .3},  // rate px
	{.1, .05, -.3, 1},  // intercept
	/*  		     {-.001, .0001, -.03, .03},  // rate px*px  */
    } ;

/*      float guess[NSP+1] ;  */

/*      float guess_delta[NSP+1] ;  */
/*    */
/*      for(int i=1 ; i <= NSP ; i++) {  */
/*  	guess[i] = parms[i-1].guess ;  */
/*  	guess_delta[i] = parms[i-1].guess_delta ;  */
/*      }  */

    double guess[NSP]  ;           /* Parameter initial conditions */
    double perror[NSP];                   /* Returned parameter errors */      
    int i;
/*      struct vars_struct v;  */
    int status;
    mp_result result;

    memset(&result,0,sizeof(result));       /* Zero results structure */
    result.xerror = perror;
/*      v.x = x;  */
/*      v.y = y;  */
/*      v.ey = ey;  */

    mp_par constraints[5] ;

    for(int i=0 ; i < 5 ; i++) {
	guess[i] = parms[i].guess ;
	constraints[i].fixed=0 ;
	constraints[i].limited[0]=1 ;
	constraints[i].limited[1]=1 ;
	constraints[i].limits[0]=parms[i].lo ;
	constraints[i].limits[1]=parms[i].hi ;
	constraints[i].parname=0 ;
	constraints[i].step=parms[i].guess_delta ;
	constraints[i].relstep=0.02 ;
	constraints[i].side=0 ;
	constraints[i].deriv_debug=0 ;
	constraints[i].deriv_reltol=0. ;
	constraints[i].deriv_abstol=0. ;
    }

    mp_config config ;
    memset(&config,0,sizeof(config));       /* Zero results structure */
    /*    config.stepfactor = 1.5 ;  */
    /*    config.epsfcn = .05 ;  */
    config.ftol = 1.e-4 ;
    config.xtol = 1.e-4 ;
    config.gtol = 1.e-4 ;


    // get overall horizon, rate and speed parameters
    g_calc_sample_error=0 ;
    status = mpfit(lm_sat_func_reduced, n_samples_to_optimize, 3, guess, constraints, &config, 0, &result);
    printf("\n*** Levenberg Marquardt SAT reduced status = %d****\n", status);

/*      guess_delta[1] = guess[1]/10. ;  */
    constraints[0].step = guess[0]/10. ;

/*      guess_delta[2] = guess[2]/10. ;  */
/*      c_lo[2] = guess[2] - guess[2]/5. ;  */
/*      c_hi[2] = guess[2] + guess[2]/5. ;  */

    constraints[1].step = guess[1]/10. ;
    constraints[1].limits[0] = guess[1] - guess[1]/5. ;
    constraints[1].limits[1] = guess[1] + guess[1]/5. ;

/*      guess_delta[3] = guess[3]/10. ;  */
/*      c_lo[3] = guess[3] - guess[3]/5. ;  */
/*      c_hi[3] = guess[3] + guess[3]/5. ;  */

    constraints[2].step = guess[2]/10. ;
    constraints[2].limits[0] = guess[2] - guess[2]/5. ;
    constraints[2].limits[1] = guess[2] + guess[2]/5. ;

    status = mpfit(lm_sat_func_full, n_samples_to_optimize, NSP, guess, constraints, &config, 0, &result);
    printf("\n*** Levenberg Marquardt SAT full model status = %d****\n", status);


    // need errors for final analysis in horizon curvature iteration
    double *dy = (double *)calloc(n_samples_to_optimize,sizeof(double)) ;
    g_calc_sample_error=1 ;
    lm_sat_func_full(NSP,n_samples_to_optimize,guess,dy,NULL,NULL) ;

    g_calc_sample_error=0 ;
    status = mpfit(lm_sat_func_full, n_samples_to_optimize, NSP, guess, constraints, &config, 0, &result);
    printf("\n*** Levenberg Marquardt SAT full model status = %d****\n", status);

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    lm_sat_func_full(NSP,n_samples_to_optimize,guess,dy,NULL,NULL) ;

    free(dy) ;

#ifdef QUADRATIC_SAT_VERSION
    coefs[0] = guess[0] ;
    coefs[1] = guess[1] ;
    coefs[2] = guess[2] ;
    coefs[3] = guess[3] ;
    coefs[4] = 0. ;
    coefs[5] = guess[4] ;
#else
    coefs[0] = guess[0] ;
    coefs[1] = guess[1] ;
    coefs[2] = guess[2] ;
    coefs[3] = guess[3] ;
    coefs[4] = guess[4] ;
#endif
    printf("SAT coefs %f %f %f %f %f\n", coefs[0], coefs[1], coefs[2], coefs[3], coefs[4]) ;
}

void get_saturation_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift)
{
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;

    float min_sat=1. ; // need to look at samples
    float max_sat=0. ; // need to look at samples
    float min_py=1.5 ;  // need to look at samples for min py_hc
    float max_py=0. ;  // need to look at samples for min py_hc

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	if(samples[i].is_outlier == 0) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	    if(max_py < py_hc) max_py = py_hc ;
	    if(min_py > py_hc) min_py = py_hc ;
	    if(min_sat > samples[i].s) min_sat= samples[i].s ;
	    if(max_sat < samples[i].s) max_sat= samples[i].s ;
	}
    }


    float (*wgt_func)(float,float) ;

    if(g_full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(g_location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }

    float best_horizon=0.1 ;
    float best_damp=0.1 ;
    float best_rate=0.1 ;
    float best_sse=1.e30 ;
    for(float horizon=0.5 ; horizon < 1.0 ; horizon+= .05) {
	struct mstat m_fit ;
	double py_coefs[2] ;
	float sse = 0. ;

	m_fit = init_reg(2) ;

	for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	    if(samples[i].is_outlier == 0) {
		float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
		double xv[2] = {1.,  1./(py_hc-horizon)} ;
		double y = 1./samples[i].s ;
		double w = 1./y ;
		sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, w*w*w)) ;
	    }
	}

	estimate_reg(&m_fit, py_coefs) ;
	for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	    if(samples[i].is_outlier == 0) {
		float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
		float s_hat = (py_hc-horizon)/(py_coefs[0]+py_coefs[1]*(py_hc-horizon)) ;
		float err = s_hat - samples[i].s ;
		sse += err*err ;
	    }
	}

	if(sse < best_sse) {
	    best_damp = py_coefs[1] ;
	    best_rate = py_coefs[0] ;
	    best_horizon = horizon ;
	    best_sse = sse ;
	}
    }

    if(best_damp < .05) best_damp=.05 ;

    lm_get_saturation_model(n_samples, coefs, pData, location_index, phase_shift, best_horizon, best_damp, best_rate) ;
    return ;

    fprintf(stderr, "============ SAT START: SAT horizon:%f, damp:%f rate:%f\n", best_horizon, best_damp, best_rate) ;
#define NSP 5
    // this little struct just makes it easier to see the parameters and their constraints
    struct valparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
    } parms[NSP] = {
	{best_horizon, .05, 0., .75}, // "horizon"
	{best_damp, .02, .001, 10.},  // 
	{best_rate, .2, .001, 10.},  // rate
	{.05, .01, -.3, .3},  // rate px
	{.1, .05, -.3, 1},  // intercept
	/*  		     {-.001, .0001, -.03, .03},  // rate px*px  */
    } ;

    float guess[NSP+1] ;
    float guess_delta[NSP+1] ;
    float c_lo[NSP+1] ;
    float c_hi[NSP+1] ;

    for(int i=1 ; i <= NSP ; i++) {
	guess[i] = parms[i-1].guess ;
	guess_delta[i] = parms[i-1].guess_delta ;
	c_lo[i] = parms[i-1].lo ;
	c_hi[i] = parms[i-1].hi ;
    }

    g_calc_sample_error=0 ;
    // get overall horizon, rate and speed parameters
    AmoebaFit(sat_optimize_function_reduced,1.e-4,200,3,guess,guess_delta,c_lo,c_hi,0) ;


    guess_delta[1] = guess[1]/10. ;

    guess_delta[2] = guess[2]/10. ;
    c_lo[2] = guess[2] - guess[2]/5. ;
    c_hi[2] = guess[2] + guess[2]/5. ;

    guess_delta[3] = guess[3]/10. ;
    c_lo[3] = guess[3] - guess[3]/5. ;
    c_hi[3] = guess[3] + guess[3]/5. ;

    AmoebaFit(sat_optimize_function,1.e-4,200,NSP,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    sat_optimize_function(guess) ;

    g_calc_sample_error=0 ;
    AmoebaFit(sat_optimize_function,1.e-4,200,NSP,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    sat_optimize_function(guess) ;

#ifdef QUADRATIC_SAT_VERSION
    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;
    coefs[4] = 0. ;
    coefs[5] = guess[5] ;
#else
    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;
    coefs[4] = guess[5] ;
#endif

}

int lm_hue_func(int m, int n, double *coefs, double *dy, double **dvec, void *vars)
{
    int i;
    float hue_coefs[4] ;

    hue_coefs[0] = coefs[0] ;
    hue_coefs[1] = coefs[1] ;
    hue_coefs[2] = coefs[2] ;
    hue_coefs[3] = coefs[3] ;
    hue_coefs[4] = coefs[4] ;

    float sse=0. ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;
    float sum_v_hat=0. ;

    wgt_func = px_wgt_single ;

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float h_hat = local_h_from_pxpy(samples[i].px, py_hc, hue_coefs) ;

	float err = samples[i].h - h_hat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;

	if(samples[i].is_outlier != 0) {
	    wgt *= .001 ;
	}

	sse += err*err*wgt ;
	dy[i] = -err*wgt ; // for levenberg marquardt divide the error by the wgt

	if(g_calc_sample_error) {
	    samples[i].abs_h_error[g_location_index] = fabsf(err) ;
	    if(fabsf(err) > 40.) {
#ifdef SATURDAY
		samples[i].is_outlier |= HUE_OUTLIER ;
#endif
	    }
	}
    }

/*      fprintf(stderr, "LM_VOF: %7.4f %7.4f %7.4f %7.4f = %f\n",  */
/*  	val_coefs[0],  */
/*  	val_coefs[1],  */
/*  	val_coefs[2],  */
/*  	val_coefs[3],  */
/*  	sse) ;  */
    return 0;
}

float hue_optimize_function(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;
    float sse=0. ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;

    if(g_full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(g_location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float h_hat = local_h_from_pxpy(samples[i].px, py_hc, coefs) ;

	float err = samples[i].h - h_hat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
	if(samples[i].is_outlier == 0) {
	    sse += err*err*wgt ;
	}
	if(g_calc_sample_error) {
	    samples[i].abs_h_error[g_location_index] = fabsf(err) ;
	    if(fabsf(err) > 40.) {
#ifdef SATURDAY
		samples[i].is_outlier |= HUE_OUTLIER ;
#endif
	    }
	}
    }

    return sse ;
}

void get_hue_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift)
{
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;

    float min_hue=360. ; // need to look at samples
    float max_hue=0. ; // need to look at samples
    float min_py=1.5 ;  // need to look at samples for min py_hc

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	if(samples[i].is_outlier == 0) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	    if(min_py > py_hc) min_py = py_hc ;
	    if(min_hue > samples[i].h) min_hue= samples[i].h ;
	    if(max_hue < samples[i].h) max_hue= samples[i].h ;
	}
    }


    float hue_spread= max_hue-min_hue ;
    if(hue_spread > 30) hue_spread=30 ;

    // this little struct just makes it easier to see the parameters and their constraints
    struct hueparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
	} parms[7] = {
		     {min_hue, 1., min_hue/2., min_hue+20},
		     {hue_spread, .1, hue_spread-20., hue_spread+20},
		     {40., .01, 5., 80.},
		     {min_py, .05, min_py-.1, 1.},
		     {.05, .01, -.2, .2},   // now center is coefs[3]+coefs[4]*px
		     } ;

    if(USE_AMOEBA) {
	float guess[8] ;
	float guess_delta[8] ;
	float c_lo[8] ;
	float c_hi[8] ;

	for(int i=1 ; i <= 4 ; i++) {
	    guess[i] = parms[i-1].guess ;
	    guess_delta[i] = parms[i-1].guess_delta ;
	    c_lo[i] = parms[i-1].lo ;
	    c_hi[i] = parms[i-1].hi ;
	}
	g_calc_sample_error=0 ;
	AmoebaFit(hue_optimize_function,1.e-4,200,4,guess,guess_delta,c_lo,c_hi,0) ;

	// need errors for final analysis in horizon curvature iteration
	g_calc_sample_error=1 ;
	hue_optimize_function(guess) ;

	g_calc_sample_error=0 ;
	AmoebaFit(hue_optimize_function,1.e-4,200,4,guess,guess_delta,c_lo,c_hi,0) ;

	// need errors for final analysis in horizon curvature iteration
	g_calc_sample_error=1 ;
	hue_optimize_function(guess) ;

	coefs[0] = guess[1] ;
	coefs[1] = guess[2] ;
	coefs[2] = guess[3] ;
	coefs[3] = guess[4] ;
    } else {
	// Levenberg Marquardt

	double guess[5]  ;           /* Parameter initial conditions */
	double perror[5];                   /* Returned parameter errors */      
	int i;
    /*      struct vars_struct v;  */
	int status;
	mp_result result;

	memset(&result,0,sizeof(result));       /* Zero results structure */
	result.xerror = perror;
    /*      v.x = x;  */
    /*      v.y = y;  */
    /*      v.ey = ey;  */

	mp_par constraints[5] ;

	for(int i=0 ; i < 5 ; i++) {
	    guess[i] = parms[i].guess ;
	    constraints[i].fixed=0 ;
	    constraints[i].limited[0]=1 ;
	    constraints[i].limited[1]=1 ;
	    constraints[i].limits[0]=parms[i].lo ;
	    constraints[i].limits[1]=parms[i].hi ;
	    constraints[i].parname=0 ;
	    constraints[i].step=parms[i].guess_delta ;
	    constraints[i].relstep=0.05 ;
	    constraints[i].side=0 ;
	    constraints[i].deriv_debug=0 ;
	    constraints[i].deriv_reltol=0. ;
	    constraints[i].deriv_abstol=0. ;
	}

	mp_config config ;
	memset(&config,0,sizeof(config));       /* Zero results structure */
	/*    config.stepfactor = 1.5 ;  */
	/*    config.epsfcn = .05 ;  */
	config.ftol = 1.e-4 ;
	config.xtol = 1.e-4 ;
	config.gtol = 1.e-4 ;


	/* Call fitting function for 10 data points and 2 parameters */
    /*      status = mpfit(lm_val_func, n_samples_to_optimize, 4, initial_guess, constraints, &config, (void *) &v, &result);  */
	status = mpfit(lm_hue_func, n_samples_to_optimize, 4, guess, constraints, &config, 0, &result);
	printf("\n*** Levenberg Marquardt HUE status = %d****\n", status);

/*  	printf("     NITER = %d\n", result.niter);  */
/*  	printf("      NFEV = %d\n", result.nfev);  */
/*  	printf("\n");  */
/*  	for (i=0; i < 5; i++) {  */
/*  	    printf("  P[%d] = %f +/- %f\n", i, guess[i], result.xerror[i]);  */
/*  	}  */

	coefs[0] = guess[0] ;
	coefs[1] = guess[1] ;
	coefs[2] = guess[2] ;
	coefs[3] = guess[3] ;
	printf("HUE coefs %f %f %f %f\n", coefs[0], coefs[1], coefs[2], coefs[3]) ;
    }
}

// Allow for fitting in separate steps by optimizing reduced models
float val_coefs[9] ;
int val_fit_mode = -1 ;
int completely_model_val=0 ;
int n_val_func_calls=0 ;

#define VAL_MODEL_IS_SINGLE

float val_optimize_function(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;
    n_val_func_calls++ ;

    val_coefs[0] = coefs[0] ;
    val_coefs[1] = coefs[1] ;
    val_coefs[2] = coefs[2] ;
    val_coefs[3] = coefs[3] ;

    float sse=0. ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;
    float sum_v_hat=0. ;

#ifdef VAL_MODEL_IS_SINGLE
    wgt_func = px_wgt_single ;
#else
    if(g_full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(g_location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }
#endif

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float v_hat ;

	v_hat = local_v_from_pxpy(samples[i].px, py_hc, val_coefs) ;
/*  	fprintf(stderr, "vh:%f px:%f py_hc:%f\n", v_hat, samples[i].px, py_hc) ;  */
	sum_v_hat += v_hat ;

	float err = v_hat - samples[i].v ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
	wgt=1. ;

	if(1 || samples[i].is_outlier == 0) {
	    sse += err*err*wgt ;
	} else {
	    printf("outlier at v:%f angle%f\n", samples[i].v, samples[i].px*pData_fit->FOV_horizontal) ;
	}

	if(g_calc_sample_error) {
	    samples[i].abs_v_error[g_location_index] = fabsf(err) ;
	}
    }

/*      fprintf(stderr, "VOF: %7.4f %7.4f %7.4f hpy:%7.4f = %f s_v_hat:%f\n",  */
/*  	val_coefs[0],  */
/*  	val_coefs[1],  */
/*  	val_coefs[2],  */
/*  	pData_fit->horizon_py,  */
/*  	sse, sum_v_hat) ;  */

    return sse ;
}

int lm_val_cie_func(int m, int n, double *coefs, double *dy, double **dvec, void *vars)
{
    int i;

    n_val_func_calls++ ;

    pData_fit->perez_A = coefs[0] ;
    pData_fit->perez_B = coefs[1] ;
    pData_fit->perez_C = coefs[2] ;
    pData_fit->perez_D = coefs[3] ;
    pData_fit->perez_E = coefs[4] ;
    pData_fit->horizon_py = coefs[5] ;
    pData_fit->sun_x = coefs[6] ;
    pData_fit->sun_py = coefs[7] ;
    set_FOV_factor() ;
    V_sun= image_relative_pixel_to_V3(sun_x_angle2px(pData_fit->sun_x), pData_fit->sun_py) ;

    float sse=0. ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;

#ifdef VAL_MODEL_IS_SINGLE
    wgt_func = px_wgt_single ;
#else
    if(g_full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(g_location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }
#endif

    float sum_v=0 ;
    float sum_v_hat=0 ;

    // recalculate value correction using sky samples and
    // current prediction coefficients
/*      sample_based_v_correction = 1.f ;  */
    float gamma, theta, cos_gamma ;

    for(i = 0 ; i < n_samples_to_optimize ; i++) {
	//float F_CIE2003(float px, float py, float *pGamma, float *pTheta, float *pCos_gamma)
	float v_hat = F_CIE2003(samples[i].px, samples[i].py, &gamma, &theta, &cos_gamma) ;

	samples[i].v_hat = v_hat ;
	sum_v += samples[i].v ;
	sum_v_hat += samples[i].v_hat ;
    }

    pData_fit->sample_based_v_correction = sum_v/sum_v_hat ;

    for(i = 0 ; i < n_samples_to_optimize ; i++) {
	float v_hat = samples[i].v_hat*pData_fit->sample_based_v_correction ;
	float err = (v_hat - samples[i].v) ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
	wgt=1. ;

	if(samples[i].is_outlier != 0) {
	    wgt *= .001 ;
	    //printf("outlier at v:%f angle%f\n", samples[i].v, samples[i].px*pData_fit->FOV_horizontal) ;
	}

	sse += err*err*wgt ;
	dy[i] = err*wgt ; // for levenberg marquardt divide by the weight

	if(g_calc_sample_error) {
	    samples[i].abs_v_error[g_location_index] = fabsf(err) ;
	}

    }

/*      for (i=0; i < 8; i++) {  */
/*  	printf("  %f", coefs[i]);  */
/*      }  */

/*      printf(" = %f\n", sse) ;  */

    return 0;
}

float get_value_model_CIE(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift, float amplitude_damp)
{
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;
    n_val_func_calls = 0 ;

#define NVP 8

    // this little struct just makes it easier to see the parameters and their constraints
    struct valparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
    } parms[NVP] = {
/*      -1.000000 -0.550000 10.000000 -3.000000 0.450000  */
/*  -1.010000 -0.200000 9.990000 -2.970000 0.410000  */
	{-0.5, .01, -8., 8}, // A
	{-.32, .01, -4, 0.}, // B
	{  0.1, .1, 0., 24.}, // C
	{-1., .1, -6., -.1}, // D
	{.45, .01, 0., 1.7}, // E
	{.5, .01, 0.0, .7}, // horizon_py
	{pData_fit->sun_x, 5, pData_fit->sun_x-20., pData_fit->sun_x+20.},
	{pData_fit->sun_py, .05, pData_fit->sun_py-.25, pData_fit->sun_py+.25}
    } ;

    // Levenberg Marquardt

    double initial_guess[NVP]  ;           /* Parameter initial conditions */
    double perror[NVP];                   /* Returned parameter errors */      
    int i;
/*      struct vars_struct v;  */
    int status;
    mp_result result;

    memset(&result,0,sizeof(result));       /* Zero results structure */
    result.xerror = perror;
/*      v.x = x;  */
/*      v.y = y;  */
/*      v.ey = ey;  */

    mp_par constraints[NVP] ;

    for(int i=0 ; i < NVP ; i++) {
	initial_guess[i] = parms[i].guess ;
	constraints[i].fixed=0 ;
	constraints[i].limited[0]=1 ;
	constraints[i].limited[1]=1 ;
	constraints[i].limits[0]=parms[i].lo ;
	constraints[i].limits[1]=parms[i].hi ;
	constraints[i].parname=0 ;
	constraints[i].step=parms[i].guess_delta ;
	constraints[i].relstep=0.05 ;
	constraints[i].side=0 ;
	constraints[i].deriv_debug=0 ;
	constraints[i].deriv_reltol=0. ;
	constraints[i].deriv_abstol=0. ;
    }

    mp_config config ;
    memset(&config,0,sizeof(config));       /* Zero results structure */
    /*    config.stepfactor = 1.5 ;  */
    /*    config.epsfcn = .05 ;  */
    config.ftol = 1.e-4 ;
    config.xtol = 1.e-4 ;
    config.gtol = 1.e-4 ;


    status = mpfit(lm_val_cie_func, n_samples_to_optimize, NVP, initial_guess, constraints, &config, 0, &result);

    printf("\n*** Levenberg Marquardt status = %d****\n", status);

    printf("     NITER = %d\n", result.niter);
    printf("      NFEV = %d\n", result.nfev);
    printf("\n");
    for (i=0; i < NVP; i++) {
	printf("  P[%d] = %f +/- %f\n", i, initial_guess[i], result.xerror[i]);
    }
    pData_fit->perez_A = initial_guess[0] ;
    pData_fit->perez_B = initial_guess[1] ;
    pData_fit->perez_C = initial_guess[2] ;
    pData_fit->perez_D = initial_guess[3] ;
    pData_fit->perez_E = initial_guess[4] ;
    pData_fit->horizon_py = initial_guess[5] ;
    pData_fit->sun_x = initial_guess[6] ;
    pData_fit->sun_py = initial_guess[7] ;

    return result.bestnorm ; // best chi squared found
}

int lm_val_func(int m, int n, double *coefs, double *dy, double **dvec, void *vars)
{
    int i;
    /*    struct vars_struct *v = (struct vars_struct *) vars;  */
    /*    double *x, *y, *ey, f;  */
    /*    double sse=0. ;  */

    /*    x = v->x;  */
    /*    y = v->y;  */
    /*    ey = v->ey;  */

    /*    for (i=0; i<m; i++) {  */
    /*      f = 1./(1.+exp(param[0] + param[1]*x[i]));  */
    /*      dy[i] = (y[i] - f) ;  */
    /*      sse += dy[i]*dy[i] ;  */
    /*      dy[i] /= ey[i] ;  */
    /*    }  */
    /*    */
    /*    fprintf(stderr, "%14g %14g=%f\n", param[0], param[1], sse) ;  */


    n_val_func_calls++ ;

    val_coefs[0] = coefs[0] ;
    val_coefs[1] = coefs[1] ;
    val_coefs[2] = coefs[2] ;
    val_coefs[3] = coefs[3] ;

    float sse=0. ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;
    float sum_v_hat=0. ;

#ifdef VAL_MODEL_IS_SINGLE
    wgt_func = px_wgt_single ;
#else
    if(g_full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(g_location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }
#endif

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float v_hat ;

	v_hat = local_v_from_pxpy(samples[i].px, py_hc, val_coefs) ;
	/*  	fprintf(stderr, "vh:%f px:%f py_hc:%f\n", v_hat, samples[i].px, py_hc) ;  */
	sum_v_hat += v_hat ;

	float err = v_hat - samples[i].v ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
	wgt=1. ;

	if(1 || samples[i].is_outlier == 0) {
	} else {
	    wgt *= .001 ;
	    printf("outlier at v:%f angle%f\n", samples[i].v, samples[i].px*pData_fit->FOV_horizontal) ;
	}

	sse += err*err*wgt ;
	dy[i] = err*wgt ; // for levenberg marquardt divide by the weight

	if(g_calc_sample_error) {
	    samples[i].abs_v_error[g_location_index] = fabsf(err) ;
	}
    }

/*      fprintf(stderr, "LM_VOF: %7.4f %7.4f %7.4f %7.4f = %f\n",  */
/*  	val_coefs[0],  */
/*  	val_coefs[1],  */
/*  	val_coefs[2],  */
/*  	val_coefs[3],  */
/*  	sse) ;  */
    return 0;
}

float get_value_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift, float amplitude_damp)
{
    // may try to use the original CIE model, now that a few bugs were fixed in the implementation
    // of image coordinates to sky angle calculations

    if(uses_CIE_model & CIE_FULL_MODEL) {
	return get_value_model_CIE(n_samples, coefs, pData, location_index, phase_shift, amplitude_damp) ;
    }

    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;
    n_val_func_calls = 0 ;

    // this little struct just makes it easier to see the parameters and their constraints
    struct valparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
    } parms[4] = {
	{0.2, .01, 0.01, 0.8},
	{0.5, .01, .1, 0.9},
	{5.5, .1, 1., 8.},
	{0.25, .05, .05, 2.},
    } ;

    if(USE_AMOEBA) {

    float guess[5] ;
    float guess_delta[5] ;
    float c_lo[5] ;
    float c_hi[5] ;

    for(int i=1 ; i <= 4 ; i++) {
	guess[i] = parms[i-1].guess ;
	guess_delta[i] = parms[i-1].guess_delta ;
	c_lo[i] = parms[i-1].lo ;
	c_hi[i] = parms[i-1].hi ;
    }
    g_calc_sample_error=0 ;
    AmoebaFit(val_optimize_function,1.e-4,200,4,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    val_optimize_function(guess) ;

    g_calc_sample_error=0 ;
    AmoebaFit(val_optimize_function,1.e-4,200,4,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    float sse = val_optimize_function(guess) ;

    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;

    fprintf(stderr, "VOF: %d func calls, %7.4f %7.4f %7.4f %7.4f hpy,s,c:%7.4f,%7.4f,%7.4f = %f\n",
	    n_val_func_calls,
	    val_coefs[0],
	    val_coefs[1],
	    val_coefs[2],
	    val_coefs[3],
	    pData_fit->horizon_py,
	    pData_fit->horizon_slope,
	    pData_fit->horizon_curvature,
	    sse) ;
    } else {
	// Levenberg Marquardt

	double initial_guess[4]  ;           /* Parameter initial conditions */
	double perror[4];                   /* Returned parameter errors */      
	int i;
    /*      struct vars_struct v;  */
	int status;
	mp_result result;

	memset(&result,0,sizeof(result));       /* Zero results structure */
	result.xerror = perror;
    /*      v.x = x;  */
    /*      v.y = y;  */
    /*      v.ey = ey;  */

	mp_par constraints[4] ;

	for(int i=0 ; i < 4 ; i++) {
	    initial_guess[i] = parms[i].guess ;
	    constraints[i].fixed=0 ;
	    constraints[i].limited[0]=1 ;
	    constraints[i].limited[1]=1 ;
	    constraints[i].limits[0]=parms[i].lo ;
	    constraints[i].limits[1]=parms[i].hi ;
	    constraints[i].parname=0 ;
	    constraints[i].step=parms[i].guess_delta ;
	    constraints[i].relstep=0.05 ;
	    constraints[i].side=0 ;
	    constraints[i].deriv_debug=0 ;
	    constraints[i].deriv_reltol=0. ;
	    constraints[i].deriv_abstol=0. ;
	}

	if(pData->fix_sky_hue) {
	    constraints[1].fixed=1 ;
	    constraints[3].fixed=1 ;
	}

	mp_config config ;
	memset(&config,0,sizeof(config));       /* Zero results structure */
	/*    config.stepfactor = 1.5 ;  */
	/*    config.epsfcn = .05 ;  */
	config.ftol = 1.e-4 ;
	config.xtol = 1.e-4 ;
	config.gtol = 1.e-4 ;


	/* Call fitting function for 10 data points and 2 parameters */
    /*      status = mpfit(lm_val_func, n_samples_to_optimize, 4, initial_guess, constraints, &config, (void *) &v, &result);  */
	status = mpfit(lm_val_func, n_samples_to_optimize, 4, initial_guess, constraints, &config, 0, &result);

/*  	printf("\n*** Levenberg Marquardt status = %d****\n", status);  */
/*    */
/*  	printf("     NITER = %d\n", result.niter);  */
/*  	printf("      NFEV = %d\n", result.nfev);  */
/*  	printf("\n");  */
/*  	for (i=0; i < 4; i++) {  */
/*  	    printf("  P[%d] = %f +/- %f\n", i, initial_guess[i], result.xerror[i]);  */
/*  	}  */

	coefs[0] = initial_guess[0] ;
	coefs[1] = initial_guess[1] ;
	coefs[2] = initial_guess[2] ;
	coefs[3] = initial_guess[3] ;
	return result.bestnorm ; // best chi squared found
    }

/*      exit(1) ;  */
}

#ifdef NEED_PHASE_SHIFT
float g_peak = 0 ;

float phase_shift_optimize_function(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;

    float sse=0. ;

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float x_angle = pData_fit->FOV_horizontal*samples[i].px ;
	float vhat = coefs[0]/(1.+abs(x_angle-g_peak)/180.) + coefs[1] + coefs[2]*py_hc ;

	float err = samples[i].v - vhat ;
	sse += err*err ;
    }

    return sse ;
}

float get_phase_shift(int n_samples, int location_index, int calc_error, SKYFILL_DATA_t *pData, int verbose)
{
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;

    // reduced, constrained to detect best location of peak, 
/*      value = b0/(1.+abs(x_angle-peak)/180.) + b1 + b2*py_hc ;  */
/*      b0 > 0.  */
/*      b2 < 0. ;  */


    float peak, best_peak=0. ;
    float best_sse=1.e30 ;

    for(peak = -180 ; peak < 180. ; peak += 5.) {
	g_peak = peak ;

	struct peakparms {
	    float guess ;
	    float guess_delta ;
	    float lo ;
	    float hi ;
	} parms[3] = {
	    {1,.1,0.1,3.},
	    {.1, .05, -1., 2.},
	    {-.6, .05, -5., 0.},
	} ;

	float guess[8] ;
	float guess_delta[8] ;
	float c_lo[8] ;
	float c_hi[8] ;

	for(int i=1 ; i <= 3 ; i++) {
	    guess[i] = parms[i-1].guess ;
	    guess_delta[i] = parms[i-1].guess_delta ;
	    c_lo[i] = parms[i-1].lo ;
	    c_hi[i] = parms[i-1].hi ;
	}

	g_calc_sample_error=0 ;
	float sse = AmoebaFit(phase_shift_optimize_function,1.e-4,200,3,guess,guess_delta,c_lo,c_hi,0) ;

	if(verbose)
	    fprintf(stderr, "Peak:%4.0f, coef %5.2f %5.2f %5.2f sse=%f\n", peak, guess[1], guess[2], guess[3], sse) ;

	if(sse < best_sse) {
	    best_sse = sse ;
	    best_peak = peak ;
	}
    }

    return -best_peak ;

}
#endif


float get_full_HSV_model_sse(int n_samples, int calc_error, SKYFILL_DATA_t *pData)
{
    //compute combined SSE from each r,g,b component
    float sse_r = 0. ;
    float sse_g = 0. ;
    float sse_b = 0. ;
    float sse_v = 0. ;

    for(int i = 0 ; i < n_samples ; i++) {
	float hsv_hat[3] ;
	predict_sky_hsv(samples[i].px, samples[i].py, hsv_hat) ;

	uint16_t rgb[3] ;
	hsv2rgb16(hsv_hat,rgb) ;

	float err_r = (float)rgb[0] - (float)samples[i].r ;
	float err_g = (float)rgb[1] - (float)samples[i].g ;
	float err_b = (float)rgb[2] - (float)samples[i].b ;
	float err_h = (float)hsv_hat[0] - (float)samples[i].h ;
	float err_s = (float)hsv_hat[1] - (float)samples[i].s ;
	float err_v = (float)hsv_hat[2] - (float)samples[i].v ;

/*  	if(calc_error) {  */
/*  	    float e = (fabsf(err_r) + fabsf(err_g) + fabsf(err_b)) ;  */
/*  	    samples[i].abs_error[location_index] = e ;  */
/*  	    samples[i].abs_h_error[location_index] = fabsf(err_h) ;  */
/*  	    samples[i].abs_s_error[location_index] = fabsf(err_s) ;  */
/*  	    samples[i].abs_v_error[location_index] = fabsf(err_v) ;  */
/*  	}  */

	sse_r += err_r*err_r ;
	sse_g += err_g*err_g ;
	sse_b += err_b*err_b ;

	sse_v += err_v*err_v ;
    }

    return sse_v ;
    return sse_r + sse_g + sse_b ;
}

void fit_sky_model_step2(SKYFILL_DATA_t *pData, int n_samples)
{
    pData->angle_factor = 2.0 ;

    // a little hack -- horizon_py has not effect on the models now, so 
    // by setting horizon_was_set to 1 it is not used in the grid search
    //pData->horizon_was_set=1 ;

    if(! pData->horizon_was_set) {
	pData->horizon_py = pData->lowest_sky_py_found-.05  ;
	if(pData->horizon_py > 1.) pData->horizon_py=1. ;
    }

    float starting_horizon_py = pData->horizon_py ;

    float best_sse = 1.e30;
    float best_horizon_py = pData->horizon_py ;
    pData->horizon_curvature = 0. ;
    pData->horizon_slope = 0. ;
    float best_horizon_slope = pData->horizon_slope ;
    float best_angle_factor = pData->angle_factor ;
    float lr_phase_shift[2] ;

#ifdef NEED_PHASE_SHIFT
    if(pData->valsat_phase_shift > -5000)
	pData->valsat_phase_shift = get_phase_shift(n_samples, 0, 0, pData, 0) ;

    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	lr_phase_shift[lr] = pData->valsat_phase_shift ;
    }
#else
    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	lr_phase_shift[lr] = 0. ;
    }
#endif

    // use only the value model to find the best slope and horizon py
    for(pData->horizon_slope=-.2 ; pData->horizon_slope < .201 ; pData->horizon_slope += .04) {

#ifdef NOT_USED_HERE
	for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	    reset_sample_errors(n_samples) ;

	    float phase_shift = lr_phase_shift[lr] ;

	    get_saturation_model(n_samples, raw_hsv_coefs[lr].s_coefs, pData, lr, phase_shift) ;
	    get_hue_model(n_samples, raw_hsv_coefs[lr].h_coefs, pData, lr, phase_shift) ;

	}
#endif

	for(pData->horizon_py = starting_horizon_py ; pData->horizon_py > .05 ; pData->horizon_py -= .05) {
	    reset_sample_errors(n_samples) ;
	    float sse = 0 ;

	    float phase_shift = lr_phase_shift[0] ;
	    sse = get_value_model(n_samples, raw_hsv_coefs[0].v_coefs, pData, 0, phase_shift, .1) ;

#ifdef VAL_MODEL_IS_SINGLE
	    for(int lr = 1 ; lr < pData->full_hsv_model_level ; lr++) {
		raw_hsv_coefs[lr].v_coefs[0] = raw_hsv_coefs[0].v_coefs[0] ;
		raw_hsv_coefs[lr].v_coefs[1] = raw_hsv_coefs[0].v_coefs[1] ;
		raw_hsv_coefs[lr].v_coefs[2] = raw_hsv_coefs[0].v_coefs[2] ;
		raw_hsv_coefs[lr].v_coefs[3] = raw_hsv_coefs[0].v_coefs[3] ;
		raw_hsv_coefs[lr].v_coefs[4] = raw_hsv_coefs[0].v_coefs[4] ;
	    }
#else
	    for(int lr = 1 ; lr < pData->full_hsv_model_level ; lr++) {
		float phase_shift = lr_phase_shift[lr] ;
		sse += get_value_model(n_samples, raw_hsv_coefs[lr].v_coefs, pData, lr, phase_shift, .1) ;
	    }
#endif

	    if(sse < best_sse) {
		best_sse = sse ;
		best_horizon_py = pData->horizon_py ;
		best_horizon_slope = pData->horizon_slope ;
	    }

	    // full CIE finds best horizon_py and sets it
	    if(uses_CIE_model & CIE_FULL_MODEL)
		break ;
	}
    }

    // final fit with best horizon_py and horizon_slope
    pData->angle_factor = best_angle_factor ;
    pData->horizon_py = best_horizon_py ;
    pData->horizon_slope = best_horizon_slope ;

    float chisq = get_value_model(n_samples, raw_hsv_coefs[0].v_coefs, pData, 0, lr_phase_shift[0], .1) ;

#ifdef VAL_MODEL_IS_SINGLE
    for(int lr = 1 ; lr < pData->full_hsv_model_level ; lr++) {
	raw_hsv_coefs[lr].v_coefs[0] = raw_hsv_coefs[0].v_coefs[0] ;
	raw_hsv_coefs[lr].v_coefs[1] = raw_hsv_coefs[0].v_coefs[1] ;
	raw_hsv_coefs[lr].v_coefs[2] = raw_hsv_coefs[0].v_coefs[2] ;
	raw_hsv_coefs[lr].v_coefs[3] = raw_hsv_coefs[0].v_coefs[3] ;
	raw_hsv_coefs[lr].v_coefs[4] = raw_hsv_coefs[0].v_coefs[4] ;
    }
#else
    for(int lr = 1 ; lr < pData->full_hsv_model_level ; lr++) {
	float phase_shift = lr_phase_shift[lr] ;
	sse += get_value_model(n_samples, raw_hsv_coefs[lr].v_coefs, pData, lr, phase_shift, .1) ;
    }
#endif

    fprintf(stderr, "*******************  CALCULATED PHASE SHIFT %3.0f*************************\n", lr_phase_shift[0]) ;
    fprintf(stderr, "FINAL hslope is %f, hpy is %f, hcrv:%f angle_factor %f, chisq= %f\n",
		    pData->horizon_slope,
		    pData->horizon_py,
		    pData->horizon_curvature,
		    pData->angle_factor, chisq) ;

    fprintf(stderr, "\n") ;

    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	float phase_shift = lr_phase_shift[lr] ;
	get_saturation_model(n_samples, raw_hsv_coefs[lr].s_coefs, pData, lr, phase_shift) ;
	get_hue_model(n_samples, raw_hsv_coefs[lr].h_coefs, pData, lr, phase_shift) ;
    }

    float sse = get_full_HSV_model_sse(n_samples,0,pData) ; // last time through, calculate errors on points
}

#define N_SKY_SECTIONS 4
double coefs_v[N_SKY_SECTIONS][5] ;
double coefs_s[N_SKY_SECTIONS][5] ;

void print_sky_section_info(float px, float py)
{
    float section_width= (float)1./(float)N_SKY_SECTIONS ;
    // interpolation of predicted normalized value from section coefs
    int sample_section = (int)((px+0.5)/section_width) ;

    // this happens at px=0.5
    if(sample_section == N_SKY_SECTIONS)
	sample_section=N_SKY_SECTIONS-1 ;

    float lower = -.5 + (float)sample_section/(float)N_SKY_SECTIONS ;
    float upper = lower + section_width ;

    // sample loc in section is 0.0 at left edge of section, 1.0 at right, 0.5 in middle
    float sample_loc_in_section = (px-lower)/(upper-lower) ;
    int section_left=-1, section_right=-1 ;

    float prop_left=0. ; // proportion of left section estimate to use
    float prop_right=0. ; // proportion of right section estimate to use

    if(sample_loc_in_section < 0.5) {
	section_right = sample_section ;
	prop_right= .5+sample_loc_in_section ;

	if(sample_section > 0) {
	    section_left = sample_section-1 ;
	    prop_left = 1.-prop_right ;
	} else {
	    prop_right=1. ;
	}
    } else {
	section_left = sample_section ;
	prop_left= 1.5-sample_loc_in_section ;

	if(sample_section < N_SKY_SECTIONS-1) {
	    section_right = sample_section+1 ;
	    prop_right = 1.-prop_left ;
	} else {
	    prop_left=1. ;
	}
    }

    fprintf(stderr, "px:%5.2f section(l,r):%2d,%2d prop (l,r):%5.2f,%5.2f\n", px+0.5, section_left, section_right, prop_left, prop_right) ;
}

float normalized_value(float px, float py, double coefs[N_SKY_SECTIONS][5])
{
    float section_width= (float)1./(float)N_SKY_SECTIONS ;
    // interpolation of predicted normalized value from section coefs
    int sample_section = (int)((px+0.5)/section_width) ;

    // this happens at px=0.5
    if(sample_section == N_SKY_SECTIONS)
	sample_section=N_SKY_SECTIONS-1 ;

    float lower = -.5 + (float)sample_section/(float)N_SKY_SECTIONS ;
    float upper = lower + section_width ;

    // sample loc in section is 0.0 at left edge of section, 1.0 at right, 0.5 in middle
    float sample_loc_in_section = (px-lower)/(upper-lower) ;
    int section_left=-1, section_right=-1 ;

    float prop_left=0. ; // proportion of left section estimate to use
    float prop_right=0. ; // proportion of right section estimate to use

    if(sample_loc_in_section < 0.5) {
	section_right = sample_section ;
	prop_right= .5+sample_loc_in_section ;

	if(sample_section > 0) {
	    section_left = sample_section-1 ;
	    prop_left = 1.-prop_right ;
	} else {
	    prop_right=1. ;
	}
    } else {
	section_left = sample_section ;
	prop_left= 1.5-sample_loc_in_section ;

	if(sample_section < N_SKY_SECTIONS-1) {
	    section_right = sample_section+1 ;
	    prop_right = 1.-prop_left ;
	} else {
	    prop_left=1. ;
	}
    }

    float v_left=0., v_right=0. ;
    if(section_left != -1) {
	double *scoefs = coefs[section_left] ;
	v_left = scoefs[0] + scoefs[1]*px
	    +scoefs[2]*log(py)
	    +scoefs[3]*px*px
	    +scoefs[4]*py
	    ;
    }
    if(section_right != -1) {
	double *scoefs = coefs[section_right] ;
	v_right = scoefs[0] + scoefs[1]*px
	    +scoefs[2]*log(py)
	    +scoefs[3]*px*px
	    +scoefs[4]*py
	    ;
    }


    // compute normalized values for all samples (makes it easier for plotting )
    return prop_left*v_left + prop_right*v_right ;
}

float get_splined_sat(float px, float py)
{
    return normalized_value(px, py, coefs_s) ;
}

void compute_normalized_s_and_v(int n_samples, SKYFILL_DATA_t *pData)
{
    // break sky up into N_SKY_SECTIONS sections, compute a normalized value
    struct mstat m_fit_v, m_fit_s ;
    float section_width= (float)1./(float)N_SKY_SECTIONS ;

    // for computing the grand total mean py ;
    float sum_py=0. ;
    float n_py=0. ;

    // compute coefs for each section
    for(int section = 0 ; section < N_SKY_SECTIONS ; section++) {

	float lower = -.5 + (float)section/(float)N_SKY_SECTIONS ;
	float upper = lower + section_width ;
	fprintf(stderr, "Compute normalized values for %f to %f\n", lower, upper) ;

	int n_found=0 ;

	while(1) {
	    float x_min=1.e30 ;
	    float x_max=-1.e30 ;
	    int n=0 ;
	    for(int i = 0 ; i < n_samples ; i++) {
		if(samples[i].px >= lower && samples[i].px < upper && samples[i].is_outlier == 0) {
		    n++ ;
		    if(x_min > samples[i].px) x_min = samples[i].px ;
		    if(x_max < samples[i].px) x_max = samples[i].px ;
		}
	    }

	    if(n < 9 || (x_max-x_min) < .01) {
		lower -= section_width ;
		upper += section_width ;
		if(lower < 0. && upper > 1.) {
		    fprintf(stderr, "FATAL: no good samples found for sky.  Adjust masking or end of sky tolerances\n") ;
		    fprintf(stderr, "       Rerun with -d3 flag to for clues on how to proceed.\n") ;
		    exit(1) ;
		}
	    } else {
		break ;
	    }
	}

	m_fit_v = init_reg(5) ;
	m_fit_s = init_reg(5) ;
	for(int i = 0 ; i < n_samples ; i++) {
	    if(samples[i].px >= lower && samples[i].px < upper && samples[i].is_outlier == 0) {
/*  		double xv[4] = {1., samples[i].px, samples[i].py, samples[i].px*samples[i].py} ;  */
		double xv[5] = {1., samples[i].px, log(samples[i].py), samples[i].px*samples[i].px, samples[i].py} ;
		sum_reg(&m_fit_v, xv, samples[i].v, 1.) ;
		sum_reg(&m_fit_s, xv, samples[i].s, 1.) ;
		sum_py += samples[i].py ;
		n_py += 1. ;
		n_found++ ;
	    }
	}

	if(n_found > 8) {
	    estimate_reg(&m_fit_v, coefs_v[section]) ;
	    estimate_reg(&m_fit_s, coefs_s[section]) ;
	} else {
	    fprintf(stderr, "No normalized coefficients found for %f to %f section\n", lower, upper) ;
	    coefs_s[section][0] = 0. ;
	    coefs_s[section][1] = 0. ;
	    coefs_s[section][2] = 0. ;
	    coefs_s[section][3] = 0. ;
	    coefs_s[section][4] = 0. ;

	    coefs_v[section][0] = 0. ;
	    coefs_v[section][1] = 0. ;
	    coefs_v[section][2] = 0. ;
	    coefs_v[section][3] = 0. ;
	    coefs_v[section][4] = 0. ;
	}
    }

    n_normalized_samples=0 ;

    float mean_py = sum_py / n_py ;

    if(0) {
	// interpolation of predicted normalized value from section coefs
	for(int i = 0 ; i < n_samples ; i++) {
	    float normalized_v = normalized_value(samples[i].px, mean_py, coefs_v) ;
	    float normalized_s = normalized_value(samples[i].px, mean_py, coefs_s) ;

	    // compute normalized values for all samples (makes it easier for plotting )
	    samples[i].normalized_v = normalized_v ;
	    samples[i].normalized_s = normalized_s ;

	    // samples are in order of X first, no need to repeat same X value
	    // so the normalized data can be very small (and quick to fit!)
	    if(n_normalized_samples > 0) {
		if(samples[i].x != normalized_samples[n_normalized_samples-1].x) {
		    fprintf(stderr, "add normalized sample at x:%d, n_normalized:%d\n", samples[i].x, n_normalized_samples) ;

		    normalized_samples[n_normalized_samples] = samples[i] ;
		    n_normalized_samples++ ;
		} else {
    /*  		fprintf(stderr, "skip normalized sample at x:%d\n", samples[i].x) ;  */
		}
	    } else {
		normalized_samples[n_normalized_samples] = samples[i] ;
		fprintf(stderr, "add normalized sample at x:%d\n", samples[i].x) ;
		n_normalized_samples++ ;
	    }
	}
    } else {


	// interpolation of predicted normalized value from section coefs
	for(int i = 0 ; i < n_samples ; i++) {
	    float normalized_v = normalized_value(samples[i].px, mean_py, coefs_v) ;
	    float normalized_s = normalized_value(samples[i].px, samples[i].py, coefs_s) ;

	    // compute normalized values for all samples (makes it easier for plotting )
	    samples[i].normalized_v = normalized_v ;
	    samples[i].normalized_s = normalized_s ;
	}

	for(float px=-0.5 ; px < 0.500001 ; px += KNOT_PWIDTH) {
	    float normalized_v = normalized_value(px, mean_py, coefs_v) ;
	    float normalized_s = normalized_value(px, mean_py, coefs_s) ;

	    int i = n_normalized_samples ;
	    if(i == N_BASE_KNOTS) {
		fprintf(stderr, "FATAL:  overran base_knots arrays\n") ;
		exit(1) ;
	    }

	    val_base_knots[i] = normalized_v ;
	    sat_base_knots[i] = normalized_s ;
	    normalized_samples[i].px = px ;
	    normalized_samples[i].x = (uint16_t)((px+0.5)*(float)IMAGE_WIDTH) ;
	    normalized_samples[i].py = 1.0 ;
	    normalized_samples[i].y = IMAGE_HEIGHT-1 ;
	    normalized_samples[i].v = normalized_v ;
	    normalized_samples[i].s = normalized_s ;
	    normalized_samples[i].normalized_v = normalized_v ;
	    normalized_samples[i].normalized_s = normalized_s ;
	    n_normalized_samples++ ;
	}

	if(fabsf(pData->reduce_spline_value_range) > .001) {
	    float sum_v=0., sum_v_reduced=0. ;
	    float reduce_factor = 1.- pData->reduce_spline_value_range ;

	    for(int i=0 ; i < N_BASE_KNOTS ; i++) {
		sum_v += val_base_knots[i] ;
		sum_v_reduced += val_base_knots[i]*reduce_factor ;
	    }

	    float reduce_intercept = (sum_v - sum_v_reduced)/(float)N_BASE_KNOTS ;

	    fprintf(stderr, "RV %f, %f\n", reduce_intercept, reduce_factor) ;

	    for(int i=0 ; i < N_BASE_KNOTS ; i++) {
		val_base_knots[i] = reduce_intercept + val_base_knots[i]*reduce_factor ;
	    }

	    for(int i=0 ; i < n_samples ; i++) {
		samples[i].normalized_v = reduce_intercept + samples[i].normalized_v*reduce_factor ;
	    }

	    for(int i=0 ; i < n_normalized_samples ; i++) {
		normalized_samples[i].normalized_v = reduce_intercept + normalized_samples[i].normalized_v*reduce_factor ;
	    }
	}

    }

    // normalize the normalize values around the mean
    float sum=0. ;
    for(int i=0 ; i < N_BASE_KNOTS ; i++) {
	sum += val_base_knots[i] ;
    }

    float factor = 1./( sum/(float)N_BASE_KNOTS ) ;

    for(int i=0 ; i < N_BASE_KNOTS ; i++) {
	val_base_knots[i] *= factor ;
    }

    for(int i=0 ; i < n_samples ; i++) {
	samples[i].normalized_v *= factor ;
    }

    fprintf(stderr, "final n_normalized:%d\n", n_normalized_samples) ;

    // check
    for(int i=0 ; i < 2 ; i++) {
	fprintf(stderr, "knot[%d] = %f\n", i, val_base_knots[i]) ;
    }

    for(float px=-.5 ; px < -0.5+KNOT_PWIDTH*1.1 ; px += KNOT_PWIDTH/2.) {
	float v ;
	PX2_base_v(px,v) ;
	fprintf(stderr, "px:%f, px2knotvalue %f\n", px, v) ;
    }

    for(float px=-.5 ; px < 0.51 ; px += 1./20.) {
	print_sky_section_info(px, 0.) ;
    }
/*      exit(1) ;  */

}

void fit_sky_model(int n_samples, SKYFILL_DATA_t *pData)
{
    fprintf(stderr, "Fit models\n") ;

    if(! pData->fix_sky_hue)
	compute_normalized_s_and_v(n_samples, pData) ;

    fit_sky_model_step2(pData, n_samples) ;

    // is the top of the sky too dark relative to mid-sky values  and the full value model used?, if so try the
    // reduced value model

    if(pData->val_model_full == 1) {
	float sum_py = 0. ;
	for(int i = 0 ; i < n_samples ; i++) {
	    sum_py += samples[i].py ;
	}
	float py_mid = sum_py/(float)n_samples ;

	// 
	int n_failed=0 ;
	for(float px = -0.5 ; px < 0.505 ; px += .01) {
	    float hsv_hat_top[3] ;
	    float hsv_hat_mid[3] ;
	    // estimate hsv, will used blended model if requested
	    predict_sky_hsv(px, 1.0, hsv_hat_top) ;
	    predict_sky_hsv(px, py_mid, hsv_hat_mid) ;


	    if(hsv_hat_top[2] < 100./255 && hsv_hat_top[2] < hsv_hat_mid[2]) n_failed++ ;
	}

	if(n_failed > 10) {
	    fprintf(stderr, "\n*******************************************************************************\n") ;
	    fprintf(stderr, "********** TOP OF SKY IS TOO DARK, changing to reduced value model  ***********\n") ;
	    fprintf(stderr, "*******************************************************************************\n\n") ;

/*  	    pData->val_model_full = 0 ;  */
/*  	    pData->val_model_is_linear = 1 ;  */
/*    */
/*  	    fit_sky_model_step2(pData, n_samples) ;  */
	}

    }
}

int sample_sky_points(int n_per_column, int n_columns,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky, SKYFILL_DATA_t *pData, int need_model_fit)
{
    int x, y ;
    float mean_column_length ;
    float sum_column_lengths=0. ;
    int n_columns_in_sum=0 ;

    // for this processing, do not consider localized sun image.
    int save_uses_CIE_model = uses_CIE_model ;

    if(uses_CIE_model & CIE_SUN_MODEL)
	uses_CIE_model = CIE_NO_MODEL ;

    V_sun= image_relative_pixel_to_V3(sun_x_angle2px(pData->sun_x), pData->sun_py) ;

    for(x=0 ; x < IMAGE_WIDTH ; x++) {
	if(pData->column_mask[x] == 1) continue ;
	if(pData->column_sample_mask[x] == 1) continue ;

	 sum_column_lengths += pData->end_of_sky[x] - pData->start_of_sky[x] + 1 ;
	 n_columns_in_sum++ ;
    }

    mean_column_length = sum_column_lengths/(float)n_columns_in_sum ;

#define MAX_SAMPLES 1100
    // try to get 1000 sample points, each will be a mean of a 5x5 area around the actual point
    // do this with a rectangular grid over the image area, refining the grid to the sky area only

/*      n_columns = 50 ;  */
/*      int n_rows = (int)(1000./(float)n_columns * (float)IMAGE_HEIGHT/mean_column_length+0.5) ;  */
    n_columns = 100 ;
    int n_rows = (int)(1000./(float)n_columns * (float)IMAGE_HEIGHT/mean_column_length+0.5) ;

    int dx = IMAGE_WIDTH/n_columns ;
    int dy = IMAGE_HEIGHT/n_rows ;
    int n_samples=0 ;

    // how many valid samples in the grid?
    for(x=0 ; x < IMAGE_WIDTH ; x += dx) {
	if(pData->column_mask[x] == 1) continue ;
	if(pData->column_sample_mask[x] == 1) continue ;

	for(y=0 ; y < IMAGE_HEIGHT ; y += dy) {

	    if(y >= pData->start_of_sky[x] && y <= pData->end_of_sky[x])
		n_samples++ ;
	}
    }

    printf("Using a grid over the sky with %d samples\n", n_samples) ;

    if(samples == NULL) {
	samples = (struct sample_point *) calloc(n_samples+100, sizeof(struct sample_point)) ;
	normalized_samples = (struct sample_point *) calloc(n_samples+100, sizeof(struct sample_point)) ;
	val_base_knots = (float *) calloc(N_BASE_KNOTS, sizeof(float)) ;
	sat_base_knots = (float *) calloc(N_BASE_KNOTS, sizeof(float)) ;
    }

    int max_samples = n_samples+100 ;

    n_samples=0 ;
    n_normalized_samples=0 ;

    float sum_val_wgts = 0. ;
    float sum_val=0. ;
    float sum_hue_wgts = 0. ;
    float sum_hue=0. ;
    float sum_sat_wgts = 0. ;
    float sum_sat=0. ;
    int n_accks=0 ;

    for(x=0 ; x < IMAGE_WIDTH ; x += dx) {
	uint16_t rgb[3] ;

	if(pData->column_mask[x] == 1) continue ; // don't attempt a mean sky sample inside a masked area
	if(pData->column_sample_mask[x] == 1) continue ;
	if(pData->fix_sky_hue && x>=pData->min_sky_hue_mask && x <= pData->max_sky_hue_mask) continue ;

	// want to get close to end of sky since the sky changes most near end of sky
	for(y = end_of_sky[x] - 6 ; y > start_of_sky[x]+3 ; y -= dy) {

/*  	for(y=0 ; y < IMAGE_HEIGHT ; y += dy) BRACKET  */
/*    */
/*  	    if(y < pData->start_of_sky[x] || y > pData->end_of_sky[x]) continue ;  */

	    int sx, sy ;
	    int n=0 ;
	    double sum_r=0. ;
	    double sum_g=0. ;
	    double sum_b=0. ;
	    double sumsq_r=0. ;
	    double sumsq_g=0. ;
	    double sumsq_b=0. ;

	    if(n_samples == max_samples) {
		fprintf(stderr, "ACCCK!  sample_sky_points(), too many samples\n") ;
		exit(0) ;
	    }

	    int x0=x-2 ;
	    int x1=x+3 ;
	    int y0=y-2 ;
	    int y1=y+3 ;
	    if(x0 < 0) x0 = 0 ;
	    if(y0 < 0) y0 = 0 ;
	    if(x1 > IMAGE_WIDTH) x1 = IMAGE_WIDTH ;
	    if(y1 > IMAGE_HEIGHT) y1 = IMAGE_HEIGHT ;

	    for(sx = x0 ; sx < x1 ; sx++) {
		if(pData->column_mask[x] == 1) continue ; // no sky sample from inside a masked area
		if(pData->column_sample_mask[x] == 1) continue ;
		if(pData->fix_sky_hue && sx>=pData->min_sky_hue_mask && sx <= pData->max_sky_hue_mask) continue ;

		for(sy = y0 ; sy < y1 ; sy++) {
		    if(sy >= pData->start_of_sky[sx] && sy <= pData->end_of_sky[sx]) {

			if(!is_in_test_mask(x,y)) {
			    tif_get3cv(image,sx,sy,rgb) ;
/*  			    float hsv[3] ;  */
/*  			    rgb2hsv16(rgb,hsv) ;  */

			    sum_r += rgb[0] ;
			    sum_g += rgb[1] ;
			    sum_b += rgb[2] ;
			    sumsq_r += (double)rgb[0]*(double)rgb[0] ;
			    sumsq_g += (double)rgb[1]*(double)rgb[1] ;
			    sumsq_b += (double)rgb[2]*(double)rgb[2] ;
			    n++ ;


			}
		    }
		}
	    }

	    // n == 0 could potentially happen near edges of detected sky area
	    // or areas with low saturation (i.e. clouds)
	    if(n == 0) continue ;

	    rgb[0] = (uint16_t)((float)sum_r/(float)n) ;
	    rgb[1] = (uint16_t)((float)sum_g/(float)n) ;
	    rgb[2] = (uint16_t)((float)sum_b/(float)n) ;

	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;
	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;

	    float hsv[3] ;
	    rgb2hsv16(rgb,hsv) ;

	    if(hsv[2] < .01 || hsv[1] < .01) {
		fprintf(stderr, "ACCK: px,py:%3.2f %3.2f s:%f v:%f r,g,b %d,%d,%d\n", px, py, hsv[1],hsv[2],rgb[0],rgb[1],rgb[2]) ;
		n_accks++ ;
	    }

	    initialize_sample(n_samples, x, y, px, py, rgb, hsv) ;

	    samples[n_samples].cv_r = cv_calc(n, sum_r, sumsq_r) ;
	    samples[n_samples].cv_g = cv_calc(n, sum_g, sumsq_g) ;
	    samples[n_samples].cv_b = cv_calc(n, sum_b, sumsq_b) ;

	    if(n_samples < 100)
		fprintf(stderr, "%d, %f %f\n", n, sum_r, sumsq_r) ;

	    float wgt = (1.-hsv[1])*PX_WGT(px) ;

	    sum_hue += hsv[0]*wgt ;
	    sum_hue_wgts += wgt ;
	    sum_sat += hsv[1] ;
	    sum_sat_wgts += 1.f ;
	    sum_val += hsv[2]*wgt ;
	    sum_val_wgts += wgt ;

	    n_samples++ ;

	}

	//fprintf(stderr, "column %d, h=%f, s=%f, v=%f\n", x, sum_h/sum_n, sum_s/sum_n, sum_v/sum_n) ;
    }

    fprintf(stderr, "N_ACCKS=%d\n", n_accks) ;

    find_sky_outliers(n_samples) ;


    pData->hue_sky = sum_hue/sum_hue_wgts ;
    pData->sat_sky = sum_sat/sum_sat_wgts ;
    pData->val_sky = sum_val/sum_val_wgts ;

    if(need_model_fit) {
	fit_sky_model(n_samples, pData) ;

	output_model_results(n_samples, pData) ;
    }

    uses_CIE_model = save_uses_CIE_model ;

    return n_samples ;
}

void output_model_results(int n_samples, SKYFILL_DATA_t *pData)
{

    fprintf(stderr, "*******************  CALCULATED PHASE SHIFT, horizon py, slope and curvature %3.0f, (%4.2f,%5.3f,%5.3f)*************************\n",
	    pData->valsat_phase_shift, pData->horizon_py,pData->horizon_slope,pData->horizon_curvature) ;

    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	fprintf(stderr, "FINAL fit_full_HSV[%d of %d]:\n", lr,pData->full_hsv_model_level) ;

	fprintf(stderr, "raw_hsv_coefs H ... half=>%3.0f spread=>%5.1f speed=>%5.2f center=>%5.2f px_center:%5.2f\n",
	        raw_hsv_coefs[lr].h_coefs[0], raw_hsv_coefs[lr].h_coefs[1],
	        raw_hsv_coefs[lr].h_coefs[2], raw_hsv_coefs[lr].h_coefs[3],
		raw_hsv_coefs[lr].h_coefs[4]) ;

#ifdef QUADRATIC_SAT_VERSION
	fprintf(stderr, "raw_hsv_coefs S .................................\n") ;
	fprintf(stderr, "                    intercept:%7.4f horizon:%7.4f,%7.4f,%7.4f damp:%7.4f rate:%7.4f\n",
										raw_hsv_coefs[lr].s_coefs[5],
										raw_hsv_coefs[lr].s_coefs[0],
										raw_hsv_coefs[lr].s_coefs[3],
										raw_hsv_coefs[lr].s_coefs[4],
										raw_hsv_coefs[lr].s_coefs[1],
										raw_hsv_coefs[lr].s_coefs[2]) ;
#else
	fprintf(stderr, "raw_hsv_coefs S .................................\n") ;
	fprintf(stderr, "                    intercept:%7.4f horizon:%7.4f,%7.4f damp:%7.4f rate:%7.4f\n",
										raw_hsv_coefs[lr].s_coefs[4],
										raw_hsv_coefs[lr].s_coefs[0],
										raw_hsv_coefs[lr].s_coefs[3],
										raw_hsv_coefs[lr].s_coefs[1],
										raw_hsv_coefs[lr].s_coefs[2]) ;
	fprintf(stderr, "                    0-4:%7.4f %7.4f %7.4f %7.4f %7.4f\n",
										raw_hsv_coefs[lr].s_coefs[0],
										raw_hsv_coefs[lr].s_coefs[1],
										raw_hsv_coefs[lr].s_coefs[2],
										raw_hsv_coefs[lr].s_coefs[3],
										raw_hsv_coefs[lr].s_coefs[4]) ;
#endif
	if(pData->val_model_is_linear) {
	    fprintf(stderr, "raw_hsv_coefs V .................................\n") ;
	    fprintf(stderr, "                    Intercept:%7.4f base_v:%7.4f\n", raw_hsv_coefs[lr].v_coefs[0],
						raw_hsv_coefs[lr].v_coefs[1]
						) ;
	    fprintf(stderr, "                    py_hc:%7.4f py_hc**2:%7.4f\n", raw_hsv_coefs[lr].v_coefs[2],
						raw_hsv_coefs[lr].v_coefs[3]
						) ;
	} else if(uses_CIE_model & CIE_FULL_MODEL) {
	    fprintf(stderr, "CIE model coefs .................................\n") ;
	    fprintf(stderr, "                a,b,c,d,e: %7.4f %7.4f %7.4f %7.4f %7.4f\n",
					     pData_fit->perez_A,
					     pData_fit->perez_B,
					     pData_fit->perez_C,
					     pData_fit->perez_D,
					     pData_fit->perez_E) ;

	} else {
	    fprintf(stderr, "raw_hsv_coefs V .................................\n") ;
	    fprintf(stderr, "                    a,b,c,d: %7.4f %7.4f %7.4f %7.4f\n", raw_hsv_coefs[lr].v_coefs[0],
					     raw_hsv_coefs[lr].v_coefs[1],
					     raw_hsv_coefs[lr].v_coefs[2],
					     raw_hsv_coefs[lr].v_coefs[3]) ;

	}

	fprintf(stderr, "\n") ;
    }

    output_sample_data(n_samples, "samples.dat") ;
    output_normalized_sample_data(n_normalized_samples, "normalized_samples.dat") ;
}

void output_normalized_sample_data(int n_normalized_samples, char filename[])
{
    if(pData_fit->output_sample_data == 0)
	return ;

    FILE *fp = fopen(filename, "w") ;

    if(fp != NULL) {
	fprintf(fp,"x y ha va px py nv ns gpsymbol\n") ;

	for(int i = 0 ; i < n_normalized_samples ; i++) {
	    float horz_angle = normalized_samples[i].px*pData_fit->FOV_horizontal ;
	    float vert_angle = normalized_samples[i].py*pData_fit->FOV_vertical ;
	    fprintf(fp, "%6d %6d", normalized_samples[i].x, normalized_samples[i].y) ;
	    fprintf(fp, " %5.1f", horz_angle) ;
	    fprintf(fp, " %5.1f", vert_angle) ;
	    fprintf(fp, " %6.4f %6.4f", normalized_samples[i].px, normalized_samples[i].py) ;
	    fprintf(fp, " %5.3f", normalized_samples[i].normalized_v) ;
	    fprintf(fp, " %5.3f", normalized_samples[i].normalized_s) ;
	    fprintf(fp, " 7") ;  // circle gnuplot symbol
	    fprintf(fp, "\n") ;

	}

	fclose(fp) ;

	fprintf(stderr, "Wrote %d normalized samples to %s\n", n_normalized_samples, filename) ;
    } else {
	fprintf(stderr, "WARNING: could not open %s for writing normalized sample data output\n", filename) ;
    }
}

void output_sample_data(int n_samples, char filename[])
{
    if(pData_fit->output_sample_data == 0)
	return ;

#ifdef HOW_CAN_THIS_WORK_WITH_360_DEGREE_PANORAMAS
	float FLmm=12. ; // focal length, mm
	float sensor_width_mm = 23.5 ;
	float n_pixels_x=6045 ;
	float image_scale=0.5 ; // scale factor of input image relative to original camera image dimensions ;
	float pixel_width = sensor_width_mm / n_pixels_x / image_scale ;
#endif

    int save_uses_CIE_model = uses_CIE_model ;

    if(uses_CIE_model & CIE_SUN_MODEL)
	uses_CIE_model = CIE_NO_MODEL ;

    FILE *fp = fopen(filename, "w") ;
    if(fp != NULL) {
	fprintf(stderr, "Writing %d samples to %s\n", n_samples, filename) ;

	fprintf(fp,"x y ha va px py hue sat val nv ns hue_hat sat_hat val_hat hue_err sat_err val_err n_close r g b cv_r cv_g cv_b  gpsymbol\n") ;

	for(int i = 0 ; i < n_samples ; i++) {
	    // estimate hsv, will used blended model if requested
	    float hsv_hat[3] ;
	    predict_sky_hsv(samples[i].px, samples[i].py, hsv_hat) ;

/*  	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;  */
/*  	    hsv_hat[2] = local_v_from_pxpy(samples[i].px, py_hc, val_coefs) ;  */

	    uint16_t rgb_hat[3] ;
	    hsv2rgb16(hsv_hat,rgb_hat) ;
/*  	    hsv_hat[1] = samples[i].normalized_s ;  */

	    float horz_angle = px_to_AOY(samples[i].px)*180./M_PI ;
	    float vert_angle = py_to_AOX(samples[i].py)*180./M_PI ;

	    //fprintf(fp, "%6.4f %6.4f", samples[i].px, samples[i].py+samples[i].px/50.) ;
	    fprintf(fp, "%6d %6d", samples[i].x, samples[i].y) ;
	    fprintf(fp, " %5.1f", horz_angle) ;
	    fprintf(fp, " %5.1f", vert_angle) ;
	    /*  	fprintf(fp, " %f", horz_dist_mm) ;  */
	    /*  	fprintf(fp, " %f", vert_dist_mm) ;  */
	    fprintf(fp, " %6.4f %6.4f", samples[i].px, samples[i].py) ;
	    fprintf(fp, " %6.2f %5.3f %5.3f", samples[i].h, samples[i].s, samples[i].v) ;
	    fprintf(fp, " %5.3f", samples[i].normalized_v) ;
	    fprintf(fp, " %5.3f", samples[i].normalized_s) ;
	    fprintf(fp, " %6.2f %5.3f %5.3f", hsv_hat[0], hsv_hat[1], hsv_hat[2]) ;
	    fprintf(fp, " %6.2f %5.3f %5.3f", hsv_hat[0]-samples[i].h, hsv_hat[1]-samples[i].s, hsv_hat[2]-samples[i].v) ;
	    fprintf(fp, " %d", samples[i].n_close_neighbors) ;
	    // normalize to blue values
	    fprintf(fp, " %7.5f %7.5f %7.5f", (float)samples[i].r/(float)samples[i].b, (float)samples[i].g/(float)samples[i].b, (float)samples[i].b/MAX16f) ;
/*  	    fprintf(fp, " %7.5f %7.5f %7.5f", (float)samples[i].r/MAX16f, (float)samples[i].g/MAX16f, (float)samples[i].b/MAX16f) ;  */
	    fprintf(fp, " %f %f %f", samples[i].cv_r, samples[i].cv_g, samples[i].cv_b) ;

	    if(samples[i].is_outlier == 0) {
		fprintf(fp, " 7") ;  // circle gnuplot symbol
	    } else {
		fprintf(fp, " 1") ;  // X gnuplot symbol
	    }
	    fprintf(fp, "\n") ;

	}

	fclose(fp) ;

	fprintf(stderr, "Wrote %d samples to %s\n", n_samples, filename) ;
    } else {
	fprintf(stderr, "WARNING: could not open %s for writing sample data output\n", filename) ;
    }

    uses_CIE_model = save_uses_CIE_model ;
}

int read_samples_and_fit_model(char *samplefile, SKYFILL_DATA_t *pData)
{

    FILE *fp = fopen(samplefile, "r") ;
    int n_samples=0 ;
    float px,py, h_hat,s_hat,v_hat,h_err,s_err,v_err,angle,cv_r,cv_g,cv_b ;
    float hsv[3] ;
    float rgb[3] ;
    int n_close ;
    pData->lowest_sky_py_found=1.  ;

    if(fp != NULL) {
	fprintf(stderr, "Reading sample data from %s\n", samplefile) ;

	IMAGE_HEIGHT = 600 ; // not used for predictions, but necessary to be set
	IMAGE_WIDTH =  800 ; // not used for predictions, but necessary to be set

	p_half_image_width = 0.5 ;
	V_zenith= P3toV3(0., 1., 0.) ;
	V_sun= P3toV3(0., 1., 0.) ; // default with sun directly overhead ;
	pData->maximum_CIE_vhat = find_maximum_CIE_vhat() ;

	if(samples == NULL) {
	    samples = (struct sample_point *) calloc(MAX_SAMPLES*2, sizeof(struct sample_point)) ;
	    normalized_samples = (struct sample_point *) calloc(MAX_SAMPLES*2, sizeof(struct sample_point)) ;
	}

	char buf[1024] ; // no line will ever be this big, really

	while( fgets(buf, sizeof(buf), fp)) {
	    if(sscanf(buf, "%f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f",
		&px,&py,&hsv[0],&hsv[1],&hsv[2],&h_hat,&s_hat,&v_hat,&h_err,&s_err,&v_err,&angle,&n_close,&rgb[0],&rgb[1],&rgb[2],&cv_r,&cv_g,&cv_b) == 19) {

		pData->FOV_horizontal = angle / px ;

		uint16_t x = (uint16_t)((px+.5)*100000.+0.5) ; // need x initialized for normalized sample logic to filter samples
		uint16_t y = (uint16_t)((py+.5)*100000.+0.5) ; // may as well do y too
		uint16_t zeros[3] = {0,0,0} ;

		initialize_sample(n_samples, x, y, px, py, zeros, hsv) ;

		samples[n_samples].r = (uint16_t)(rgb[0]*MAX16f+0.5) ;
		samples[n_samples].g = (uint16_t)(rgb[1]*MAX16f+0.5) ;
		samples[n_samples].b = (uint16_t)(rgb[2]*MAX16f+0.5) ;
		samples[n_samples].cv_r = cv_r ;
		samples[n_samples].cv_g = cv_g ;
		samples[n_samples].cv_b = cv_b ;

		if(pData->lowest_sky_py_found > py)
		    pData->lowest_sky_py_found = py ;
		n_samples++ ;
	    }
	}
	fprintf(stderr, "Read %d samples from %s\n", n_samples, samplefile) ;

	fclose(fp) ;
	fprintf(stderr, "FOV %f\n", pData->FOV_horizontal) ;
	set_FOV_factor() ;

	find_sky_outliers(n_samples) ;
	fit_sky_model(n_samples, pData) ;

	output_model_results(n_samples, pData) ;
	free(normalized_samples) ;
	free(samples) ;
	free(val_base_knots) ;
	free(sat_base_knots) ;

    } else {
	fprintf(stderr, "Could not open %s\n", samplefile) ;
    }

    return 0 ;

}

void find_sky_outliers(int n_samples)
{

    fprintf(stderr, "Finding outliers...\n") ;

    // outlier detection -- compute mean distances to all other points
    for(int i = 0 ; i < n_samples-1 ; i++) {
	for(int j=i+1 ; j < n_samples ; j++) {
	    float ds = samples[i].s - samples[j].s ;
	    float dv = samples[i].v - samples[j].v ;
	    float dh = samples[i].h - samples[j].h ;
	    float dy = samples[i].py - samples[j].py ;
	    float dx = samples[i].px - samples[j].px ;

	    float dist_sx = ds*ds + dx*dx ;
	    float dist_sv = ds*ds + dv*dv ;
	    float dist_vy = dy*dy + dv*dv ;
	    float dist_hv = dh*dh + dv*dv ;

	    if(dist_sv < .005 && dist_sx < .01 && dist_hv < 36) {
		samples[i].n_close_neighbors++ ;
		samples[j].n_close_neighbors++ ;
	    }
	}
    }

    for(int i = 0 ; i < n_samples ; i++) {
	if(samples[i].n_close_neighbors < 2)
	    samples[i].is_outlier = 255 ;
    }

    // recompute mean v,s
    float sum_v=0., sum_s=0., n=0. ;
    for(int i = 0 ; i < n_samples ; i++) {
	if(samples[i].is_outlier == 0) {
	    sum_v += samples[i].v ;
	    sum_s += samples[i].s ;
	    n += 1. ;
	}
    }

    float mean_v = sum_v / n ;
    float mean_s = sum_s / n ;

    // look for suspicious samples with low s and high v, relative to mean values
    for(int i = 0 ; i < n_samples ; i++) {
	if(samples[i].is_outlier == 0) {
	    float relative_v = samples[i].v / mean_v ;
	    float relative_s = samples[i].s / mean_s ;

	    if(relative_s < .05 && relative_v > .95) {
		samples[i].is_outlier = 255 ;
	    }

	    // normal sky hues
	    if(samples[i].h < 150 || samples[i].h > 250) {
		samples[i].is_outlier = 255 ;
	    }
	}
    }
}

void initialize_sample(int i, uint16_t x, uint16_t y, float px, float py, uint16_t *rgb, float *hsv)
{
    samples[i].x = x ;
    samples[i].y = y ;
    samples[i].px = px ;
    samples[i].py = py ;
    samples[i].r = rgb[0] ;
    samples[i].g = rgb[1] ;
    samples[i].b = rgb[2] ;
    samples[i].abs_error[0] = 1. ;
    samples[i].abs_error[1] = 1. ;
    samples[i].is_outlier = 0 ;

    samples[i].h = hsv[0] ;
    samples[i].s = hsv[1] ;
    samples[i].v = hsv[2] ;
    samples[i].cv_r=0. ;
    samples[i].cv_g=0. ;
    samples[i].cv_b=0. ;
    samples[i].n_close_neighbors=0 ;
}
