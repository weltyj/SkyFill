#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "tiffio.h"

#include "skyfill_tif.h"
#include "mstat.h"
#include "colorspace_conversions.h"
#include "pixel_tests.h"
#include "sample_and_fit_sky_model.h"
#include "amoeba_06.h"

void output_model_results(int n_samples, SKYFILL_DATA_t *pData) ;
void initialize_sample(int n_samples, uint16_t x, uint16_t y, float px, float py, uint16_t r, uint16_t g, uint16_t b, float h, float s, float v) ;
void find_sky_outliers(int n_samples) ;

struct sample_point *samples = NULL ;
int n_samples_to_optimize = 0 ;

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

    float sse=0. ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;

    wgt_func = px_wgt_single ;

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float s_hat = local_s_from_pxpy(samples[i].px, py_hc, coefs_all) ;

	float err = samples[i].s - s_hat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;

	if(samples[i].is_outlier == 0) {
	    sse += err*err*wgt ;
	}

	if(g_calc_sample_error) {
	    samples[i].abs_s_error[g_location_index] = fabs(err) ;
/*  	    if(fabs(err) > 0.1) {  */
/*  		samples[i].is_outlier |= SAT_OUTLIER ;  */
/*  	    }  */
	}
    }

    return sse ;
}
float sat_optimize_function(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;
    float coefs_all[7] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = coefs[3] ;
    coefs_all[4] = coefs[4] ;
    coefs_all[5] = coefs[5] ;
    coefs_all[6] = coefs[6] ;
    coefs_all[7] = coefs[7] ;

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
	    samples[i].abs_s_error[g_location_index] = fabs(err) ;
/*  	    if(fabs(err) > 0.1) {  */
/*  		samples[i].is_outlier |= SAT_OUTLIER ;  */
/*  	    }  */
	}
    }

    return sse ;
}

void get_saturation_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift)
{
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;


/*      float x_mapped = (pData_fit->FOV_horizontal*px+pData_fit->valsat_phase_shift) ;  */
/*      return -cos(   pow((x_mapped)/180.,coefs[0])*M_PI  ) * coefs[1] + (coefs[2] + coefs[3]*py_hc) ;  */

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


#ifndef TEST_SAT2
    float F=1. ;  // flatness
    F=0.75 ;  // flatness
    float sat_range = max_sat-min_sat - (max_py-min_py)/2. ;
    sat_range = max_sat-min_sat ;
    fprintf(stderr, "============ SAT START0: min_sat:%f  max_sat:%f F:%f  sat_range:%f\n", min_sat, max_sat, F, sat_range) ;
    if(sat_range < .1) sat_range=.1 ;
    if(sat_range > .3) sat_range=.3 ;
    fprintf(stderr, "============ SAT START1: min_sat:%f  max_sat:%f F:%f  sat_range:%f\n", min_sat, max_sat, F, sat_range) ;

#define NSP 7
    // this little struct just makes it easier to see the parameters and their constraints
    struct valparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
	} parms[NSP] = {
		     {F, -.2, 0.45, 3.2}, // flatness at the bottom of the "pit"
		     {sat_range*1.1, .05, .01, sat_range*2.},  // scale
		     {min_sat,  .05, -3., 3.}, // min 
		     {0.5,  .05, -3.0, 8.}, // coef on py
		     {-0.1,  .001, -8.0, 8.}, // coef on log(py)
		     {1,  .5, -5.0, 5.}, // phase shift adjust
		     {1.5,  .1, 0.8, 2.0}, // angle factor
		     } ;
#else

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
    fprintf(stderr, "============ SAT START: horizon:%f, damp:%f rate:%f\n", best_horizon, best_damp, best_rate) ;
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
		     {-.001, .0001, -.03, .03},  // rate px*px
		     } ;
#endif

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
if(1) {

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

#ifndef TEST_SAT2
    min_sat=1. ; // need to look at samples
    max_sat=0. ; // need to look at samples
    min_py=1.5 ;  // need to look at samples for min py_hc
    max_py=0. ;  // need to look at samples for min py_hc

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	if(samples[i].is_outlier == 0) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	    if(max_py < py_hc) max_py = py_hc ;
	    if(min_py > py_hc) min_py = py_hc ;
	    if(min_sat > samples[i].s) min_sat= samples[i].s ;
	    if(max_sat < samples[i].s) max_sat= samples[i].s ;
	} else {
/*  	    fprintf(stderr, "S:%f, HA:%f is an outlier\n", samples[i].s, samples[i].px*pData_fit->FOV_horizontal) ;  */
	}
    }

    sat_range = max_sat-min_sat - (max_py-min_py)/2. ;
    sat_range = max_sat-min_sat ;
    fprintf(stderr, "============ SAT START2: min_sat:%f  max_sat:%f F:%f  sat_range:%f\n", min_sat, max_sat, F, sat_range) ;
    if(sat_range < .1) sat_range=.1 ;
    if(sat_range > .3) sat_range=.3 ;
    fprintf(stderr, "============ SAT START3: min_sat:%f  max_sat:%f F:%f  sat_range:%f\n", min_sat, max_sat, F, sat_range) ;

    guess[1] = 1.5 ;
    guess[2] = sat_range*1.1 ;
    c_hi[2] = sat_range*2. ;
    guess[3] = min_sat ;
#endif

    g_calc_sample_error=0 ;
    AmoebaFit(sat_optimize_function,1.e-4,200,NSP,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    sat_optimize_function(guess) ;

    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;
    coefs[4] = guess[5] ;
} else {
    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = 0. ;
    coefs[4] = 0. ;
}
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
	    samples[i].abs_h_error[g_location_index] = fabs(err) ;
/*  	    if(fabs(err) > 40.) {  */
/*  		samples[i].is_outlier |= HUE_OUTLIER ;  */
/*  	    }  */
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
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	if(min_py > py_hc) min_py = py_hc ;
	if(min_hue > samples[i].h) min_hue= samples[i].h ;
	if(max_hue < samples[i].h) max_hue= samples[i].h ;
    }

    float hue_spread= max_hue-min_hue ;

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
}

#ifndef LINEAR_VAL_MODEL
#ifdef TEST_VAL_MODEL
float val_optimize_function(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;
    float coefs_all[7] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = coefs[3] ;
    coefs_all[4] = coefs[4] ;
    coefs_all[5] = coefs[5] ;
    coefs_all[6] = coefs[6] ;

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
	float v_hat = local_v_from_pxpy(samples[i].px, py_hc, coefs_all) ;

	float err = samples[i].v - v_hat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;

	if(samples[i].is_outlier == 0) {
	    sse += err*err*wgt ;
	}

	if(g_calc_sample_error) {
	    samples[i].abs_v_error[g_location_index] = fabs(err) ;
/*  	    if(fabs(err) > 0.2) {  */
/*  		samples[i].is_outlier |= VAL_OUTLIER ;  */
/*  	    }  */
	}
    }

    return sse ;
}

void get_value_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift)
{
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;


/*      float x_mapped = (pData_fit->FOV_horizontal*px*(coefs[1]+coefs[4]*py_hc)+90+pData_fit->valsat_phase_shift)*M_PI/180. ;  */
/*      float Z = 1./(1. + abs(cos(x_mapped))**(coefs[2]+coefs[5]*py_hc) ;  */
/*      return coefs[0]*(Z-.5) + coefs[3] ;  */

    float min_val=1. ; // need to look at samples
    float max_val=0. ; // need to look at samples
    float min_py=1.5 ;  // need to look at samples for min py_hc
    float max_py=0. ;  // need to look at samples for min py_hc

    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	if(max_py < py_hc) max_py = py_hc ;
	if(min_py > py_hc) min_py = py_hc ;
	if(min_val > samples[i].v) min_val= samples[i].v ;
	if(max_val < samples[i].v) max_val= samples[i].v ;
    }


    float val_range = max_val-min_val ;
    fprintf(stderr, "============ VAL START0: min_val:%f  max_val:%f val_range:%f\n", min_val, max_val, val_range) ;
    if(val_range < .1) val_range=.1 ;
    if(val_range > .8) val_range=.8 ;
    fprintf(stderr, "============ VAL START1: min_val:%f  max_val:%f val_range:%f\n", min_val, max_val, val_range) ;

    // this little struct just makes it easier to see the parameters and their constraints
    struct valparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
	} parms[7] = {
		     {val_range*2., .1, val_range*0.5, 4.},  // (max-min)*2
		     {.8,  .01, 0., 1.6}, // S0
		     {2.5,  .05, 2., 3.}, // F0
		     {(1.+min_val)/2.,  .001, min_val, 1.}, // min value, M0
		     {-.1,  .02, -0.2, 0.}, // S1
		     {-4., .05, -8., 0.}, // F1
		     {-.1, .05, -.2, 0.}, // M1
		     } ;

    float guess[8] ;
    float guess_delta[8] ;
    float c_lo[8] ;
    float c_hi[8] ;

    for(int i=1 ; i <= 7 ; i++) {
	guess[i] = parms[i-1].guess ;
	guess_delta[i] = parms[i-1].guess_delta ;
	c_lo[i] = parms[i-1].lo ;
	c_hi[i] = parms[i-1].hi ;
    }

    g_calc_sample_error=0 ;
    AmoebaFit(val_optimize_function,1.e-4,200,7,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    val_optimize_function(guess) ;

    g_calc_sample_error=0 ;
    AmoebaFit(val_optimize_function,1.e-4,200,7,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    val_optimize_function(guess) ;

    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;
    coefs[4] = guess[5] ;
    coefs[5] = guess[6] ;
    coefs[6] = guess[7] ;
}
#else

float value_optimize_function(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;
    float coefs_all[9] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = coefs[3] ;
    coefs_all[4] = coefs[4] ;
    coefs_all[5] = coefs[5] ;
    coefs_all[6] = 0. ;
    coefs_all[7] = 0. ;
    coefs_all[8] = 0. ;

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
	float v_hat = local_v_from_pxpy(samples[i].px, py_hc, coefs_all) ;

	float err = samples[i].v - v_hat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
	if(samples[i].is_outlier == 0) {
	    sse += err*err*wgt ;
	}
	if(g_calc_sample_error) {
	    samples[i].abs_v_error[g_location_index] = fabs(err) ;
/*  	    if(fabs(err) > 0.2) {  */
/*  		samples[i].is_outlier |= VAL_OUTLIER ;  */
/*  	    }  */
	}
    }

    return sse ;
}

void get_value_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift)
{
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;

    g_val_model_full = pData->val_model_full ;
    double py_coefs[3] = {0,0,0} ;

    if(1) {
	float (*wgt_func)(float,float) ;

	if(g_full_hsv_model_level == 1) {
	    wgt_func = px_wgt_single ;
	} else if(g_location_index == 0) {
	    wgt_func = px_wgt_left ;
	} else {
	    wgt_func = px_wgt_right ;
	}

	struct mstat m_fit ;

	m_fit = init_reg(2+g_val_model_full) ;

	for(int i = 0 ; i < n_samples ; i++) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
/*  	    double xv[3] = {1.,  py_hc, py_hc*py_hc} ;  */
	    double xv[3] = {1.,  py_hc, log(py_hc)} ;
	    double y = samples[i].v ;
	    double w = 1. ;
	    sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, w*w*w)) ;
	}

	estimate_reg(&m_fit, py_coefs) ;
	fprintf(stderr, "***************** VAL QUAD COEFS  %f,%f,%f\n", py_coefs[0], py_coefs[1], py_coefs[2]) ;
    }

    // x_mapped = (FOV*px-B0)*M_PI/180.
    // model form is B1/(cos(x_mapped)*B2+1.35)+B3+B4*py_hc + B5*py_hc*py_hc ;
    // constraint B2 < 1.35 ;

    // this little struct just makes it easier to see the parameters and their constraints
    struct valparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
	} parms[6] = {
		     {2.0, .05, 1.9, 2.1}, // aka angle_factor intercept
		     {-.001, .00001, 0.01, 0.}, // slope coef on angle_factor
		     {.33, .05, 0.1, 1.},  // scale
		     {1.,  .05, 0.3, VAL_INV_ADD-.01}, // min to max scale -- max constraint very important
		     {.64, .05, -1., 1.}, // simple intercept coef
		     {-0.3, .01, -1., -.1}, // coef on py_hc
		     } ;

    float guess[10] ;
    float guess_delta[10] ;
    float c_lo[10] ;
    float c_hi[10] ;

    for(int i=1 ; i <= 6 ; i++) {
	guess[i] = parms[i-1].guess ;
	guess_delta[i] = parms[i-1].guess_delta ;
	c_lo[i] = parms[i-1].lo ;
	c_hi[i] = parms[i-1].hi ;
    }

    g_calc_sample_error=0 ;
    AmoebaFit(value_optimize_function,1.e-4,200,6+g_val_model_full,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    value_optimize_function(guess) ;

    g_calc_sample_error=0 ;
    AmoebaFit(value_optimize_function,1.e-4,200,6+g_val_model_full,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    value_optimize_function(guess) ;

    if(1) {
	float (*wgt_func)(float,float) ;

	// these are what we will find, but must set to 0. first so
	// local_v_from_pxpy uses unmodified prediction in this step
	guess[7] = 0. ;
	guess[8] = 0. ;
	guess[9] = 0. ;

	if(g_full_hsv_model_level == 1) {
	    wgt_func = px_wgt_single ;
	} else if(g_location_index == 0) {
	    wgt_func = px_wgt_left ;
	} else {
	    wgt_func = px_wgt_right ;
	}

	struct mstat m_fit ;

	m_fit = init_reg(3) ;

	for(int i = 0 ; i < n_samples ; i++) {

	    if(samples[i].is_outlier == 0) {
		float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
		double xv[3] = {1.,  py_hc, py_hc*py_hc} ;
		float v_hat = local_v_from_pxpy(samples[i].px, py_hc, &guess[1]) ;
		double y = samples[i].v - v_hat ;
		double w = 1. ;
		sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, w*w*w)) ;
	    }

	}

	estimate_reg(&m_fit, py_coefs) ;
	fprintf(stderr, "***************** VAL ERROR COEFS  %f,%f,%f\n", py_coefs[0], py_coefs[1], py_coefs[2]) ;
    }

    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;
    coefs[4] = guess[5] ;
    coefs[5] = guess[6] ;
    coefs[6] = py_coefs[0] ;
    coefs[7] = py_coefs[1] ;
    coefs[8] = py_coefs[2] ;

    if(g_val_model_full) {
    } else {
	fprintf(stderr, "FATAL: settting coefs[5] to 0.!\n") ;
	exit(1) ;
	coefs[5] = 0. ;
    }
}
#endif

// for ifndef LINEAR_VAL_MODEL
#endif

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


    float peak, best_peak ;
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


float fit_full_HSV_model(int n_samples, int location_index, int calc_error, SKYFILL_DATA_t *pData)
{
    int save_uses_CIE_model = uses_CIE_model ;

    // make sure CIE model is not used in this step
    uses_CIE_model = 0 ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;

    if(pData->full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }

    struct HSV_model_coefs *p = &raw_hsv_coefs[location_index] ;

    // Full model of HSV


    // will be filled by estimate_reg
    double coefs[7] ;
    struct mstat m_fit ;


/*      Note -- this function is 0 at at 0.5, so subtracting .5 from py_hc ...  */
/*      hue = 800 * (py_hc-.5) / (1 + abs(py_hc-.5)*20) + 180 ;  */
    // Hue
    get_hue_model(n_samples, p->h_coefs, pData, location_index, pData->valsat_phase_shift) ;

    // Sat
    // No longer done in this module, have switched to a nonlinear function

    // Val

#ifdef LINEAR_VAL_MODEL
    // Val
    if(pData->val_model_full) {
	m_fit = init_reg(5) ;
	for(int i = 0 ; i < n_samples ; i++) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	    float x_mapped = pData->FOV_horizontal*samples[i].px*M_PI/180.*pData->angle_factor ;
	    double xv[5] = {1.,  cos(x_mapped),sin(x_mapped),py_hc,py_hc*py_hc} ;
	    double y = samples[i].v ;
	    sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, samples[i].abs_v_error[location_index])) ;
	}

	estimate_reg(&m_fit, coefs) ;
	p->v_coefs[0] = coefs[0] ;
	p->v_coefs[1] = coefs[1] ;
	p->v_coefs[2] = coefs[2] ;
	p->v_coefs[3] = coefs[3] ;
	p->v_coefs[4] = coefs[4] ;
    } else {
	m_fit = init_reg(4) ;
	for(int i = 0 ; i < n_samples ; i++) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	    float x_mapped = pData->FOV_horizontal*samples[i].px*M_PI/180.*pData->angle_factor ;
	    double xv[4] = {1.,  cos(x_mapped),sin(x_mapped),py_hc} ;
	    double y = samples[i].v ;
	    sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, samples[i].abs_v_error[location_index])) ;
	}

	estimate_reg(&m_fit, coefs) ;
	p->v_coefs[0] = coefs[0] ;
	p->v_coefs[1] = coefs[1] ;
	p->v_coefs[2] = coefs[2] ;
	p->v_coefs[3] = coefs[3] ;
	p->v_coefs[4] = 0. ;
    }
#endif

    //compute combined SSE from each r,g,b component
    float sse_r = 0. ;
    float sse_g = 0. ;
    float sse_b = 0. ;

    for(int i = 0 ; i < n_samples ; i++) {
	float h,s,v ;
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	h = local_h_from_pxpy(samples[i].px, py_hc, raw_hsv_coefs[location_index].h_coefs) ;
	s = local_s_from_pxpy(samples[i].px, py_hc, raw_hsv_coefs[location_index].s_coefs) ;
	v = local_v_from_pxpy(samples[i].px, py_hc, raw_hsv_coefs[location_index].v_coefs) ;

/*  	float fpsat = 1.0 - (float)(samples[i].y)/(float)pData->start_of_sky[samples[i].x] ;  */
/*  	if(fpsat < 0.0) fpsat = 0.0 ;  */
/*  	s *= 1.0 - sqrt(fpsat)*(1.-pData->final_saturation_factor) ;  */

	uint16_t r,g,b ;
	hsv2rgb16(h,s,v,&r,&g,&b) ;

	float err_r = (float)r - (float)samples[i].r ;
	float err_g = (float)g - (float)samples[i].g ;
	float err_b = (float)b - (float)samples[i].b ;
	float err_h = (float)h - (float)samples[i].h ;
	float err_s = (float)s - (float)samples[i].s ;
	float err_v = (float)v - (float)samples[i].v ;

	if(calc_error) {
	    float e = (fabs(err_r) + fabs(err_g) + fabs(err_b)) ;
	    samples[i].abs_error[location_index] = e ;
	    samples[i].abs_h_error[location_index] = fabs(err_h) ;
	    samples[i].abs_s_error[location_index] = fabs(err_s) ;
	    samples[i].abs_v_error[location_index] = fabs(err_v) ;
	}

	sse_r += err_r*err_r ;
	sse_g += err_g*err_g ;
	sse_b += err_b*err_b ;
    }

    uses_CIE_model = save_uses_CIE_model ;

    return sse_r + sse_g + sse_b ;
}

void fit_sky_model(SKYFILL_DATA_t *pData, int n_samples)
{
    pData->angle_factor = 2.0 ;

    // a little hack -- horizon_py has not effect on the models now, so 
    // by setting horizon_was_set to 1 it is not used in the grid search
    pData->horizon_was_set=1 ;

    if(! pData->horizon_was_set) {
	pData->horizon_py = pData->lowest_sky_py_found-.1  ;
    }

    float starting_horizon_py = pData->horizon_py ;

    float best_sse = 1.e30;
    float best_horizon_py = pData->horizon_py ;
    pData->horizon_curvature = 0. ;
    float best_horizon_curvature = pData->horizon_curvature ;
    float best_angle_factor = pData->angle_factor ;
    float lr_phase_shift[2] ;

    if(pData->valsat_phase_shift > -5000)
	pData->valsat_phase_shift = get_phase_shift(n_samples, 0, 0, pData, 0) ;

    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	lr_phase_shift[lr] = pData->valsat_phase_shift ;
    }

    for(pData->horizon_curvature=-.1 ; pData->horizon_curvature < .101 ; pData->horizon_curvature += .02) {
    //for(pData->horizon_curvature=0. ; pData->horizon_curvature < .01 ; pData->horizon_curvature += .02) {

	for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	    reset_sample_errors(n_samples) ;

	    float phase_shift = lr_phase_shift[lr] ;
#ifndef LINEAR_VAL_MODEL
	    get_value_model(n_samples, raw_hsv_coefs[lr].v_coefs, pData, lr, phase_shift) ;
#endif
	    get_saturation_model(n_samples, raw_hsv_coefs[lr].s_coefs, pData, lr, phase_shift) ;

	}

	if(! pData->horizon_was_set) {

	    for(pData->horizon_py = starting_horizon_py ; pData->horizon_py > .05 ; pData->horizon_py -= .05) {
		reset_sample_errors(n_samples) ;
		float best_sse_this=1.e30 ;

/*  #define N_LS_ITERATIONS 3  */
#define N_LS_ITERATIONS 1
		// iteratively reweighted least squares
		for(int iteration=0 ; iteration < N_LS_ITERATIONS ; iteration++) {
		    float sse=0. ;
		    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
			sse += fit_full_HSV_model(n_samples,lr,1,pData) ;
		    }
		    if(sse < best_sse) {
			best_sse = sse ;
			best_horizon_py = pData->horizon_py ;
			best_horizon_curvature = pData->horizon_curvature ;
		    }
		    if(sse < best_sse_this) {
			best_sse_this = sse ;
		    }
		}

		fprintf(stderr, "GRID fit_full_HSV  hpy,hcurv(%f,%f) = %f\n", pData->horizon_py,pData->horizon_curvature,best_sse_this) ;
	    }

	} else {

	    for(pData->angle_factor=1. ; pData->angle_factor < 3.1 ; pData->angle_factor += .5) {
		reset_sample_errors(n_samples) ;
		float best_sse_this=1.e30 ;

		// iteratively reweighted least squares
		for(int iteration=0 ; iteration < N_LS_ITERATIONS ; iteration++) {
		    float sse=0. ;
		    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
			sse += fit_full_HSV_model(n_samples,lr,1,pData) ;
		    }
		    if(sse < best_sse) {
			best_sse = sse ;
			best_horizon_curvature = pData->horizon_curvature ;
			best_angle_factor = pData->angle_factor ;
		    }
		    if(sse < best_sse_this) {
			best_sse_this = sse ;
		    }
		}

		fprintf(stderr, "GRID fit_full_HSV  vfactor,hcurv(%f,%f) = %f\n", pData->angle_factor,pData->horizon_curvature,best_sse_this) ;
	    }

	}
    }

    // final fit with best horizon_py and horizon_curvature
    pData->angle_factor = best_angle_factor ;
    pData->horizon_py = best_horizon_py ;
    pData->horizon_curvature = best_horizon_curvature ;


    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	float phase_shift = lr_phase_shift[lr] ;
#ifndef LINEAR_VAL_MODEL
	get_value_model(n_samples, raw_hsv_coefs[lr].v_coefs, pData, lr, phase_shift) ;
#endif
	get_saturation_model(n_samples, raw_hsv_coefs[lr].s_coefs, pData, lr, phase_shift) ;

	float sse = fit_full_HSV_model(n_samples,lr,0,pData) ; // last time through, calculate errors on points

	fprintf(stderr, "*******************  CALCULATED PHASE SHIFT %3.0f*************************\n", phase_shift) ;

	fprintf(stderr, "FINAL fit_full_HSV[%d of %d]:  hcurv is %f, angle_factor %f, sse= %f\n",
			lr,pData->full_hsv_model_level,pData->horizon_curvature,pData->angle_factor, sse) ;

	fprintf(stderr, "\n") ;
    }
}

void fit_sky_model_step1(int n_samples, SKYFILL_DATA_t *pData)
{
    fprintf(stderr, "Fit models\n") ;
    // for this processing, do not consider localized sun image.
    int save_uses_CIE_model = uses_CIE_model ;
    uses_CIE_model = 0 ;

    fit_sky_model(pData, n_samples) ;

    // is the top of the sky too dark, and the full value model used?, if so try the
    // reduced value model

#ifdef NEED_ANOTHER_IDEA
    if(pData->val_model_full == 1) {
	int n_failed=0 ;
	for(float px = -0.5 ; px < 0.505 ; px += .01) {
	    int x = IMAGE_RELATIVE_TO_PIXEL_X(px) ;
	    float py_mid = IMAGE_PIXEL_Y_TO_RELATIVE(end_of_sky[x]/2.) ;
	    float hhat_top, shat_top, vhat_top ;
	    float hhat_mid, shat_mid, vhat_mid ;
	    // estimate hsv, will used blended model if requested
	    predict_sky_hsv(px, 1.0, &hhat_top, &shat_top, &vhat_top) ;
	    predict_sky_hsv(px, py_mid, &hhat_mid, &shat_mid, &vhat_mid) ;


	    if(vhat_top < 100./255 && vhat_top < vhat_mid) n_failed++ ;
	}

	if(n_failed > 10) {
	    fprintf(stderr, "\n*******************************************************************************\n") ;
	    fprintf(stderr, "********** TOP OF SKY IS TOO DARK, changing to reduced value model  ***********\n") ;
	    fprintf(stderr, "*******************************************************************************\n\n") ;

	    pData->val_model_full = 0 ;
	    fit_sky_model(pData, n_samples) ;
	}
    }
#endif

    uses_CIE_model = save_uses_CIE_model ;
}

int sample_sky_points(int n_per_column, int n_columns,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky, SKYFILL_DATA_t *pData)
{
    int x, y ;
    float mean_column_length ;
    float sum_column_lengths=0. ;
    int n_columns_in_sum=0 ;

    V_sun= image_relative_pixel_to_V3_noclip(sun_x_angle2px(pData->sun_x), pData->sun_py) ;

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
    }

    int max_samples = n_samples+100 ;

    n_samples=0 ;

    float sum_val_wgts = 0. ;
    float sum_val=0. ;
    float sum_hue_wgts = 0. ;
    float sum_hue=0. ;
    float sum_sat_wgts = 0. ;
    float sum_sat=0. ;

    for(x=0 ; x < IMAGE_WIDTH ; x += dx) {
	uint16_t r,g,b ;

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
	    float sum_r=0. ;
	    float sum_g=0. ;
	    float sum_b=0. ;

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
			    tif_get3c(image,x,y,r,g,b) ;
			    float h,s,v ;
			    rgb2hsv16(r,g,b,&h,&s,&v) ;

			    sum_r += r ;
			    sum_g += g ;
			    sum_b += b ;
			    n++ ;


			}
		    }
		}
	    }

	    // n == 0 could potentially happen near edges of detected sky area
	    // or areas with low saturation (i.e. clouds)
	    if(n == 0) continue ;

	    r = (uint16_t)(sum_r/(float)n) ;
	    g = (uint16_t)(sum_g/(float)n) ;
	    b = (uint16_t)(sum_b/(float)n) ;

	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;
	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    initialize_sample(n_samples, x, y, px, py, r, g, b, h, s, v) ;

	    float wgt = (1.-s)*PX_WGT(px) ;

	    sum_hue += h*wgt ;
	    sum_hue_wgts += wgt ;
	    sum_sat += s ;
	    sum_sat_wgts += 1.f ;
	    sum_val += v*wgt ;
	    sum_val_wgts += wgt ;

	    n_samples++ ;

	}

	//fprintf(stderr, "column %d, h=%f, s=%f, v=%f\n", x, sum_h/sum_n, sum_s/sum_n, sum_v/sum_n) ;
    }

    find_sky_outliers(n_samples) ;


    pData->hue_sky = sum_hue/sum_hue_wgts ;
    pData->sat_sky = sum_sat/sum_sat_wgts ;
    pData->val_sky = sum_val/sum_val_wgts ;

    fit_sky_model_step1(n_samples, pData) ;

    output_model_results(n_samples, pData) ;

    return n_samples ;
}

void output_model_results(int n_samples, SKYFILL_DATA_t *pData)
{

    fprintf(stderr, "*******************  CALCULATED PHASE SHIFT, horizon curvature %3.0f, %4.2f*************************\n",
	    pData->valsat_phase_shift, pData->horizon_curvature) ;

    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	fprintf(stderr, "FINAL fit_full_HSV[%d of %d]:\n", lr,pData->full_hsv_model_level) ;

	fprintf(stderr, "raw_hsv_coefs H ... half=>%3.0f spread=>%5.1f speed=>%5.2f center=>%5.2f px_center:%5.2f\n",
	        raw_hsv_coefs[lr].h_coefs[0], raw_hsv_coefs[lr].h_coefs[1],
	        raw_hsv_coefs[lr].h_coefs[2], raw_hsv_coefs[lr].h_coefs[3],
		raw_hsv_coefs[lr].h_coefs[4]) ;

#ifdef TEST_SAT2
	fprintf(stderr, "raw_hsv_coefs S .................................\n") ;
	fprintf(stderr, "                    horizon:%7.4f,%7.4f,%7.4f damp:%7.4f rate:%7.4f\n", raw_hsv_coefs[lr].s_coefs[0],
										raw_hsv_coefs[lr].s_coefs[3],
										raw_hsv_coefs[lr].s_coefs[4],
										raw_hsv_coefs[lr].s_coefs[1],
										raw_hsv_coefs[lr].s_coefs[2]) ;
#elif TEST_SAT_MODEL
	fprintf(stderr, "raw_hsv_coefs S .................................\n") ;
	fprintf(stderr, "                    flatness:%7.4f range:%7.4f\n", raw_hsv_coefs[lr].s_coefs[0], raw_hsv_coefs[lr].s_coefs[1]) ;
	fprintf(stderr, "                    min_min:%7.4f max_min:%7.4f, %7.4f\n", raw_hsv_coefs[lr].s_coefs[2],
										raw_hsv_coefs[lr].s_coefs[3], raw_hsv_coefs[lr].s_coefs[4]) ;
	fprintf(stderr, "                    phase tune:%7.4f angle_factor:%7.4f\n", raw_hsv_coefs[lr].s_coefs[5],raw_hsv_coefs[lr].s_coefs[6]) ;
#endif

#ifdef LINEAR_VAL_MODEL
	fprintf(stderr, "raw_hsv_coefs V .................................\n") ;
	fprintf(stderr, "                    Intercept:%7.4f cos:%7.4f sin:%7.4f\n", raw_hsv_coefs[lr].v_coefs[0],
						raw_hsv_coefs[lr].v_coefs[1],
						raw_hsv_coefs[lr].v_coefs[2]
						) ;
	fprintf(stderr, "                    py_hc:%7.4f py_hc**2:%7.4f\n", raw_hsv_coefs[lr].v_coefs[3],
						raw_hsv_coefs[lr].v_coefs[4]
						) ;
#else
#ifdef TEST_VAL_MODEL
	fprintf(stderr, "raw_hsv_coefs V .................................\n") ;
	fprintf(stderr, "                    B0:%7.4f S0:%7.4f S1:%7.4f\n", raw_hsv_coefs[lr].v_coefs[0],
						raw_hsv_coefs[lr].v_coefs[1],
						raw_hsv_coefs[lr].v_coefs[4]
						) ;

	fprintf(stderr, "                    F0:%7.4f F1:%7.4f\n", raw_hsv_coefs[lr].v_coefs[2], raw_hsv_coefs[lr].v_coefs[5]) ;
	fprintf(stderr, "                    M0:%7.4f M1:%7.4f\n", raw_hsv_coefs[lr].v_coefs[3], raw_hsv_coefs[lr].v_coefs[6]) ;
#else
	fprintf(stderr, "raw_hsv_coefs V .................................\n") ;
	fprintf(stderr, "                    angle_factor:%5.2f angle_factor_py_hc:%5.2f\n", raw_hsv_coefs[lr].v_coefs[0], raw_hsv_coefs[lr].v_coefs[1]) ;
	fprintf(stderr, "                    scale:%5.2f peak2peakscale:%5.2f\n", raw_hsv_coefs[lr].v_coefs[1], raw_hsv_coefs[lr].v_coefs[2]) ;
	fprintf(stderr, "                    QUAD(py_hc) %5.2f,%5.2f,%5.2f\n",
					    raw_hsv_coefs[lr].v_coefs[3], raw_hsv_coefs[lr].v_coefs[4], raw_hsv_coefs[lr].v_coefs[5]) ;
	fprintf(stderr, "                    ERROR_CORRECT(py_hc) %5.2f,%5.2f,%5.2f\n",
					    raw_hsv_coefs[lr].v_coefs[6], raw_hsv_coefs[lr].v_coefs[7], raw_hsv_coefs[lr].v_coefs[8]) ;
#endif
#endif

	fprintf(stderr, "\n") ;
    }

    FILE *fp = fopen("samples.dat", "w") ;
    for(int i = 0 ; i < n_samples ; i++) {
	float h_hat, s_hat, v_hat ;
	// estimate hsv, will used blended model if requested
	uint16_t r_hat, g_hat, b_hat ;
	predict_sky_hsv(samples[i].px, samples[i].py, &h_hat, &s_hat, &v_hat) ;
	hsv2rgb16(h_hat,s_hat,v_hat,&r_hat,&g_hat,&b_hat) ;
	float angle = pData->FOV_horizontal*samples[i].px ;

	//fprintf(fp, "%6.4f %6.4f", samples[i].px, samples[i].py+samples[i].px/50.) ;
	fprintf(fp, "%6.4f %6.4f", samples[i].px, samples[i].py) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", samples[i].h, samples[i].s, samples[i].v) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", h_hat, s_hat, v_hat) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", h_hat-samples[i].h, s_hat-samples[i].s, v_hat-samples[i].v) ;
	fprintf(fp, " %5.1f", angle) ;
	fprintf(fp, " %f %f", samples[i].sum_dist2_sv, samples[i].sum_dist2_vy) ;
	fprintf(fp, "\n") ;

    }

    fclose(fp) ;
}

int read_samples_and_fit_model(char *samplefile, SKYFILL_DATA_t *pData)
{

    FILE *fp = fopen(samplefile, "r") ;
    int n_samples=0 ;
    float px,py,h,s,v,h_hat,s_hat,v_hat,h_err,s_err,v_err,angle,n_close,n_close2 ;

    if(fp != NULL) {
	fprintf(stderr, "Reading sample data from %s\n", samplefile) ;

	if(samples == NULL) {
	    samples = (struct sample_point *) calloc(MAX_SAMPLES*2, sizeof(struct sample_point)) ;
	}

	while( fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f",
	    &px,&py,&h,&s,&v,&h_hat,&s_hat,&v_hat,&h_err,&s_err,&v_err,&angle,&n_close,&n_close2) == 14) {

	    pData->FOV_horizontal = angle / px ;

	    initialize_sample(n_samples, 0, 0, px, py, 0, 0, 0, h, s, v) ;
	    n_samples++ ;
	}
	fprintf(stderr, "Read %d samples from %s\n", n_samples, samplefile) ;

	fclose(fp) ;

	fit_sky_model_step1(n_samples, pData) ;

	output_model_results(n_samples, pData) ;
    } else {
	fprintf(stderr, "Could not open %s\n", samplefile) ;
    }

}

void find_sky_outliers(int n_samples)
{

    fprintf(stderr, "Finding outliers...\n") ;

    // outlier detection -- compute mean distances to all other points
    for(int i = 0 ; i < n_samples-1 ; i++) {
	for(int j=i+1 ; j < n_samples ; j++) {
	    float ds = samples[i].s - samples[j].s ;
	    float dv = samples[i].v - samples[j].v ;
	    float dy = samples[i].py - samples[j].py ;
	    float dist_sv = ds*ds + dv*dv ;
	    float dist_vy = dy*dy + dv*dv ;
	    if(dist_sv < .1 && dist_vy < .1) {
/*  		samples[i].sum_dist2_sv += ds*ds + dv*dv ;  */
/*  		samples[i].sum_dist2_vy += dy*dy + dv*dv ;  */
/*  		samples[j].sum_dist2_sv += ds*ds + dv*dv ;  */
/*  		samples[j].sum_dist2_vy += dy*dy + dv*dv ;  */
		samples[i].n_close_neighbors++ ;
		samples[j].n_close_neighbors++ ;
	    }
	}
    }

/*      float sum_close=0. ;  */
/*      float sum_close2=0. ;  */
/*    */
/*      for(int i = 0 ; i < n_samples ; i++) {  */
/*  	sum_close += (float)samples[i].n_close_neighbors ;  */
/*  	sum_close2 += (float)samples[i].n_close_neighbors * (float)samples[i].n_close_neighbors ;  */
/*      }  */
/*    */
/*      float N=n_samples ;  */
/*    */
/*      float mean_close = sum_close/N ;  */
/*    */
/*      float sd_close = sqrt( 1./(N-1.)*sum_close2 - 1./(N*(N-1.)) * sum_close*sum_close ) ;  */

    // put classic Z score in sum_dist2
    for(int i = 0 ; i < n_samples ; i++) {
/*  	samples[i].sum_dist2_sv = mean_close-(float)samples[i].n_close_neighbors ;  */
/*  	samples[i].sum_dist2_sv /= sd_close ;  */

	// let sum_dist2_sv hold the n_close_neighbors value so it will available to make a plot with gnuplot
	samples[i].sum_dist2_sv = (float)samples[i].n_close_neighbors ;

	samples[i].sum_dist2_vy = samples[i].sum_dist2_sv ;

	if(samples[i].n_close_neighbors < 125)
	    samples[i].is_outlier = 255 ;
    }
}

void initialize_sample(int n_samples, uint16_t x, uint16_t y, float px, float py, uint16_t r, uint16_t g, uint16_t b, float h, float s, float v)
{
	    samples[n_samples].x = x ;
	    samples[n_samples].y = y ;
	    samples[n_samples].px = px ;
	    samples[n_samples].py = py ;
	    samples[n_samples].r = r ;
	    samples[n_samples].g = g ;
	    samples[n_samples].b = b ;
	    samples[n_samples].abs_error[0] = 1. ;
	    samples[n_samples].abs_error[1] = 1. ;
	    samples[n_samples].is_outlier = 0 ;

	    samples[n_samples].h = h ;
	    samples[n_samples].s = s ;
	    samples[n_samples].v = v ;
	    samples[n_samples].sum_dist2_sv=0. ;
	    samples[n_samples].sum_dist2_vy=0. ;
	    samples[n_samples].n_close_neighbors=0 ;
}
