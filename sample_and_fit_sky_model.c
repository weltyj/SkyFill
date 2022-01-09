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

float sat_optimize_function(float one_based_coefs[])
{
    float sse=0. ;
    float *coefs = &one_based_coefs[1] ;

    //static inline float (*wgt_func)(float) ;
    float (*wgt_func)(float,float) ;

    if(g_full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(g_location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }

#define SAT_COEF7 (1)
    for(int i = 0 ; i < n_samples_to_optimize ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;

	float x_adjust=0. ;

	if(g_sat_model_full) {
	    float angle_factor = coefs[5] ;
	    float x_mapped = ( (pData_fit->FOV_horizontal*samples[i].px+pData_fit->valsat_phase_shift)*angle_factor-180.)
	                      *M_PI/180. ;
	    x_adjust= (coefs[6]-1./(cos(x_mapped)+ 1.+1./coefs[6]))/10. ;
	}

	float shat = coefs[0]/(1.+exp(-(py_hc-coefs[1])*coefs[2]))+coefs[3]+coefs[4]*(py_hc-coefs[1]) + x_adjust ;
	float err = samples[i].s - shat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
	sse += err*err*wgt ;
	if(g_calc_sample_error) {
	    samples[i].abs_s_error[g_location_index] = fabs(err) ;
	}
    }

/*      fprintf(stderr, "SOF: SSE(%f) = %f,%f,%f,%f\n", sse, coefs[1],coefs[2],coefs[3],coefs[4]) ;  */

    return sse ;
}

void get_saturation_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift)
{
    float max_sat=0. ;
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;

    for(int i = 0 ; i < n_samples ; i++) {
	if(max_sat < samples[i].s) max_sat = samples[i].s ;
    }

    // model form is B0/(1.+exp(-(py_hc-B1)*B2))+B3+B4*py_hc ;

    // this little struct just makes it easier to see the parameters and their constraints
    struct satparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
	} parms[8] = {
		     {0.,0.,0.,0.},
		     {(max_sat+1.)/2., .1, max_sat, 1.},  // saturation asymptot (sorta)
		     {0.1, .05, -.5, .5}, // adjustment to horizon ... py_hc -= coefs[2]
		     {10., 1., .1, 100.}, // factor on py_hc ... log(1. + py_hc*coefs[3]
		     {-.1, .05, -1., 1.}, // simple intercept term
		     {.01, .001, -.1, .1}, // linear on py_hc
		     {1.5, .2, 1.0, 3.0}, // aka angle_factor

/*  		     {phase_shift, 1, phase_shift-10, phase_shift+10.}, // centering for x_adjust  */

		     {SAT_COEF7, .2, 1, 10} // factor on x_adjust
		     } ;

    float guess[8] ;
    float guess_delta[8] ;
    float c_lo[8] ;
    float c_hi[8] ;

    for(int i=1 ; i <= 7 ; i++) {
	guess[i] = parms[i].guess ;
	guess_delta[i] = parms[i].guess_delta ;
	c_lo[i] = parms[i].lo ;
	c_hi[i] = parms[i].hi ;
    }

    g_calc_sample_error=0 ;
    AmoebaFit(sat_optimize_function,1.e-4,200,5+g_sat_model_full*2,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    sat_optimize_function(guess) ;

    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;
    coefs[4] = guess[5] ;

    if(g_sat_model_full) {
	coefs[5] = guess[6] ;
	coefs[6] = guess[7] ;
    } else {
#ifdef TEST_SAT_MODEL
	coefs[5] = SAT_COEF6 ;
	coefs[6] = SAT_COEF7 ;
#else
	coefs[5] = 0. ;
	coefs[6] = 0. ;
#endif
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
	sse += err*err*wgt ;
	if(g_calc_sample_error) {
	    samples[i].abs_h_error[g_location_index] = fabs(err) ;
	}
    }

    return sse ;
}

void get_hue_model(int n_samples, float coefs[], SKYFILL_DATA_t *pData, int location_index, float phase_shift)
{
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;

/*      float half=coefs[0] ;  */
/*      float spread=coefs[1] ;  */
/*      float speed=coefs[2] ;  */
/*      float center=coefs[3] ;  */

/*      return spread*(px-center)*speed/(1.+abs(px-center)*speed) + half ;  */
    float min_hue=180. ; // need to look at samples
    float max_hue=220. ; // need to look at samples
    float center=0.5 ;  // need to look at samples for min py_hc

    float hue_spread= max_hue-min_hue ;

    // this little struct just makes it easier to see the parameters and their constraints
    struct hueparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
	} parms[7] = {
		     {0.,0.,0.,0.},
		     {min_hue, 1., min_hue/2., min_hue+20},
		     {hue_spread, .1, hue_spread-20., hue_spread+20},
		     {40., .01, 5., 80.},
		     {center, .05, center-.1, 1.}
		     } ;

    float guess[8] ;
    float guess_delta[8] ;
    float c_lo[8] ;
    float c_hi[8] ;

    for(int i=1 ; i <= 4 ; i++) {
	guess[i] = parms[i].guess ;
	guess_delta[i] = parms[i].guess_delta ;
	c_lo[i] = parms[i].lo ;
	c_hi[i] = parms[i].hi ;
    }

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

float value_optimize_function(float one_based_coefs[])
{
    float *coefs = &one_based_coefs[1] ;
    float coefs_all[6] ;
    coefs_all[0] = coefs[0] ;
    coefs_all[1] = coefs[1] ;
    coefs_all[2] = coefs[2] ;
    coefs_all[3] = coefs[3] ;
    coefs_all[4] = coefs[4] ;
    if(g_val_model_full) {
	coefs_all[5] = coefs[5] ;
    } else {
	coefs_all[5] = 0. ;
    }
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
	float vhat = local_v_from_pxpy(samples[i].px, py_hc, coefs_all) ;

	float err = samples[i].v - vhat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
	sse += err*err*wgt ;
	if(g_calc_sample_error) {
	    samples[i].abs_v_error[g_location_index] = fabs(err) ;
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

    // x_mapped = (FOV*px-B0)*M_PI/180.
    // model form is B1/(cos(x_mapped)*B2+1.35)+B3+B4*py_hc + B5*py_hc*py_hc ;
    // constraint B2 < 1.35 ;

    // this little struct just makes it easier to see the parameters and their constraints
    struct valparms {
	float guess ;
	float guess_delta ;
	float lo ;
	float hi ;
	} parms[7] = {
		     {0.,0.,0.,0.},
		     {1.5, .05, 0.3, 2.2}, // aka angle_factor intercept

/*  		     {phase_shift, 10., phase_shift-40, phase_shift+40}, // centering for x_adjust  */

		     {.33, .05, 0.1, 1.},  // scale
		     {1.,  .05, 0.3, VAL_INV_ADD-.01}, // min to max scale -- max constraint very important
		     {.64, .05, -1., 1.}, // simple intercept coef
		     {-0.3, .01, -1., -.1}, // coef on py_hc
		     {-1, .1, -10., -.1}, // slope coef on angle_factor
		     } ;

    float guess[8] ;
    float guess_delta[8] ;
    float c_lo[8] ;
    float c_hi[8] ;

    for(int i=1 ; i <= 6 ; i++) {
	guess[i] = parms[i].guess ;
	guess_delta[i] = parms[i].guess_delta ;
	c_lo[i] = parms[i].lo ;
	c_hi[i] = parms[i].hi ;
    }

    g_calc_sample_error=0 ;
    AmoebaFit(value_optimize_function,1.e-4,200,5+g_val_model_full*1,guess,guess_delta,c_lo,c_hi,0) ;

    // need errors for final analysis in horizon curvature iteration
    g_calc_sample_error=1 ;
    value_optimize_function(guess) ;

    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;
    coefs[4] = guess[5] ;

    if(g_val_model_full) {
	coefs[5] = guess[6] ;
    } else {
	fprintf(stderr, "FATAL: settting coefs[5] to 0.!\n") ;
	exit(1) ;
	coefs[5] = 0. ;
    }
}

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

	struct valparms {
	    float guess ;
	    float guess_delta ;
	    float lo ;
	    float hi ;
	} parms[7] = {
	    {0.,0.,0.,0.},
	    {1,.1,0.1,3.},
	    {.1, .05, -1., 2.},
	    {-.6, .05, -5., 0.},
	} ;

	float guess[8] ;
	float guess_delta[8] ;
	float c_lo[8] ;
	float c_hi[8] ;

	for(int i=1 ; i <= 3 ; i++) {
	    guess[i] = parms[i].guess ;
	    guess_delta[i] = parms[i].guess_delta ;
	    c_lo[i] = parms[i].lo ;
	    c_hi[i] = parms[i].hi ;
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

float get_phase_shift_v2(int n_samples, int location_index, int calc_error, SKYFILL_DATA_t *pData, int verbose)
{

    float (*wgt_func)(float,float) ;

    if(pData->full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }

    double coefs[8] ;
    struct mstat m_fit ;

    float best_angle_factor = 1. ;

    m_fit = init_reg(5) ;

    for(int i = 0 ; i < n_samples ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float x_angle = pData->FOV_horizontal*samples[i].px ;
	double xv[5] = {1.,  x_angle,x_angle*x_angle,x_angle*x_angle*x_angle,py_hc} ;
	double y = samples[i].v ;
	double w = 1. ;
	sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, w*w*w)) ;
    }

    estimate_reg(&m_fit, coefs) ;

    float max_v_hat=-1.e30 ;
    float max_angle = 0 ;
    for(float angle = -180 ; angle < 181 ; angle +=10.) {
	float v_hat = 
			coefs[0] +
			coefs[1]*angle +
			coefs[2]*angle*angle +
			coefs[3]*angle*angle*angle +
			coefs[4]*0.5 ;

    if(verbose) {
	fprintf(stderr, "get_phase_shift: angle:%3.1f v_hat:%5.2f\n", angle, v_hat) ;
    }
	if(v_hat > max_v_hat) {
	    max_v_hat = v_hat ;
	    max_angle = angle ;
	}
    }

    float phase_shift = -max_angle ;

    if(verbose) {
	fprintf(stderr, "get_phase_shift: AF:%3.1f PS:%3.0f\n", pData_fit->angle_factor, phase_shift) ;
    }


    pData_fit->angle_factor = best_angle_factor ;
    return phase_shift ;
}

float get_phase_shift_v1(int n_samples, int location_index, int calc_error, SKYFILL_DATA_t *pData, int verbose)
{

    float (*wgt_func)(float,float) ;

    if(pData->full_hsv_model_level == 1) {
	wgt_func = px_wgt_single ;
    } else if(location_index == 0) {
	wgt_func = px_wgt_left ;
    } else {
	wgt_func = px_wgt_right ;
    }

    double coefs[8] ;
    struct mstat m_fit ;

    float best_phase_shift = 1. ;
    float best_angle_factor = 1. ;
    float best_sse = 1.e30 ;

    for(pData_fit->angle_factor = 1. ; pData_fit->angle_factor < 8.0 ; pData_fit->angle_factor += .1) {
/*      for(pData_fit->angle_factor = 1.8 ; pData_fit->angle_factor < 1.81 ; pData_fit->angle_factor += .1) {  */
	if(pData->val_model_full) {
	    m_fit = init_reg(5) ;
	    for(int i = 0 ; i < n_samples ; i++) {
		float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
		float x_mapped = pData->FOV_horizontal*samples[i].px*M_PI/180.*pData_fit->angle_factor ;
		double xv[5] = {1.,  cos(x_mapped),sin(x_mapped),py_hc,py_hc*py_hc} ;
		double y = samples[i].v ;
/*  		sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, samples[i].abs_v_error[location_index])) ;  */
		double w = samples[i].py ;
		sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, w*w*w)) ;
	    }

	    estimate_reg(&m_fit, coefs) ;
	} else {
	    m_fit = init_reg(4) ;
	    for(int i = 0 ; i < n_samples ; i++) {
		float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
		float x_mapped = pData->FOV_horizontal*samples[i].px*M_PI/180.*pData_fit->angle_factor ;
		double xv[4] = {1.,  cos(x_mapped),sin(x_mapped),py_hc} ;
		double y = samples[i].v ;
/*  		sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, samples[i].abs_v_error[location_index])) ;  */
		double w = samples[i].py ;
		sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, w*w*w)) ;
	    }

	    estimate_reg(&m_fit, coefs) ;
	    coefs[4] = 0. ;
	}

	float phase_shift = atan(coefs[1]/coefs[2]) * 180./M_PI ;
	float sse=0 ;
	for(int i = 0 ; i < n_samples ; i++) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	    float x_mapped = pData->FOV_horizontal*samples[i].px*M_PI/180.*pData_fit->angle_factor ;
	    double xv[5] = {1.,  cos(x_mapped),sin(x_mapped),py_hc,py_hc*py_hc} ;
	    float hat = coefs[0] +coefs[1]*xv[1] +coefs[2]*xv[2] +coefs[3]*xv[3] +coefs[4]*xv[4] ;
	    float err = hat - samples[i].v ;
	    double w = samples[i].py ;
	    sse += err*err * (*wgt_func)(samples[i].px, w*w*w) ;
	}

	if(verbose) {
	    fprintf(stderr, "get_phase_shift: AF:%3.1f PS:%3.0f SSE:%f\n", pData_fit->angle_factor, phase_shift, sse) ;
	}

	if(sse < best_sse) {
	    best_sse = sse ;
	    best_angle_factor = pData_fit->angle_factor ;
	    best_phase_shift = phase_shift ;
	}

    }


    pData_fit->angle_factor = best_angle_factor ;
    return best_phase_shift ;
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
    if(0) {
	m_fit = init_reg(3) ;
	for(int i = 0 ; i < n_samples ; i++) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	    double xv[3] = {1.,samples[i].px,log(0.+py_hc)} ;
	    double y = samples[i].h ;
	    sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, samples[i].abs_h_error[location_index])) ;
	}

	estimate_reg(&m_fit, coefs) ;
	p->h_coefs[0] = coefs[0] ;
	p->h_coefs[1] = coefs[1] ;
	p->h_coefs[2] = coefs[2] ;
	p->h_coefs[3] = 0. ;
    } else if (0) {
	m_fit = init_reg(2) ;
	for(int i = 0 ; i < n_samples ; i++) {
	    float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	    double xv[2] = {1.,1./(.001+py_hc)} ;
	    double y = log(samples[i].h) ;
	    sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, samples[i].abs_h_error[location_index])) ;
	}

	estimate_reg(&m_fit, coefs) ;
	p->h_coefs[0] = exp(coefs[0]) ;
	p->h_coefs[1] = coefs[1] ;
	p->h_coefs[2] = 0. ;
	p->h_coefs[3] = 0. ;
    } else {
	get_hue_model(n_samples, p->h_coefs, pData, location_index, pData->valsat_phase_shift) ;
    }

    // Sat
    // No longer done in this module, have switched to a nonlinear function

    // Val

#ifdef OLD_CODE
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
	p->coefs[2][0] = coefs[0] ;
	p->coefs[2][1] = coefs[1] ;
	p->coefs[2][2] = coefs[2] ;
	p->coefs[2][3] = coefs[3] ;
	p->coefs[2][4] = coefs[4] ;
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
	p->coefs[2][0] = coefs[0] ;
	p->coefs[2][1] = coefs[1] ;
	p->coefs[2][2] = coefs[2] ;
	p->coefs[2][3] = coefs[3] ;
	p->coefs[2][4] = 0. ;
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
	s = local_s_from_pxpy(samples[i].px, py_hc, &raw_hsv_coefs[location_index]) ;
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
    float lr_phase_shift[2] ;

    if(pData->valsat_phase_shift > -5000)
	pData->valsat_phase_shift = get_phase_shift(n_samples, 0, 0, pData, 0) ;

    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	lr_phase_shift[lr] = pData->valsat_phase_shift ;
    }

    //for(pData->horizon_curvature=-.1 ; pData->horizon_curvature < .101 ; pData->horizon_curvature += .02) {
    for(pData->horizon_curvature=0. ; pData->horizon_curvature < .01 ; pData->horizon_curvature += .02) {

	for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	    reset_sample_errors(n_samples) ;

	    float phase_shift = lr_phase_shift[lr] ;
	    get_value_model(n_samples, raw_hsv_coefs[lr].v_coefs, pData, lr, phase_shift) ;
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
		}
		if(sse < best_sse_this) {
		    best_sse_this = sse ;
		}
	    }

	    fprintf(stderr, "GRID fit_full_HSV  vfactor,hcurv(%f,%f) = %f\n", pData->angle_factor,pData->horizon_curvature,best_sse_this) ;

	}
    }

    // final fit with best horizon_py and horizon_curvature
    pData->horizon_py = best_horizon_py ;
    pData->horizon_curvature = best_horizon_curvature ;


    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	float phase_shift = lr_phase_shift[lr] ;
	get_value_model(n_samples, raw_hsv_coefs[lr].v_coefs, pData, lr, phase_shift) ;
	get_saturation_model(n_samples, raw_hsv_coefs[lr].s_coefs, pData, lr, phase_shift) ;

	float sse = fit_full_HSV_model(n_samples,lr,0,pData) ; // last time through, calculate errors on points

	fprintf(stderr, "*******************  CALCULATED PHASE SHIFT, angle_factor are %3.0f %3.1f*************************\n", phase_shift,pData_fit->angle_factor) ;

	fprintf(stderr, "FINAL fit_full_HSV[%d of %d]:  vfactor,hcurv from is (%f,%f) = %f\n",
			lr,pData->full_hsv_model_level,pData->angle_factor,pData->horizon_curvature,sse) ;

	fprintf(stderr, "raw_hsv_coefs H = half:%f spread:%f speed:%f center:%f\n",
	        raw_hsv_coefs[lr].h_coefs[0], raw_hsv_coefs[lr].h_coefs[1],
	        raw_hsv_coefs[lr].h_coefs[2], raw_hsv_coefs[lr].h_coefs[3]
		) ;

	fprintf(stderr, "raw_hsv_coefs S = %f %f %f %f %f\n\t%f %f\n",
	        raw_hsv_coefs[lr].s_coefs[0], raw_hsv_coefs[lr].s_coefs[1], raw_hsv_coefs[lr].s_coefs[2],
		raw_hsv_coefs[lr].s_coefs[3], raw_hsv_coefs[lr].s_coefs[4],
		raw_hsv_coefs[lr].s_coefs[5], raw_hsv_coefs[lr].s_coefs[6]
		) ;

	fprintf(stderr, "raw_hsv_coefs V = %f %f %f %f %f %f\n",
		raw_hsv_coefs[lr].v_coefs[0], raw_hsv_coefs[lr].v_coefs[1],
		raw_hsv_coefs[lr].v_coefs[2], raw_hsv_coefs[lr].v_coefs[3],
		raw_hsv_coefs[lr].v_coefs[4], raw_hsv_coefs[lr].v_coefs[5]
		) ;

	fprintf(stderr, "\n") ;
    }
}

int sample_sky_points(int n_per_column, int n_columns,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky, SKYFILL_DATA_t *pData)
{
    int x, y ;
    float mean_column_length ;
    float sum_column_lengths=0. ;
    int n_columns_in_sum=0 ;

    // for this processing, do not consider localized sun image.
    int save_uses_CIE_model = uses_CIE_model ;
    uses_CIE_model = 0 ;

    V_sun= image_relative_pixel_to_V3_noclip(sun_x_angle2px(pData->sun_x), pData->sun_py) ;

    for(x=0 ; x < IMAGE_WIDTH ; x++) {
	if(pData->column_mask[x] == 1) continue ;
	if(pData->column_sample_mask[x] == 1) continue ;

	 sum_column_lengths += pData->end_of_sky[x] - pData->start_of_sky[x] + 1 ;
	 n_columns_in_sum++ ;
    }

    mean_column_length = sum_column_lengths/(float)n_columns_in_sum ;

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
#ifdef PRINTOUT_HUE_BY_ROW
    int n_max_rows= IMAGE_HEIGHT/dy ;

    float *sum_h_by_row = (float *)calloc(n_max_rows, sizeof(float)) ;
    int *n_h_by_row = (int *)calloc(n_max_rows, sizeof(int)) ;

    for(int i = 0 ; i < n_max_rows ; i++) {
	sum_h_by_row[i] = 0. ;
	n_h_by_row[i] = 0 ;
    }
#endif


    for(x=0 ; x < IMAGE_WIDTH ; x += dx) {
	uint16_t r,g,b ;

	if(pData->column_mask[x] == 1) continue ; // don't attempt a mean sky sample inside a masked area
	if(pData->column_sample_mask[x] == 1) continue ;
	if(pData->fix_sky_hue && x>=pData->min_sky_hue_mask && x <= pData->max_sky_hue_mask) continue ;

#ifdef PRINTOUT_HUE_BY_ROW
	int row=0 ;
#endif

	// want to get close to end of sky since the sky changes most near end of sky
	for(y = end_of_sky[x] - 6 ; y > start_of_sky[x]+3 ; y -= dy) {

/*  	for(y=0 ; y < IMAGE_HEIGHT ; y += dy) {  */
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

	    samples[n_samples].x = x ;
	    samples[n_samples].y = y ;
	    samples[n_samples].px = px ;
	    samples[n_samples].py = py ;
	    samples[n_samples].r = r ;
	    samples[n_samples].g = g ;
	    samples[n_samples].b = b ;
	    samples[n_samples].abs_error[0] = 1. ;
	    samples[n_samples].abs_error[1] = 1. ;

	    samples[n_samples].h = h ;
	    samples[n_samples].s = s ;
	    samples[n_samples].v = v ;

	    float wgt = (1.-s)*PX_WGT(px) ;

	    sum_hue += h*wgt ;
	    sum_hue_wgts += wgt ;
	    sum_sat += s ;
	    sum_sat_wgts += 1.f ;
	    sum_val += v*wgt ;
	    sum_val_wgts += wgt ;

	    n_samples++ ;

#ifdef PRINTOUT_HUE_BY_ROW
	    sum_h_by_row[row] += h ;
	    n_h_by_row[row]++ ;
	    row++ ; // end of y for loop now
#endif
	}

	//fprintf(stderr, "column %d, h=%f, s=%f, v=%f\n", x, sum_h/sum_n, sum_s/sum_n, sum_v/sum_n) ;
    }

    fprintf(stderr, "\n") ;
#ifdef PRINTOUT_HUE_BY_ROW
    for(int i = 0 ; i < n_max_rows ; i++) {
	if(n_h_by_row[i] > 0) {
	    float mean = sum_h_by_row[i]/(float)n_h_by_row[i] ;
	    fprintf(stderr, "INFO:  Mean sky hue at y=%d is %f, n=%d\n", i*dy, mean, n_h_by_row[i]) ;
	}
    }
    fprintf(stderr, "\n") ;

    free(sum_h_by_row) ;
    free(n_h_by_row) ;
#endif

    pData->hue_sky = sum_hue/sum_hue_wgts ;
    pData->sat_sky = sum_sat/sum_sat_wgts ;
    pData->val_sky = sum_val/sum_val_wgts ;

    fit_sky_model(pData, n_samples) ;

    // is the top of the sky too dark, and the full value model used?, if so try the
    // reduced value model

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

/*  	    printf("CTOS: px:%f x:%d t:%f m:%f\n", px, x, vhat_top, vhat_mid) ;  */

	    if(vhat_top < 100./255 && vhat_top < vhat_mid) n_failed++ ;
	}

	if(n_failed > 10) {
	    fprintf(stderr, "\n*******************************************************************************\n") ;
	    fprintf(stderr, "********** TOP OF SKY IS TOO DARK, changing to reduced value model  ***********\n") ;
	    fprintf(stderr, "*******************************************************************************\n\n") ;

/*  	    pData->val_model_full = 0 ;  */
/*  	    fit_sky_model(pData, n_samples) ;  */
	}
    }

    FILE *fp = fopen("samples.dat", "w") ;
    for(int i = 0 ; i < n_samples ; i++) {
	float h_hat, s_hat, v_hat ;
	// estimate hsv, will used blended model if requested
	uint16_t r_hat, g_hat, b_hat ;
	predict_sky_hsv(samples[i].px, samples[i].py, &h_hat, &s_hat, &v_hat) ;
	hsv2rgb16(h_hat,s_hat,v_hat,&r_hat,&g_hat,&b_hat) ;
	float angle = pData->FOV_horizontal*samples[i].px ;

	// temporarily replace v with intensity
/*  	samples[i].v = (float)(samples[i].r+samples[i].b+samples[i].g)/3./(float)MAX16 ;  */
/*  	v_hat = (float)(r_hat+b_hat+g_hat)/3./(float)MAX16 ;  */

	fprintf(fp, "%6.4f %6.4f", samples[i].px, samples[i].py+samples[i].px/50.) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", samples[i].h, samples[i].s, samples[i].v) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", h_hat, s_hat, v_hat) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", h_hat-samples[i].h, s_hat-samples[i].s, v_hat-samples[i].v) ;
	fprintf(fp, " %5.1f %5.2lf\n", angle, log(1.+log(samples[i].py+.0001))) ;

    }

    fclose(fp) ;
    uses_CIE_model = save_uses_CIE_model ;

/*      for(int i = 0 ; i < n_samples ; i++) {  */
/*  	float hhat, shat, vhat ;  */
/*  	predict_sky_hsv(samples[i].px, samples[i].py, &hhat, &shat, &vhat) ;  */
/*  	float angle = pData->FOV_horizontal*samples[i].px ;  */
/*    */
/*  	fprintf(stderr, "%6.4f %6.4f %5.3f\n", samples[i].px, samples[i].py+samples[i].px/50., vhat) ;  */
/*    */
/*      }  */
/*    */
/*      fprintf(stderr, "sample points, sky hsv=%f,%f,%f\n", pData->hue_sky, pData->sat_sky, pData->val_sky) ;  */
/*      exit(1) ;  */

    return n_samples ;
}
