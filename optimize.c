#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "tiffio.h"


// in case Windows MSVC, or 
#define _USE_MATH_DEFINES
#include <math.h>

#include "skyfill_tif.h"
#include "amoeba_06.h"
#include "string.h"
#include "mstat.h"
#include "find_sky.h"
#include "repair_sky.h"
#include "colorspace_conversions.h"
#include "sample_and_fit_sky_model.h"
#include "feather_factor.h"
#include "estimate_sky.h"
#include "optimize.h"

// these are the "standard" CIE parameters for various types of skies.
#define N_CIE_STD_PARMS 16
struct CIE_STD_PARMS cie_parms[N_CIE_STD_PARMS] = {
	{ 1,1,   4.0, -0.70,  0, -1.0, 0.00, "CIE standard overcast sky, alternative form.  Steep lum gradation towars zenith, azimuthal uniformity"},
	{ 2,1,   4.0, -0.70,  2, -1.5, 0.15, "Overcast, with steep lum gradation and slight brightening towards sun"},
	{ 3,2,   1.1, -0.80,  0, -1.0, 0.00, "Overcast, moderately graded with azimuthal uniformity"},
	{ 4,2,   1.1, -0.80,  2, -1.5, 0.15, "Overcast, moderately graded with slight brightening towards sun"},
	{ 5,3,   0.0, -1.00,  0, -1.0, 0.00, "Sky of uniform luminance"},
	{ 6,3,   0.0, -1.00,  2, -1.5, 0.15, "Partly cloudy sky, no gradation towards zenith, slight brigtening towards sun"},
	{ 7,3,   0.0, -1.00,  5, -2.5, 0.35, "Partly cloudy sky, no gradation towards zenith, brighter circumsolar region"},
	{ 8,3,   0.0, -1.00, 10, -3.0, 0.45, "Partly cloudy sky, no gradation towards zenith, distinct solar corona"},
	{ 9,4,  -1.0, -0.55,  2, -1.5, 0.15, "Partly cloudy, with the obscured sun"},
	{10,4,  -1.0, -0.55,  5, -2.5, 0.30, "Partly cloudy, with brighter circumsolar region"},
	{11,4,  -1.0, -0.55, 10, -3.0, 0.45, "White-blue sky with distinct solar corona"},
	{12,5,  -1.0, -0.32, 10, -3.0, 0.45, "CIE Standard Clear Sky, low illuminance"},
	{13,5,  -1.0, -0.32, 16, -3.0, 0.30, "CIE Standard Clear Sky, polluted atmosphere"},
	{14,6,  -1.0, -0.15, 16, -3.0, 0.30, "Cloudless turbid sky with broad solar corona"},
	{15,6,  -1.0, -0.15, 24, -2.8, 0.15, "White-blue turbid sky with broad solar corona"},
	{16,3,  -1.0, -0.30,  0, -1.0, 0.00, "Sky of uniform luminance, maximum horizon glow"},
	} ;

void set_xy_constraints(struct OPT_PARM opt_parms[])
{
    // sun_x, MAX,MIN in case a tighter constraint already
    opt_parms[0].lo = MAX(-180.f, opt_parms[0].lo) ;
    opt_parms[0].lo = MIN(180.f, opt_parms[0].lo) ;

    opt_parms[0].hi = MIN(180.f, opt_parms[0].hi) ;
    opt_parms[0].hi = MAX(-180.f, opt_parms[0].hi) ;
}


int optimize_type = OPTIMIZE_V_from_HSV ;
int calln=0 ;
int autoset_lumf_flag=1 ;

float samples_error(void)
{
    int i ;
    float sum_err_sq = 0. ;
    //optimize_type = FULL_RGB ;

    set_FOV_factor() ;


    V_sun= image_relative_pixel_to_V3_noclip(sun_x_angle2px(pData_fit->sun_x), pData_fit->sun_py) ;

    pData_fit->maximum_CIE_vhat = find_maximum_CIE_vhat() ;

    if(isinf(pData_fit->maximum_CIE_vhat)) {
	fprintf(stderr, "maximum_CIE_vhat is NAN in grid optimize: sun_x:%f sun_py%f\n", pData_fit->sun_x, pData_fit->sun_py) ;
	exit(1) ;
    }

    float sum_v=0 ;
    float sum_v_hat=0 ;

    // recalculate value correction using sky samples and
    // current prediction coefficients
/*      sample_based_v_correction = 1.f ;  */

    for(i = 0 ; i < n_samples_to_optimize ; i++) {
	float h,s,v ;
	predict_sky_hsv(samples[i].px, samples[i].py, &h, &s, &v) ;
	samples[i].v_hat = v ;
	sum_v += samples[i].v ;
	sum_v_hat += samples[i].v_hat ;
    }

/*      sample_based_v_correction = sum_v/sum_v_hat ;  */

    // XXXX
    // clamp this -- gets way to high if the sun is in the image and the CIE parms allow for a sun
/*      if(sample_based_v_correction > 4.0) sample_based_v_correction = 4.0 ;   */

    if(optimize_type == OPTIMIZE_V_from_HSV) {

	float sum_vhat =0. ;

	for(i = 0 ; i < n_samples_to_optimize ; i++) {
	    float v_hat = samples[i].v_hat ;
	    float err_v = (samples[i].v - v_hat) ;

	    if(isinf(err_v)) {
		fprintf(stderr, "err_v is inf in grid optimize: v:%f px:%f py%f\n", v_hat, samples[i].px, samples[i].py) ;
		exit(1) ;
	    }

	    sum_err_sq += err_v*err_v*PX_WGT(samples[i].px) ; // a test for local weighting
	    //sum_err_sq += err_v*err_v ;
	    sum_vhat += v_hat ;
	}

	if(!(calln % 1000)) fprintf(stderr, "mean vhat=%f\n", sum_vhat/(float)n_samples_to_optimize) ;
    } else {
	for(i = 0 ; i < n_samples_to_optimize ; i++) {
	    uint16_t rhat,ghat,bhat ;
	    predict_sky_color(samples[i].px, samples[i].py, &rhat, &ghat, &bhat) ;
	    float err_r = (samples[i].r - rhat) ;
	    float err_g = (samples[i].g - ghat) ;
	    float err_b = (samples[i].b - bhat) ;
	    sum_err_sq += (err_r*err_r + err_g*err_g + err_b*err_b)*PX_WGT(samples[i].px) ; // a test for local weighting
	    //sum_err_sq += (err_r*err_r + err_g*err_g + err_b*err_b) ;
	}
    }

    calln++ ;


    return sum_err_sq ;
}

void set_cie_parms(int CIE_index)
{
    pData_fit->perez_A = cie_parms[CIE_index].A ;
    pData_fit->perez_B = cie_parms[CIE_index].B ;
    pData_fit->perez_C = cie_parms[CIE_index].C ;
    pData_fit->perez_D = cie_parms[CIE_index].D ;
    pData_fit->perez_E = cie_parms[CIE_index].E ;
}

float optimize_grid_function(int CIE_index)
{

    set_cie_parms(CIE_index) ;

    return samples_error() ;

}

float optimize_grid(int n_samples, int verbose, int allowed_sky_type, int CIE_sky_index,struct OPT_PARM opt_parms[])
{
    int i ;
    float best_feval = 1.e30 ;
    int best_i = -1 ;
    float feval ;
    optimize_type = OPTIMIZE_V_from_HSV ;

    autoset_lumf_flag=1 ;
    for(i=0 ; i < MAX_OPT_PARMS ; i++) {
	if(opt_parms[i].grid_optimize_flag == 0) {
	    if(!strcmp(opt_parms[i].name, "sun_lum")) {
		autoset_lumf_flag=0 ;
	    }
	}
    }

    n_samples_to_optimize=n_samples ;

    if(CIE_sky_index > 0) {
	i=CIE_sky_index-1 ;

	feval = optimize_grid_function(i) ;
	if(feval < best_feval) {
	    best_i = i ;
	    best_feval = feval ;
	}
    } else {
	for(i = 0 ; i < N_CIE_STD_PARMS ; i++) {
	    int is_uniform_sky=0 ;
	    if(cie_parms[i].C < .01 && cie_parms[i].D < .01) is_uniform_sky=1;

	    if(is_uniform_sky && allowed_sky_type == NONUNIFORM_CIE_SKIES) continue ;
	    if(!is_uniform_sky && allowed_sky_type == UNIFORM_CIE_SKIES) continue ;

	    feval = optimize_grid_function(i) ;
	    if(feval < best_feval) {
		best_i = i ;
		best_feval = feval ;
	    }
	}
	feval = optimize_grid_function(best_i) ;
    }


    if(verbose) printf("GRID feval:%g CIE index %d, %s\n", feval, cie_parms[best_i].type, cie_parms[best_i].description) ; ;

    if(verbose) dump_parameters_to_stdout() ;

    return feval ;
}

struct OPT_PARM *amoeba_opt_parms ;

float smart_optimize_function(float *pParams)
{
    int i ;

    int j=1 ;

    for(i=0 ; i < MAX_OPT_PARMS ; i++) {
	if(amoeba_opt_parms[i].optimize_flag == 1) {
	    *(amoeba_opt_parms[i].variable_address) = pParams[j] ;
	    j++ ;
	}
    }

    return samples_error() ;
}

void smart_optimize(int n_samples, SKYFILL_DATA_t *pData,struct OPT_PARM opt_parms[])
{
    float guess[MAX_OPT_PARMS+1] = {0,pData->sun_x,pData->sun_py,pData->sun_lum,pData->horizon_py,pData->perez_A,
				pData->perez_B,pData->perez_C,pData->perez_D,pData->perez_E,pData->perez_F,pData->perez_G,pData->FOV_horizontal} ;
    float guess_delta[MAX_OPT_PARMS+1], c_lo[MAX_OPT_PARMS+1], c_hi[MAX_OPT_PARMS+1] ;
    int i ;
    int nparms_to_optimize=0 ;

    // set the global pointer
    pData_fit = pData ;

    n_samples_to_optimize=n_samples ;

    int j=1 ;

    for(i=0 ; i < MAX_OPT_PARMS ; i++) {
	if(opt_parms[i].optimize_flag == 1) {
	    guess[j] = *opt_parms[i].variable_address ;
	    c_lo[j] = opt_parms[i].lo ;
	    c_hi[j] = opt_parms[i].hi ;
	    if(fabs(guess[j]) > 1.e-5) {
		guess_delta[j] = fabs(guess[j])/10. ;
	    } else {
		guess_delta[j] = fabs(opt_parms[i].hi-opt_parms[i].lo)/10. ;
	    }

	    if(!strcmp(opt_parms[i].name, "sun_x")) {
		// special case for sun_x_angle
		guess_delta[j] = 5. ;
		c_lo[j] = guess[j]-20. ;
		c_hi[j] = guess[j]+20. ;
	    }

	    if(!strcmp(opt_parms[i].name, "horizon_py")) {
		// special case for horizon_py
		guess_delta[j] = .01 ;
		c_lo[j] = 0. ;
		c_hi[j] = ((float)(IMAGE_HEIGHT)-(float)pData->max_end_of_sky)/(float)IMAGE_HEIGHT ;
	    }

	    j++ ;
	    printf("OPTIMIZE->%s\n", opt_parms[i].name) ;
	    nparms_to_optimize++ ;
	}
    }

    amoeba_opt_parms = opt_parms ;

    float feval = smart_optimize_function(guess) ;

    printf("Initial guess feval:%g\n", feval) ;

    dump_parameters_to_stdout() ;

    feval = AmoebaFit(smart_optimize_function,1.e-7,1000,nparms_to_optimize,guess,guess_delta,c_lo,c_hi,0) ;

    feval = smart_optimize_function(guess) ;

    printf("Optimized feval:%g\n", feval) ;

    dump_parameters_to_stdout() ;

}

void grid_search_sun_xy(SKYFILL_DATA_t *pData, struct OPT_PARM opt_parms[], int n_samples)
{
    // default case, do a grid search for sun position
    float best_s_x = 0. ;
    float best_s_y = 1. ;
    float best_err=FMAX ;
    int n_tot = 0 ;
    float sun_dx = 1.e30 ;
    float sun_dy = 1.e30 ;

    opt_parms[0].lo = -180. ;
    opt_parms[0].hi = 180. ;
    opt_parms[1].lo = 0. ;
    opt_parms[1].hi = 1.5 ;

    fprintf(stderr, "1  Sunx,suny,horz_y=%f,%f,%f\n", pData->sun_x, pData->sun_py, pData->horizon_py) ;

    if(opt_parms[0].grid_optimize_flag == 1) {
	pData->sun_x=opt_parms[0].lo ;
	sun_dx = 10. ;
	fprintf(stderr, "Grid optimize sun_x\n") ;
    }
    if(opt_parms[1].grid_optimize_flag == 1) {
	pData->sun_py=opt_parms[1].lo ;
	sun_dy = .5 ;
	fprintf(stderr, "Grid optimize sun_py\n") ;
    }

    float start_sun_x = pData->sun_x ;
    float start_sun_py = pData->sun_py ;

    fprintf(stderr, "Grid optimize sun_x from %f to %f by %f\n", start_sun_x, opt_parms[0].hi, sun_dx) ;
    fprintf(stderr, "Grid optimize sun_py from %f to %f by %f\n", start_sun_py, opt_parms[1].hi, sun_dy) ;

#define FTOL .0001
    for(pData->sun_x=start_sun_x; pData->sun_x < opt_parms[0].hi+FTOL ; pData->sun_x += sun_dx) {
	for(pData->sun_py=start_sun_py; pData->sun_py < opt_parms[1].hi+FTOL ; pData->sun_py += sun_dy) {
	    //fprintf(stderr, "grid evaluations %d, sx:%f sy:%f hy:%f\n", n_tot, sun_x, sun_py, horizon_py) ;
	    n_tot++ ;
	}
    }

    int i_grid_search_no=0 ;
    fprintf(stderr, "%d grid evaluations will be done\n", n_tot) ;

    for(pData->sun_x=start_sun_x; pData->sun_x < opt_parms[0].hi+FTOL ; pData->sun_x += sun_dx) {
	for(pData->sun_py=start_sun_py; pData->sun_py < opt_parms[1].hi+FTOL ; pData->sun_py += sun_dy) {

		float err = optimize_grid(n_samples,0,pData->allowed_sky_type, pData->CIE_sky_index,opt_parms) ;
		i_grid_search_no++ ;
		printf("Completed %d of %d grid evaluations for sun position\n", i_grid_search_no, n_tot) ;

		if(err < best_err) {
		    best_err = err ;
		    best_s_x = pData->sun_x ;
		    best_s_y = pData->sun_py ;
		}

	}
    }

    /* reset for best values */
    pData->sun_x = best_s_x ;
    pData->sun_py = best_s_y ;
    optimize_grid(n_samples,1,pData->allowed_sky_type, pData->CIE_sky_index,opt_parms) ;
}

void fit_sky_modelv2(int force_grid_optimization, int n_custom_optimizations, int optimization_var[10][MAX_OPT_PARMS+1], int n_samples, SKYFILL_DATA_t *pData,struct OPT_PARM opt_parms[])
{
    pData->model_is_being_fit=1 ;

    // set the global pointer
    pData_fit = pData ;

    if(force_grid_optimization || n_custom_optimizations == 0) {
	if(uses_CIE_model)
	    grid_search_sun_xy(pData, opt_parms, n_samples) ;


	if(n_custom_optimizations == 0) {

    //-O sx hy B C sl sg -O sx sy sl sg A D E

	    char *opt0[6] ;
	    char *opt1[7] ;

	    // note -- order of elements is important to compare final optimations
	    char *opt0_a[4] = {"hy","sl","sd","sx"} ;
	    char *opt1_a[5] = {"sl","sd","sx","sy","D"} ;

	    int i,j ;

	    int np0, np1 ;

	    np0 = uses_CIE_model ? 4 : 3 ;
	    np1 = uses_CIE_model ? 5 : 2 ;
	    for(j = 0 ; j < np0 ; j++) {
		opt0[j] = opt0_a[j] ;
	    }
	    for(j = 0 ; j < np1 ; j++) {
		opt1[j] = opt1_a[j] ;
	    }


	    int n_vars=0 ;
	    optimization_var[n_custom_optimizations][0] = -1 ;

	    for(j = 0 ; j < np0 ; j++) {
		for(i = 0 ; i < MAX_OPT_PARMS ; i++) {
		    if(!strcmp(opt0[j], opt_parms[i].abreviation)) {
			if(opt_parms[i].optimize_flag == 1) {
			    optimization_var[n_custom_optimizations][n_vars] = i ;
			    n_vars++ ;
			    optimization_var[n_custom_optimizations][n_vars] = -1 ;
			}
		    }
		}
	    }
	    n_custom_optimizations++ ;

	    n_vars=0 ;
	    optimization_var[n_custom_optimizations][0] = -1 ;

	    for(j = 0 ; j < np1 ; j++) {
		for(i = 0 ; i < MAX_OPT_PARMS ; i++) {
		    if(!strcmp(opt1[j], opt_parms[i].abreviation)) {
			if(opt_parms[i].optimize_flag == 1) {
			    optimization_var[n_custom_optimizations][n_vars] = i ;
			    n_vars++ ;
			    optimization_var[n_custom_optimizations][n_vars] = -1 ;
			}
		    }
		}
	    }
	    n_custom_optimizations++ ;
	    optimization_var[n_custom_optimizations][0] = -1 ;

	    // are the two optimizations the same after they have been filtered for command line settings?
	    // if so reduce to one optimization

	    fprintf(stderr, "ARE THE AUTOMATIC OPTIMIZATIONS THE SAME after filtering?\n") ;

	    int same=1 ;
	    for(int var = 0 ; var < MAX_OPT_PARMS ; var++) {
		if(optimization_var[0][var] != optimization_var[1][var]) {
		    same=0 ;
		}

		fprintf(stderr, "\t%d(%d) => %d %d\n", var, same, optimization_var[0][var], optimization_var[1][var]) ;

		if(optimization_var[0][var] < 0 || optimization_var[1][var] < 0) {
		    break ;
		}

	    }

	    if(same == 1) {
		// eliminate the redundant optimization
		fprintf(stderr, " -- YES, delete the second optimization.\n") ;
		n_custom_optimizations-- ;
		optimization_var[n_custom_optimizations][0] = -1 ;
	    }

	    fprintf(stderr, " -- perform %d prebuilt optimizations.\n\n", n_custom_optimizations) ;
	}

    }

    {
	// perform custom optimizations either from defaults or specified on the command line
	int i,o ;

	for(o = 0 ; o < n_custom_optimizations ; o++) {
	    optimize_type = OPTIMIZE_V_from_HSV ;

	    for(i = 0 ; i < MAX_OPT_PARMS ; i++) {
		opt_parms[i].optimize_flag=0 ;
	    }

	    for(i = 0 ; i < MAX_OPT_PARMS ; i++) {
		int j = optimization_var[o][i] ;
		if( j < 0) break ;
		opt_parms[j].optimize_flag=1 ;
		if(!strcmp(opt_parms[j].name, "sat_depth")) optimize_type = OPTIMIZE_FULL_RGB ;
		if(!strcmp(opt_parms[j].name, "hue_sky")) optimize_type = OPTIMIZE_FULL_RGB ;
		if(!strcmp(opt_parms[j].name, "hue_horizon")) optimize_type = OPTIMIZE_FULL_RGB ;
	    }

	    smart_optimize(n_samples, pData, opt_parms) ;

	}
    }
}
