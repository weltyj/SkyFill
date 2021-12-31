#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
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

#define FULL_RAW_HSV_MODEL


// coefs to model sky HSV from input file sky area
// model is H|S|V = B0 + B1*px + B2*py + B3*px*py ;
struct HSV_model_coefs raw_hsv_coefs[2] ;

float sat_optimize_function(float coefs[])
{
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
	float shat = coefs[1]/(1.+exp(-(py_hc-coefs[2])*coefs[3]))+coefs[4]+coefs[5]*py_hc ;
	float err = samples[i].s - shat ;
	float wgt = (*wgt_func)(samples[i].px, 1.) ;
	sse += err*err*wgt ;
    }

/*      fprintf(stderr, "SOF: SSE(%f) = %f,%f,%f,%f\n", sse, coefs[1],coefs[2],coefs[3],coefs[4]) ;  */

    return sse ;
}

void get_saturation_model(int n_samples, double coefs[], SKYFILL_DATA_t *pData, int location_index)
{
    float max_sat=0. ;
    n_samples_to_optimize = n_samples ;
    g_full_hsv_model_level = pData->full_hsv_model_level ;
    g_location_index = location_index ;

    for(int i = 0 ; i < n_samples ; i++) {
	if(max_sat < samples[i].s) max_sat = samples[i].s ;
    }

    // model form is B0/(1.+exp(-(py_hc-B1)*B2))+B3+B4*py_hc ;

    float guess[6] = {0,(max_sat+1.)/2.,0.1,10.,-.1,.01} ;
    float guess_delta[6] = {0,.1,.05,1,.05,.001} ;
    float c_lo[6] = {0,max_sat,-.5,.1,-1.,-.1} ;
    float c_hi[6] = {0,     1., .5,100.,1.,.1};

    float feval = sat_optimize_function(guess) ;

    feval = AmoebaFit(sat_optimize_function,1.e-4,200,5,guess,guess_delta,c_lo,c_hi,0) ;

    feval = sat_optimize_function(guess) ;

    coefs[0] = guess[1] ;
    coefs[1] = guess[2] ;
    coefs[2] = guess[3] ;
    coefs[3] = guess[4] ;
    coefs[4] = guess[5] ;
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
    double coefs[5] ;
    struct mstat m_fit ;

    // Hue
    m_fit = init_reg(3) ;
    for(int i = 0 ; i < n_samples ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	double xv[3] = {1.,samples[i].px,log(1.+py_hc)} ;
	double y = samples[i].h ;
	sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, samples[i].abs_h_error[location_index])) ;
    }

    estimate_reg(&m_fit, coefs) ;
    p->raw_HSV_coefs[0][0] = coefs[0] ;
    p->raw_HSV_coefs[0][1] = coefs[1] ;
    p->raw_HSV_coefs[0][2] = coefs[2] ;
    p->raw_HSV_coefs[0][3] = 0. ;

    // Sat
    // No longer done in this module, have switch to a nonlinear function

    // Val
    m_fit = init_reg(4) ;
    for(int i = 0 ; i < n_samples ; i++) {
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	float horz_angle = pData->FOV_horizontal*samples[i].px*M_PI/180. ;
/*  	double xv[4] = {1.,  samples[i].px,py_hc,py_hc*py_hc} ;  */
	double xv[4] = {1.,  cos(2.*horz_angle),sin(2.*horz_angle),py_hc} ;
	double y = samples[i].v ;
	sum_reg(&m_fit, xv, y, (*wgt_func)(samples[i].px, samples[i].abs_v_error[location_index])) ;
    }

    estimate_reg(&m_fit, coefs) ;
    p->raw_HSV_coefs[2][0] = coefs[0] ;
    p->raw_HSV_coefs[2][1] = coefs[1] ;
    p->raw_HSV_coefs[2][2] = coefs[2] ;
    p->raw_HSV_coefs[2][3] = coefs[3] ;

    //compute combined SSE from each r,g,b component
    float sse_r = 0. ;
    float sse_g = 0. ;
    float sse_b = 0. ;

    for(int i = 0 ; i < n_samples ; i++) {
	float h,s,v ;
	float py_hc = py_adjusted_for_horizon_curvature(samples[i].px, samples[i].py) ;
	h = h_from_pxpy(samples[i].px, py_hc, pData->full_hsv_model_level) ;
	s = s_from_pxpy(samples[i].px, py_hc, pData->full_hsv_model_level, pData->horizon_py) ;
	v = v_from_pxpy(samples[i].px, py_hc, pData->full_hsv_model_level) ;

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

    // try to get 1000 sample points, each will be a mean of a 5x5 area around the actual point
    // do this with a rectangular grid over the image area, refining the grid to the sky area only

    n_columns = 50 ;
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
    int n_max_rows= IMAGE_HEIGHT/dy ;

    float *sum_h_by_row = (float *)calloc(n_max_rows, sizeof(float)) ;
    int *n_h_by_row = (int *)calloc(n_max_rows, sizeof(int)) ;

    for(int i = 0 ; i < n_max_rows ; i++) {
	sum_h_by_row[i] = 0. ;
	n_h_by_row[i] = 0 ;
    }


    for(x=0 ; x < IMAGE_WIDTH ; x += dx) {
	uint16_t r,g,b ;

	if(pData->column_mask[x] == 1) continue ; // don't attempt a mean sky sample inside a masked area
	if(pData->column_sample_mask[x] == 1) continue ;
	if(pData->fix_sky_hue && x>=pData->min_sky_hue_mask && x <= pData->max_sky_hue_mask) continue ;

	int row=0 ;

	for(y=0 ; y < IMAGE_HEIGHT ; y += dy) {

	    if(y < pData->start_of_sky[x] || y > pData->end_of_sky[x]) continue ;

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
			tif_get3c(image,x,y,r,g,b) ;
			float h,s,v ;
			rgb2hsv16(r,g,b,&h,&s,&v) ;

/*  			if(s > 0.25) {  */
/*  			    // only use sample points with saturation > 55%  */

			    sum_r += r ;
			    sum_g += g ;
			    sum_b += b ;
			    n++ ;

/*  			}  */

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
	    sum_h_by_row[row] += h ;
	    n_h_by_row[row]++ ;

	    float wgt = (1.-s)*PX_WGT(px) ;

	    sum_hue += h*wgt ;
	    sum_hue_wgts += wgt ;
	    sum_sat += s ;
	    sum_sat_wgts += 1.f ;
	    sum_val += v*wgt ;
	    sum_val_wgts += wgt ;

	    n_samples++ ;

	    row++ ; // end of y for loop now
	}

	//fprintf(stderr, "column %d, h=%f, s=%f, v=%f\n", x, sum_h/sum_n, sum_s/sum_n, sum_v/sum_n) ;
    }

    fprintf(stderr, "\n") ;
    for(int i = 0 ; i < n_max_rows ; i++) {
	if(n_h_by_row[i] > 0) {
	    float mean = sum_h_by_row[i]/(float)n_h_by_row[i] ;
	    fprintf(stderr, "INFO:  Mean sky hue at y=%d is %f, n=%d\n", i*dy, mean, n_h_by_row[i]) ;
	}
    }
    fprintf(stderr, "\n") ;

    free(sum_h_by_row) ;
    free(n_h_by_row) ;

    pData->hue_sky = sum_hue/sum_hue_wgts ;
    pData->sat_sky = sum_sat/sum_sat_wgts ;
    pData->val_sky = sum_val/sum_val_wgts ;

    float best_sse = 1.e30;
    float best_horizon_py = pData->horizon_py ;
    float best_horizon_curvature = pData->horizon_py ;

    double coefs[5] ;

#define START_S_COEF 1.
    for(int lr = 0 ; lr < pData->full_hsv_model_level ; lr++) {
	for(pData->horizon_curvature=-.1 ; pData->horizon_curvature < .101 ; pData->horizon_curvature += .02) {
	    get_saturation_model(n_samples, coefs, pData, lr) ;

	    raw_hsv_coefs[lr].raw_HSV_coefs[1][0] = coefs[0] ;
	    raw_hsv_coefs[lr].raw_HSV_coefs[1][1] = coefs[1] ;
	    raw_hsv_coefs[lr].raw_HSV_coefs[1][2] = coefs[2] ;
	    raw_hsv_coefs[lr].raw_HSV_coefs[1][3] = coefs[3] ;
	    raw_hsv_coefs[lr].raw_HSV_coefs[1][4] = coefs[4] ;

	    if(! pData->horizon_was_set) {
		for(pData->horizon_py = pData->horizon_py-.05 ; pData->horizon_py > .25 ; pData->horizon_py -= .05) {
		    reset_sample_errors(n_samples) ;
		    for(int iteration=0 ; iteration < 8 ; iteration++) {
			// iteratively reweighted least squares
			float sse = fit_full_HSV_model(n_samples,lr,1,pData) ;
			if(sse < best_sse) {
			    best_sse = sse ;
			    best_horizon_py = pData->horizon_py ;
			    best_horizon_curvature = pData->horizon_curvature ;
			}
		    }
		}
	    } else {
		reset_sample_errors(n_samples) ;
		for(int iteration=0 ; iteration < 8 ; iteration++) {
		    // iteratively reweighted least squares
		    float sse = fit_full_HSV_model(n_samples,lr,1,pData) ;
		    if(sse < best_sse) {
			best_sse = sse ;
			best_horizon_curvature = pData->horizon_curvature ;
		    }
		}
	    }
	}
	// final fit
	pData->horizon_py = best_horizon_py ;
	pData->horizon_curvature = best_horizon_curvature ;

	get_saturation_model(n_samples, coefs, pData, lr) ;

	raw_hsv_coefs[lr].raw_HSV_coefs[1][0] = coefs[0] ;
	raw_hsv_coefs[lr].raw_HSV_coefs[1][1] = coefs[1] ;
	raw_hsv_coefs[lr].raw_HSV_coefs[1][2] = coefs[2] ;
	raw_hsv_coefs[lr].raw_HSV_coefs[1][3] = coefs[3] ;
	raw_hsv_coefs[lr].raw_HSV_coefs[1][4] = coefs[4] ;

	float sse = fit_full_HSV_model(n_samples,lr,0,pData) ; // last time through, calculate errors on points

	fprintf(stderr, "FINAL fit_full_HSV[%d of %d]:  hpy,Icoef SSE from is (%f,%f) = %f\n", lr,pData->full_hsv_model_level,pData->horizon_py,raw_hsv_coefs[lr].S_inv_coef,sse) ;
	printf("raw_hsv_coefs S = %f %f %f %f %f\n", raw_hsv_coefs[lr].raw_HSV_coefs[1][0], raw_hsv_coefs[lr].raw_HSV_coefs[1][1],
		raw_hsv_coefs[lr].raw_HSV_coefs[1][2], raw_hsv_coefs[lr].raw_HSV_coefs[1][3], raw_hsv_coefs[lr].raw_HSV_coefs[1][4]) ;
	printf("raw_hsv_coefs V = %f %f %f %f\n", raw_hsv_coefs[lr].raw_HSV_coefs[2][0], raw_hsv_coefs[lr].raw_HSV_coefs[2][1],
		raw_hsv_coefs[lr].raw_HSV_coefs[2][2], raw_hsv_coefs[lr].raw_HSV_coefs[2][3]) ;
    }

    int save_uses_CIE_model = uses_CIE_model ;
    uses_CIE_model = 0 ;
    FILE *fp = fopen("samples.dat", "w") ;
    for(int i = 0 ; i < n_samples ; i++) {
	float hhat, shat, vhat ;
	// estimate hsv, will used blended model if requested
	predict_sky_hsv(samples[i].px, samples[i].py, &hhat, &shat, &vhat) ;
	float angle = pData->FOV_horizontal*samples[i].px ;

	fprintf(fp, "%6.4f %6.4f", samples[i].px, samples[i].py+samples[i].px/50.) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", samples[i].h, samples[i].s, samples[i].v) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", hhat, shat, vhat) ;
	fprintf(fp, " %6.2f %5.3f %5.3f", hhat-samples[i].h, shat-samples[i].s, vhat-samples[i].v) ;
	fprintf(fp, " %5.1f\n", angle) ;

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

#ifdef OLD_CODE
// weighting for all points when fitting model for single function for entire sky
static inline float px_wgt_single(float px, float abs_error)
{
    return 1./(1.+abs_error)  ;
}

// weighting for all points when fitting model for left half of sky
static inline float px_wgt_left(float px, float abs_error)
{
    px *= 1.175 ;
    px -= .25 ;
    if(px < 0.) px = 0. ;
    return (1.-3.*px*px+2.*px*px*px) / (1.+abs_error) ;
}

// weighting for all points when fitting model for right half of sky
static inline float px_wgt_right(float px, float abs_error)
{
    px = (1.-px)*1.175 -.25;
    if(px < 0.) px = 0. ;
    return (1.-3.*px*px+2.*px*px*px) / (1.+abs_error) ;
}


static inline float local_h_from_pxpy(float px, float py, struct HSV_model_coefs *p)
{
    return p->raw_HSV_coefs[0][0] + p->raw_HSV_coefs[0][1]*px + p->raw_HSV_coefs[0][2]*py + p->raw_HSV_coefs[0][3]*1./(py-horizon_py+.001) ;
}

static inline float local_ty_from_pxpy(float px, float py, struct HSV_model_coefs *p)
{
    float py_hc = py_adjusted_for_horizon_curvature(px, py) ;
    py_hc -= horizon_py ;

    return 1./(1.+exp(-py_hc*p->S_inv_coef)) ; // really a simple logistic centered on horizon_py
}

static inline float local_s_from_pxpy(float px, float py, struct HSV_model_coefs *p)
{
    float ty = local_ty_from_pxpy(px,py, p) ;
    return p->raw_HSV_coefs[1][0] + p->raw_HSV_coefs[1][1]*px + p->raw_HSV_coefs[1][2]*ty + p->raw_HSV_coefs[1][3]*px*ty ;
}

static inline float local_v_from_pxpy(float px, float py, struct HSV_model_coefs *p)
{
    return p->raw_HSV_coefs[2][0] + p->raw_HSV_coefs[2][1]*px + p->raw_HSV_coefs[2][2]*py + p->raw_HSV_coefs[2][3]*px*px ;
}

#define hsv_blended(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,&raw_hsv_coefs[0]) ;\
    float rvalue = f(px,py,&raw_hsv_coefs[1]) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }


static inline float h_from_pxpy(float px, float py)
{
    if(full_hsv_model_level==1)
	return local_h_from_pxpy(px,py,&raw_hsv_coefs[0]) ;

    float h ;
    hsv_blended(local_h_from_pxpy, px, py, h) ;
    return h ;
}

static inline float s_from_pxpy(float px, float py)
{
    if(full_hsv_model_level==1)
	return local_s_from_pxpy(px,py,&raw_hsv_coefs[0]) ;

    float s ;
    hsv_blended(local_s_from_pxpy, px, py, s) ;
    return s ;
}

static inline float v_from_pxpy(float px, float py)
{
    if(full_hsv_model_level==1)
	return local_v_from_pxpy(px,py,&raw_hsv_coefs[0]) ;

    float v ;
    hsv_blended(local_v_from_pxpy, px, py, v) ;
    return v ;
}
#endif
