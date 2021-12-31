#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>

#include "skyfill_tif.h"
#include "colorspace_conversions.h"
#include "find_sky.h"
#include "feather_factor.h"
#include "pixel_tests.h"

// return luminance for image data at pixel x,y
float get_xy_L(tdata_t *image, int x, int y)
{
    if(y < 0) y = 0 ;
    if(y > IMAGE_HEIGHT-1) y = IMAGE_HEIGHT-1 ;
    uint16_t r,g,b ;
    tif_get3c(image,x,y,r,g,b) ;

    float fr = (float)r/MAX16f ;
    float fg = (float)g/MAX16f ;
    float fb = (float)b/MAX16f ;
    return 0.21126*fr + 0.7152*fg + 0.0722*fb ;
}

void estimate_sky(int x0,int x1,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky,int16_t *final_end_of_sky,
		    float depth_of_sample, int h, int w, float feather_factor, int debug, float depth_of_fill, int extra_feather_length, SKYFILL_DATA_t *pData)
{
    int x, y ;
    fprintf(stderr, "estimate sky x:%d to %d, feather_factor=%f, nonlinear_feather=%f\n", x0,x1,feather_factor, pData->nonlinear_feather) ;
    dump_parameters_to_stdout() ;

    float sun_px = sun_x_angle2px(pData->sun_x) ;

    if(x0 < 0) x0=0 ;
    if(x1 > w-1) x1=w-1 ;

    V_sun = image_relative_pixel_to_V3(sun_x_angle2px(pData->sun_x), pData->sun_py) ;

    pData->maximum_CIE_vhat = find_maximum_CIE_vhat() ;

    float *vcorrect_at_x, *hcorrect_at_x, *scorrect_at_x ;
    float *avg_vcorrect_at_x, *avg_hcorrect_at_x, *avg_scorrect_at_x ;

    // actual h,s,v divided by predicted at column x
    hcorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;
    scorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;
    vcorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;

    // moving average actual h,s,v divided by predicted at column x
    avg_hcorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;
    avg_scorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;
    avg_vcorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;

#ifdef OLD_CODE
    // is the sky saturation prediction not right, saturation should increases from max start of sky to y=0 ;
    {
	struct mstat m_fit = init_reg(2) ; // will be a simple linear in y
	for(x = x0 ; x <= x1 ; x += (x1-x0)/20) {
	    if(column_mask[x] == 1) continue ;
	    for(y = 0 ; y < max_start_of_sky ; y++) {
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

		float hhat,shat,vhat ;
		predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;
		double xv[2] = {1.,py} ;
		sum_reg(&m_fit, xv, shat, 1.) ;
	    }
	}

	double coefs[2] ;
	estimate_reg(&m_fit, coefs) ;

	if(coefs[1] < 0.0) {
	    fprintf(stderr, "*********************************************************************************\n") ;
	    fprintf(stderr, "*********************************************************************************\n") ;
	    fprintf(stderr, "*****************  WARNING WARNING WARNING **************************************\n") ;
	    fprintf(stderr, "It appears that the saturation prediction is not correct\n") ;

	    if(sat_prediction_method == 1) {
		fprintf(stderr, " --reverting to a simpler\n") ;
		fprintf(stderr, "prediction method -- but a better fit could be had by adding \"-spm 0\" to the command line\n") ;
	    } else {
		fprintf(stderr, " results may look unnatural.\n") ;
	    }

	    fprintf(stderr, "*********************************************************************************\n") ;
	    fprintf(stderr, "*********************************************************************************\n") ;

	    S_from_V_coefs[0] = sat_sky ;
	    S_from_V_coefs[1] = 0. ;
	    S_from_V_coefs[2] = 0. ;
	} else {
	    fprintf(stderr, "*********************************************************************************\n") ;
	    fprintf(stderr, "Modeled sky saturation coefs are b0:%f b1:%f\n", coefs[0], coefs[1]) ;
	    fprintf(stderr, "*********************************************************************************\n") ;
	}

    }
#endif


    if(pData->show_raw_prediction || pData->estimate_only) {
	for(x = x0 ; x <= x1 ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    int ymax = pData->show_raw_prediction ? pData->start_of_sky[x]+1 : pData->end_of_sky[x]  ;

	    for(y = 0 ; y < ymax ; y++) {
		// completely model sky color
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

		float hhat,shat,vhat ;
		uint16_t rhat, ghat, bhat ;

		predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;

/*  		shat *= final_saturation_factor ;  */
/*  		predict_sky_huesat_from_val(vhat, &hhat, &shat, &vhat, px, py) ;  */

		hsv2rgb16(hhat,shat,vhat,&rhat,&ghat,&bhat) ;

		tif_set4c(image,x,y, (uint16_t)(rhat+0.5), (uint16_t)(ghat+0.5), (uint16_t)(bhat+0.5), MAX16);
	    }

	}

	goto estimate_sky_cleanup ;
    }

    if(pData->full_sky_replacement || pData->show_raw_error) {
	if(pData->show_raw_error) pData->full_sky_replacement=0 ; // show_raw_error overrides 
	for(x = x0 ; x <= x1 ; x++) {
/*  	    if(pData->column_mask[x] == 1) continue ;  */
	    int ymax = pData->end_of_sky[x]+10  ;

	    if(ymax > IMAGE_HEIGHT-1)
		ymax = IMAGE_HEIGHT-1 ;

	    ymax = IMAGE_HEIGHT-1 ;

	    int feather_end_y ;
	    int feather_length = compute_feather_length(x, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length, pData) ;

	    for(y = 0 ; y < ymax ; y++) {
		// completely model sky color
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

		if(py_adjusted_for_horizon_curvature(px,py) < pData->horizon_py+.005)
		    break ;

		float hhat,shat,vhat ;
		uint16_t rhat, ghat, bhat ;
		float r_err=0., g_err=0., b_err=0. ;
		float max_err=0. ;
		float p_hat = 1. ;
		float sv_err = 0. ;
		float prob_sky=1. ;

		predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;
		hsv2rgb16(hhat,shat,vhat,&rhat,&ghat,&bhat) ;

		if(y >= pData->start_of_sky[x]) {
		    uint16_t r, g, b ;
		    float feather_p_hat ;

		    if(y <= feather_end_y)
			feather_p_hat = compute_feather_at_y(x, y, feather_length, b, feather_factor, pData) ;
		    else
			feather_p_hat = 0.0 ;

		    tif_get3c(image,x,y,r,g,b) ;
		    float h,s,v ;
		    rgb2hsv16(r,g,b,&h,&s,&v) ;

		    sv_err=(fabs(h-hhat)/36.*s+fabs(s-shat)+fabs(v-vhat))*HALF16 ;

		    if(sv_err < 9.*256.) {
			p_hat = 1. ;
		    } else {
			p_hat = 1. - (sv_err-9.*256)/(30.*256.) ;
		    }

		    if(rhat > r)
			r_err = rhat-r ;
		    else
			r_err = r-rhat ;

		    if(ghat > g)
			g_err = ghat-g ;
		    else
			g_err = g-ghat ;

		    if(bhat > b)
			b_err = bhat-b ;
		    else
			b_err = b-bhat ;

		    max_err = r_err ;

		    if(max_err < g_err) max_err = g_err ;
		    if(max_err < b_err) max_err = b_err ;
/*  		    p_hat = 1.-2.*max_err/(float)(HALF16) ;  */

		    // compute hsv of r_err,g_err,b_err
		    // looks like if saturation is low, it is probably
		    // a localized error of the full sky model
		    float h_of_e,s_of_e,v_of_e ;
		    rgb2hsv16(r_err,g_err,b_err,&h_of_e,&s_of_e,&v_of_e) ;

		    prob_sky = (1.-s_of_e)*(1.-v_of_e) ;
		    if(s_of_e > .01) s_of_e = .01 ;
		    prob_sky = (1.-v_of_e)*(1.-s_of_e*s_of_e) ;

#define PS_THRESH 0.95
		    if(prob_sky > PS_THRESH) {
			p_hat = 1. ;
		    } else {
			p_hat = 1. - (PS_THRESH-prob_sky)/(PS_THRESH-.85) ;
		    }


		    if(p_hat > 1.) p_hat = 1. ;
		    if(p_hat < 0.) p_hat = 0. ;
		    if(! pData->show_raw_error && p_hat < feather_p_hat) p_hat = feather_p_hat ;

		    // cloud detector, doesn't seem to work
/*  		    if(y < pData->end_of_sky[x] && shat/.8 > s && vhat*.8 < v) {  */
/*  			p_hat=0. ;  */
/*  			max_err=0. ;  */
/*  		    }  */

		    if(1) {
			rhat = (uint16_t)(p_hat*(float)rhat + (1.-p_hat)*(float)r +.5) ;
			ghat = (uint16_t)(p_hat*(float)ghat + (1.-p_hat)*(float)g +.5) ;
			bhat = (uint16_t)(p_hat*(float)bhat + (1.-p_hat)*(float)b +.5) ;
		    } else {
			// doesn't seem to work any better, and is slower
/*  			hhat = (p_hat*(float)hhat + (1.-p_hat)*(float)h) ;  */
/*  			shat = (p_hat*(float)shat + (1.-p_hat)*(float)s) ;  */
/*  			vhat = (p_hat*(float)vhat + (1.-p_hat)*(float)v) ;  */
/*    */
/*  			hsv2rgb16(hhat,shat,vhat,&rhat,&ghat,&bhat) ;  */
		    }
		}


		if(pData->full_sky_replacement) {
		    tif_set4c(image,x,y, rhat, ghat, bhat, MAX16);
		} else {
		    uint16_t val = (uint16_t)(p_hat * (float)MAX16) ;
		    tif_set4c(image,x,y, val, val, val, MAX16);
/*  		    tif_set4c(image,x,y, r_err, g_err, b_err, MAX16);  */
/*  		    tif_set4c(image,x,y, r_err, g_err, b_err, MAX16);  */
/*  		    tif_set4c(image,x,y, sv_err, sv_err, sv_err, MAX16);  */
		}
	    }

	}

	goto estimate_sky_cleanup ;
    }

    // need to generate correction factors by column

    for(x = x0 ; x <= x1 ; x++) {
	if(pData->column_mask[x] == 1) continue ;


	// to collect correction data for this column only look 85% into the column
	float sky_height = pData->end_of_sky[x] - pData->start_of_sky[x] ;

/*  	int feather_end_y = pData->end_of_sky[x] ;  */
	int feather_end_y = pData->start_of_sky[x] + (int)(sky_height*.25+0.5) ;

	if(feather_end_y > IMAGE_HEIGHT-1) feather_end_y = IMAGE_HEIGHT-1 ;
	if(feather_end_y > end_of_sky[x]) feather_end_y = pData->end_of_sky[x] ;
	int feather_length = feather_end_y-pData->start_of_sky[x]+1 ;

	float sum_h=0., sum_s=0., sum_v=0. ;
	float sum_hhat=0., sum_shat=0., sum_vhat=0. ;

	float sum_wgts=0. ;

	// original
	int y0 = pData->start_of_sky[x] ;
	int y1 = feather_end_y ;

	for(y = y0 ; y < y1 ; y++) {
	    float hhat,shat,vhat ;
	    float fp0 ;

	    fp0 = compute_feather_at_y(x, y, feather_length, HALF16, feather_factor, pData) ;

	    if(xy_has_nonblack_pixel(image, x, y) == 0) continue ;
	    sum_wgts=0. ;

	    // completely model sky color
	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;
	    shat *= 1.0 - sqrt(fp0)*(1.-pData->final_saturation_factor) ;

/*  	    predict_sky_huesat_from_val(vhat, &hhat, &shat, &vhat, px, py) ;  */

	    uint16_t r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    sum_h += h*s ;
	    sum_s += s ;
	    sum_v += v*s ;

	    sum_hhat += hhat*s ;
	    sum_shat += shat ;
	    sum_vhat += vhat*s ;
	    sum_wgts += s ;
	}

#define CORRECTION_TYPE_IS_FACTOR
#ifdef CORRECTION_TYPE_IS_FACTOR
	if(sum_wgts == 0) {
	    // can happen at edges
	    hcorrect_at_x[x] = 1.f ;
	    scorrect_at_x[x] = 1.f ;
	    vcorrect_at_x[x] = 1.f ;
	} else {
	    hcorrect_at_x[x] = sum_h/sum_hhat ;
	    scorrect_at_x[x] = sum_s/sum_shat ;
	    vcorrect_at_x[x] = sum_v/sum_vhat ;
	}
#else
	if(sum_wgts == 0) {
	    // can happen at edges
	    hcorrect_at_x[x] = 0.f ;
	    scorrect_at_x[x] = 0.f ;
	    vcorrect_at_x[x] = 0.f ;
	} else {
	    hcorrect_at_x[x] = (sum_h-sum_hhat)/sum_wgts ;
	    scorrect_at_x[x] = (sum_s-sum_shat)/sum_wgts ;
	    vcorrect_at_x[x] = (sum_v-sum_vhat)/sum_wgts ;
	}
#endif
    }

    if(0 && pData->allowed_sky_type == UNIFORM_CIE_SKIES) {
	float sum_h=0., sum_s=0., sum_v=0. ;
	int n_correct = 0 ;

	for(x = x0 ; x <= x1 ; x++) {
	    if(pData->column_mask[x0] == 1) continue ;

	    // vcorrect is 0. if there were no samples to calculate it
	    if(vcorrect_at_x[x0] > 1.e-10) {
		sum_h += hcorrect_at_x[x] ;
		sum_s += scorrect_at_x[x] ;
		sum_v += vcorrect_at_x[x] ;
		n_correct++ ;
	    }
	}

	if(n_correct > 0) {
	    for(x = x0 ; x <= x1 ; x++) {
		hcorrect_at_x[x] = (hcorrect_at_x[x] + sum_h / (float)n_correct)/2. ;
		scorrect_at_x[x] = (scorrect_at_x[x] + sum_s / (float)n_correct)/2. ;
		vcorrect_at_x[x] = (vcorrect_at_x[x] + sum_v / (float)n_correct)/2. ;
	    }
	}
    }

    for(x = x0 ; x <= x1 ; x++) {
	if(pData->column_mask[x] == 1) continue ;
	if(pData->show_raw_prediction) continue ;

	float sum_h=0., sum_s=0., sum_v=0. ;
	int n_correct = 0 ;

#ifdef CORRECTION_TYPE_IS_FACTOR
	avg_hcorrect_at_x[x] = 1. ;
	avg_scorrect_at_x[x] = 1. ;
	avg_vcorrect_at_x[x] = 1. ;
#else
	avg_hcorrect_at_x[x] = 0. ;
	avg_scorrect_at_x[x] = 0. ;
	avg_vcorrect_at_x[x] = 0. ;
#endif

	for(int dx=-126 ; dx < 127 ; dx++) {
	    if(x+dx < 0 || x+dx > IMAGE_WIDTH-1) continue ;
	    if(pData->column_mask[x+dx] == 1) continue ;

	    // vcorrect is 0. if there were no samples to calculate it
	    if(vcorrect_at_x[x+dx] > 1.e-10) {
		sum_h += hcorrect_at_x[x+dx] ;
		sum_s += scorrect_at_x[x+dx] ;
		sum_v += vcorrect_at_x[x+dx] ;
		n_correct++ ;
	    }
	}

	if(feather_factor > .01 && n_correct > 0) {
	    avg_hcorrect_at_x[x] = sum_h / (float)n_correct ;
	    avg_scorrect_at_x[x] = sum_s / (float)n_correct ;
	    avg_vcorrect_at_x[x] = sum_v / (float)n_correct ;
	}
    }



    for(x = x0 ; x <= x1 ; x++) {
	if(pData->column_mask[x] == 1) continue ;

	int feather_end_y ;
	int feather_length = compute_feather_length(x, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length, pData) ;


	for(y = 0 ; y < feather_end_y ; y++) {
	    uint16_t a = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+3] ;

	    if(a < HALF16 && y > pData->raw_start_of_sky[x]) continue ; // do not adjust transparent pixels below original start of sky

	    uint16_t rhat, ghat, bhat ;

	    // current pixel value
	    uint16_t r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    // completely model sky color
	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    float fp0 = compute_feather_at_y(x, y, feather_length, b, feather_factor, pData) ;
	    float fp1 = 1.-fp0 ;

	    float hhat,shat,vhat ;
	    predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;

/*  	    shat *= final_saturation_factor ;  */
/*  	    predict_sky_huesat_from_val(vhat, &hhat, &shat, &vhat, px, py) ;  */

	    float fpsat = 1.0 - (float)(y)/(float)pData->start_of_sky[x] ;
	    shat *= 1.0 - fpsat*(1.-pData->final_saturation_factor) ;


	    // are we at or above the sun?
/*  	    if(sun_py < 1. && sun_px > -0.5 && sun_px < 0.5 && py >= sun_py && fabs(px-sun_px) < py-sun_py) {  */
	    if(0) {
		float px_dist = (py-pData->sun_py) ;
/*  		float fdist_from_center = fabs(px-sun_px)/(py-sun_py) ;  */

		// x0, x1 are left and right extents where we want to
		// modify values more directly 
		int x0 = IMAGE_RELATIVE_TO_PIXEL_X(sun_px-px_dist) ;
		int x1 = IMAGE_RELATIVE_TO_PIXEL_X(sun_px+px_dist) ;
		float px1=(float)(x-x0)/(float)(x1-x0) ;
		float px0=1.-px1 ;

/*  		if(fdist_from_center > .99) {  */
/*  		    fprintf(stderr, "sun est x,y:%d,%d, x0,x1:%d,%d px0, px1 %f,%f\n", x,y, x0,x1, px0, px1) ;  */
/*  		}  */

#ifdef CORRECTION_TYPE_IS_FACTOR
		hhat *= px0*avg_hcorrect_at_x[x0] + px1*avg_hcorrect_at_x[x1] ;
		shat *= px0*avg_scorrect_at_x[x0] + px1*avg_scorrect_at_x[x1] ;
		vhat *= px0*avg_vcorrect_at_x[x0] + px1*avg_vcorrect_at_x[x1] ;
#else
		hhat += px0*avg_hcorrect_at_x[x0] + px1*avg_hcorrect_at_x[x1] ;
		shat += px0*avg_scorrect_at_x[x0] + px1*avg_scorrect_at_x[x1] ;
		vhat += px0*avg_vcorrect_at_x[x0] + px1*avg_vcorrect_at_x[x1] ;
#endif

		fp0 = 1. ;
		fp1 = 1.-fp0 ;
	    } else {
#ifdef CORRECTION_TYPE_IS_FACTOR
/*  		if(full_sky_replacement) {  */
/*  		    hhat *= 1 + fp1*(1.-avg_hcorrect_at_x[x]) ;  */
/*  		    shat *= 1 + fp1*(1.-avg_scorrect_at_x[x]) ;  */
/*  		    vhat *= 1 + fp1*(1.-avg_vcorrect_at_x[x]) ;  */
/*  		} else {  */
		    hhat *= avg_hcorrect_at_x[x] ;
		    shat *= avg_scorrect_at_x[x] ;
		    vhat *= avg_vcorrect_at_x[x] ;
/*  		}  */
#else
		hhat += avg_hcorrect_at_x[x] ;
		shat += avg_scorrect_at_x[x] ;
		vhat += avg_vcorrect_at_x[x] ;
#endif
	    }

	    if(vhat < 0.f) vhat=0.f ;
	    if(vhat > 1.f) vhat=1.f ;
	    if(shat < 0.f) shat=0.f ;
	    if(shat > 1.f) shat=1.f ;
	    hsv2rgb16(hhat,shat,vhat,&rhat,&ghat,&bhat) ;



	    r = (uint16_t)(fp0*(float)rhat + fp1*(float)r+0.5) ;
	    g = (uint16_t)(fp0*(float)ghat + fp1*(float)g+0.5) ;
	    b = (uint16_t)(fp0*(float)bhat + fp1*(float)b+0.5) ;
	    tif_set4c(image,x,y, r,g,b,MAX16) ;
	}

	// did end of sky not go far enough?   See if there is an increase in L
	// at end of sky +1 followed by a decrease in L at end of sky +2

	int yeos ;
	pData->final_end_of_sky[x] = pData->end_of_sky[x] ;

	for(yeos = pData->end_of_sky[x]-2 ; yeos < pData->end_of_sky[x]+2 ; yeos++) {

	    float L_eos0 = get_xy_L(image, x, yeos) ;
	    float L_eos1 = get_xy_L(image, x, yeos+1) ;
	    float L_eos2 = get_xy_L(image, x, yeos+2) ;

	    if(L_eos1 > L_eos0 && L_eos1 > L_eos2) {
		float fac = L_eos0/L_eos1 ; // scale rgb values back to L found at end of sky
		y = yeos+1 ;
		pData->final_end_of_sky[x] = y ;

		uint16_t r,g,b,a ;
		tif_get4c(image,x,y,r,g,b,a) ;

		if(a < HALF16 && y > pData->raw_start_of_sky[x]) continue ; // do not adjust transparent pixels below original start of sky

		tif_set3c(image,x,y, (uint16_t)((float)r*fac+0.5), (uint16_t)((float)g*fac+0.5), (uint16_t)((float)b*fac+0.5));
	    }
	}
    }

    //for(x = x0 ; x <= x1 ; x++) {
//	float sky_height = end_of_sky[x] - start_of_sky[x] ;
//	int feather_end_y = start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) ;
 //   }

estimate_sky_cleanup:

    free(hcorrect_at_x) ;
    free(scorrect_at_x) ;
    free(vcorrect_at_x) ;

    free(avg_hcorrect_at_x) ;
    free(avg_scorrect_at_x) ;
    free(avg_vcorrect_at_x) ;

    fprintf(stderr, "done with estimate sky, column=%d\n", x) ;

}
