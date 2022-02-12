#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>

#include "skyfill_tif.h"
#include "colorspace_conversions.h"
#include "find_sky.h"
#include "feather_factor.h"
#include "pixel_tests.h"
#include "estimate_sky.h"
#include "mstat.h"

// return luminance for image data at pixel x,y
float get_xy_L(tdata_t *image, int x, int y)
{
    if(y < 0) y = 0 ;
    if(y > IMAGE_HEIGHT-1) y = IMAGE_HEIGHT-1 ;
    uint16_t rgb[3] ;
    tif_get3cv(image,x,y,rgb) ;

    float fr = (float)rgb[0]/MAX16f ;
    float fg = (float)rgb[1]/MAX16f ;
    float fb = (float)rgb[2]/MAX16f ;
    return 0.21126*fr + 0.7152*fg + 0.0722*fb ;
}

int get_end_of_sky(SKYFILL_DATA_t *pData, int x)
{
    if(pData->column_mask[x] == 0) return pData->end_of_sky[x] ;

    // column is masked, search for nearest left, nearest right
    int left_eos=-1 ;
    int right_eos=-1 ;
    int left_x, right_x ;

    for(left_x = x-1 ; left_x > -1 ; left_x--) {
	if(pData->column_mask[left_x] == 0) {
	    left_eos = pData->end_of_sky[left_x] ;
	    break ;
	}
    }

    for(right_x = x+1 ; right_x < IMAGE_WIDTH ; right_x++) {
	if(pData->column_mask[right_x] == 0) {
	    right_eos = pData->end_of_sky[right_x] ;
	    break ;
	}
    }

    if(left_x == -1) return right_eos ;
    if(right_x == -1) return left_eos ;

    float p_left = 1.-(float)(x-left_x)/(float)(right_x-left_x) ;
    return (int)(p_left*(float)left_eos + (1.-p_left)*(float)right_eos + 0.5) ;
}

float get_prob_sky(tdata_t *image, int x,int y, SKYFILL_DATA_t *pData, struct sky_pixel_data *p)
{

    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

    p->rgb_err[0]=0., p->rgb_err[1]=0., p->rgb_err[2]=0. ;
    p->hsv_of_e[0]=0.,p->hsv_of_e[1]=0.,p->hsv_of_e[2]=0. ;
    p->hsv[0]=0.,p->hsv[1]=0.,p->hsv[2]=0. ;
    p->max_err=0. ;
    p->p_hat = 1. ;
    p->hsv_err = 0. ;
    p->prob_sky=1. ;

    predict_sky_hsv(px, py, p->hsv_hat) ;
    hsv2rgb16(p->hsv_hat,p->rgb_hat) ;

    if(y >= pData->start_of_sky[x]) {
	tif_get4c(image,x,y,p->rgb[0],p->rgb[1],p->rgb[2], p->a) ;
	rgb2hsv16(p->rgb,p->hsv) ;

	// intensity
	float I_hat = (float)(p->rgb_hat[0]+p->rgb_hat[1]+p->rgb_hat[2])/3. ;
	float I = (float)(p->rgb[0]+p->rgb[1]+p->rgb[2])/3. ;

	p->rgb_err[0] = fabsf((float)p->rgb_hat[0] - (float)p->rgb[0]) ;
	p->rgb_err[1] = fabsf((float)p->rgb_hat[1] - (float)p->rgb[1]) ;
	p->rgb_err[2] = fabsf((float)p->rgb_hat[2] - (float)p->rgb[2]) ;

/*  	if(0) {  */
/*    */
/*  	if(p->r_hat > p->r)  */
/*  	    p->r_err = p->r_hat-p->r ;  */
/*  	else  */
/*  	    p->r_err = p->r-p->r_hat ;  */
/*    */
/*  	if(p->g_hat > p->g)  */
/*  	    p->g_err = p->g_hat-p->g ;  */
/*  	else  */
/*  	    p->g_err = p->g-p->g_hat ;  */
/*    */
/*  	if(p->b_hat > p->b)  */
/*  	    p->b_err = p->b_hat-p->b ;  */
/*  	else  */
/*  	    p->b_err = p->b-p->b_hat ;  */
/*  	}  */

	// a few attempts at computing probablity this  pixel is a clear sky pixel
	if(0) {
/*  	    float I_err = fabsf(1.-I/I_hat) ;  */
/*  	    float s_err = fabsf(1.-p->hsv[1]/p->hsv_hat[1]) ;  */
/*  	    float h_err = fabsf(1.-p->hsv[0]/p->hsv_hat[0]) ;  */
/*  	    float wgt_I = 5. ;  */
/*  	    float wgt_s = 3. ;  */
/*  	    float wgt_h = 10.*p->hsv_hat[1] ; // low saturation means low weight for hue  */
/*  	    float hsI_err = (wgt_I*I_err + wgt_s*s_err + wgt_h*h_err)/(wgt_I+wgt_s+wgt_h) ;  */
/*    */
/*  	    if(hsI_err > 1.)  */
/*  		hsI_err = 1. ;  */
/*    */
/*  	    p->prob_sky = 1.-hsI_err ;  */
	} else if(0) {

	    // compute hsv of r_err,g_err,b_err
	    // looks like if saturation is low, it is probably
	    // a localized error of the full sky model

/*  	    rgb2hsv16(p->rgb_err,p->hsv_of_e) ;  */

	    // v_of_e will be the max of (r_err,g_err,b_err) now

/*  	    if(p->hsv_of_e[1] > .01) p->hsv_of_e[1] = .01 ;  */
/*    */
/*  	    p->prob_sky = (1.-p->hsv_of_e[2])*(1.-p->hsv_of_e[1]*p->hsv_of_e[1]) ;  */
	} else if(1) {
	    float wgt_r = p->rgb_hat[0] ;
	    float wgt_g = p->rgb_hat[1] ;
	    float wgt_b = p->rgb_hat[2] ;

	    float wgted_err = (
		wgt_r*(float)p->rgb_err[0]/(float)MAX16 +
		wgt_g*(float)p->rgb_err[1]/(float)MAX16 +
		wgt_b*(float)p->rgb_err[2]/(float)MAX16) / (wgt_r + wgt_g + wgt_b) ;

	    p->prob_sky= 1. - wgted_err ;
	} else {

	    p->prob_sky= 1. - (float)p->rgb_err[2]/(float)MAX16 ;
	}

	p->p_hat = 1./(1+exp(-(p->prob_sky-pData->full_sky_replacement_thresh)*pData->full_sky_replacement_ramp)) ;

	if(p->hsv[1] < p->hsv_hat[1]) {
	    if(I > I_hat && p->hsv[1] < p->hsv_hat[1]) {
		// if the actual is less saturated, with higher intensity, it could be a cloud
		p->p_hat *= I_hat/I * p->hsv[1]/p->hsv_hat[1] ;
	    } else {
		// if the actual is less saturated, with lower intensity, it could also be a cloud
		p->p_hat *= I/I_hat * p->hsv[1]/p->hsv_hat[1] ;
	    }
	}




	if(p->p_hat > 1.) p->p_hat = 1. ;
	if(p->p_hat < 0.) p->p_hat = 0. ;
    }

    return p->prob_sky ;
}

void estimate_sky(int x0,int x1,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky,int16_t *final_end_of_sky,
		    float depth_of_sample, int h, int w, float feather_factor, int debug, float depth_of_fill, int extra_feather_length, SKYFILL_DATA_t *pData)
{
/*      int x, y ;  */
    fprintf(stderr, "estimate sky x:%d to %d, depth_of_fill:%f feather_factor=%f, nonlinear_feather=%f\n", x0,x1,depth_of_fill,feather_factor, pData->nonlinear_feather) ;
    dump_parameters_to_stdout() ;

    pData->n_estimate_sky_clipped = 0 ; // number of pixels clipped at MAX16 output by estimate_sky

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


#define set4cv_clip_check(image,x,y,rgb,a,pData) {\
    if( ((rgb[0]) == MAX16) || ((rgb[1]) == MAX16) || ((rgb[2]) == MAX16) ) pData->n_estimate_sky_clipped++;\
    tif_set4c(image,x,y, rgb[0], rgb[1], rgb[2], a);\
    }

    if(pData->show_raw_prediction || pData->estimate_only) {
	for(int x = x0 ; x <= x1 ; x++) {
	    if(!pData->estimate_only && pData->column_mask[x] == 1) continue ;
	    int ymax = pData->show_raw_prediction ? pData->start_of_sky[x]+1 : pData->end_of_sky[x]  ;
/*  	    int ymax = pData->show_raw_prediction ? pData->start_of_sky[x]+1 : IMAGE_HEIGHT  ;  */

	    for(int y = 0 ; y < ymax ; y++) {
		// completely model sky color
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

		float hsv_hat[3] ;
		uint16_t rgb_hat[3] ;

		predict_sky_hsv(px, py, hsv_hat) ;

/*  		s_hat *= final_saturation_factor ;  */

		hsv2rgb16(hsv_hat,rgb_hat) ;

		set4cv_clip_check(image,x,y, rgb_hat, MAX16, pData);
	    }

	}

	goto estimate_sky_cleanup ;
    }

    if(uses_CIE_model)
	fprintf(stderr, "****************** Blending in sun model *****************\n") ;

    if(pData->full_sky_replacement || pData->show_raw_error || pData->show_sky_prediction) {
	if(pData->show_raw_error) pData->full_sky_replacement=0 ; // show_raw_error overrides 
	if(pData->show_raw_error) pData->show_sky_prediction=0 ; // show_raw_error overrides 
	if(pData->show_sky_prediction) pData->full_sky_replacement=0 ; // show_sky_prediction overrides 



	for(int x = x0 ; x <= x1 ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    int ymax = get_end_of_sky(pData,x)  ;
	    if(ymax < pData->manual_end_of_sky[x]) ymax = pData->manual_end_of_sky[x] ;
	    if(ymax > pData->manual_end_of_sky[x] && pData->manual_end_of_sky[x] != -1) ymax = pData->manual_end_of_sky[x] ;

	    int feather_end_y ;
	    int feather_length = compute_feather_length_with_eos(x, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length, pData, ymax) ;

/*  	    ymax = IMAGE_HEIGHT-1 ;  */

	    int one_more_applied=0 ;

	    for(int y = 0 ; y < ymax ; y++) {
		// completely model sky color
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

/*  		if(py_adjusted_for_horizon_curvature(px,py) < pData->horizon_py+.005)  */
/*  		    break ;  */

		// if we are near the sun, but below the start of sky, incorporation
		// more of the base image so things like sun flare can come through
		float sun_prob = 0. ;
		float py_dist = -1. ;
		float px_dist = -1. ;

		if(uses_CIE_model != CIE_NO_MODEL) {
		    py_dist = fabs(py-pData->sun_py)*(float)IMAGE_HEIGHT/(float)IMAGE_WIDTH ;
		    px_dist = fabs(px-sun_px) ;

		    // weight of perez_F at sun center, falls off with cube of distance
/*  		    sun_prob = pData->perez_F/(1.+pData->perez_G*(px_dist*px_dist*px_dist+py_dist*py_dist*py_dist)) ;  */

		    float pdist = sqrtf(px_dist*px_dist+py_dist*py_dist) ;
		    float r = pdist/pData->perez_G ;
		    float damp ;

		    if(r < 1.) {
			damp=1. ;
		    } else {
			r -= 1. ;
			// over the next 50% of perez_G dist, take damp from 1.0 to 0.0
/*  			damp = 1. - r*2. ;  */
/*  			if(damp < 0.) damp=0. ;  */

			damp = 1.-1./(1.+exp(5.-10*r)) ;
		    }
		    sun_prob = pData->perez_F*sqrtf(damp) ;
/*  		    sun_prob=0. ;  */
/*  		    if(damp > 0.) sun_prob=1. ;  */

		}

		struct sky_pixel_data sp ;
		float prob_sky = get_prob_sky(image, x, y, pData, &sp) ;
		float feather_prob = 1. ;

		if(y >= pData->start_of_sky[x]) {

		    if(y <= feather_end_y)
			feather_prob = compute_feather_at_y(x, y, feather_length, MAX16, feather_factor, pData) ;
		    else
			feather_prob = 0.0 ;

		    if(pData->column_mask[x] == 1) {
			feather_prob = 0.0 ;
		    } else {
			if(y < start_of_sky[x]) {
			    feather_prob = 1. ;
			} else {
/*  			    feather_prob = 1. - (float)(y-start_of_sky[x])/((float)(ymax-(float)start_of_sky[x]) ;  */
/*  			    feather_prob = powf(feather_prob,pData->nonlinear_feather) ;  */
			}
		    }

		    // convert linear blend to logistic blend
#define L2L_B0 .5
#define L2L_B1 5
		    float l2l_offset = 1./(1.+exp(-(0.-L2L_B0)*L2L_B1)) ;
		    float l2l_scale = 1./(1.-l2l_offset*2.) ;
		    feather_prob = l2l_scale*(1./(1.+exp(-(feather_prob-L2L_B0)*L2L_B1))-l2l_offset) ;

		    // Now we have 3 factors to consider of how much of the predicted pixel color to use
		    // 1. p_hat -- the amount of predicted value to use only with consideration if the pixel looks like a sky pixel
		    // 2. feather_prob -- the amount of predicted value to use 100% at the top of the original sky, to 0% at the depth of fill level
		    // 3. image_wgt_sun -- override for how much predicted value to use given the nearness to the sun

		    if(pData->column_nonsky_repair_mask[x] == 0) {
			if(sp.p_hat < feather_prob) sp.p_hat = feather_prob ;

			// beyond depth of fill, make sure p_hat goes to 0 at ymax (this happens when manual end of sky > detected end of sky
			if(y > feather_end_y && pData->noclip_end_of_sky) {
			    sp.p_hat *= 1. - (float)(y-feather_end_y)/(ymax-feather_end_y) ;
			}

			// increase p_hat as x,y gets more near to the sun
			if(sp.p_hat < sun_prob) sp.p_hat = sun_prob ;

		    } else {
			// user has requested that pixels in the sky that are not close to sky color be unchanged
			// if p_hat is low, do not change it
			float raw_p_hat = sp.p_hat ;
			if(sp.p_hat < feather_prob) sp.p_hat = feather_prob ;
			if(sp.p_hat < sun_prob) sp.p_hat = sun_prob ;

			// now correct based on raw_p_hat
			// if raw_p_hat is 1., use all of the new value
			// if raw_p_hat is 0., use all of the initial value
			sp.p_hat = raw_p_hat * (sp.p_hat) + (1.-raw_p_hat)*raw_p_hat ;
		    }

		    if(sp.a < MAX16) {
			// uh-oh, somehow a transparent pixel slipped through, make it all sky
			sp.p_hat=1.0 ;
		    }

		    sp.rgb_hat[0] = (uint16_t)(sp.p_hat*(float)sp.rgb_hat[0] + (1.-sp.p_hat)*(float)sp.rgb[0] +.5) ;
		    sp.rgb_hat[1] = (uint16_t)(sp.p_hat*(float)sp.rgb_hat[1] + (1.-sp.p_hat)*(float)sp.rgb[1] +.5) ;
		    sp.rgb_hat[2] = (uint16_t)(sp.p_hat*(float)sp.rgb_hat[2] + (1.-sp.p_hat)*(float)sp.rgb[2] +.5) ;
		} // end of if y > start of sky block

/*  #define DEBUG_COL 1515  */
#ifdef DEBUG_COL
		if(x == DEBUG_COL && !(y%10)) {
		    printf("SUN: x:%d y:%d px_dist:%4.2f py_dist:%4.2f IWS:%4.2f p_hat:%4.2f fp_hat:%4.2f prob_sky:%4.2f\n",
		    x, y, px_dist, py_dist, image_wgt_sun, sp.p_hat, feather_prob, prob_sky) ;
		    printf("\t hsv hat %3.0f %5.2f %5.2f\n", sp.h_hat, sp.s_hat, sp.v_hat) ;
		    printf("\t b,b_hat %5d %5d\n", sp.b, sp.b_hat) ;
		}
#endif


		if(pData->full_sky_replacement) {
/*  		    set4cv_clip_check(image,x,y, sp.rgb_hat, MAX16, pData);  */
		    tif_set4c(image,x,y, sp.rgb_hat[0], sp.rgb_hat[1], sp.rgb_hat[2], MAX16) ;
		} else if(pData->show_sky_prediction) {
		    uint16_t val = (uint16_t)(sp.p_hat * (float)MAX16) ;
		    uint16_t vec[3] = {val,val,val} ;
/*  		    uint16_t val = (uint16_t)(feather_prob * (float)MAX16) ;  */
		    set4cv_clip_check(image,x,y, vec, MAX16, pData);
/*  		    uint16_t val = (uint16_t)(hsv_err*MAX16) ;  */
/*  		    set4c_clip_check(image,x,y, val, val, val, MAX16, pData);  */
		} else {
/*  		    r_err = (uint16_t)(sqrt((float)r_err/(float)MAX16)*(float)MAX16+0.5) ;  */
/*  		    g_err = (uint16_t)(sqrt((float)g_err/(float)MAX16)*(float)MAX16+0.5) ;  */
/*  		    b_err = (uint16_t)(sqrt((float)b_err/(float)MAX16)*(float)MAX16+0.5) ;  */
		    set4cv_clip_check(image,x,y, sp.rgb_err, MAX16, pData);
		}

		if(! one_more_applied) {
		    // remove possible halo like pixel at next y, if ymax wasn't large enough
		    if(y == ymax-1 && ymax-1 < IMAGE_HEIGHT-3) {
			float c_max0=sp.rgb_hat[0] ;
			if(c_max0 < sp.rgb_hat[1]) c_max0 = sp.rgb_hat[1] ;
			if(c_max0 < sp.rgb_hat[2]) c_max0 = sp.rgb_hat[2] ;

			uint16_t r,g,b ;
			tif_get3c(image,x,y+1,r,g,b) ;

			float c_max1=r ;
			if(c_max1 < g) c_max1 = g ;
			if(c_max1 < b) c_max1 = b ;

			tif_get3c(image,x,y+2,r,g,b) ;

			float c_max2=r ;
			if(c_max2 < g) c_max2 = g ;
			if(c_max2 < b) c_max2 = b ;

			if(c_max1 > c_max0 && c_max1 > c_max2) {
			    ymax++ ;
			    one_more_applied=1 ;
			}
		    }
		}
	    } // end of outer y loop

	} // end of outer x loop

	goto estimate_sky_cleanup ;
    }

    // need to generate correction factors by column

    for(int x = x0 ; x <= x1 ; x++) {
	if(pData->column_mask[x] == 1) continue ;


	// to collect correction data for this column only look 25% into the column
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

	for(int y = y0 ; y < y1 ; y++) {
	    float hsv_hat[3] ;
/*  	    float fp0 ;  */
/*    */
/*  	    fp0 = compute_feather_at_y(x, y, feather_length, HALF16, feather_factor, pData) ;  */

	    if(xy_has_nonblack_pixel(image, x, y) == 0) continue ;
	    sum_wgts=0. ;

	    // completely model sky color
	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    predict_sky_hsv(px, py, hsv_hat) ;
/*  	    shat *= 1.0 - sqrt(fp0)*(1.-pData->final_saturation_factor) ;  */

/*  	    predict_sky_huesat_from_val(vhat, &hhat, &shat, &vhat, px, py) ;  */

	    uint16_t rgb[3] ;
	    tif_get3cv(image,x,y,rgb) ;

	    float hsv[3] ;
	    rgb2hsv16(rgb,hsv) ;

	    sum_h += hsv[0]*hsv[1] ;
	    sum_s += hsv[1] ;
	    sum_v += hsv[2]*hsv[1] ;

	    sum_hhat += hsv_hat[0]*hsv[1] ;
	    sum_shat += hsv_hat[1] ;
	    sum_vhat += hsv_hat[2]*hsv[1] ;
	    sum_wgts += hsv[1] ;
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

	for(int x = x0 ; x <= x1 ; x++) {
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
	    for(int x = x0 ; x <= x1 ; x++) {
		hcorrect_at_x[x] = (hcorrect_at_x[x] + sum_h / (float)n_correct)/2. ;
		scorrect_at_x[x] = (scorrect_at_x[x] + sum_s / (float)n_correct)/2. ;
		vcorrect_at_x[x] = (vcorrect_at_x[x] + sum_v / (float)n_correct)/2. ;
	    }
	}
    }

    for(int x = x0 ; x <= x1 ; x++) {
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



    for(int x = x0 ; x <= x1 ; x++) {
	if(pData->column_mask[x] == 1) continue ;

	int feather_end_y ;
	int feather_length = compute_feather_length(x, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length, pData) ;


	for(int y = 0 ; y < feather_end_y ; y++) {
	    uint16_t a = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+3] ;

	    if(a < HALF16 && y > pData->raw_start_of_sky[x]) continue ; // do not adjust transparent pixels below original start of sky

	    uint16_t rgb_hat[3] ;

	    // current pixel value
	    uint16_t rgb[3] ;
	    tif_get3cv(image,x,y,rgb) ;


	    // completely model sky color
	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    float fp0 = compute_feather_at_y(x, y, feather_length, rgb[2], feather_factor, pData) ;
	    float fp1 = 1.-fp0 ;

	    float hsv_hat[3] ;
	    predict_sky_hsv(px, py, hsv_hat) ;

	    if(y < pData->start_of_sky[x]) {
		float fpsat = 1.0 - (float)(y)/(float)(pData->start_of_sky[x]) ;
		hsv_hat[1] *= 1.0 - fpsat*(1.-pData->final_saturation_factor) ;
	    }
#define noDEBUG_COL 997
#ifdef DEBUG_COL
	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;
	    if(x == DEBUG_COL) {
		fprintf(stderr, "ES: x,y %5d,%5d SOS %5d EOS %5d fp0:%5.3f fpsat:%5.3f hsv %5.3f %5.3f %5.3f hsv_hat %5.3f %5.3f %5.3f hsv_factors %5.3f %5.3f %5.3f\n",
			x,y, pData->start_of_sky[x], pData->end_of_sky[x], fp0, fpsat,
			h,s,v, hhat,shat,vhat,
			avg_hcorrect_at_x[x], avg_scorrect_at_x[x], avg_vcorrect_at_x[x]) ;
	    }
#endif


#ifdef CORRECTION_TYPE_IS_FACTOR
	    hsv_hat[0] *= avg_hcorrect_at_x[x] ;
	    hsv_hat[1] *= avg_scorrect_at_x[x] ;
	    hsv_hat[2] *= avg_vcorrect_at_x[x] ;
#else
	    hsv_hat[0] += avg_hcorrect_at_x[x] ;
	    hsv_hat[1] += avg_scorrect_at_x[x] ;
	    hsv_hat[2] += avg_vcorrect_at_x[x] ;
#endif

	    // clamp value between 0 and 1
	    if(hsv_hat[2] < 0.f) hsv_hat[2]=0.f ;
	    if(hsv_hat[2] > 1.f) hsv_hat[2]=1.f ;

	    // clamp sat between 0 and 1
	    if(hsv_hat[1] < 0.f) hsv_hat[1]=0.f ;
	    if(hsv_hat[1] > 1.f) hsv_hat[1]=1.f ;

	    hsv2rgb16(hsv_hat,rgb_hat) ;

/*  	    if(x == 100) {  */
/*  		fprintf(stderr, "x:%d y:%d fp_hat:%f fp_act:%f\n", x, y, fp0, fp1) ;  */
/*  	    }  */

	    rgb[0] = (uint16_t)(fp0*(float)rgb_hat[0] + fp1*(float)rgb[0]+0.5) ;
	    rgb[1] = (uint16_t)(fp0*(float)rgb_hat[1] + fp1*(float)rgb[1]+0.5) ;
	    rgb[2] = (uint16_t)(fp0*(float)rgb_hat[2] + fp1*(float)rgb[2]+0.5) ;
	    set4cv_clip_check(image,x,y, rgb,MAX16, pData) ;
	}

	// did end of sky not go far enough?   See if there is an increase in L
	// at end of sky +1 followed by a decrease in L at end of sky +2

	pData->final_end_of_sky[x] = pData->end_of_sky[x] ;
	if(feather_end_y >= pData->end_of_sky[x]) {
	    int yeos ;

	    for(yeos = pData->end_of_sky[x]-2 ; yeos < pData->end_of_sky[x]+2 ; yeos++) {

		float L_eos0 = get_xy_L(image, x, yeos) ;
		float L_eos1 = get_xy_L(image, x, yeos+1) ;
		float L_eos2 = get_xy_L(image, x, yeos+2) ;

		if(L_eos1 > L_eos0 && L_eos1 > L_eos2) {
		    float fac = L_eos0/L_eos1 ; // scale rgb values back to L found at end of sky
		    int y = yeos+1 ;
		    pData->final_end_of_sky[x] = y ;

		    uint16_t r,g,b,a ;
		    tif_get4c(image,x,y,r,g,b,a) ;

		    if(a < HALF16 && y > pData->raw_start_of_sky[x]) continue ; // do not adjust transparent pixels below original start of sky

		    tif_set3c(image,x,y, (uint16_t)((float)r*fac+0.5), (uint16_t)((float)g*fac+0.5), (uint16_t)((float)b*fac+0.5));
		}
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

    fprintf(stderr, "done with estimate sky, column=%d\n", x1) ;

}
