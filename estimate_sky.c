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

struct sky_pixel_data {
    float h_hat,s_hat,v_hat ;
    uint16_t r_hat, g_hat, b_hat ;
    float r_err, g_err, b_err ;
    float h_of_e,s_of_e,v_of_e ;
    float h,s,v ;
    float max_err ;
    float p_hat ;
    float hsv_err ;
    float prob_sky ;
    uint16_t r, g, b ;
} ;

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

    p->r_err=0., p->g_err=0., p->b_err=0. ;
    p->h_of_e=0.,p->s_of_e=0.,p->v_of_e=0. ;
    p->h=0.,p->s=0.,p->v=0. ;
    p->max_err=0. ;
    p->p_hat = 1. ;
    p->hsv_err = 0. ;
    p->prob_sky=1. ;

    predict_sky_hsv(px, py, &p->h_hat, &p->s_hat, &p->v_hat) ;
    hsv2rgb16(p->h_hat,p->s_hat,p->v_hat,&p->r_hat,&p->g_hat,&p->b_hat) ;

    if(y >= pData->start_of_sky[x]) {
	tif_get3c(image,x,y,p->r,p->g,p->b) ;
	rgb2hsv16(p->r,p->g,p->b,&p->h,&p->s,&p->v) ;

	// intensity
	float I_hat = (float)(p->r_hat+p->g_hat+p->b_hat)/3. ;
	float I = (float)(p->r+p->g+p->b)/3. ;

	p->r_err = fabs((float)p->r_hat - (float)p->r) ;
	p->g_err = fabs((float)p->g_hat - (float)p->g) ;
	p->b_err = fabs((float)p->b_hat - (float)p->b) ;

	if(0) {

	if(p->r_hat > p->r)
	    p->r_err = p->r_hat-p->r ;
	else
	    p->r_err = p->r-p->r_hat ;

	if(p->g_hat > p->g)
	    p->g_err = p->g_hat-p->g ;
	else
	    p->g_err = p->g-p->g_hat ;

	if(p->b_hat > p->b)
	    p->b_err = p->b_hat-p->b ;
	else
	    p->b_err = p->b-p->b_hat ;
	}

	// a few attempts at computing probablity this  pixel is a clear sky pixel
	if(0) {
	    float I_err = fabs(1.-I/I_hat) ;
	    float s_err = fabs(1.-p->s/p->s_hat) ;
	    float h_err = fabs(1.-p->h/p->h_hat) ;
	    float wgt_I = 5. ;
	    float wgt_s = 3. ;
	    float wgt_h = 10.*p->s_hat ; // low saturation means low weight for hue
	    float hsI_err = (wgt_I*I_err + wgt_s*s_err + wgt_h*h_err)/(wgt_I+wgt_s+wgt_h) ;

	    if(hsI_err > 1.)
		hsI_err = 1. ;

	    p->prob_sky = 1.-hsI_err ;
	} else if(0) {

	    // compute hsv of r_err,g_err,b_err
	    // looks like if saturation is low, it is probably
	    // a localized error of the full sky model

	    rgb2hsv16(p->r_err,p->g_err,p->b_err,&p->h_of_e,&p->s_of_e,&p->v_of_e) ;
	    // v_of_e will be the max of (r_err,g_err,b_err) now

	    if(p->s_of_e > .01) p->s_of_e = .01 ;

	    p->prob_sky = (1.-p->v_of_e)*(1.-p->s_of_e*p->s_of_e) ;
	} else if(1) {
	    float wgt_r = p->r_hat ;
	    float wgt_g = p->g_hat ;
	    float wgt_b = p->b_hat ;

	    float wgted_err = (
		wgt_r*(float)p->r_err/(float)MAX16 +
		wgt_g*(float)p->g_err/(float)MAX16 +
		wgt_b*(float)p->b_err/(float)MAX16) / (wgt_r + wgt_g + wgt_b) ;

	    p->prob_sky= 1. - wgted_err ;
	} else {

	    p->prob_sky= 1. - (float)p->b_err/(float)MAX16 ;
	}

	p->p_hat = 1./(1+exp(-(p->prob_sky-pData->full_sky_replacement_thresh)*pData->full_sky_replacement_ramp)) ;

	if(p->s < p->s_hat) {
	    if(I > I_hat && p->s < p->s_hat) {
		// if the actual is less saturated, with higher intensity, it could be a cloud
		p->p_hat *= I_hat/I * p->s/p->s_hat ;
	    } else {
		// if the actual is less saturated, with lower intensity, it could also be a cloud
		p->p_hat *= I/I_hat * p->s/p->s_hat ;
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
	    if(!pData->estimate_only && pData->column_mask[x] == 1) continue ;
	    int ymax = pData->show_raw_prediction ? pData->start_of_sky[x]+1 : pData->end_of_sky[x]  ;
/*  	    int ymax = pData->show_raw_prediction ? pData->start_of_sky[x]+1 : IMAGE_HEIGHT  ;  */

	    for(y = 0 ; y < ymax ; y++) {
		// completely model sky color
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

		float h_hat,s_hat,v_hat ;
		uint16_t r_hat, g_hat, b_hat ;

		predict_sky_hsv(px, py, &h_hat, &s_hat, &v_hat) ;

/*  		s_hat *= final_saturation_factor ;  */

		hsv2rgb16(h_hat,s_hat,v_hat,&r_hat,&g_hat,&b_hat) ;

		tif_set4c(image,x,y, r_hat, g_hat, b_hat, MAX16);
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



	for(x = x0 ; x <= x1 ; x++) {
/*  	    if(pData->column_mask[x] == 1) continue ;  */
	    int ymax = get_end_of_sky(pData,x)  ;
	    if(ymax < pData->manual_end_of_sky[x]) ymax = pData->manual_end_of_sky[x] ;
	    if(ymax > pData->manual_end_of_sky[x] && pData->manual_end_of_sky[x] != -1) ymax = pData->manual_end_of_sky[x] ;

	    int feather_end_y ;
	    int feather_length = compute_feather_length_with_eos(x, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length, pData, ymax) ;

/*  	    ymax = IMAGE_HEIGHT-1 ;  */

	    for(y = 0 ; y < ymax ; y++) {
		// completely model sky color
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

		if(py_adjusted_for_horizon_curvature(px,py) < pData->horizon_py+.005)
		    break ;

		// if we are near the sun, but below the start of sky, incorporation
		// more of the base image so things like sun flare can come through
		float sun_prob = 0. ;
		float py_dist = -1. ;
		float px_dist = -1. ;

		if(uses_CIE_model) {
		    py_dist = (py-pData->sun_py) ;
		    px_dist = (px-sun_px) ;

		    // weight of perez_F at sun center, falls off with square of distance
		    sun_prob = pData->perez_F/(1.+pData->perez_G*(px_dist*px_dist+py_dist*py_dist)) ;
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

			// beyond depth of fill, make sure p_hat goes to 0 at end of sky
			if(y > feather_end_y && pData->noclip_end_of_sky)
			    sp.p_hat *= 1. - (float)(y-feather_end_y)/(ymax-feather_end_y) ;

			// increase p_hat as x,y gets more near to the sun
			if(sp.p_hat < sun_prob) sp.p_hat = sun_prob ;

		    } else {
			// if p_hat is low, do not change it
			float raw_p_hat = sp.p_hat ;
			if(sp.p_hat < feather_prob) sp.p_hat = feather_prob ;
			if(sp.p_hat < sun_prob) sp.p_hat = sun_prob ;

			// now correct based on raw_p_hat
			// if raw_p_hat is 1., use all of the new value
			// if raw_p_hat is 0., use all of the initial value
			sp.p_hat = raw_p_hat * (sp.p_hat) + (1.-raw_p_hat)*raw_p_hat ;
		    }

		    sp.r_hat = (uint16_t)(sp.p_hat*(float)sp.r_hat + (1.-sp.p_hat)*(float)sp.r +.5) ;
		    sp.g_hat = (uint16_t)(sp.p_hat*(float)sp.g_hat + (1.-sp.p_hat)*(float)sp.g +.5) ;
		    sp.b_hat = (uint16_t)(sp.p_hat*(float)sp.b_hat + (1.-sp.p_hat)*(float)sp.b +.5) ;
		}

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
		    tif_set4c(image,x,y, sp.r_hat, sp.g_hat, sp.b_hat, MAX16);
		} else if(pData->show_sky_prediction) {
		    uint16_t val = (uint16_t)(sp.p_hat * (float)MAX16) ;
/*  		    uint16_t val = (uint16_t)(feather_prob * (float)MAX16) ;  */
		    tif_set4c(image,x,y, val, val, val, MAX16);
/*  		    uint16_t val = (uint16_t)(hsv_err*MAX16) ;  */
/*  		    tif_set4c(image,x,y, val, val, val, MAX16);  */
		} else {
/*  		    r_err = (uint16_t)(sqrt((float)r_err/(float)MAX16)*(float)MAX16+0.5) ;  */
/*  		    g_err = (uint16_t)(sqrt((float)g_err/(float)MAX16)*(float)MAX16+0.5) ;  */
/*  		    b_err = (uint16_t)(sqrt((float)b_err/(float)MAX16)*(float)MAX16+0.5) ;  */
		    tif_set4c(image,x,y, sp.r_err, sp.g_err, sp.b_err, MAX16);
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
