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

void repair_sky_hue(tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky, SKYFILL_DATA_t *pData)
{
    int x, y ;

    fprintf(stderr, "repair sky hue\n") ;

    // need some correction factors
    double s_coefs_left[2] ;
    double h_coefs_left[2] ;
    double s_coefs_right[2] ;
    double h_coefs_right[2] ;

    float sum_s=0. , sum_h=0. ;
    float sum_shat=0. , sum_hhat=0. ;
    int n=0, n_l=0, n_r=0 ;

    struct mstat m_fit_s_l = init_reg(2) ; // will be a simple linear in y
    struct mstat m_fit_h_l = init_reg(2) ; // will be a simple linear in y
    struct mstat m_fit_s_r = init_reg(2) ; // will be a simple linear in y
    struct mstat m_fit_h_r = init_reg(2) ; // will be a simple linear in y

    fprintf(stderr, "Repair sky hue from %d to %d\n", pData->min_sky_hue_mask, pData->max_sky_hue_mask) ;

    for(x=pData->min_sky_hue_mask-10 ; x < pData->min_sky_hue_mask ; x ++) {
	if(x < 0) continue ;
	if(x > IMAGE_WIDTH-1) continue ;
	if(pData->column_mask[x] == 1) continue ;

	for(y = start_of_sky[x] ; y <= end_of_sky[x] ; y++) {
	    if(y > IMAGE_HEIGHT-1) break ;

	    uint16_t rgb[3] ;
	    tif_get3cv(image,x,y,rgb) ;

	    float hsv[3] ;
	    rgb2hsv16(rgb,hsv) ;

	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    float hsv_hat[3] ;
	    predict_sky_modelled_huesat_from_val(hsv[2], hsv_hat, px, py) ;

	    if(hsv_hat[1] > .1) {
		sum_s += hsv[1] ;
		sum_h += hsv[0] ;
		double xv[2] = {1.,y} ;
		sum_reg(&m_fit_s_l, xv, hsv[1]/hsv_hat[1], 1.) ;
		sum_reg(&m_fit_h_l, xv, hsv[0]/hsv_hat[0], 1.) ;
		sum_shat += hsv_hat[1] ;
		sum_hhat += hsv_hat[0] ;

		n_l++ ;

		n++ ;
	    }
	}
    }

    for(x=pData->max_sky_hue_mask ; x < pData->max_sky_hue_mask+10 ; x ++) {
	if(x < 0) continue ;
	if(x > IMAGE_WIDTH-1) break ;
	if(pData->column_mask[x] == 1) continue ;

	for(y = start_of_sky[x] ; y <= end_of_sky[x] ; y++) {
	    if(y > IMAGE_HEIGHT-1) break ;

	    uint16_t rgb[3] ;
	    tif_get3cv(image,x,y,rgb) ;

	    float hsv[3] ;
	    rgb2hsv16(rgb,hsv) ;

	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    float hsv_hat[3] ;
	    predict_sky_modelled_huesat_from_val(hsv[2], hsv_hat, px, py) ;

	    if(hsv_hat[1] > .1) {
		sum_s += hsv[1] ;
		sum_h += hsv[0] ;
		double xv[2] = {1.,y} ;
		sum_reg(&m_fit_s_r, xv, hsv[1]/hsv_hat[1], 1.) ;
		sum_reg(&m_fit_h_r, xv, hsv[0]/hsv_hat[0], 1.) ;
		sum_shat += hsv_hat[1] ;
		sum_hhat += hsv_hat[0] ;

		n_r++ ;

		n++ ;
	    }
	}
    }

/*      print_reg_matrix(&m_fit_s_l, "m_fit_h_l") ;  */
/*      print_reg_matrix(&m_fit_s_l, "m_fit_s_l") ;  */
/*      print_reg_matrix(&m_fit_s_r, "m_fit_s_r") ;  */

    float s_correct=1. ;
/*      float h_correct=1. ;  */

    if(n > 0) {
	s_correct = sum_s/sum_shat ;
/*  	h_correct = sum_h/sum_hhat ;  */
    }

    if(n_l > 2) {
	estimate_reg(&m_fit_s_l, s_coefs_left) ;
	estimate_reg(&m_fit_h_l, h_coefs_left) ;
    } else {
	s_coefs_left[0] = s_correct ;
	s_coefs_left[1] = 0.0 ;
    }

    if(n_r > 2) {
	estimate_reg(&m_fit_s_r, s_coefs_right) ;
	estimate_reg(&m_fit_h_r, h_coefs_right) ;
    } else {
	s_coefs_right[0] = s_correct ;
	s_coefs_right[1] = 0.0 ;
    }

    if(n_l > 2 && n_r <= 2) {
	s_coefs_right[0] = s_coefs_left[0] ;
	s_coefs_right[1] = s_coefs_left[1] ;
	h_coefs_right[0] = h_coefs_left[0] ;
	h_coefs_right[1] = h_coefs_left[1] ;
    }

    if(n_r > 2 && n_l <= 2) {
	s_coefs_left[0] = s_coefs_right[0] ;
	s_coefs_left[1] = s_coefs_right[1] ;
	h_coefs_left[0] = h_coefs_right[0] ;
	h_coefs_left[1] = h_coefs_right[1] ;
    }

/*      fprintf(stderr, "repair sky hue, h,s correct = %f, %f\n", h_correct, s_correct) ;  */
/*      fprintf(stderr, "repair sky hue, h_coefs(l,r) = (%f, %f) (%f,%f)\n", h_coefs_left[0], h_coefs_left[1], h_coefs_right[0], h_coefs_right[1]) ;  */
/*      fprintf(stderr, "repair sky hue, s_coefs(l,r) = (%f, %f) (%f,%f)\n", s_coefs_left[0], s_coefs_left[1], s_coefs_right[0], s_coefs_right[1]) ;  */



    pData->n_sky_hue_fixed_clipped = 0 ;

    for(x=pData->min_sky_hue_mask ; x < pData->max_sky_hue_mask ; x ++) {
	uint16_t rgb[3] ;
	if(pData->column_mask[x] == 1) continue ;

	float pleft = 1.-(float)(x-pData->min_sky_hue_mask)/(float)(pData->max_sky_hue_mask-pData->min_sky_hue_mask) ;


	for(y = start_of_sky[x] ; y <= end_of_sky[x] ; y++) {
	    if(y >= IMAGE_HEIGHT) break ;
	    tif_get3cv(image,x,y,rgb) ;
	    uint16_t r0 = rgb[0], g0=rgb[1], b0=rgb[2] ;

	    float hsv[3] ;
	    rgb2hsv16(rgb,hsv) ;

	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    predict_sky_modelled_huesat_from_val(hsv[2], hsv, px, py) ;

	    //compute correction factors on left and right side
	    float h_correct_l = h_coefs_left[0] + h_coefs_left[1]*(float)y ;
	    float s_correct_l = s_coefs_left[0] + s_coefs_left[1]*(float)y ;
	    float h_correct_r = h_coefs_right[0] + h_coefs_right[1]*(float)y ;
	    float s_correct_r = s_coefs_right[0] + s_coefs_right[1]*(float)y ;

	    // interpolate correction factors
	    float pixel_h_correct = pleft*h_correct_l + (1.-pleft)*h_correct_r ;
	    float pixel_s_correct = pleft*s_correct_l + (1.-pleft)*s_correct_r ;


/*  	    if(x == pData->min_sky_hue_mask) {  */
/*  		fprintf(stderr, "repair sky hue, y=%d, h,s correct = %f, %f\n", y, pixel_h_correct, pixel_s_correct) ;  */
/*  	    }  */

	    hsv[0] *= pixel_h_correct ;
	    hsv[1] *= pixel_s_correct ;

	    hsv2rgb16(hsv,rgb) ;
	    int dr = (int)rgb[0]-(int)r0 ;
	    int dg = (int)rgb[1]-(int)g0 ;
	    int db = (int)rgb[2]-(int)b0 ;
	    float v_adjust=1. ;
	    // find the component that dropped the most
	    // and bring it back to the original value, bring other components up by the same factor
	    // reasoning is the sky was blown out, and a component with a nonclipped value > 1.0, was
	    // clipped to 1.0, changing the sky hue.  Using this logic will allow that component to
	    // be restored to its original value
	    // note -- if the image values are already too high, lowering the exposure with the -ef flag
	    // may be needed so value are not clipped to 1.0 again...
	    if(dr < dg) {
		if(dr < db)
		    v_adjust = (float)r0/(float)rgb[0] ;
		else
		    v_adjust = (float)b0/(float)rgb[2] ;
	    } else {
		if(dg < db)
		    v_adjust = (float)g0/(float)rgb[1] ;
		else
		    v_adjust = (float)b0/(float)rgb[2] ;
	    }

/*  	    v_adjust =( (float)r0/(float)r+ (float)g0/(float)g+ (float)b0/(float)b) / 3.0 ;  */

	    int32_t r32 = (uint32_t)((float)rgb[0] * v_adjust + 0.5) ;
	    int32_t g32 = (uint32_t)((float)rgb[1] * v_adjust + 0.5) ;
	    int32_t b32 = (uint32_t)((float)rgb[2] * v_adjust + 0.5) ;
	    int clipped=0 ;

	    if(r32 > MAX16) {
		r32=MAX16 ;
		clipped=1 ;
	    }
	    if(g32 > MAX16) {
		g32=MAX16 ;
		clipped=1 ;
	    }
	    if(g32 > MAX16) {
		g32=MAX16 ;
		clipped=1 ;
	    }

	    pData->n_sky_hue_fixed_clipped += clipped ;

	    tif_set3c(image,x,y,r32,g32,b32) ;

	}
    }
}

void repair_alpha(tdata_t *image, SKYFILL_DATA_t *pData)
{
    // hugin is creating pixels near the edges of the data with low r,g,b values
    // presumeably because they didn't represent a full pixel from an original
    // image but are part of a pixel.  These often also have really bad colors,
    // the max 16 bit alpha is correct at 65535
    // out of hugin, so we need to roll through the entire image and set the alpha
    // and r,g,b of those suspect pixels to 0,0,0,0
    int x,y ;

    for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	for(y = 0 ; y < IMAGE_HEIGHT ; y++)
/*  	for(y = IMAGE_HEIGHT-1 ; y > -1 ; y--)  */
	{
	    uint16_t r,g,b,a ;
	    tif_get4c(image,x,y,r,g,b,a) ;
	    if(a < MAX16) {
		tif_set4c(image,x,y,0,0,0,0) ;
		continue ;
	    }

	    // nearest neighbor search +/- 2 pixels to find max r,g,b in area
	    uint16_t max_r=0 ;
	    uint16_t max_g=0 ;
	    uint16_t max_b=0 ;
	    float sum_r=0., sum_g=0., sum_b=0. ;
	    float n_neighbor ;
	    int x0 = x-2 ; if(x0 < 0) x0 = 0 ;
	    int y0 = y-2 ; if(y0 < 0) y0 = 0 ;
	    int x1 = x+2 ; if(x1 > IMAGE_WIDTH-1) x1 = IMAGE_WIDTH-1 ;
	    int y1 = y+2 ; if(y1 > IMAGE_HEIGHT-1) y1 = IMAGE_HEIGHT-1 ;
	    int found_edge=0;  // are we near the edge of image data (i.e. there is a transparent pixel in the neighborhood?) ;

	    for(int dx=x0 ; dx<=x1 ; dx++) {
		for(int dy=y0 ; dy<=y1 ; dy++) {
		    uint16_t r_neighbor,g_neighbor,b_neighbor,a_neighbor ;
		    tif_get4c(image,dx,dy,r_neighbor,g_neighbor,b_neighbor,a_neighbor) ;

		    if(a_neighbor > 0) {
			if(max_r < r_neighbor) max_r = r_neighbor ;
			if(max_g < g_neighbor) max_g = g_neighbor ;
			if(max_b < b_neighbor) max_b = b_neighbor ;
			sum_r += r_neighbor ;
			sum_g += g_neighbor ;
			sum_b += b_neighbor ;
			n_neighbor += 1. ;
		    } else {
			found_edge=1 ;
		    }
		}
	    }

	    if(found_edge && (r < max_r/8 || g < max_g/8 || b < max_b/8)) {
		// this is a suspect pixel on the edge, with at least one very low color component
		tif_set4c(image,x,y,0,0,0,0) ;
/*  		tif_set4c(image,x,y,sum_r/n_neighbor,sum_r/n_neighbor,sum_r/n_neighbor,MAX16) ;  */
	    }
	}
    }

}

void set_minmax_sky_values(SKYFILL_DATA_t *pData)
{
    pData->max_end_of_sky=0 ;
    pData->min_end_of_sky=IMAGE_HEIGHT-1 ;
    pData->max_start_of_sky=0 ;
    pData->min_start_of_sky=IMAGE_HEIGHT-1 ;

    for(int x = 0 ; x < IMAGE_WIDTH ; x++) {
	if(pData->column_mask[x] == 1) continue ;
	if(pData->max_end_of_sky < pData->end_of_sky[x]) pData->max_end_of_sky = pData->end_of_sky[x] ;
	if(pData->max_start_of_sky < pData->start_of_sky[x]) pData->max_start_of_sky = pData->start_of_sky[x] ;
	if(pData->min_end_of_sky > pData->end_of_sky[x]) pData->min_end_of_sky = pData->end_of_sky[x] ;
	if(pData->min_start_of_sky > pData->start_of_sky[x]) pData->min_start_of_sky = pData->start_of_sky[x] ;
    }
}

#ifdef GMR_DEBUG
int gmr_debug = 0 ;
#endif

int get_mean_rgb(tdata_t *image, int xc,int yc,double rcoef[],double gcoef[], double bcoef[],int search_width, int clip_eos_flag, int repair_bad_points_flag, int recursion_level, SKYFILL_DATA_t *pData)
{
    int x0 = xc-search_width ;
    int x1 = xc+search_width ;

    if(x0 < 0) x0 = 0 ;
    if(x1 < 0) x1 = 0 ;
    if(x0 > IMAGE_WIDTH-1) x0 = IMAGE_WIDTH-1 ;
    if(x1 > IMAGE_WIDTH-1) x1 = IMAGE_WIDTH-1 ;

    rcoef[0] =  rcoef[1] =  rcoef[2] =  rcoef[3] = 0. ;
    gcoef[0] =  gcoef[1] =  gcoef[2] =  gcoef[3] = 0. ;
    bcoef[0] =  bcoef[1] =  bcoef[2] =  bcoef[3] = 0. ;

    int x,y ;
    int increment=2 ;

#ifdef GMR_DEBUG
    if(gmr_debug) {
	fprintf(stderr, "GMR: x,y:%d,%d\n", xc,yc) ;
    }
#endif

    // first find the coordinates of valid point and store them
    int max_points = (2*search_width)*(2*search_width)/2+16 ;

    if(repair_bad_points_flag) {
	max_points *= 4 ;
	increment=1 ;
    }

    uint16_t *x_possible_coords, *y_possible_coords, n_possible_points ;
    uint16_t *x_coords, *y_coords, n_points ;
    uint8_t *point_status ;

    n_points=0 ;
    point_status = (uint8_t *)calloc(max_points, sizeof(uint8_t)) ;
    x_coords = (uint16_t *)calloc(max_points, sizeof(uint16_t)) ;
    y_coords = (uint16_t *)calloc(max_points, sizeof(uint16_t)) ;

    n_possible_points=0 ;
    x_possible_coords = (uint16_t *)calloc(max_points, sizeof(uint16_t)) ;
    y_possible_coords = (uint16_t *)calloc(max_points, sizeof(uint16_t)) ;

    // look through all possible pixels
    for(x = x0 ; x <= x1 ; x += increment) {
	if(pData->column_mask[x] == 1) continue ;

	int y1 = yc+search_width ;
	if(y1 < pData->start_of_sky[x]) continue ;

	int y0 = yc-search_width ;
	if(y0 < pData->start_of_sky[x]) y0 = pData->start_of_sky[x] ;

	if(clip_eos_flag && y1 > pData->end_of_sky[x]) y1 = pData->end_of_sky[x] ;

#ifdef GMR_DEBUG
	if(gmr_debug) {
	    fprintf(stderr, "GMR: col search x:%d, y0,y1:%d,%d\n", x,y0,y1) ;
	}
#endif

	for(y = y0 ; y <= y1 ; y += increment) {
	    /* if alpha channel is small then hugin indicates no image data here */
#ifdef GMR_DEBUG
	    if(gmr_debug) {
		uint16_t r,g,b,a ;
		tif_get4c(image,x,y,r,g,b,a) ;
		fprintf(stderr, "GMR: pixel x,y:%d,%d rgba=%5d,%5d,%5d,%5d\n", x,y,r,g,b,a) ;
	    }
#endif
/*  	    XXX  */
	    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+3] < MAX16)
		continue ;

/*  	    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+0] <= BLK16) continue ;  */
/*  	    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+1] <= BLK16) continue ;  */
/*  	    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+2] <= BLK16) continue ;  */

	    if(n_possible_points == max_points) {
		fprintf(stderr, "FATAL: n_possible_points > max_points(%d) in get_mean_rgb(), search_width=%d\n", n_possible_points, search_width) ;
		exit(1) ;
	    }

#ifdef GMR_DEBUG
	    if(gmr_debug) {
		fprintf(stderr, "GMR: ADDED pixel x,y:%d,%d\n", x,y) ;
	    }
#endif

	    x_possible_coords[n_possible_points] = x ;
	    y_possible_coords[n_possible_points] = y ;
	    point_status[n_possible_points] = 1 ;
	    n_possible_points++ ;
	}
    }

    if(n_possible_points <  16) {
/*  	fprintf(stderr, "WARNING: n_possible_points is < 16 in get_mean_rgb(), x,y=%d,%d search_width=%d\n", xc, yc, search_width) ;  */

	free(point_status) ;
	free(x_possible_coords) ;
	free(y_possible_coords) ;
	free(x_coords) ;
	free(y_coords) ;

	recursion_level-- ;

	if(recursion_level <= 0)
	    return -1 ;

	search_width *= 2 ;
/*  	fprintf(stderr, "         increasing search_width to %d\n", search_width) ;  */

	return get_mean_rgb(image, xc,yc,rcoef,gcoef, bcoef,search_width, clip_eos_flag, repair_bad_points_flag, recursion_level, pData) ;
    }

    int c ;

    // compute wegihted mean r,g,b to help look for outliers
    // weight is distance from center pixel
    float sum_wgts = 0. ;
    float sum[3] = {0.,0.,0.} ;
    float mean[3] = {0.,0.,0.} ;

    for(int i=0 ; i < n_possible_points ; i++) {
	x = x_possible_coords[i] ;
	y = y_possible_coords[i] ;
	float dx = fabsf(x-xc) ;
	float dy = fabsf(y-yc) ;
	float wgt = 1./(dx+dy+1) ;
	sum[0] += wgt*(float)((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+0] ; // R
	sum[1] += wgt*(float)((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+1] ; // G
	sum[2] += wgt*(float)((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+2] ; // B

	sum_wgts += wgt ;
    }

    mean[0] = sum[0] / sum_wgts ;
    mean[1] = sum[1] / sum_wgts ;
    mean[2] = sum[2] / sum_wgts ;

#ifdef GMR_DEBUG
    if(gmr_debug) {
	fprintf(stderr, "\nGMR: mean rgb %5d,%5d,%5d\n", (int)mean[0],(int)mean[1],(int)mean[2]) ;
    }
#endif

    int min_x=IMAGE_WIDTH-1 ;
    int max_x=0 ;
    int min_y=IMAGE_HEIGHT-1 ;
    int max_y=0 ;

    // look at differences from the mean to throw out bad colored pixels (and mark them needing repair
    for(int i=0 ; i < n_possible_points ; i++) {
	x = x_possible_coords[i] ;
	y = y_possible_coords[i] ;

	float r,g,b ;
	tif_get3c(image,x,y,r,g,b) ;

	point_status[i] = 0 ; // this point might need repair

	float r_err = fabsf( (r-mean[0])/mean[0]) ;
	float g_err = fabsf( (g-mean[1])/mean[1]) ;
	float b_err = fabsf( (b-mean[2])/mean[2]) ;
#ifdef GMR_DEBUG
	fprintf(stderr, "GMR:      err x,y:%d,%d %f %f %f\n", x,y,r_err,g_err,b_err) ;
#endif

	// throw out point of different by more than desired tolerance from mean ;
	if(r_err > pData->repair_tos_thresh) continue ;
	if(g_err > pData->repair_tos_thresh) continue ;
	if(b_err > pData->repair_tos_thresh) continue ;

	point_status[i] = 1 ; // mark this point as not needing repair, if requested

	x_coords[n_points] = x ;
	y_coords[n_points] = y ;
	n_points++ ;

	if(min_x > x) min_x=x ;
	if(max_x < x) max_x=x ;
	if(min_y > y) min_y=y ;
	if(max_y < y) max_y=y ;
    }

#ifdef GMR_DEBUG
    if(gmr_debug) {
	fprintf(stderr, "\nGMR: n_possible,n %d,%d\n", n_possible_points, n_points) ;
	if(xc == 1705 && yc == 91)
	    exit(1) ;
    }
#endif

    if(n_points < 8) {
	free(point_status) ;
	free(x_possible_coords) ;
	free(y_possible_coords) ;
	free(x_coords) ;
	free(y_coords) ;

	recursion_level-- ;

	if(recursion_level <= 0) {
	    if(pData->verbose) fprintf(stderr, "WARNING: number of fitting points is < 8 in get_mean_rgb, coefs are all 0.0!\n") ;
	    return -1 ;
	}
	if(pData->verbose) fprintf(stderr, "WARNING: number of fitting points is < 8 in get_mean_rgb, doubling search width\n") ;

	search_width *= 2 ;
/*  	fprintf(stderr, "         increasing search_width to %d\n", search_width) ;  */

	return get_mean_rgb(image, xc,yc,rcoef,gcoef, bcoef,search_width, clip_eos_flag, repair_bad_points_flag, recursion_level, pData) ;
    }

// fit types, FIT_X_Y2 is color=B0 + B1*x + B2*y + B3*y*y ;
// fit types, FIT_X_Y1 is color=B0 + B1*x + B2*y + 0 *y*y ;
// fit types, FIT_0_Y2 is color=B0 + 0 *x + B2*y + B3*y*y ;
// fit types, FIT_0_Y1 is color=B0 + 0 *x + B2*y + 0 *y*y ;
// fit types, FIT_X_0  is color=B0 + B1*x + 0 *y + 0 *y*y ;
// fit types, FIT_0_0  is color=B0 + 0 *x + 0 *y + 0 *y*y ;

// note the first 2 places are the order of the "X" coord used, places 3&4 are the order of the "Y" coord used
#define FIT_X_Y2 0x0102
#define FIT_X_Y1 0x0101
#define FIT_X_Y0 0x0101
#define FIT_0_Y2 0x0002
#define FIT_0_Y1 0x0001
#define FIT_0_Y0 0x0000

    int fit_type=0x0000 ;

    if((max_x-min_x) > 0) fit_type |= 0x0100 ;

    if((max_y-min_y) > 1)
	fit_type |= 0x0002 ;
    else if((max_y-min_y) > 0)
	fit_type |= 0x0001 ;

    int order_x = (fit_type&0xff00) >> 8 ;
    int order_y = (fit_type&0x00ff) ;

    int total_vars = order_x + order_y ;

    struct mstat fit_struct[3] ;

    for(c=0 ; c < 3 ; c++) {
	fit_struct[c] = init_reg(total_vars+1) ; // +1 for intercept
    }

    for(int i=0 ; i < n_points ; i++) {
	x = x_coords[i] ;
	y = y_coords[i] ;

	double xv[4] ;
	int j=0 ;
	{ xv[j] = 1. ; j++; }
	if(order_x) { xv[j] = x ; j++; }
	if(order_y) { xv[j] = y ; j++; }
	if(order_y > 1) { xv[j] = y*y ; j++; }

	uint16_t pixel[3] ;
	tif_get3c(image,x,y,pixel[0],pixel[1],pixel[2]) ;

	sum_reg((&fit_struct[0]), xv, (double)pixel[0], 1.) ;
	sum_reg((&fit_struct[1]), xv, (double)pixel[1], 1.) ;
	sum_reg((&fit_struct[2]), xv, (double)pixel[2], 1.) ;
    }

    double fit_coefs[4] ;
    double *coef_dest[3] = {rcoef,gcoef,bcoef} ;

    for(c=0 ; c < 3 ; c++) {
	double *coefs = coef_dest[c] ;
	estimate_reg((&fit_struct[c]), fit_coefs) ;

	int j=0 ;
	int k=0 ;
	coefs[0] = fit_coefs[0] ; //B0 is always used
	j++ ;
	k++ ;

	if(order_x) {  coefs[j] = fit_coefs[k] ; k++ ; j++; } else { coefs[j] = 0.0 ; j++;}
	if(order_y) {  coefs[j] = fit_coefs[k] ; k++ ; j++; } else { coefs[j] = 0.0 ; j++;}
	if(order_y > 1) {  coefs[j] = fit_coefs[k] ; k++ ; j++; } else { coefs[j] = 0.0 ; j++;}
    }

#ifdef GMR_DEBUG
    if(gmr_debug) {
	fprintf(stderr, "GMR: rcoef %f,%f,%f,%f\n", rcoef[0], rcoef[1], rcoef[2], rcoef[3]) ;
	fprintf(stderr, "GMR: gcoef %f,%f,%f,%f\n", gcoef[0], gcoef[1], gcoef[2], gcoef[3]) ;
	fprintf(stderr, "GMR: bcoef %f,%f,%f,%f\n", bcoef[0], bcoef[1], bcoef[2], bcoef[3]) ;
    }
#endif

    int n_repaired=0 ;
    if(repair_bad_points_flag == 1) {

	for(int i=0 ; i < n_possible_points ; i++) {

	    if(point_status[i] == 1)
		continue ;

	    x = x_possible_coords[i] ;

	    // pixels in column masks have already been eliminated in the very first loop
	    // now see if the repair mask is set for this column
	    if(pData->column_repair_mask[x] == 1) continue ;
	   
	    y = y_possible_coords[i] ;

	    float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
	    float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
	    float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;
#ifdef GMR_DEBUG
	    if(gmr_debug) {
		fprintf(stderr, "   Set pixel %d,%d to  %5d,%5d,%5d\n", x,y,(int)rhat,(int)bhat,(int)ghat) ;
	    }
#endif

	    if(x == 1721 && y == 54) {
		fprintf(stderr, "********************************************************\n") ;
		fprintf(stderr, "   Set pixel %d,%d to  %5d,%5d,%5d\n", x,y,(int)rhat,(int)bhat,(int)ghat) ;
		fprintf(stderr, "********************************************************\n") ;
	    }

	    tif_set4c(image,x,y,rhat,ghat,bhat,MAX16) ;

	    n_repaired++ ;
	}
    }

    free(point_status) ;
    free(x_possible_coords) ;
    free(y_possible_coords) ;
    free(x_coords) ;
    free(y_coords) ;
#ifdef GMR_DEBUG
    if(gmr_debug) {
	if(n_possible_points == 121 && n_points == 0) exit(1) ;
    }
#endif

    return n_repaired ;
}

// adjust this as the number of phases change!!
#define MAX_PHASES 10

void repair_top_of_sky(tdata_t *image, int end_of_sky_is_known_flag, int phase, SKYFILL_DATA_t *pData)
{
    int x, y ;
    fprintf(stderr, "RTOS::: sky hsv (%3.0f,%4.2f,%4.2f)\n", pData->hue_sky, pData->sat_sky, pData->val_sky) ;

    if(phase == 1) {
	// check "edges of found sky for issues:
	// 1. are there black pixels (with alpha near MAX16 ?
	fprintf(stderr, "RTOS eliminate black pixels at top of sky\n") ;

	for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    int n_ok = 0 ;
	    for(y = pData->raw_start_of_sky[x]-1 ; y > 0 ; y--) {
		if(xy_is_black_pixel(image,x,y) == 1) {
		    fprintf(stderr, "RTOS found black pixel at %d,%d\n",x,y) ;
		    ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+3] = 0 ;
		    n_ok=0 ;
		} else if(n_ok > 10) {
		    break ;
		} else {
		    n_ok++ ;
		}
	    }
	}

	fprintf(stderr, "RTOS eliminate black pixels on vertical edges of sky\n") ;
	set_minmax_sky_values(pData) ;

	for(y = pData->min_start_of_sky ; y < pData->max_start_of_sky+1 ; y++) {

	    // left to right
	    int max_x=IMAGE_WIDTH-1 ;
	    for(x = 0 ; x < max_x ; x++) {
		if(xy_is_transparent(image,x,y) && xy_is_black_pixel(image,x+1,y) == 1) {
		    if(pData->column_mask[x+1] == 1) continue ;
		    if(pData->verbose) fprintf(stderr, "RTOS LR fix x,y %d,%d\n", x,y) ;
		    ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*(x+1)+3] = 0 ;
		    max_x = x+10 ;
		    if(max_x > IMAGE_WIDTH-1) max_x = IMAGE_WIDTH-1 ;
		}
	    }

	    // right to left
	    int min_x =0 ;
	    for(x = IMAGE_WIDTH-1 ; x > min_x ; x--) {
		if(xy_is_transparent(image,x,y) && xy_is_black_pixel(image,x-1,y) == 1) {
		    if(pData->column_mask[x-1] == 1) continue ;
		    if(pData->verbose) fprintf(stderr, "RTOS RL fix x,y %d,%d\n", x,y) ;
		    ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*(x-1)+3] = 0 ;
		    min_x = x-10 ;
		    if(min_x < 0) min_x = 0 ;
		}
	    }

	}
    }

    if(phase == 3) {
	double rcoef[4], gcoef[4], bcoef[4] ;

	fprintf(stderr, "RTOS clean up pixels on vertical edges of sky, EOS is known %d\n", end_of_sky_is_known_flag) ;
	set_minmax_sky_values(pData) ;

	for(y = pData->min_start_of_sky ; y < pData->max_start_of_sky+1 ; y++) {

	    // left to right
	    int max_x=IMAGE_WIDTH-1 ;
	    for(x = 0 ; x < max_x ; x++) {
		if(xy_is_transparent(image,x,y) && xy_is_transparent(image,x+1,y) == 0) {
		    if(pData->column_mask[x+1] == 1) continue ;
		    if(pData->verbose) fprintf(stderr, "RTOS LR clean up x,y %d,%d\n", x+1,y) ;
		    int search_width=2 ;
		    int n_repaired = get_mean_rgb(image,x,y,rcoef,gcoef,bcoef,search_width,0,1,0, pData) ;

		    if(n_repaired > 0)
			if(pData->verbose) fprintf(stderr, "Bad pixels near top of sky: Repaired %d pixel colors near %d,%d\n", n_repaired, x+1, y) ;
		}
	    }

	    // right to left
	    int min_x =0 ;
	    for(x = IMAGE_WIDTH-1 ; x > min_x ; x--) {
		if(xy_is_transparent(image,x,y) && xy_is_transparent(image,x-1,y) == 0) {
		    if(pData->column_mask[x-1] == 1) continue ;
		    if(pData->verbose) fprintf(stderr, "RTOS RL cleanup x,y %d,%d\n", x,y) ;
		    int search_width=2 ;
		    int n_repaired = get_mean_rgb(image,x-1,y,rcoef,gcoef,bcoef,search_width,0,1,1, pData) ;

		    if(n_repaired > 0)
			if(pData->verbose) fprintf(stderr, "Bad pixels near top of sky: Repaired %d pixel colors near %d,%d\n", n_repaired, x-1, y) ;

		}
	    }

	}
    }

    if(0 && phase == 4) {
	double rcoef[4], gcoef[4], bcoef[4] ;

	fprintf(stderr, "RTOS first scan for bad pixels in top 5 pixels of sky\n") ;


	// are there wrong pixels in the top 5 of the sky?

	int xc = 0 ;
	while(1) {
	    while(pData->column_mask[xc] == 1 && xc < IMAGE_WIDTH)
		xc++ ;

	    if(xc >= IMAGE_WIDTH)
		break ;

	    int search_width = 5 ;
	    int yc = pData->start_of_sky[xc]+5 ;

	    int n_repaired = get_mean_rgb(image,xc,yc,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,1,1, pData) ;

	    if(n_repaired > 0)
		if(pData->verbose) fprintf(stderr, "Bad pixels near top of sky: Repaired %d pixel colors near %d,%d\n", n_repaired, xc, yc) ;
	    xc++ ;
	}

    }

    if(phase == 5) {
	// fill in vertical gaps where alpha is 0, happens
	// when there is a mask in the image that hugin used that
	// appears on the edge

	// first, if there are gaps at the edges of the image,
	// these can extend below the actual sky, so they need
	// to be detected.
	int left_edge_max_x=-1 ;

	for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	    int found=0 ;
	    for(y = pData->raw_start_of_sky[x] ; y <(IMAGE_HEIGHT*3)/4 ; y++) {
		if(xy_is_opaque_pixel(image,x,y) == 0) {
		    left_edge_max_x=x ;
		    found=1 ;
		    break ;
		}
	    }
	    if(found == 0)
		break ;
	}

	int right_edge_min_x=IMAGE_WIDTH ;

	for(x = IMAGE_WIDTH-1 ; x > -1 ; x--) {
	    int found=0 ;
	    for(y = pData->raw_start_of_sky[x] ; y <(IMAGE_HEIGHT*3)/4 ; y++) {
		if(xy_is_opaque_pixel(image,x,y) == 0) {
		    right_edge_min_x=x ;
		    found=1 ;
		    break ;
		}
	    }
	    if(found == 0)
		break ;
	}

	fprintf(stderr, "RTOS fill vertical gap,  leftedgemask:%d,  rightedgemask,%d\n", left_edge_max_x, right_edge_min_x) ;

	for(x = left_edge_max_x+1 ; x < right_edge_min_x-1 ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    int maxy = end_of_sky_is_known_flag ? pData->end_of_sky[x] : (IMAGE_HEIGHT*3)/4 ;
	    maxy = (IMAGE_HEIGHT*3)/4 ;

	    for(y = pData->raw_start_of_sky[x] ; y < maxy ; y++) {

		if(xy_is_opaque_pixel(image,x,y) == 0) {
		    int ty = y ;
		    int by ;

		    for(by = y+1 ; by < maxy ; by++) {
			if(xy_is_opaque_pixel(image,x,by) == 1)
			    break ;
		    }

		    if(by == maxy)
			continue ; // can't fix transparent pixels at bottom of image

		    fprintf(stderr, "\nRTOS fill vertical gap at x:%d,  y %d to %d\n", x, ty, by) ;

		    int search_width = 50 ;

		    // the top and bottom pixels may be badly colored pixels from edges
		    // of images that hugin output, correct those while getting prediction
		    // coefs
		    double rcoef0[4], gcoef0[4], bcoef0[4] ;
		    double rcoef1[4], gcoef1[4], bcoef1[4] ;
		    int nt = get_mean_rgb(image,x,ty,rcoef0,gcoef0,bcoef0,search_width,end_of_sky_is_known_flag,1,2, pData) ;
		    if(nt > 0)
			if(pData->verbose) fprintf(stderr, "Repaired %d pixels near the top of the gap\n",nt) ;
		    int nb = get_mean_rgb(image,x,by,rcoef1,gcoef1,bcoef1,search_width,end_of_sky_is_known_flag,1,2, pData) ;
		    if(nb > 0)
			if(pData->verbose) fprintf(stderr, "Repaired %d pixels near the top of the gap\n",nb) ;


		    for(y = ty ; y <= by ; y++) {
			float p1 = (float)(y-ty)/(by-ty) ;
			float p0 = 1.f - p1 ;

			float rhat0 = rcoef0[0] + rcoef0[1]*(float)x + rcoef0[2]*(float)y + rcoef0[3]*(float)y*(float)y ;
			float ghat0 = gcoef0[0] + gcoef0[1]*(float)x + gcoef0[2]*(float)y + gcoef0[3]*(float)y*(float)y ;
			float bhat0 = bcoef0[0] + bcoef0[1]*(float)x + bcoef0[2]*(float)y + bcoef0[3]*(float)y*(float)y ;

			float rhat1 = rcoef1[0] + rcoef1[1]*(float)x + rcoef1[2]*(float)y + rcoef1[3]*(float)y*(float)y ;
			float ghat1 = gcoef1[0] + gcoef1[1]*(float)x + gcoef1[2]*(float)y + gcoef1[3]*(float)y*(float)y ;
			float bhat1 = bcoef1[0] + bcoef1[1]*(float)x + bcoef1[2]*(float)y + bcoef1[3]*(float)y*(float)y ;

			float rhat = p0*rhat0 + p1*rhat1 ;
			float ghat = p0*ghat0 + p1*ghat1 ;
			float bhat = p0*bhat0 + p1*bhat1 ;

			tif_set4c(image,x,y, (uint16_t)rhat, (uint16_t)ghat, (uint16_t)bhat, MAX16 );
		    }
		}

		// the search continues for the next greater values of y
	    }
	}
    }

    if(0 && phase == 6) {
	// divot repair has issues, can create bad colors when interpolating/extraapolating long distances

	set_minmax_sky_values(pData) ;

	fprintf(stderr, "RTOS fill divots\n") ;

	// look for "divots" at the top of the sky to fill in
	// this is essentially filling horizontal gaps, where a gap is black pixels or alpha == 0 ;

	for(y = pData->max_start_of_sky ; y >= pData->min_start_of_sky ; y--) {
	    continue ;

	    int lx, rx ;

	    for(x = 0 ; x < IMAGE_WIDTH ;) {
		if(pData->column_mask[x] == 1) {
		    x++ ;
		    continue ;
		}

		// find lx, which is next x where SOS[x] <= y and SOS[x+1] > y
		lx=-1 ;

		while(x < IMAGE_WIDTH-2) {
		    if(pData->column_mask[x] == 1) {
			x++ ;
			continue ;
		    }
		    if(xy_has_nonblack_pixel(image, x, y) == 1)
			lx=x ;

		    if(pData->start_of_sky[x+1] > y) {
			break ;
		    }

		    x++ ;
		}

		if(lx > -1) {
		    // find rx, which is next x where SOS[x] <= y ;
		    rx=-1 ;
		    x=lx+1 ;
		    while(x < IMAGE_WIDTH-1) {
			if(pData->column_mask[x] == 1) {
			    x++ ;
			    continue ;
			}
			if(xy_has_nonblack_pixel(image, x, y) == 1) {
			    rx=x ;
			    break ;
			}
			x++ ;
		    }
		}

		if(lx > -1 && rx > -1) {
		    float width = rx-lx ;

		    if(pData->column_mask[rx]) {
			fprintf(stderr, "RTOS: Fatal rx=%d is masked column!\n", rx) ;
			exit(1) ;
		    }

		    if(pData->column_mask[lx]) {
			fprintf(stderr, "RTOS: Fatal lx=%d is masked column!\n", lx) ;
			exit(1) ;
		    }

		    if(width < IMAGE_WIDTH) {
			fprintf(stderr, "RTOS fill horizontal gap,  y:%d, l:%d, r:%d\n", y, lx, rx) ;

			double rcoef0[4], gcoef0[4], bcoef0[4] ;
			double rcoef1[4], gcoef1[4], bcoef1[4] ;

			int search_width=5 ;

			get_mean_rgb(image,lx,y,rcoef0,gcoef0,bcoef0,search_width,end_of_sky_is_known_flag,1,2, pData) ;
			get_mean_rgb(image,rx,y,rcoef1,gcoef1,bcoef1,search_width,end_of_sky_is_known_flag,1,2, pData) ;

			for(x = lx+1 ; x < rx ; x++) {
			    if(pData->column_mask[x] == 1) {
				continue ;
			    }

			    float p1 = (float)(x-lx)/width ;
			    float p0 = 1. - p1 ;

			    float rhat0 = rcoef0[0] + rcoef0[1]*(float)x + rcoef0[2]*(float)y + rcoef0[3]*(float)y*(float)y ;
			    float ghat0 = gcoef0[0] + gcoef0[1]*(float)x + gcoef0[2]*(float)y + gcoef0[3]*(float)y*(float)y ;
			    float bhat0 = bcoef0[0] + bcoef0[1]*(float)x + bcoef0[2]*(float)y + bcoef0[3]*(float)y*(float)y ;

			    float rhat1 = rcoef1[0] + rcoef1[1]*(float)x + rcoef1[2]*(float)y + rcoef1[3]*(float)y*(float)y ;
			    float ghat1 = gcoef1[0] + gcoef1[1]*(float)x + gcoef1[2]*(float)y + gcoef1[3]*(float)y*(float)y ;
			    float bhat1 = bcoef1[0] + bcoef1[1]*(float)x + bcoef1[2]*(float)y + bcoef1[3]*(float)y*(float)y ;

			    float rhat = p0*rhat0 + p1*rhat1 ;
			    float ghat = p0*ghat0 + p1*ghat1 ;
			    float bhat = p0*bhat0 + p1*bhat1 ;

			    tif_set4c(image,x,y, (uint16_t)rhat, (uint16_t)ghat, (uint16_t)bhat, MAX16 );

			    pData->start_of_sky[x] = y ;
			}
		    } else {
			x=rx ;
		    }

		} else {
		    x++ ;
		}
	    }

	}

    }

    if(phase == 7) {
	fprintf(stderr, "RTOS phase 7: fill side edges\n") ;
	// do the side edges need to be fixed?
	int edge_width = IMAGE_WIDTH/40 ;
	int min_start_y=IMAGE_HEIGHT ;
	int max_start_y=0 ;

	// LEFT side
	for(x = 0 ; x < edge_width ; x++) {
	    if(min_start_y > pData->start_of_sky[x]) min_start_y = pData->start_of_sky[x] ;
	    if(max_start_y < pData->start_of_sky[x]) max_start_y = pData->start_of_sky[x] ;
	}

	for(int y = max_start_y ; y >= min_start_y ; y--) {
	    int xgood ;
	    for(xgood = 0 ; xy_has_nonblack_pixel(image, xgood, y) == 0 && xgood < edge_width ; xgood++) ;
	    if(xgood == 0) continue ; // left edge pixel is filled

	    if(xgood == edge_width) {
		fprintf(stderr, "WARNING: could not find good start pixel on left edge gap at y=%d\n", y) ;
		continue ;
	    }

	    double rcoef[4], gcoef[4], bcoef[4] ;
	    int search_width = (xgood-0)*10 ;
	    if(search_width > 20) search_width=20 ;
	    get_mean_rgb(image,xgood,y,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,0,2, pData) ;

	    fprintf(stderr, "RTOS fix sky left edge from x,y %d,%d to %d,%d\n", 0, y, xgood-1, y) ;
	    for(x = 0 ; x < xgood ; x++) {
		if(pData->column_mask[x] == 1) continue ;
		float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
		float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
		float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

		tif_set4c(image,x,y, (uint16_t)rhat, (uint16_t)ghat, (uint16_t)bhat, MAX16 );
		pData->start_of_sky[x] = y ;
	    }

	}

	// RIGHT side
	min_start_y=IMAGE_HEIGHT ;
	max_start_y=0 ;
	int inside_x = IMAGE_WIDTH-1-edge_width ;

	for(x = inside_x ; x < IMAGE_WIDTH ; x++) {
	    if(min_start_y > pData->start_of_sky[x]) min_start_y = pData->start_of_sky[x] ;
	    if(max_start_y < pData->start_of_sky[x]) max_start_y = pData->start_of_sky[x] ;
	}

	for(int y = max_start_y ; y >= min_start_y ; y--) {
	    int xgood ;
	    for(xgood = IMAGE_WIDTH-1 ; xy_has_nonblack_pixel(image, xgood, y) == 0 && xgood > inside_x ; xgood--) ;
	    if(xgood == IMAGE_WIDTH-1) continue ; // right edge pixel is filled

	    if(xgood == inside_x) {
		fprintf(stderr, "WARNING: could not find good start pixel on right edge gap at y=%d\n", y) ;
		continue ;
	    }

	    double rcoef[4], gcoef[4], bcoef[4] ;
	    int search_width = (IMAGE_WIDTH-xgood)*10 ;
	    if(search_width > 20) search_width=20 ;
	    get_mean_rgb(image,xgood,y,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,0,2, pData) ;
/*  	    XXX  */

	    fprintf(stderr, "RTOS fix sky right edge from x,y %d,%d to %d,%d\n", xgood+1, y, IMAGE_WIDTH-1,y) ;

	    for(x = xgood+1 ; x < IMAGE_WIDTH ; x++) {
		if(pData->column_mask[x] == 1) continue ;
		float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
		float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
		float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

		//fprintf(stderr, "    x:%d rgb:%d,%d,%d\n", x, (int)rhat,(int)ghat,(int)bhat) ;

		tif_set4c(image,x,y, (uint16_t)rhat, (uint16_t)ghat, (uint16_t)bhat, MAX16 );
		pData->start_of_sky[x] = y ;
	    }

	}



    }

    if(phase == 8) {
#ifdef GMR_DEBUG
	gmr_debug = 1000 ;
#endif
	fprintf(stderr, "RTOS second scan for bad pixels in top 50%% of sky\n") ;

	double rcoef[4], gcoef[4], bcoef[4] ;
	// are there wrong pixels in the top 50 % of the sky?


	int search_width = 5 ;
	int n_repaired_gt=0 ; // grand total repaired ;

	for(int xc=search_width ; xc < IMAGE_WIDTH-5 ; xc += search_width) {
	    int y0 = 0 ;
	    int y1 = IMAGE_HEIGHT-1 ;
	    int x0 = xc-search_width ;
	    int x1 = xc+search_width ;
	    if(x1 > IMAGE_WIDTH-1) x1 = IMAGE_WIDTH-1 ;
	    int n_unmasked_columns=0 ;

	    for(int x = x0 ; x <= x1 ; x++) {
		if(pData->column_mask[x] == 1) continue ;
		if(pData->column_repair_mask[x] == 1) continue ;
		if(y0 < pData->start_of_sky[x]) y0 = pData->start_of_sky[x] ;
		if(y1 > pData->end_of_sky[x]) y1 = pData->end_of_sky[x] ;
		n_unmasked_columns++ ;
	    }

	    if(n_unmasked_columns == 0) continue ;

	    y1 = (y0 + y1)/2. ;
	    for(int yc=y1-search_width ; yc > y0 ; yc -= search_width) {
		int n_repaired = get_mean_rgb(image,xc,yc,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,1,2, pData) ;
		n_repaired_gt += n_repaired ;
/*  		if(n_repaired > 0)  */
/*  		    fprintf(stderr, "Bad pixels in top 50%% of sky: clip EOS=%d,repaired %d pixel colors near %d,%d\n", end_of_sky_is_known_flag, n_repaired, xc, yc) ;  */

#ifdef GMR_DEBUG
		gmr_debug-- ;
		if(gmr_debug == 0)
		    exit(1) ;
#endif
	    }
	}

	fprintf(stderr, "Bad pixels in top 50%% of sky: repaired %d pixel colors\n", n_repaired_gt) ;
    }

    // are there sky slopes so steep they will affect the result?
    if(phase == 9 && pData->repair_sky_slope_correct==1) {

        // run back and forth until there is no slope too steep
	int needs_more=1 ;

        while(needs_more) {
	    needs_more=0 ;
	    fprintf(stderr, "RTOS slope correct\n") ;
	    int search_width = 15 ;

	    double rcoef[4], gcoef[4], bcoef[4] ;

	    for(x = 0 ; x < IMAGE_WIDTH-4 ; x++) {
		int x2 = x+3 ;
		// cannot go across masked columns
		if(pData->column_mask[x] == 1) continue ;
		if(pData->column_mask[x+1] == 1) continue ;
		if(pData->column_mask[x+2] == 1) continue ;
		if(pData->column_mask[x+3] == 1) continue ;

		if(pData->start_of_sky[x] > pData->start_of_sky[x2]+1) {
		    y=pData->start_of_sky[x]-1 ;
		    int rval = get_mean_rgb(image,x,pData->start_of_sky[x]+search_width/2.,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,0,2, pData) ;

		    if(rval >= 0) {
			float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
			float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
			float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

			tif_set4c(image,x,y, (uint16_t)(rhat+0.5), (uint16_t)(ghat+0.5), (uint16_t)(bhat+0.5), MAX16 );

			pData->start_of_sky[x]-- ;
			needs_more=1 ;
		    }
		}
		if(pData->start_of_sky[x]+1 < pData->start_of_sky[x2]) {
		    y=pData->start_of_sky[x2]-1 ;
		    int rval = get_mean_rgb(image,x2,pData->start_of_sky[x2]+search_width/2.,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,0,2, pData) ;

		    if(rval >= 0) {
			float rhat = rcoef[0] + rcoef[1]*(float)x2 + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
			float ghat = gcoef[0] + gcoef[1]*(float)x2 + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
			float bhat = bcoef[0] + bcoef[1]*(float)x2 + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

			tif_set4c(image,x2,y, (uint16_t)(rhat+0.5), (uint16_t)(ghat+0.5), (uint16_t)(bhat+0.5), MAX16 );
			pData->start_of_sky[x2]-- ;
			needs_more=1 ;
		    }
		}
	    }
	}
    }

    if(pData->fill_top_of_sky == 1 && phase == 10) {
	fprintf(stderr, "RTOS fill gaps to min start of sky\n") ;

	for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	    if(pData->column_mask[x] == 1) continue ;

	    if(pData->start_of_sky[x] > pData->min_start_of_sky) {
		int search_width = 30 ; // need a robust estimate, so wide search width
		double rcoef[4], gcoef[4], bcoef[4] ;
		get_mean_rgb(image,x,pData->start_of_sky[x],rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,1,2, pData) ;

		for(y = pData->start_of_sky[x] ; y > pData->min_start_of_sky ; y--) {
		    float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
		    float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
		    float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

		    tif_set4c(image,x,y, (uint16_t)(rhat+0.5), (uint16_t)(ghat+0.5), (uint16_t)(bhat+0.5), MAX16 );
		}

		pData->start_of_sky[x] = pData->min_start_of_sky ;


	    }
	}
    }

    if(!pData->full_sky_replacement && phase == 11) {
	fprintf(stderr, "RTOS fill holes under opaque sections\n") ;

	int *opaque_x = (int *)calloc(IMAGE_HEIGHT, sizeof(int)) ;

	for(x = 0 ; x < IMAGE_WIDTH-1 ; x++) {
	    if(pData->column_mask[x] == 1) continue ;

	    // looking for start of sky changes > 1
	    int mid_x = -1, mid_y, y0, y1 ;

	    if(pData->start_of_sky[x+1] - pData->start_of_sky[x] > 1) {
		y0 = pData->start_of_sky[x] ;
		y1 = pData->start_of_sky[x+1] ;
		fprintf(stderr, "check alpha gap at x:%d, y %d to %d\n", x, y0, y1) ;
		mid_x = IMAGE_WIDTH ;
		mid_y = 0 ;

		for(int y=y0+1 ; y < y1 ; y++) {
		    int xs ;
		    for(xs=x ; xs > 0 ; xs--) {
			if(xy_has_nonblack_pixel(image, xs, y)) {
			    opaque_x[y] = xs ;
			    if(xs < mid_x) {
				mid_x = xs ;
				mid_y = y ;
			    }
			    break  ;
			}
		    }

		    if(xs < 0) {
			// hit the edge of the image, this is not a fixable problem
			mid_x=-1 ;
			break ;
		    }
		}

		if(mid_x == x) {
		    // no transparent pixels found beneath neigboring start of sky, nothing to be done
		    mid_x=-1 ;
		}
	    } else if(pData->start_of_sky[x] - pData->start_of_sky[x+1] > 1) {
		y0 = pData->start_of_sky[x+1] ;
		y1 = pData->start_of_sky[x] ;
		fprintf(stderr, "check alpha gap at x:%d, y %d to %d\n", x, y0, y1) ;
		mid_x = -1 ;
		mid_y = 0 ;

		for(int y=y0+1 ; y < y1 ; y++) {
		    int xs ;
		    for(xs=x+1 ; xs < IMAGE_WIDTH ; xs++) {
			if(xy_has_nonblack_pixel(image, xs, y)) {
			    opaque_x[y] = xs ;
			    if(xs > mid_x) {
				mid_x = xs ;
				mid_y = y ;
			    }
			    break  ;
			}
		    }

		    if(xs >= IMAGE_WIDTH) {
			// hit the edge of the image, this is not a fixable problem
			mid_x=-1 ;
			break ;
		    }
		}

		if(mid_x == x+1) {
		    // no transparent pixels found beneath neigboring start of sky, nothing to be done
		    mid_x=-1 ;
		}
	    } else {
		continue ;
	    }

	    if(mid_x > -1) 
		fprintf(stderr, "mid_x:%d, mid_y %d\n", mid_x, mid_y) ;

	    if(mid_x >= 0 && mid_x < IMAGE_WIDTH) {
		int x0, x1 ;
		int pts[3][2] ;
		// get three sets of coefs for prediction, at (min_x, min_y), (x,y0) and (x+1,y1)
		pts[0][0] = mid_x ;
		pts[0][1] = mid_y ;
		if(mid_x < x+1) {
		    x0=mid_x+1 ;
		    x1=x+1 ;
		    pts[1][0] = x ;
		    pts[1][1] = y0 ;
		    pts[2][0] = x+1 ;
		    pts[2][1] = y1 ;
		} else {
		    x0=x+1 ;
		    x1=mid_x ;
		    pts[1][0] = x ;
		    pts[1][1] = y1 ;
		    pts[2][0] = x+1 ;
		    pts[2][1] = y0 ;
		}

		double rcoef[3][4], gcoef[3][4], bcoef[3][4] ;

		for(int i=0 ; i < 3 ; i++) {
		    int search_width = 30 ; // need a robust estimate, so wide search width
/*  int get_mean_rgb(tdata_t *image, int xc,int yc,double rcoef[],double gcoef[], double bcoef[],int search_width, int clip_eos_flag, int repair_bad_points_flag, int recursion_level, SKYFILL_DATA_t *pData)  */
		    get_mean_rgb(image,pts[i][0],pts[i][1],rcoef[i],gcoef[i],bcoef[i],search_width,0,1,2, pData) ;
		}

#define V1 0
#define V2 1
#define V3 2
		// get denominator for barycentric weighting
		float denom = (pts[V2][1]-pts[V3][1])*(pts[V1][0]-pts[V3][0]) + (pts[V3][0]-pts[V2][0])*(pts[V1][1]-pts[V3][1]) ;

		// now fill in the missing values
		for(int xs=x0 ; xs < x1 ; xs++) {
		    for(y=y0+1 ; y < y1 ; y++) {
			if(xy_has_nonblack_pixel(image, xs, y) == 0) {
			    float rgb_hat[3][3] ;
			    for(int i=0 ; i < 3 ; i++) {
				rgb_hat[i][0] = rcoef[i][0] + rcoef[i][1]*(float)xs + rcoef[i][2]*(float)y + rcoef[i][3]*(float)y*(float)y ;
				rgb_hat[i][1] = gcoef[i][0] + gcoef[i][1]*(float)xs + gcoef[i][2]*(float)y + gcoef[i][3]*(float)y*(float)y ;
				rgb_hat[i][2] = bcoef[i][0] + bcoef[i][1]*(float)xs + bcoef[i][2]*(float)y + bcoef[i][3]*(float)y*(float)y ;
			    }
			    float num1 = (pts[V2][1]-pts[V3][1])*(xs-pts[V3][0]) + (pts[V3][0]-pts[V2][0])*(y-pts[V3][1]) ;
			    float w1 = num1/denom ;
			    float num2 = (pts[V3][1]-pts[V1][1])*(xs-pts[V3][0]) + (pts[V1][0]-pts[V3][0])*(y-pts[V3][1]) ;
			    float w2 = num2/denom ;
			    float w3 = 1. - w1 - w2 ;

			    float rhat = w1*rgb_hat[0][0] + w2*rgb_hat[1][0] + w3*rgb_hat[2][0] ;
			    float ghat = w1*rgb_hat[0][1] + w2*rgb_hat[1][1] + w3*rgb_hat[2][1] ;
			    float bhat = w1*rgb_hat[0][2] + w2*rgb_hat[1][2] + w3*rgb_hat[2][2] ;

			    tif_set4c(image,xs,y, (uint16_t)(rhat+0.5), (uint16_t)(ghat+0.5), (uint16_t)(bhat+0.5), MAX16 );
			}
		    }
		}

	    }


	}

	free(opaque_x) ;
    }

/*      fprintf(stderr, "Finished RTOS::: sky hsv (%3.0f,%4.2f,%4.2f)\n", pData->hue_sky, pData->sat_sky, pData->val_sky) ;  */

}
