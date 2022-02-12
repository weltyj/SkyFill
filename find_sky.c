#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>

#include "skyfill_tif.h"
#include "colorspace_conversions.h"
#include "pixel_tests.h"
#include "find_sky.h"
#include "mstat.h"

// some shared variables in this module
static float *sky_hue, *sky_val, *sky_sat  ;
static uint8_t *is_clear_sky ;

void simple_find_start_of_sky(int16_t *start_of_sky,uint16_t w,uint16_t h,tdata_t *image, SKYFILL_DATA_t *pData)
{
    int32_t x, y ;

    for(x = 0 ; x < w ; x++) {
	pData->start_of_sky[x] = 0 ;

	if(pData->column_mask[x] == 1) {
	    pData->end_of_sky[x] = 0 ;
	    continue ;
	}

	pData->end_of_sky[x] = IMAGE_HEIGHT-1 ;

	for(y = 0 ; y < h ; y++) {
	    if(xy_has_nonblack_pixel(image, x, y) == 1) {
		pData->start_of_sky[x] = y ;
		break ;
	    }
	}
    }

    fprintf(stderr, "Finished simple find start of sky\n") ;
}

void find_start_of_sky(int16_t *start_of_sky,uint16_t w,uint16_t h,tdata_t *image,int fix_SOS_edges, SKYFILL_DATA_t *pData)
{
    simple_find_start_of_sky(start_of_sky,w,h,image,pData) ;
    return ;

    // the old code tried to correct problems with the image here.
    // the new method uses the repair_top_of_sky() function to find problems and adjust the start_of_sky in a series of phases

#ifdef OLD_CODE
    /* look for vertical gaps in sky (alpha=0), which can happen after lens correction */

    for(x = 0 ; x < w ; x++) {
	int32_t gap_start=-1, gap_end=-1 ;

	for(y = start_of_sky[x] ; y < h ; y++) {
	    if(xy_has_nonblack_pixel(image, x, y) == 0) {
		gap_start = y ;
		break ;
	    }
	}

	if(gap_start > -1) {
	    for(y = gap_start ; y < h ; y++) {
		/* if alpha channel is zero then hugin has placed image data here */
		if(xy_has_nonblack_pixel(image, x, y) == 1) {
		    gap_end = y-1 ;
		    break ;
		}
	    }
	}

	if(gap_end> -1 && gap_start > -1) {
	    /* interpolate gap values from start to end */
	    float start_r,start_g,start_b ;
	    float end_r,end_g,end_b ;
	    tif_get3c(image,x,gap_start-1,start_r,start_g,start_b) ;
	    tif_get3c(image,x,gap_end+1,start_r,start_g,start_b) ;

	    for(y = gap_start ; y < gap_end+1 ; y++) {
		float pStart = 1.f - ((float)y-((float)gap_start-1.f))/((float)gap_end-(float)gap_start+3.f) ;
		float pEnd = 1.f - pStart ;

		tif_set4c(image,x,y, (uint16_t)(pStart*start_r + pEnd*end_r), (uint16_t)(pStart*start_g + pEnd*end_g), (uint16_t)(pStart*start_b + pEnd*end_b), MAX16) ; 
	    }
	}

    }

    for(x = 0 ; x < w ; x++) {
	if(pData->column_mask[x] == 1) start_of_sky[x]=0 ;
    }


    if(fix_SOS_edges == 1) {
	fprintf(stderr, "Fixing start of sky edges\n") ;

	// eliminate long vertical edges in the sky end, which can create
	// visible edges in the feathering result
	for(x = 0 ; x < w-1 ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    if(pData->column_mask[x+1] == 1) continue ;

	    if( (start_of_sky[x] - start_of_sky[x+1]) > 2) {
		start_of_sky[x+1] = start_of_sky[x]-2 ;
	    }
	}

	for(x = w-1 ; x > 0 ; x--) {
	    if(pData->column_mask[x] == 1) continue ;
	    if(pData->column_mask[x-1] == 1) continue ;

	    if( (start_of_sky[x] - start_of_sky[x-1]) > 2) {
		start_of_sky[x-1] = start_of_sky[x]-2 ;
	    }
	}

    }
#endif
}


int get_sky_mean_var(tdata_t *image, int x, int y0, int n_want, float mean[], float var[])
{
    float sum[3] = {0., 0., 0.} ;  // for r,g,b
    float sumsq[3] = {0., 0., 0.} ;  // for r,g,b
    int y ;
    int n=0 ;

    for(y = y0 ; y < y0+n_want ; y++) {

	if(y < 0) continue ;
	if(y >= IMAGE_HEIGHT) break ;
	if(x < 0 || x >= IMAGE_WIDTH) break ;

	float r,g,b ;
	tif_get3c(image,x,y,r,g,b) ;
	r /= MAX16f ;
	g /= MAX16f ;
	b /= MAX16f ;

	sum[0] += r ;
	sum[1] += g ;
	sum[2] += b ;
	sumsq[0] += r*r ;
	sumsq[1] += g*g ;
	sumsq[2] += b*b ;
	n++ ;
    }

    if(n > 1) {
	mean[0] = sum[0]/(float)n ;
	mean[1] = sum[1]/(float)n ;
	mean[2] = sum[2]/(float)n ;

	var[0] = sumsq[0] - sum[0]*mean[0] ; // equiv to sumsq - sum*sum/n ;
	var[1] = sumsq[1] - sum[1]*mean[1] ; // equiv to sumsq - sum*sum/n ;
	var[2] = sumsq[2] - sum[2]*mean[2] ; // equiv to sumsq - sum*sum/n ;
    } else {
	mean[0] = 0. ;
	mean[1] = 0. ;
	mean[2] = 0. ;

	var[0] = 10. ;
	var[1] = 10. ;
	var[2] = 10. ;
    }

    return n ;
}


int is_end_of_sky(uint16_t y, int n_samples, float *hdiff, float *sdiff, float *vdiff, SKYFILL_DATA_t *pData)
{
    float m_hue0=0. ;
    float m_val0=0. ;
    float m_sat0=0. ;
    float m_hue1=0. ;
    float m_val1=0. ;
    float m_sat1=0. ;
    int iy ;

    if(1.f-(float)y/((float)IMAGE_HEIGHT-1.) > pData->min_sky_end_p) return 0 ;  // is it greater than the min_sky_end_p requested?
    float sum_wgts=0. ;

    for(iy=y-n_samples+1 ; iy <= y ; iy++) {
	if(iy < 0 || iy > IMAGE_HEIGHT-1) {
	    fprintf(stderr, "iy=%d, y=%d n_samples=%d\n", (int)iy, (int)y, n_samples) ;
	    exit(1) ;
	}
	m_hue0 += sky_hue[iy]*is_clear_sky[iy] ;
	m_sat0 += sky_sat[iy]*is_clear_sky[iy] ;
	m_val0 += sky_val[iy]*is_clear_sky[iy] ;
	sum_wgts += is_clear_sky[iy] ;
    }

    m_hue0 /= sum_wgts ;
    m_sat0 /= sum_wgts ;
    m_val0 /= sum_wgts ;

    sum_wgts=0. ;
    for(iy=y+1 ; iy < y+1+n_samples/2. ; iy++) {
	if(sky_hue[iy] < -0.01) {
	    fprintf(stderr, "FATAL: sky_hue not computed in is_end_of_sky at y=%d\n", iy) ;
	}
	m_hue1 += sky_hue[iy]*is_clear_sky[iy] ;
	m_sat1 += sky_sat[iy]*is_clear_sky[iy] ;
	m_val1 += sky_val[iy]*is_clear_sky[iy] ;
	sum_wgts += is_clear_sky[iy] ;
    }

    m_hue1 /= sum_wgts ;
    m_sat1 /= sum_wgts ;
    m_val1 /= sum_wgts ;

    float huediff = fabsf(m_hue0 - m_hue1) ;
    float satdiff = fabsf(m_sat0 - m_sat1) ;
    float valdiff = fabsf(m_val0 - m_val1) ;

    huediff *= (m_sat0+m_sat1)/2. ;   // adjust hue difference by amount of saturation
    *hdiff = huediff ;
    *sdiff = satdiff ;
    *vdiff = valdiff ;


    if(huediff > pData->sky_hue_tolerance) {
	//fprintf(stderr, "ieos %d %3.0f %0.3f %0.3f\n", y, huediff, satdiff, valdiff) ;
	return 1 ;
    }

    if(valdiff > pData->sky_val_tolerance) {
	//fprintf(stderr, "ieos %d %3.0f %0.3f %0.3f\n", y, huediff, satdiff, valdiff) ;
	return 1 ;
    }

    if(satdiff > pData->sky_sat_tolerance) {
	//fprintf(stderr, "ieos %d %3.0f %0.3f %0.3f\n", y, huediff, satdiff, valdiff) ;
	return 1 ;
    }

    return 0 ;
}

float sobel(tdata_t *image, int xc, int yc)
{
    uint16_t r[3][3] ;
    uint16_t g[3][3] ;
    uint16_t b[3][3] ;

    for(int x = xc-1 ; x <= xc+1 ; x++) {
	int i = x - (xc-1) ;
	for(int y = yc-1 ; y <= yc+1 ; y++) {
	    int j = y - (yc-1) ;
	    tif_get3c(image,x,y,r[i][j],g[i][j],b[i][j]) ;
	}
    }

    float Gx_r = 
	-  r[0][0] +   r[0][2]
	-2*r[1][0] + 2*r[0][2]
	-  r[2][0] +   r[2][2] ;
    float Gx_g = 
	-  g[0][0] +   g[0][2]
	-2*g[1][0] + 2*g[0][2]
	-  g[2][0] +   g[2][2] ;
    float Gx_b = 
	-  b[0][0] +   b[0][2]
	-2*b[1][0] + 2*b[0][2]
	-  b[2][0] +   b[2][2] ;

    float Gy_r = 
	    r[0][0] + 2*r[0][1] + r[0][2]
	   -r[2][0] - 2*r[2][1] - r[2][2] ;
    float Gy_g = 
	    g[0][0] + 2*g[0][1] + g[0][2]
	   -g[2][0] - 2*g[2][1] - g[2][2] ;
    float Gy_b = 
	    b[0][0] + 2*b[0][1] + b[0][2]
	   -b[2][0] - 2*b[2][1] - b[2][2] ;

    return fabsf(Gx_r) + fabsf(Gx_g) + fabsf(Gx_b) + (fabsf(Gy_r) + fabsf(Gy_g) + fabsf(Gy_b))/2. ;
}


void print_sobel(tdata_t *image, int x, SKYFILL_DATA_t *pData)
{
    int16_t y ;

    for(y = pData->end_of_sky[x]-5 ; y < pData->end_of_sky[x]+5 ; y++) {
	if( y == pData->end_of_sky[x]) {
	    fprintf(stderr, "SOBEL col %d,%d * %f\n", x, y, sobel(image,x,y)) ;
	} else {
	    fprintf(stderr, "SOBEL col %d,%d   %f\n", x, y, sobel(image,x,y)) ;
	}
    }
    exit(1) ;
}

int get_hires_eos(tdata_t *image, int x, SKYFILL_DATA_t *pData)
{
    int y ;

    // trying out new end of sky detection
    float obuf[20] ; 
    //current eos-10 will be obuf index 0 ;
    int offset = pData->end_of_sky[x]-10 ;
    int y0 = pData->end_of_sky[x]-10 ;
    int y1 = y0+19 ;

    int sample_length=3 ;
    for(y = y0 ; y <= y1 ; y++) {
	float r0,g0,b0,a0 ;
	float r1,g1,b1,a1 ;
	r0=g0=b0=a0=0. ;
	r1=g1=b1=a1=0. ;
	float r,g,b ;

	if(y-sample_length > 0 && y+sample_length <= IMAGE_HEIGHT-1) {
	    for(int z = y-sample_length ; z < y ; z++) {
		tif_get3c(image,x,z,r,g,b) ;
		r0 += r ;
		g0 += g ;
		b0 += b ;
	    }

	    for(int z = y+1 ; z <= y+sample_length ; z++) {
		tif_get3c(image,x,z,r,g,b) ;
		r1 += r ;
		g1 += g ;
		b1 += b ;
	    }

	    r0 /= (float)sample_length ;
	    r1 /= (float)sample_length ;
	    g0 /= (float)sample_length ;
	    g1 /= (float)sample_length ;
	    b0 /= (float)sample_length ;
	    b1 /= (float)sample_length ;

	    float dr = fabsf(r1-r0) ;
	    float dg = fabsf(g1-g0) ;
	    float db = fabsf(b1-b0) ;

	    obuf[y-offset] = dr+dg+db ;
	} else {
	    obuf[y-offset] = 0. ;
	}

    }

    int eos = 0 ;
    float val_tol = pData->sky_val_tolerance/.015* 4000. ;
    for(y = y0 ; y <= y1  ; y++) {
	if(obuf[y-offset] > val_tol) {
	    eos = y+sample_length-1 ;
	    return eos ;
	}
    }

    return pData->end_of_sky[x] ;

}

int new_get_hires_eos(tdata_t *image, int x, SKYFILL_DATA_t *pData)
{
    int y ;

    int sample_length=20 ;
    int y0 = pData->end_of_sky[x]-sample_length ;
    int y1 = y0+2*sample_length-1 ;
    float best_metric=0. ;
    int best_y=pData->end_of_sky[x] ;

    for(y = y0 ; y <= y1 ; y++) {
	float mean_above[3], var_above[3] ;
	float mean_below[3], var_below[3] ;

	int n_above = get_sky_mean_var(image, x, y-sample_length-1,  sample_length, mean_above, var_above) ;
	int n_below = get_sky_mean_var(image, x, y+1,  sample_length, mean_below, var_below) ;
	if(n_above < 2 || n_below < 2) continue ;

	float sum_diff=0. ;
	float sum_var_above=0. ;

	for(int c=0 ; c < 3 ; c++) {
	    sum_diff += fabsf(mean_above[c]-mean_below[c]) ;
	    sum_var_above += var_above[c] ;
	}

	float metric = 1./(sum_var_above+1.)+sum_diff ;
	if(x == 1417) {
	    printf("NGHE: y:%d sum_var_above:%f sum_diff:%f metric:%f\n", y, sum_var_above, sum_diff, metric) ;
	}
	if(metric > best_metric) {
	    best_metric=metric ;
	    best_y=y ;
	}
    }

    return best_y ;

}

//int16_t find_eos_at_column(int16_t x, int16_t starting_y, int16_t *end_of_sky,int16_t *start_of_sky,int w,int h,tdata_t *image, SKYFILL_DATA_t *pData)
int16_t find_eos_at_column(int16_t x, int16_t starting_y, tdata_t *image, SKYFILL_DATA_t *pData)
{
    int16_t y ;

    for(y=0; y < IMAGE_HEIGHT ; y++) {
	is_clear_sky[y]=1 ;
    }

    for(int y = 0 ; y < IMAGE_HEIGHT ; y++) {
	sky_hue[y] = -1. ;
    }

    int first_y ;

    if(starting_y > 0) {

	int first_valid_y = max_test_mask(x, 0) ;

	if(first_valid_y >= IMAGE_HEIGHT) {
	    pData->end_of_sky[x] = 0 ;
	    return 0 ;
	}

	first_y=starting_y ;
    } else {

	int first_valid_y = max_test_mask(x, 0) ;

	if(first_valid_y >= IMAGE_HEIGHT) {
	    pData->end_of_sky[x] = 0 ;
	    return 0 ;
	}

	// start looking at maximum start of sky in local area
	first_y = first_valid_y+1 ;
	int x0 = x-20 ;
	if(x0 < 0) x0 = 0 ;
	int x1 = x+20 ;
	if(x1 > IMAGE_WIDTH-1) x1 = IMAGE_WIDTH-1 ;

	for(int dx=x0 ; dx <= x1 ; dx++) {
	    if(pData->column_mask[dx] == 1) continue ;

	    if(first_y < pData->raw_start_of_sky[dx])
		first_y = pData->raw_start_of_sky[dx]  ;
	}

	// start looking 10 pixels below max raw start of sky in local area
	first_y += 10 ;
    }

    float hdiff, sdiff, vdiff ;


#ifdef EOS_STEP1
    //first_y = end_of_sky[x]-6 ;  // temporarily see how this looks 
    float top_mean[3], top_var[3] ;
    float bot_mean[3], bot_var[3] ;

    //fprintf(stderr, "\nNEW EOS x:%5d eos:%5d\n", x, end_of_sky[x]) ;

    for(y = first_y ; y < IMAGE_HEIGHT-10 ; y++) {
	// compare the mean value of the previous 10 pixels,
	// to the  mean value of the next 5 pixels

	// FIXME -- can't use raw_start_of_sky when start_y was > 0 ;
	int n_prev = y - pData->raw_start_of_sky[x] ;
	if(n_prev > 20) n_prev=20 ;

	int n_above = get_sky_mean_var(image, x, y-n_prev, n_prev, top_mean, top_var) ;

	int n_below = get_sky_mean_var(image, x, y, 5, bot_mean, bot_var) ;
	float rr = fabs(top_mean[0]/bot_mean[0]-1.) ;
	float rg = fabs(top_mean[1]/bot_mean[1]-1.) ;
	float rb = fabs(top_mean[2]/bot_mean[2]-1.) ;
#ifdef DEBUG_COL
	if(x == DEBUG_COL) {
	    fprintf(stderr, "x:%d, y:%5d n_prev:%2d rr,rg,rb = %5.2f %5.2f %5.2f\n", x, y, n_prev, rr, rg, rb) ;
	}
#endif

	// do we have the same color above and below the test area?
	if(rr > .1 || rg > .1 || rb > .1) {
	    // no, colors are too different.  But mark sky as clear
	    // for this and all pixels below so they can be used
	    // in the next step for calculating sky hsv values
	    for(; y < IMAGE_HEIGHT-10 ; y++) {
		is_clear_sky[y]=1 ;
	    }
	    // break out of the loop -- last_y
	    break ;
	}

	is_clear_sky[y]=1 ;

	// if transparent or black, must be at bottom of image
	uint16_t r,g,b,a ;
	tif_get4c(image,x,y+10,r,g,b,a) ;

	if(a < HALF16) break ;
	if( r < BLK16 && g < BLK16 && b < BLK16 ) break ;



	// try to skip small highlights (jet trails) in sky
	if(0) {
	    int y1 = y+5 ;
	    for(int y0 = y ; y0 < y1 ; y0++) {
		tif_get3c(image,x,y0,r,g,b) ;

		float rtop=(float)r/top_mean[0] ;
		float gtop=(float)g/top_mean[1] ;
		float btop=(float)b/top_mean[2] ;
		float rbot=(float)r/bot_mean[0] ;
		float gbot=(float)g/bot_mean[1] ;
		float bbot=(float)b/bot_mean[2] ;

		if( rtop > 1. && gtop > 1. && btop > 1. && rbot > 1. && gbot > 1. && bbot > 1.) {
		    // pixel at y0 is brighter than pixels above and below it
		    // mark it as not clear, and set y to y0 to prevent algorithm
		    // from marking is_clear_sky[y] to 1 in the next iteration of the outer loop
		    is_clear_sky[y0] = 0 ;
		    y = y0 ;
		} else {
		    is_clear_sky[y0] = 1 ;
		}
	    }
	}

#ifdef DEBUG_COL
	if(x == DEBUG_COL) printf("x:%d y:%d max_test_mask:%d\n", x, y, max_test_mask(x,y)) ;
#endif
	y = max_test_mask(x, y)+1 ;

    }

#ifdef DEBUG_COL
    if(x == DEBUG_COL)
	printf("\n") ;
#endif

    // STEP1
#endif


    // load up arrays

    for(y = first_y ; y < first_y+20 ; y++) {
	uint16_t rgb[3] ;
	tif_get3cv(image,x,y,rgb) ;

	float hsv[3] ;
	rgb2hsv16(rgb,hsv) ;
	float l = (.299*(float)rgb[0]+.587*(float)rgb[1]+.114*(float)rgb[2])/MAX16f ;

	// any areas to skip because of a test mask get the hsl of the previous valid y
	int next_valid_y = max_test_mask(x, y)+1 ;
	for(int y0 = y ; y0 < next_valid_y ; y0++) {
	    if(y0 == IMAGE_HEIGHT) break ;
	    sky_hue[y0] = hsv[0] ;
	    sky_sat[y0] = hsv[1] ;
	    sky_val[y0] = l ;

#ifdef DEBUG_COL
	    if(x == DEBUG_COL) {
		fprintf(stderr, "x:%d, y:%5d PRELOAD load hsl of  = %5.2f %5.2f %5.2f\n", x, y0, sky_hue[y0],sky_sat[y0],sky_val[y0]) ;
	    }
#endif
	}

	y = next_valid_y-1 ;
    }

    for(y = first_y+4 ; y < IMAGE_HEIGHT-10 ; y++) {
	pData->end_of_sky[x] = y ;
	int end_flag = is_end_of_sky(y,5,&hdiff,&sdiff,&vdiff,pData) ;

	if(end_flag) {
	    int y0 = y-4 ;

	    if(y0 < pData->raw_start_of_sky[x])
		y0 = pData->raw_start_of_sky[x] ;

	    for(int testy = y0 ; testy <= y ; testy++) {
		if(is_in_test_mask(x,testy)) {
		    // end of sky cannot happen in masked area
		    end_flag=0 ;
		}
	    }
	}

#ifdef DEBUG_COL
	if(x == DEBUG_COL) {
	    fprintf(stderr, "x:%d, y:%5d end_flag:%2d hdiff,sdiff,vdiff = %5.2f %5.2f %5.2f\n", x, y, end_flag, hdiff,sdiff,vdiff) ;
	}
#endif

	if(end_flag) {
	    break ;
	}

	// need to get the hsl for y+10 ;
	if(is_in_test_mask(x,y+10)) {
	    sky_hue[y+10] = sky_hue[y+9] ;
	    sky_sat[y+10] = sky_sat[y+9] ;
	    sky_val[y+10] = sky_val[y+9] ;
#ifdef DEBUG_COL
	    if(x == DEBUG_COL) {
		fprintf(stderr, "x:%d, y:%5d MASK LOAD hsl of  = %5.2f %5.2f %5.2f\n", x, y+10, sky_hue[y+10],sky_sat[y+10],sky_val[y+10]) ;
	    }
#endif
	} else {
	    uint16_t rgb[3] ;
	    tif_get3cv(image,x,y+10,rgb) ;

	    float hsv[3] ;
	    rgb2hsv16(rgb,hsv) ;
	    float l = (.299*(float)rgb[0]+.587*(float)rgb[1]+.114*(float)rgb[2])/MAX16f ;
	    sky_hue[y+10] = hsv[0] ;
	    sky_sat[y+10] = hsv[1] ;
	    sky_val[y+10] = l ;
#ifdef DEBUG_COL
	    if(x == DEBUG_COL) {
		fprintf(stderr, "x:%d, y:%5d INC LOAD hsl of  = %5.2f %5.2f %5.2f\n", x, y+10, sky_hue[y+10],sky_sat[y+10],sky_val[y+10]) ;
	    }
#endif
	}

    }

    // fine tune the detection -- increase y as long as vdiff becomes larger
    float vdiff_prev=vdiff ;

    for(y = pData->end_of_sky[x]+1 ; y < IMAGE_HEIGHT-10 ; y++) {

	if(is_in_test_mask(x,y+5)) {
	    sky_hue[y+5] = sky_hue[y+4] ;
	    sky_sat[y+5] = sky_sat[y+4] ;
	    sky_val[y+5] = sky_val[y+4] ;
	} else {
	    uint16_t rgb[3] ;
	    tif_get3cv(image,x,y+5,rgb) ;

	    float hsv[3] ;
	    rgb2hsv16(rgb,hsv) ;

	    float l = (.299*(float)rgb[0]+.587*(float)rgb[1]+.114*(float)rgb[2])/MAX16f ;
	    sky_hue[y+5] = hsv[0] ;
	    sky_sat[y+5] = hsv[1] ;
	    sky_val[y+5] = l ;
	}

	is_end_of_sky(y,5,&hdiff,&sdiff,&vdiff,pData) ;

#ifdef DEBUG_COL
	if(x == DEBUG_COL) {
	    fprintf(stderr, "Fine tune x:%d, y:%5d vdiff, vdiff_prev = %5.2f %5.2f\n", x, y, vdiff,vdiff_prev) ;
	}
#endif

	if(vdiff <= vdiff_prev)
	    break ;

	vdiff_prev = vdiff ;

	pData->end_of_sky[x]++ ;
    }

    // set up for a more refined estimate of end of sky
    pData->end_of_sky[x]-=2 ;
    pData->end_of_sky[x] = get_hires_eos(image,x,pData) ;

    return pData->end_of_sky[x] ;
}

#define DEFAULT_SKY_SLIVER_REPAIR 0x01
#define TEST_SKY_SLIVER_REPAIR 0x02
#define ED_SKY_SLIVER_REPAIR 0x04

void fix_sky_slivers(tdata_t *image, SKYFILL_DATA_t *pData, int method, int debug_col)
{

    // check for 1 pixel in sky that failed threshold test
    if(method == TEST_SKY_SLIVER_REPAIR) {
	fprintf(stderr, "Fixing sky slivers, testing method\n") ;
	for(int x = 1 ; x < IMAGE_WIDTH-1 ; x++) {
	    if(pData->end_of_sky[x] == 0) continue ;
	    int dleft = pData->end_of_sky[x-1] - pData->end_of_sky[x] ;
	    int dright = pData->end_of_sky[x+1] - pData->end_of_sky[x] ;
	    int done=0 ;
	    int current_eos = pData->end_of_sky[x] ;
	    int iter=0 ;
	    while(!done && (dleft > 10 || dright > 10)) {
		fprintf(stderr, "EOS sliver: x:%d current:%d dleft:%d dright:%d\n", x, current_eos, dleft, dright) ;
		int new_eos = find_eos_at_column(x, current_eos+8, image, pData) ;
		fprintf(stderr, "       new:%d\n", new_eos) ;
		if(new_eos-current_eos < 10) done=1 ;
		dleft = pData->end_of_sky[x-1] - pData->end_of_sky[x] ;
		dright = pData->end_of_sky[x+1] - pData->end_of_sky[x] ;
		iter++ ;
		if(iter > 4) break ;
	    }
	}
    } else if(method == ED_SKY_SLIVER_REPAIR) {
	fprintf(stderr, "Fixing sky slivers, ED t-test method\n") ;
	float *ed = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;

	for(int x = 1 ; x < IMAGE_WIDTH-1 ; x++) {
	    if(pData->end_of_sky[x] == 0) continue ;
	    int current_eos = pData->end_of_sky[x] ;
	    int dleft = pData->end_of_sky[x-1] - current_eos ;
	    int dright = pData->end_of_sky[x+1] - current_eos ;
	    int compare_x = -1 ;

	    if(dleft > 10 || dright > 10) {
		compare_x = dleft > dright ? x-1 : x+1 ;
	    }

	    int compare_eos = pData->end_of_sky[compare_x] ;
	    int stepsize=10 ;

	    if(compare_x > -1) {
		// simple case check, just do paired t-tests in blocks of 10 pixels until it fails
		// mean == 0 test
		while(current_eos < compare_eos) {
		    float sum_d[3] = {0.,0.,0.} ;
		    float sum_d2[3] = {0.,0.,0.} ;
		    float n=0. ;
		    for(int y = current_eos ; y < current_eos+stepsize ; y++) {
			float rgb[3] ;
			float rgb_neighbor[3] ;
			if(y >= compare_eos) break ;
			tif_get3c(image, x, y, rgb[0], rgb[1], rgb[2]) ;
			tif_get3c(image, compare_x, y, rgb_neighbor[0], rgb_neighbor[1], rgb_neighbor[2]) ;

/*  			if(x == debug_col)  */
/*  			    fprintf(stderr, "EDSLIV x:%d y:%d {rgb}", x, y) ;  */

			for(int c=0 ; c < 3 ; c++) {
			    float d = (rgb[c] - rgb_neighbor[c])/MAX16f ;
			    sum_d[c] += d ;
			    sum_d2[c] += d*d ;
/*  			    if(x == debug_col)  */
/*  				fprintf(stderr, " %f,%f,%f,%f", rgb[c]/MAX16f,rgb_neighbor[c]/MAX16f,d,sum_d[c]) ;  */
			}

/*  			if(x == debug_col)  */
/*  			    fprintf(stderr, "\n") ;  */
			n++ ;
		    }

		    int failed=0 ;

		    if(x == debug_col)
			fprintf(stderr, "EDSLIV x:%d eos:%d neighbor eos:%d d,t{rgb}", x, current_eos, compare_eos) ;

		    for(int c=0 ; c < 3 ; c++) {
			float sd = sqrt( (n*sum_d2[c]-sum_d[c]*sum_d[c])/(n-1.)) ;
			float t = (sum_d[c]/n) / sd ;

			if(x == debug_col)
			    fprintf(stderr, " %f,%f", sum_d[c],t) ;

			if(fabsf(t) > 1.833 && sum_d[c]/n > .02) {// two tailed, 90% probability, df=(10-1), and the mean diff has to be practically meaningful
			    if(x == debug_col)
				fprintf(stderr, "EDSLIV t-test fail, x:%d eos:%d c:%d, mean_d:%6.3f,t:%5.2f", x, current_eos, c, sum_d[c]/n, t) ;
			    failed=1 ;
			}
		    }

		    if(x == debug_col)
			fprintf(stderr, "\n") ;

		    if(failed) {
			pData->end_of_sky[x] = current_eos-stepsize ;
			break ;
		    } else {
			pData->end_of_sky[x] = current_eos ;
			current_eos += stepsize ;
		    }
		}

		if(pData->end_of_sky[x] > compare_eos)
		    pData->end_of_sky[x] = compare_eos ;
	    }

	} // for(x
    } else if(method == DEFAULT_SKY_SLIVER_REPAIR) {
	//look for and fix slivers
	int fix_slivers=1 ;
	while(fix_slivers >0) {
	    fprintf(stderr, "Fixing sky slivers\n") ;

	    for(int x = 1 ; x < IMAGE_WIDTH-3 ; x++) {
		if(pData->column_mask[x] == 1) continue ;
		if(pData->column_mask[x-1] == 1) continue ;

		if( (pData->end_of_sky[x-1] - pData->end_of_sky[x]) > 10 && (pData->end_of_sky[x+1] - pData->end_of_sky[x]) > 10) {
		    if(pData->column_mask[x+1] == 1) continue ;
		    //fprintf(stderr,"fix sliver at x=%d\n", x) ;
		    pData->end_of_sky[x] = (int)((pData->end_of_sky[x-1]+pData->end_of_sky[x+1])/2+0.5) ;

		} else if( (pData->end_of_sky[x-1] - pData->end_of_sky[x]) > 10 && (pData->end_of_sky[x+2] - pData->end_of_sky[x]) > 10) {
		    if(pData->column_mask[x+2] == 1) continue ;
		    //fprintf(stderr,"fix 2-wide sliver at x=%d\n", x) ;
		    pData->end_of_sky[x] = (int)((pData->end_of_sky[x-1]+pData->end_of_sky[x+2])/2+0.5) ;
		    pData->end_of_sky[x+1] = (int)((pData->end_of_sky[x-1]+pData->end_of_sky[x+2])/2+0.5) ;

		} else if( (pData->end_of_sky[x-1] - pData->end_of_sky[x]) > 10 && (pData->end_of_sky[x+3] - pData->end_of_sky[x]) > 10) {
		    if(pData->column_mask[x+3] == 1) continue ;
		    //fprintf(stderr,"fix 3-wide sliver at x=%d\n", x) ;
		    pData->end_of_sky[x] = (int)((pData->end_of_sky[x-1]+pData->end_of_sky[x+3])/2+0.5) ;
		    pData->end_of_sky[x+1] = (int)((pData->end_of_sky[x-1]+pData->end_of_sky[x+3])/2+0.5) ;
		    pData->end_of_sky[x+2] = (int)((pData->end_of_sky[x-1]+pData->end_of_sky[x+3])/2+0.5) ;
		}
	    }
	    fix_slivers-- ;
	}
    }
}



int find_end_of_sky_method_1(int16_t *end_of_sky,int16_t *start_of_sky,int w,int h,tdata_t *image,int fix_edges,int fix_slivers_flag, SKYFILL_DATA_t *pData)
{

    int16_t x,y ;

    // is_clear_sky[y] is used as a weight for calculating sky colors in 
    // is_end_of_sky(), the weight is 0. or 1.
    is_clear_sky = (uint8_t *)calloc(IMAGE_HEIGHT, sizeof(uint8_t)) ;

    sky_hue = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_val = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_sat = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;

// if defined, will output some debug info to consol
/*  #define DEBUG_COL 913  */

    // look for end of sky now

    for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	find_eos_at_column(x, -1, image, pData) ;
    }

/*      print_sobel(image,w/2) ;  */

    free(is_clear_sky) ;
    free(sky_hue) ;
    free(sky_sat) ;
    free(sky_val) ;

    return 0 ;
}

int find_end_of_sky_method_6(int16_t *end_of_sky,int16_t *start_of_sky,int w,int h,tdata_t *image,int fix_edges,int fix_slivers_flag, SKYFILL_DATA_t *pData, int debug_col_want, float sky_ratio_threshold)
{
    fprintf(stderr, "Find end of sky method 6\n") ;
    int *sos=(int *)calloc(IMAGE_WIDTH, sizeof(int)) ;
    int *next_row_check=(int *)calloc(IMAGE_WIDTH, sizeof(int)) ;
    int *eos_status=(int *)calloc(IMAGE_WIDTH, sizeof(int)) ;
    int *eos=(int *)calloc(IMAGE_WIDTH, sizeof(int)) ;

    int window_size=4 ; // fixing for now, ultimately creates a 9x9 square window for the statistical regression and analysis
    int full_window_size=window_size*2+1 ;

#define EOS_UNTESTED 1
#define EOS_MOVED 2
#define EOS_FINAL 3
#define BAD_FIT_PIXEL 0x10

    for(int x = 0 ; x < IMAGE_WIDTH ; x++) {
	eos_status[x] = EOS_UNTESTED ;
	eos[x] = IMAGE_HEIGHT-1 ;

	int W=4 ;
	int edge_W=window_size ;

	if(x < full_window_size || x > IMAGE_WIDTH-1-(full_window_size) ) {
	    // where x is 0 or IMAGE_WIDTH-1, no change
	    float p=1. ;

	    if(x < full_window_size) 
		p = (float)x/(float)full_window_size ;
	    else
		p = (float)(IMAGE_WIDTH-1-x)/(float)full_window_size ;

	    int w = (int)((float)W*p+0.5) ;

	    sos[x] = pData->raw_start_of_sky[x]+w ;
	} else {
	    int x0 = x-edge_W ;
	    int x1 = x+edge_W ;
	    if(x0 < 0) x0=0 ;
	    if(x1 > IMAGE_WIDTH-1) x1=IMAGE_WIDTH-1 ;
	    sos[x] = pData->raw_start_of_sky[x1] ;

	    // first check to make sure "edges" near other lower start of sky are moved down to 
	    // match the lowest neighbor
	    for(int x_test = x0 ; x_test < x1 ; x_test++) {
		int test = pData->raw_start_of_sky[x_test] ;
		if(sos[x] < test) sos[x] = test ;
	    }

	    // now lower the start of sky by W pixels
	    sos[x] += W ;
	}

	if(x < 20)
	    fprintf(stderr, "x:%d, rsos:%d sos:%d\n", x, pData->raw_start_of_sky[x], sos[x]) ;
    }

    for(int x = 0 ; x < IMAGE_WIDTH ; x++) {
	pData->raw_start_of_sky[x] = sos[x] ;
    }

    int delta_xc = 2*window_size+1 ;

    struct mstat m_avg_err[3] ;
    m_avg_err[0]=init_reg(2) ;
    m_avg_err[1]=init_reg(2) ;
    m_avg_err[2]=init_reg(2) ;

    int total_window_size = 2*window_size+1 ;

    float **wgts = (float **)calloc(total_window_size, sizeof(float *)) ;
    for(int i=0 ; i < total_window_size ; i++) {
	wgts[i] = (float *)calloc(total_window_size, sizeof(float)) ;
    }

    for(int xc = window_size ; xc+window_size < IMAGE_WIDTH ; xc += delta_xc) {
	int y_start = sos[xc]+window_size ;
	int yc = y_start ;
	int x0=xc-window_size ;
	int x1=xc+window_size ;
	if(x0 < 0) x0=0 ;
	if(x1 > IMAGE_WIDTH-1) x1=IMAGE_WIDTH-1 ;

	int debug_window=0 ;
	if(debug_col_want >= x0 && debug_col_want <= x1) debug_window=1 ;

/*  	if(xc < window_size*2) debug_window = 1 ;  */

	double highest_ratio[3] = {0.,0.,0.} ; // for debugging, highest err2/avg_err2 found on debug_col_want


#define C0 0
#define C3 3
#define BAD_FIT_PIXEL_WGT 0.01

#define DEBUG_COL 81

	while(yc < IMAGE_HEIGHT-1-window_size) {
	    struct mstat m[3] ;
	    m[0]=init_reg(4) ;
	    m[1]=init_reg(4) ;
	    m[2]=init_reg(4) ;
	    int y0=yc-window_size ;
	    int y1=yc+window_size ;

	    int low_wgt=0 ;

	    for(int x=x0 ; x <= x1 ; x++) {
		for(int y=y0 ; y <= y1 ; y++) {
/*  		    if(y >= sos[x] && y <= eos[x])  */
		    if(y >= sos[x])
		    {
			double rgb[3] ;
			tif_get3cv(image,x,y,rgb) ;
/*  			double rgb[4] ;  */
/*  			tif_get4cv(image,x,y,rgb) ;  */
/*  			if(rgb[3] < HALF16) continue ;  */
			double wgt=1. ;
			if(y > eos[x] || is_in_test_mask(x,y)) wgt=.05 ; // need to keep pixels in test masks, but with lower weights
			if((pData->pixel_class[y][x] & BAD_FIT_PIXEL)) {
			    wgt = BAD_FIT_PIXEL_WGT ;
			}
#ifdef DEBUG_COL
			wgts[y-y0][x-x0] = wgt ;
			if(wgt < 1. && debug_window) {
			    low_wgt=1 ;
			}
#endif

			double xv[4] = {1., (double)x, (double)y, (double)y*(double)y} ;
			for(int c=C0 ; c < C3 ; c++) {
			    sum_reg(&m[c], xv, rgb[c]/MAX16f, wgt) ;
			}
		    }
		}
	    }
#ifdef DEBUG_COL
	    if(debug_window || low_wgt) {
		fprintf(stderr, "      ") ;
		for(int x=x0 ; x <= x1 ; x++) {
		    fprintf(stderr, " %5d", x) ;
		}
		fprintf(stderr, "\n") ;

		for(int y=y0 ; y <= y1 ; y++) {
		    fprintf(stderr, " %5d", y) ;

		    for(int x=x0 ; x <= x1 ; x++) {
			fprintf(stderr, " %5.3f", wgts[y-y0][x-x0]) ;
		    }
		    fprintf(stderr, "\n") ;
		}
		fprintf(stderr, "\n") ;
	    }
#endif

	    double B[3][5] ;

	    int n_removed=0 ;


	    double avg_err2[3] = {.0004, .0004, .0004} ;
	    double act_avg_err2[3] = {.0004, .0004, .0004} ;
	    int iter=1 ;

	    do {

		for(int c=C0 ; c < C3 ; c++) {
		    estimate_reg(&m[c], B[c]) ;
		    avg_err2[c] = m[c].sse/(double)m[c].nobs ;
		    act_avg_err2[c] = m[c].sse/(double)m[c].nobs ;
		}

/*  		if(avg_err2[0] < .0004*49.) avg_err2[0]=.0004*49. ;  */
/*  		if(avg_err2[1] < .0004*36.) avg_err2[1]=.0004*36. ;  */
/*  		if(avg_err2[2] < .0004) avg_err2[2]=.0004 ;  */
/*  		if(avg_err2[0] < .0004) avg_err2[0]=.0004 ;  */
/*  		if(avg_err2[1] < .0004) avg_err2[1]=.0004 ;  */
/*  		if(avg_err2[2] < .0004) avg_err2[2]=.0004 ;  */

		avg_err2[0] += pData->sky_min_err_squared[0] ;
		avg_err2[1] += pData->sky_min_err_squared[1] ;
		avg_err2[2] += pData->sky_min_err_squared[2] ;


		n_removed=0 ;

#ifdef DEBUG_COL
		if(debug_window) {
		    fprintf(stderr, "AVG errors, yc:%d nobs:%d iter %d", yc, m[0].nobs, iter) ;
		   
		    for(int c=C0 ; c < C3 ; c++) {
			fprintf(stderr, " avgerr%d:%5.3f", c, sqrt(act_avg_err2[c])) ;
		    }
		    fprintf(stderr, "\n") ;
		}
#endif

		for(int x=x0 ; x <= x1 ; x++) {
		    for(int y=y0 ; y <= y1 ; y++) {
			if((pData->pixel_class[y][x] & BAD_FIT_PIXEL) == 0 ) {
/*  			    if(y >= sos[x] && y <= eos[x])  */
			    if(y >= sos[x])
			    {
				if(is_in_test_mask(x,y)) continue ; // ignore pixels in test masks

				double rgb[3] ;
				tif_get3cv(image,x,y,rgb) ;
/*  				double rgb[4] ;  */
/*  				tif_get4cv(image,x,y,rgb) ;  */
/*  				if(rgb[3] < HALF16) continue ;  */

				double prod_ratios=1. ;

#ifdef DEBUG_COL
				if(x == debug_col_want) {
				    fprintf(stderr, "iter:%d x:%d, y:%d nobs:%d", iter, x, y, m[0].nobs) ;
				}
#endif
				for(int c=C0 ; c < C3 ; c++) {
				    double y_hat = B[c][0] + B[c][1]*x + B[c][2]*y + B[c][3]*y*y ;
				    double err = y_hat - rgb[c]/MAX16f ;
				    double err2 = err*err ;
				    double ratio=err2/avg_err2[c] ;
				    prod_ratios *= 1.+ratio ;

#ifdef DEBUG_COL
				    if(x == debug_col_want) {
					if(x == debug_col_want && ratio > highest_ratio[c]) highest_ratio[c] = ratio ;
					fprintf(stderr, " h,a[%3d,%3d] ratio[%d]:%f", (int)(y_hat*255.), (int)(255.*rgb[c]/MAX16f), c, ratio) ;
				    }
#endif
				}
#ifdef DEBUG_COL
				if(x == debug_col_want) {
				    fprintf(stderr, "\n") ;
				}
#endif

				if(prod_ratios > sky_ratio_threshold) {
#ifdef DEBUG_COL
				    if(debug_window) {
					fprintf(stderr, "iter:%d nobs:%d, REMOVE POINT x:%d, y:%d, prod_ratios:%f\n", iter, m[0].nobs, x, y, prod_ratios) ;
				    }

#endif

				    n_removed++ ;

				    // remove this pixel from each color component regression
				    // and add it back in as a bad pixel with low weight
				    double xv[4] = {1., x, y, y*y} ;
				    for(int d=C0 ; d < C3 ; d++) {
					sum_reg(&m[d], xv, rgb[d]/MAX16f, -wgts[y-y0][x-x0]) ;
					sum_reg(&m[d], xv, rgb[d]/MAX16f, BAD_FIT_PIXEL_WGT) ;
				    }

				    pData->pixel_class[y][x] |= BAD_FIT_PIXEL ;
				}

			    }
			}
		    }
		}

		if(m[C0].nobs <= 0) {
		    printf("Died at find eos, x0=%d, yc=%d, iter=%d\n", x0, yc, iter) ;
		    exit(1) ;
		}

		iter++ ;
	    } while(iter < 3 && n_removed > 0) ;

	    // collect some statistics about actual average SSE found in first few windows near top of sky
	    if(yc - y_start < 10 ) {
		for(int c=C0 ; c < C3 ; c++) {
		    double xv[2] = {1., xc};
		    sum_reg(&m_avg_err[c], xv, act_avg_err2[c], 1.) ;
		}
	    }

	    // check the next row of pixels beyond the window, to see if they look like sky pixels
	    int more=0 ;
	    for(int x=x0 ; x <= x1 ; x++) {
		next_row_check[x] = 0 ;
	    }

	    for(int x=x0 ; x <= x1 ; x++) {
		for(int y=y1+1 ; y <= y1+2 ; y++) {
		    if(y >= sos[x] && y <= eos[x]) {

			if(is_in_test_mask(x,y)) {
			    more++ ;
			    continue ;
			}

/*  			double rgb[4] ;  */
/*  			tif_get4cv(image,x,y,rgb) ;  */
/*  			if(rgb[3] < HALF16) {  */
/*  			    more++ ;  */
/*  			    continue ;  */
/*  			}  */

			double rgb[3] ;
			tif_get3cv(image,x,y,rgb) ;

#ifdef DEBUG_COL
			double err2_debug[3] ;
			double y_hat_debug[3] ;
#endif

			double prod_ratios=1. ;
			for(int c=C0 ; c < C3 ; c++) {
			    double y_hat = B[c][0] + B[c][1]*x + B[c][2]*y + B[c][3]*y*y ;
			    double err = y_hat - rgb[c]/MAX16f ;
			    double err2 = err*err ;

#ifdef DEBUG_COL
			    err2_debug[c] = err2 ;
			    y_hat_debug[c] = y_hat ;
#endif
			    double ratio=err2/avg_err2[c] ;
			    prod_ratios *= 1.+ratio ;
			}


#ifdef DEBUG_COL
			int failed = (prod_ratios < sky_ratio_threshold) ? 0 : 1 ;
			if(x == debug_col_want || (debug_window && failed)) {

			    fprintf(stderr, "more %s x:%d, y:%d eos:%d", failed ? "FALSE" : "TRUE", x, y, eos[x]) ;

			    for(int c=C0 ; c < C3 ; c++) {
				fprintf(stderr, " ratio%d:%5.3f", c, err2_debug[c]/avg_err2[c]) ;
			    }
			    fprintf(stderr, " ratioAll:%5.3f", prod_ratios) ;

			    for(int c=C0 ; c < C3 ; c++) {
				fprintf(stderr, " highest_ratio%d:%4.2f", c, highest_ratio[c]) ;
			    }

			    for(int c=C0 ; c < C3 ; c++) {
				fprintf(stderr, " avgerr%d:%5.3f", c, sqrt(act_avg_err2[c])) ;
			    }
			    fprintf(stderr, "\n") ;

			    if(failed) {
				fprintf(stderr, "\tyhat, actual ") ;
				for(int c=C0 ; c < C3 ; c++) {
				    fprintf(stderr, ": %5.3f, %5.3f :", y_hat_debug[c], rgb[c]/MAX16f) ;
				}
				fprintf(stderr, "\n\n") ;
			    }
			}
#endif

			if(prod_ratios > sky_ratio_threshold) {
			    next_row_check[x]++ ;
			    pData->pixel_class[y][x] |= BAD_FIT_PIXEL ;
			}
		    }
		}
	    }

	    more=0 ;
	    for(int x=x0 ; x <= x1 ; x++) {
		if(next_row_check[x] == 2) {
		    // this column has two out of bounds pixels in a row, mark it as eos
		    if(y1+1 < eos[x]) eos[x]=y1+1 ;
		} else {
		    if(eos[x] == IMAGE_HEIGHT-1) {
			// only increment more (for this column) if eos has not be detected in this column
			more++ ;
		    }
		}
	    }

	    if(more > 0)
	    {
		/*  			fprintf(stderr, "at xc:%d, more\n", xc) ;  */
		yc++ ;
	    } else {
		for(int x=x0 ; x <= x1 ; x++) {
		    if(y1+1 < eos[x]) eos[x]=y1+1 ;
		}
		break ;
	    }
	}

	// do this or it will loop endlessly on last window, which is set in the next step
	if(x1 == IMAGE_WIDTH-1)
	    break ;

	// set up for last window
	if(xc+delta_xc+window_size > IMAGE_WIDTH-1)
	    xc = IMAGE_WIDTH-1-delta_xc-window_size ;

    }

    for(int i=0 ; i < total_window_size ; i++) {
	free(wgts[i]) ;
    }
    free(wgts) ;

    double B[3][2] ;
    for(int c=C0 ; c < C3 ; c++) {
	estimate_reg(&m_avg_err[c], B[c]) ;
    }

    for(int x=0 ; x < IMAGE_WIDTH-1 ; x += IMAGE_WIDTH/4-1) {
	fprintf(stderr, "AVERAGE ERR2 near top of sky, nobs:%d, x:%5d, %8g, %8g, %8g\n", m_avg_err[0].nobs, x, 
	B[0][0]+B[0][1]*(double)x, B[1][0]+B[1][1]*(double)x, B[2][0]+B[2][1]*(double)x) ;
    }

#ifdef DEBUG_COL
    {
	int x = debug_col_want ;
	fprintf(stderr, "End of algorithm, end_of_sky[%d] = %d, eos[%d]=%d\n", x, pData->end_of_sky[x], x, eos[x]) ;
    }
#endif

    for(int x = 0 ; x < IMAGE_WIDTH ; x++) {
	pData->end_of_sky[x] = eos[x] ;
    }

    free(sos) ;
    free(next_row_check) ;
    free(eos_status) ;
    free(eos) ;

    return 1 ;
}

int find_end_of_sky(int16_t *end_of_sky,int16_t *start_of_sky,int w,int h,tdata_t *image,int fix_edges,int fix_slivers_flag, SKYFILL_DATA_t *pData)
{


    int rval = test_find_end_of_sky(end_of_sky,start_of_sky,w,h,image,fix_edges,fix_slivers_flag, pData) ;

    if(rval == 0) {
	// no test method was done, do default method
	if(pData->end_of_sky_method == 6) {
	    find_end_of_sky_method_6(end_of_sky,start_of_sky,w,h,image,fix_edges,fix_slivers_flag, pData, -1, pData->sky_ratio_threshold) ;
	} else {
	    find_end_of_sky_method_1(end_of_sky,start_of_sky,w,h,image,fix_edges,fix_slivers_flag, pData) ;
	}

	// and continue on to fix sky slivers etc
    } else {
	// test method was performed, return 1 to signal calling routine
	fprintf(stderr, "End of sky test method %d, no more processing...\n", rval) ;
	return 1 ;
    }

    if(fix_slivers_flag == 1) {
	fix_sky_slivers(image, pData, DEFAULT_SKY_SLIVER_REPAIR,-1) ;
    }

/*      if(0 && ! pData->full_sky_replacement && fix_edges == 1)  */
    {
	fprintf(stderr, "Fixing end of sky edges\n") ;

	// eliminate long vertical edges in the sky end, which can create
	// visible edges in the feathering result
	for(int x = 0 ; x < w-1 ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    if(pData->column_mask[x+1] == 1) continue ;

	    if( (pData->end_of_sky[x+1] - pData->end_of_sky[x]) > 2) {
		//fprintf(stderr, "set EOS %d to %d\n", x+1, pData->end_of_sky[x]+2) ;
		pData->end_of_sky[x+1] = pData->end_of_sky[x]+2 ;

		if(pData->end_of_sky[x+1] < pData->start_of_sky[x+1])
		    pData->end_of_sky[x+1] = pData->start_of_sky[x+1] ;
	    }
	}

	for(int x = w-1 ; x > 0 ; x--) {
	    if(pData->column_mask[x] == 1) continue ;
	    if(pData->column_mask[x-1] == 1) continue ;

	    if( (pData->end_of_sky[x-1] - pData->end_of_sky[x]) > 2) {
		pData->end_of_sky[x-1] = pData->end_of_sky[x]+2 ;

		if(pData->end_of_sky[x-1] < pData->start_of_sky[x-1])
		    pData->end_of_sky[x-1] = pData->start_of_sky[x-1] ;
	    }
	}
    }

/*      ZZZZ  */
    if(0 && pData->full_sky_replacement && fix_edges == 1) {
	fprintf(stderr, "Fixing DOWNWARD sky edges\n") ;

	// eliminate long vertical edges in the sky end, which can create
	// visible edges in the feathering result
	// but only correct "downward" slivers
	for(int x = 0 ; x < w-1 ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    if(pData->column_mask[x+1] == 1) continue ;

	    if( (pData->end_of_sky[x+1] - pData->end_of_sky[x]) > 2) {
		// x+1 is a potential downward sliver, look for next "good" x
		int xgood ;
		for(xgood = x+2 ; xgood < IMAGE_WIDTH ; xgood++) {
		    float slope = (float)(pData->end_of_sky[xgood] - pData->end_of_sky[x]) / (float)(xgood-x) ;
		    if(fabs(slope) < 4.) {
			for(int x0=x+1 ; x0 < xgood ; x0++) {
			    pData->end_of_sky[x0] = pData->end_of_sky[x]+slope*(float)(x0-x) ;
			}
			x = xgood-1 ;
			break ;
		    }
		}
	    }
	}
    }

    // set lowest_sky_py_found
    pData->lowest_sky_py_found = 1 ;
    for(int x = 0 ; x < w ; x++) {
	if(pData->column_mask[x] == 1) continue ;
	float py = 1.-(float)pData->end_of_sky[x]/(float)(h-1) ;

	if(pData->lowest_sky_py_found > py) {
	    pData->lowest_sky_py_found = py ;
	}
    }

    return 0 ;
}

float guassian(float x, float sigma)
{
    return expf(-((x*x))/(2 * (sigma*sigma))) / (2. * M_PI * (sigma*sigma));
}

float guassian_presq(float x_sq, float sigma_sq)
{
    return expf(-((x_sq))/(2 * (sigma_sq))) / (2. * M_PI * (sigma_sq));
}

void bl_set_output_grey_value(tdata_t *image, float *output, int x, int ystart, int ymax, float wr, float wg, float wb, int width, float alpha, float beta)
{
    float r, g, b;
    float r0, g0, b0; // rgb of center pixel
    float alpha_sq = alpha*alpha ;
    float beta_sq = beta*beta ;

    for(int y = ystart ; y < ymax ; y++) {
	output[y] = 0. ;
	float sum_wgts=0. ;
	float sum=0. ;
	tif_get3c(image,x,y,r0,g0,b0) ;
	r0 *= wr/MAX16f ;
	g0 *= wg/MAX16f ;
	b0 *= wb/MAX16f ;

	for(int sx=x-width ; sx <= x+width ; sx++) {
	    if(sx < 0 || sx > IMAGE_WIDTH-1) continue ;

	    for(int sy = y-width ; sy <= y+width ; sy++) {
		if(sy < 0 || sy > IMAGE_HEIGHT-1) continue ;

		if(sy == y && sx == x) {
		    sum += r0 + g0 + b0 ;
		    sum_wgts += 1. ;
		} else {
		    tif_get3c(image,sx,sy,r,g,b) ;
		    r *= wr/MAX16f ;
		    g *= wg/MAX16f ;
		    b *= wb/MAX16f ;
		    int dx = x-sx ;
		    int dy = y-sy ;
		    int dist_sq= dx*dx + dy*dy ;
		    float dist_wgt = guassian_presq(dist_sq,alpha_sq) ;
		    float dr = r-r0 ;
		    float dg = g-g0 ;
		    float db = b-b0 ;
		    float dist_sq_color = dr*dr + dg*dg + db*db ;
		    float color_wgt = guassian_presq(dist_sq_color,beta_sq) ;
/*  		    float color_wgt = 1. ;  */
		    float wgt = dist_wgt * color_wgt ; 
		    sum += wgt *(r+g+b) ;
		    sum_wgts += wgt ;
		}

	    }
	}

	output[y] = sum/sum_wgts ;
    }
}

void set_output_grey_value(tdata_t *image, float *output, int x, int ystart, int ymax, float wr, float wg, float wb)
{
    float r, g, b;
    for(int y=ystart ; y < ymax ; y++) {
	tif_get3c(image,x,y,r,g,b) ;
	output[y] = (r*wr+g*wg+b*wb)/MAX16f ;
    }
}

void emboss1d(float *in, float *out, int ystart, int ymax)
{
    for(int y=ystart ; y < ymax ; y++) {
	out[y] =in[y-2] ; ;
	out[y]+=in[y-1] ;

	out[y]-=in[y+1] ;
	out[y]-=in[y+2] ;
	out[y] += 0.5 ;
    }
}

void bilateral_filter1d(float *in, float *out,int x,int ystart,int ymax,int width,float alpha,float beta)
{
    float g_d[width+1] ;
    for(int d=0 ; d < width+1 ; d++) {
	g_d[d] = guassian((float)d, alpha) ;
    }

    float g_r[100] ;
    for(int r=0 ; r < 100 ; r++) {
	g_r[r] = guassian((float)r/99., beta) ;
    }

    // 1 dimensional bilateral filter
    for(int y = ystart+5 ; y < ymax-5 ; y++) {
	float sum_wgts=0. ;
	float new=0. ;
	for(int sy = y-5 ; sy <= y+5 ; sy++) {
	    int dist=sy-y ;
	    if(dist < 0 ) dist = -dist ;
	    float wgt_dist = g_d[dist] ;

	    float value_diff = fabsf(in[y]-in[sy]) ;
	    int vdist = (int)(value_diff*99.+0.5) ;
	    if(vdist > 99) vdist = 99 ;
	    float wgt_value = g_r[vdist] ;

	    float wgt = wgt_dist * wgt_value ;
	    new += in[sy]*wgt ;
	    sum_wgts += wgt ;
	}

	out[y] = new/sum_wgts ;
    }
}

int add_to_rgb_sums(tdata_t *image, int x, int y, float *sum, float *sumsq)
{
    if(is_in_test_mask(x,y)) return 0 ;

    for(int c=0 ; c < 3 ; c++) {
	float value = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+c] ;
	value /= MAX16f ;
	sum[c] += value ;
	sumsq[c] += value*value ;
    }

    return 1 ;
}

int subtract_from_rgb_sums(tdata_t *image, int x, int y, float *sum, float *sumsq)
{
    if(is_in_test_mask(x,y)) return 0 ;

    for(int c=0 ; c < 3 ; c++) {
	float value = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+c] ;
	value /= MAX16f ;
	sum[c] -= value ;
	sumsq[c] -= value*value ;
    }

    return 1 ;
}

int get_rgb_sums(tdata_t *image, int x, int y0, int y1, float *sum, float *sumsq)
{
    for(int c=0 ; c < 3 ; c++) {
	sum[c]=0. ;
	sumsq[c]=0. ;
    }

    int n=0 ;

    for(int y=y0 ; y < y1 ; y++) {

	if(is_in_test_mask(x,y)) continue ;
	n++ ;

	for(int c=0 ; c < 3 ; c++) {
	    float value = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+c] ;
	    value /= MAX16f ;
	    sum[c] += value ;
	    sumsq[c] += value*value ;
	}
    }

    return n ;
}

int get_rgb_sum(tdata_t *image, int x, int y0, int y1, float *sum)
{
    for(int c=0 ; c < 3 ; c++) {
	sum[c]=0. ;
    }

    int n=0 ;

    for(int y=y0 ; y < y1 ; y++) {

	if(is_in_test_mask(x,y)) continue ;
	n++ ;

	for(int c=0 ; c < 3 ; c++) {
	    float value = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+c] ;
	    value /= MAX16f ;
	    sum[c] += value ;
	}
    }

    return n ;
}

int test_find_eos_at_column(tdata_t *image, SKYFILL_DATA_t *pData, int x, float alpha, float beta, FILE *fp)
{
    float ed_hat=0. ;
    float cv_diff_hat=0. ;
    float ed_min=1.e30 ;
    float ed_max=0. ;
    int N_ABOVE=20 ;
    int N_NEXT=5 ;
    int first_y = pData->raw_start_of_sky[x]+N_ABOVE ;
    float sum1[3], sumsq1[3] ;
    float sum2[3], sumsq2[3] ;
    float *ed_by_y = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_hue = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_sat = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_val = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    is_clear_sky = (uint8_t *)calloc(IMAGE_HEIGHT, sizeof(uint8_t)) ;

    float save_hue_tol = pData->sky_hue_tolerance ;
    float save_sat_tol = pData->sky_sat_tolerance ;
    float save_val_tol = pData->sky_val_tolerance ;


    for(int y = 0 ; y < IMAGE_HEIGHT ; y++) {
	ed_by_y[y] = -1. ;
	is_clear_sky[y] = 1 ;
	if(y < pData->raw_start_of_sky[x]) {
	    sky_hue[y] = 0. ;
	    sky_sat[y] = 0. ;
	    sky_val[y] = 0. ;
	} else {
	    uint16_t rgb[3] ;
	    float hsv[3] ;
	    tif_get3cv(image,x,y,rgb) ;
	    rgb2hsv16(rgb,hsv) ;
	    sky_hue[y] = hsv[0] ;
	    sky_sat[y] = hsv[1] ;
	    sky_val[y] = hsv[2] ;
	}
    }

    if(fp != NULL) {
	fprintf(fp, "y ed arma hratio sratio vratio\n") ;
    }



    int n1 = get_rgb_sums(image, x, first_y-N_ABOVE, first_y-1, sum1, sumsq1) ;
    int n2 = get_rgb_sums(image, x, first_y, first_y+N_NEXT-1, sum2, sumsq2) ;

    int n_rows_processed=0 ;
    int y0 = first_y+N_ABOVE ;
    int analysis_y = y0 ;

    for(int y = y0 ; y < IMAGE_HEIGHT-4-N_NEXT ; y++) {
	n1 += add_to_rgb_sums(image, x, y, sum1, sumsq1) ;
	n2 += add_to_rgb_sums(image, x, y+N_NEXT-1, sum2, sumsq2) ;

	// compare the mean value of the previous 10 pixels,
	// to the  mean value of the next 5 pixels
	/*  		get_rgb_sums(image, x, y-N_ABOVE, y, sum1, sumsq1) ;  */
	/*  		get_rgb_sums(image, x, y, y+N_NEXT, sum2, sumsq2) ;  */
	float mean_cv_diff=0. ;
	float euclidean_dist=0. ;
	float mean_cv1 = 0. ;

	for(int c=0 ; c < 3 ; c++) {
	    float mean1 = sum1[c] / (float)n1 ;
	    float mean2 = sum2[c] / (float)n2 ;
	    float dist = (mean2-mean1)/mean1 ;
	    euclidean_dist += dist*dist ;
	    float cv1, cv2 ;
	    cv1 = cv_calc(n1, sum1[c], sumsq1[c]) ;
	    cv2 = cv_calc(n2, sum2[c], sumsq2[c]) ;
	    mean_cv_diff += fabsf(cv2-cv1) ;
	    mean_cv1 += cv1 ;
	}

	euclidean_dist = sqrt(euclidean_dist) ;
	ed_by_y[y] = euclidean_dist ;

	if(n_rows_processed < 10) {
	    if(ed_min > euclidean_dist) ed_min = euclidean_dist ;
	    if(ed_max < euclidean_dist) ed_max = euclidean_dist ;
	    if(cv_diff_hat < mean_cv_diff) cv_diff_hat = mean_cv_diff ;
	} else {

	    if(y - analysis_y < 20) continue ;

	    float spread_hat = ed_max - ed_min ;
	    float spread_this = euclidean_dist - ed_min ;

	    float sum=0., n=0. ;
	    for(int y0 = analysis_y-20 ; y0 < analysis_y+20 ; y0++) {
		if(y0 < first_y) continue ;
		if(y0 > IMAGE_HEIGHT-1) continue ;
		if(ed_by_y[y0] < 0.) continue ;
		sum += ed_by_y[y0] ;
		n += 1. ;
	    }

	    float arma_ed = sum / n ;
	    float hdiff, sdiff, vdiff ;
	    int end_flag = is_end_of_sky(analysis_y,5,&hdiff,&sdiff,&vdiff,pData) ;
	    
	    float hratio = hdiff / pData->sky_hue_tolerance ;
	    float sratio = sdiff / pData->sky_sat_tolerance ;
	    float vratio = vdiff / pData->sky_val_tolerance ;


	    if(sratio < 1.) {
		float reduce = MAX(sratio, alpha) ;
		pData->sky_sat_tolerance *= reduce ;
	    }

	    if(vratio < 1.) {
		float reduce = MAX(vratio, alpha) ;
		pData->sky_val_tolerance *= reduce ;
	    }

	    if(fp != NULL) {
		fprintf(fp, "%d %f %f %f %f %f\n", analysis_y, ed_by_y[analysis_y], arma_ed, hratio, sratio, vratio) ;
	    }

	    if(fp == NULL) {
		if(analysis_y > IMAGE_HEIGHT/2. || hratio > 1. || sratio > 1. || vratio > 1.) {
		    pData->sky_hue_tolerance = save_hue_tol ;
		    pData->sky_sat_tolerance = save_sat_tol ;
		    pData->sky_val_tolerance = save_val_tol ;
		    free(sky_hue) ;
		    free(sky_sat) ;
		    free(sky_val) ;
		    free(is_clear_sky) ;
		    free(ed_by_y) ;
		    pData->end_of_sky[x] = analysis_y+1 ;
		    return pData->end_of_sky[x] ;
		}
	    } else {

		if(euclidean_dist < ed_min) {
		    ed_min = MAX(ed_min*alpha, euclidean_dist) ;
		} else if(euclidean_dist > ed_max) {
		    ed_max = MIN(ed_max/alpha, euclidean_dist) ;
/*  		    ed_min = ed_min/alpha ;  */
		} else {
		    ed_max = ed_max*(1.-(1.-alpha)/2.) ;
		}
	    }

	    analysis_y++ ;
	}

	n_rows_processed++ ;


	n1 -= subtract_from_rgb_sums(image, x, y-N_ABOVE, sum1, sumsq1) ;
	n2 -= subtract_from_rgb_sums(image, x, y, sum2, sumsq2) ;

    }

    pData->sky_hue_tolerance = save_hue_tol ;
    pData->sky_sat_tolerance = save_sat_tol ;
    pData->sky_val_tolerance = save_val_tol ;
    free(sky_hue) ;
    free(sky_sat) ;
    free(sky_val) ;
    free(is_clear_sky) ;
    free(ed_by_y) ;
    return IMAGE_HEIGHT-1 ;
}


int test_find_end_of_sky(int16_t *end_of_sky,int16_t *start_of_sky,int w,int h,tdata_t *image,int fix_edges,int fix_slivers_flag, SKYFILL_DATA_t *pData)
{
    int x0=0 ;
    int x1=w-1 ;

    FILE *fp = fopen("eos_test.prms", "r") ;

    if(fp != NULL) {
	int eos_method ;
	int debug_col_want ;
	int window_size = 20 ;
	float alpha, beta ;

	if(fscanf(fp, "%d %d %f %f", &eos_method, &debug_col_want, &alpha, &beta) != 4) {
	    fprintf(stderr, "something weird happened reading eos_test.prms\n") ;
	    exit(1) ;
	} else {
	    int tmp ;
	    if(fscanf(fp, "%d", &tmp) == 1) {
		window_size=tmp ;
	    }
	}

	fclose(fp) ;
	fprintf(stderr, "EOS test start processing, method=%d,window size=%d, alpha=%f, beta=%f\n", eos_method, window_size, alpha, beta) ;

	if(eos_method == 6) {
	    find_end_of_sky_method_6(end_of_sky,start_of_sky,w,h,image,fix_edges,fix_slivers_flag, pData, debug_col_want, alpha) ;
	    return 6 ;
	}


	if(eos_method == 5) {
	    // autoregression method,  with overlapping windows for first two samples, non-overlapping last sample
	    // removal of observations when n=10, so the regression can adapt to changing sky hue/sat/val as the
	    // process moves from the top most sky pixels towards the horizon
	    if(0) {
		FILE *fp_out = fopen("find_sky.dat", "w") ;
		fprintf(fp_out, "y r g b\n") ;
		float sum[3] ;

		int first_y = pData->raw_start_of_sky[debug_col_want]+window_size/2 ;

		for(int y = first_y ; y < IMAGE_HEIGHT-1-window_size ; y += window_size/2) {
		    get_rgb_sum(image, debug_col_want, y-window_size/2, y+window_size/2, sum) ;
		    fprintf(fp_out, "%d %f %f %f\n", y, sum[0], sum[1], sum[2]) ;
		}

		fclose(fp_out) ;
		exit(1) ;
	    } else {
		float *sums = (float *)calloc(IMAGE_HEIGHT*3, sizeof(float)) ;

		for(int x = 0 ; x < IMAGE_WIDTH-1 ; x++) {
		    float sum[3] ;
		    struct mstat m[3] ;
		    int have_left_x = x > 0 ? 1 : 0 ;
		    int have_right_x = x < IMAGE_WIDTH-1 ? 1 : 0 ;

		    for(int y = 0 ; y < IMAGE_HEIGHT-1 ; y++) {
			sums[3*y+0] = -1. ;
			sums[3*y+1] = -1. ;
			sums[3*y+2] = -1. ;
		    }

		    for(int c=0 ; c < 3 ; c++)
			m[c] = init_reg(2) ;

		    int first_y = pData->raw_start_of_sky[x]+window_size/2 ;

/*  		    for(int y = first_y ; y < IMAGE_HEIGHT-1-window_size ; y += window_size/2)  */
		    for(int y = first_y ; y < IMAGE_HEIGHT-1-window_size ; y ++ )
		    {
			int n = get_rgb_sum(image, x, y-window_size/2, y+window_size/2, sum) ;
			sums[3*y+0] = sum[0] ;
			sums[3*y+1] = sum[1] ;
			sums[3*y+2] = sum[2] ;

			if(have_left_x) {
			    n += get_rgb_sum(image, x-1, y-window_size/2, y+window_size/2, sum) ;
			    sums[3*y+0] += sum[0] ;
			    sums[3*y+1] += sum[1] ;
			    sums[3*y+2] += sum[2] ;
			}

			if(have_right_x) {
			    n += get_rgb_sum(image, x+1, y-window_size/2, y+window_size/2, sum) ;
			    sums[3*y+0] += sum[0] ;
			    sums[3*y+1] += sum[1] ;
			    sums[3*y+2] += sum[2] ;
			}

			if(n > 0) {
			    sums[3*y+0] /= (float)n ;
			    sums[3*y+1] /= (float)n ;
			    sums[3*y+2] /= (float)n ;
			}

			if(y < first_y+window_size+window_size/2) continue ;

			float err[3] ;
			float err2[3] ;
			float mse[3] ;

			if(m[0].nobs > 2) {
			    double B[2] ;
			    for(int c=0 ; c < 3 ; c++) {
				estimate_reg(&m[c], B) ;
#define AR_Y0Y1(y,y0,y1) { y1=(y)-window_size ; y0=y1-window_size/2 ; }
				int y0,y1 ;
				AR_Y0Y1(y,y0,y1) ;
				float s0 = sums[3*y0+c] ;
				float s1 = sums[3*y1+c] ;
				float s_this = sums[3*y+c] ;
				if(s0 < 0. || s1 < 0.) {
				    fprintf(stderr, "WOH: y:%d y0:%d y1:%d c:%d s0:%f s1%f\n", y, y0, y1, c, s0, s1) ;
				    exit(1) ;
				}

				float hat = B[0]*s0 + B[1]*s1 ;
				err[c] = fabsf(hat - s_this) ;
				err2[c] = err[c]*err[c] ;
				mse[c] = m[c].sse/(float)m[c].sum_wgts ;
				// scale mse until nobs gets to 4 to avoid unreasonably low mse
				mse[c] *= (m[c].sum_wgts < 3.99 ? pow((4./(float)m[c].sum_wgts),10.) : 1.) ;
			    }

			    float mean_err = (err[0]+err[1]+err[2])/3. ;
			    float mean_err2 = (err2[0]+err2[1]+err2[2])/3. ;
			    float mean_mse= (mse[0]+mse[1]+mse[2])/3. ;

			    if(x == debug_col_want)
				printf("x:%d y:%d sum_w:%f sqrtmean_err2:%10g, mse:%10g ratio %f\n", x, y, m[0].sum_wgts, sqrt(mean_err2), mean_mse, mean_err2/mean_mse) ;

/*  			    if( mean_err > beta)  */
			    if( mean_err2/mean_mse > beta && mean_err2 > alpha*alpha)
			    {
/*  				printf("EOS tol limit x:%d y:%d sum_w:%f sqrtmean_err2:%10g, mse:%10g ratio %f\n", x, y, m[0].sum_wgts, sqrt(mean_err2), mean_mse, mean_err2/mean_mse) ;  */
				pData->end_of_sky[x] = y ;
				break ;
			    }
			}

			for(int c=0 ; c < 3 ; c++) {
			    int y0,y1 ;
			    AR_Y0Y1(y,y0,y1) ;
			    double xv[2] = {sums[3*y0+c], sums[3*y1+c]} ;

/*  			    if(x == debug_col_want)  */
/*  				printf("sum_reg y:%d c:%d x:%f %f dep:%f\n", y, c, xv[0], xv[1], sums[3*y+c]) ;  */

			    sum_reg(&m[c], xv, (double)sums[3*y+c], 1.) ;

			    int y_delete = y-5*window_size ;
			    if(y_delete >= first_y+window_size) {
				// as y increases, remove observations higher in the image from the regression
				int y0=y_delete - window_size ;
				int y1=y_delete - window_size/2 ;
				double xv_delete[2] = {sums[3*y0+c], sums[3*y1+c]} ;
				double c_delete = sums[3*y_delete+c] ;

				sum_reg(&m[c], xv_delete, c_delete, -1.) ;
			    }
			}
		    }

		    pData->end_of_sky[x] -= window_size ;

		    if(pData->end_of_sky[x] < pData->raw_start_of_sky[x]+20)
			pData->end_of_sky[x] = pData->raw_start_of_sky[x]+20 ;


		    if(x == debug_col_want)
			printf("x:%d eos:%d\n", x, pData->end_of_sky[x]) ;

		    pData->end_of_sky[x] = new_get_hires_eos(image, x, pData) ;

		    if(x == debug_col_want)
			printf("x:%d new hires eos:%d\n", x, pData->end_of_sky[x]) ;
		}
		free(sums) ;

		if(fix_slivers_flag == 1) {
		    //fix_sky_slivers(image, pData, DEFAULT_SKY_SLIVER_REPAIR) ;
/*  		    fix_sky_slivers(image, pData, ED_SKY_SLIVER_REPAIR,debug_col_want) ;  */
		}

		return 5 ;
	    }

	}

	if(eos_method == 4) {
	    if(0) {
		FILE *fp_out = fopen("find_sky.dat", "w") ;
		int eos = test_find_eos_at_column(image, pData, debug_col_want, alpha, beta, fp_out) ;
		fclose(fp_out) ;
		return 4 ;
	    } else {
		for(int x = 0 ; x < IMAGE_WIDTH-1 ; x++) {
		    int eos = test_find_eos_at_column(image, pData, x, alpha, beta, NULL) ;
		    pData->end_of_sky[x] = new_get_hires_eos(image, x, pData) ;
/*  		    tif_set3c(image,x,eos+2,0,MAX16,0) ;  */
		}

		if(fix_slivers_flag == 1) {
		    fix_sky_slivers(image, pData, DEFAULT_SKY_SLIVER_REPAIR,debug_col_want) ;
		}

		return 4 ;
	    }

	}

	if(eos_method == 3) {
	    int x = (int)(alpha+0.5) ;
	    int N_ABOVE=10 ;
	    int N_NEXT=5 ;
	    int first_y = pData->raw_start_of_sky[x]+N_ABOVE ;
	    FILE *fp_out = fopen("find_sky.dat", "w") ;
	    fprintf(fp_out, "x y p_r p_g p_b pdiff ed cv_diff\n") ;
	    float sum1[3], sumsq1[3] ;
	    float sum2[3], sumsq2[3] ;

	    int n1 = get_rgb_sums(image, x, first_y-N_ABOVE, first_y-1, sum1, sumsq1) ;
	    int n2 = get_rgb_sums(image, x, first_y, first_y+N_NEXT-1, sum2, sumsq2) ;

	    for(int y = first_y+N_ABOVE ; y < IMAGE_HEIGHT-1-N_NEXT ; y++) {
		n1 += add_to_rgb_sums(image, x, y, sum1, sumsq1) ;
		n2 += add_to_rgb_sums(image, x, y+N_NEXT-1, sum2, sumsq2) ;

		fprintf(fp_out, "%d %d", x, y) ;

		// compare the mean value of the previous 10 pixels,
		// to the  mean value of the next 5 pixels
/*  		get_rgb_sums(image, x, y-N_ABOVE, y, sum1, sumsq1) ;  */
/*  		get_rgb_sums(image, x, y, y+N_NEXT, sum2, sumsq2) ;  */
		int n1 = N_ABOVE ;
		int n2 = N_NEXT ;
		float mean_p=0. ;
		float mean_diff=0. ;
		float mean_cv_diff=0. ;
		float euclidean_dist=0. ;

		for(int c=0 ; c < 3 ; c++) {
		    float t = means_t_test(n1, sum1[c],  sumsq1[c], n2, sum2[c], sumsq2[c], 0) ;
		    float mean1 = sum1[c] / (float)n1 ;
		    float mean2 = sum2[c] / (float)n2 ;
		    float dist = (mean2-mean1)/mean1 ;
		    mean_diff += fabsf(mean1-mean2) ;
		    euclidean_dist += dist*dist ;
		    float cv1, cv2 ;
		    cv1 = cv_calc(n1, sum1[c], sumsq1[c]) ;
		    cv2 = cv_calc(n2, sum2[c], sumsq2[c]) ;
		    mean_cv_diff += fabsf(cv2-cv1) ;

		    float p = Student_t_probability((double)(fabsf(t)), N_ABOVE-1) ;
		    fprintf(fp_out, " %7.6f", p) ;
		    mean_p += p ;
		    //fprintf(fp_out, " (%5.3f:%5.3f) t:%7.1f", mean1,mean2, t) ;
/*  		    fprintf(fp_out, " :%f %f", sum1[c], sum2[c]) ;  */
		}

		fprintf(fp_out, "   %5.4f  %5.4f %f\n", mean_p/3.* mean_diff/3., sqrt(euclidean_dist), mean_cv_diff/3.) ;
		n1 -= subtract_from_rgb_sums(image, x, y-N_ABOVE, sum1, sumsq1) ;
		n1 -= subtract_from_rgb_sums(image, x, y, sum2, sumsq2) ;

	    }

	    fclose(fp_out) ;
	    exit(1) ;
	}


	if(eos_method == 2) {
	    int delta=10 ;
	    int ymax = IMAGE_HEIGHT-delta ;  ;
	    // embossing technique
	    float *output = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
	    float *bl_output = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
	    float *exp_py = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;

	    for(int y = 0 ; y < IMAGE_HEIGHT ; y++) {
		float py = (float)y/(float)IMAGE_HEIGHT ;
		exp_py[y] = exp(py) ;
	    }

	    for(int x = x0 ; x <= x1 ; x++) {
		int ystart=pData->raw_start_of_sky[x]+5 ;
		bl_set_output_grey_value(image,output,x,ystart-5,IMAGE_HEIGHT,.2,.2,.6,2,alpha,beta) ;
		bilateral_filter1d(output,bl_output,x,ystart,ymax,5,alpha,beta) ;
		emboss1d(bl_output,output,ystart+5,IMAGE_HEIGHT-5) ;

		int done=0 ;

		for(int y = pData->raw_start_of_sky[x]+delta ; !done && y < ymax-delta ; y++) {
		    if(y < pData->raw_start_of_sky[x]+delta) continue ;

		    int try_prediction_error=1 ;

		    if(!try_prediction_error) {
			float value = output[y]*MAX16f ;
			tif_set3c(image,x,y,value,value,value) ;
		    } else {
			struct mstat m_fit ;
			m_fit = init_reg(3) ;
			double B[3] ;

			// regress previous delta samples
			for(int sample_y = y-delta ; sample_y < y ; sample_y++) {
			    float py = (float)sample_y/(float)IMAGE_HEIGHT ;
			    double xv[3] = {1.,  py, exp_py[sample_y]} ;
			    double y = bl_output[sample_y] ;
			    sum_reg(&m_fit, xv, y, 1.) ;

			}

			estimate_reg(&m_fit, B) ;

			// what is the expected RMSE of the samples used in the fit;
			float sum2_err=0. ;
			for(int sample_y = y-delta ; sample_y < y ; sample_y++) {
			    float py = (float)sample_y/(float)IMAGE_HEIGHT ;
			    float output_hat = B[0] + B[1]*py +B[2]*exp_py[sample_y] ;
			    float err = output_hat - bl_output[sample_y] ;
			    sum2_err += err*err ;
			}

			float RMSE = sqrt(sum2_err/(float)delta) ;

			sum2_err=0. ;

			// what is the combined prediction error of the next 2 samples (y, y+1) ;
			for(int sample_y = y ; sample_y < y+2 ; sample_y++) {
			    float py = (float)sample_y/(float)IMAGE_HEIGHT ;
			    float output_hat = B[0] + B[1]*py +B[2]*exp_py[sample_y] ;
			    float err = output_hat - bl_output[sample_y] ;
			    sum2_err += err*err ;
			}

			float RMSE_next = sqrt(sum2_err/2.) ;

			float p_err = RMSE_next / RMSE ;
			p_err /= 500. ;

			if(p_err > 1.) p_err=1. ;

			if(fabsf(p_err) > 20.) {
			    if(y < 220)
				fprintf(stderr, "EOS test %d, %d: output:%f RMSE_next:%f RMSE:%f p_err:%f\n", x, y, output[y], RMSE_next, RMSE, p_err) ;
/*  			    tif_set3c(image,x,y+2,HALF16,MAX16,HALF16) ;  */
			    done=1 ;
			}

			float v=p_err * MAX16f ;
			tif_set3c(image,x,y,v,v,v) ;


		    }

		}
	    }

	    free(bl_output) ;
	    free(output) ;
	    free(exp_py) ;
	} else {

	    int delta=10 ;
	    int ymax = IMAGE_HEIGHT-delta ;  ;
	    float *tvals_r = (float *)calloc(ymax, sizeof(float)) ;
	    float *tvals_g = (float *)calloc(ymax, sizeof(float)) ;
	    float *tvals_b = (float *)calloc(ymax, sizeof(float)) ;

	    for(int x = x0 ; x <= x1 ; x++) {

		for(int y = delta ; y < ymax ; y++) {
		    if(y < pData->raw_start_of_sky[x]+delta) continue ;


		    float sum1_r=0., sum1_g=0., sum1_b=0 ;
		    float sum1sq_r=0., sum1sq_g=0., sum1sq_b=0 ;
		    int n1=0 ;

		    float sum2_r=0., sum2_g=0., sum2_b=0 ;
		    float sum2sq_r=0., sum2sq_g=0., sum2sq_b=0 ;
		    int n2=0 ;

		    for(int sample_y = y-delta ; sample_y < y+1 ; sample_y++) {
			float r, g, b ;

			r = ((uint16_t *)(image[sample_y]))[IMAGE_NSAMPLES*x+0] ;
			g = ((uint16_t *)(image[sample_y]))[IMAGE_NSAMPLES*x+1] ;
			b = ((uint16_t *)(image[sample_y]))[IMAGE_NSAMPLES*x+2] ;

			if(sample_y < y) {
			    sum1_r += (float)r ;
			    sum1_g += (float)g ;
			    sum1_b += (float)b ;
			    sum1sq_r += (float)r*(float)r ;
			    sum1sq_g += (float)g*(float)g ;
			    sum1sq_b += (float)b*(float)b ;
			    n1++ ;
			} else {
			    sum2_r += (float)r ;
			    sum2_g += (float)g ;
			    sum2_b += (float)b ;
			    sum2sq_r += (float)r*(float)r ;
			    sum2sq_g += (float)g*(float)g ;
			    sum1sq_b += (float)b*(float)b ;
			    n2++ ;
			}
		    }

		    float dist_r = fabsf((sum1_r/(float)n1-sum2_r/(float)n2)) ;
		    float dist_g = fabsf((sum1_g/(float)n1-sum2_g/(float)n2)) ;
		    float dist_b = fabsf((sum1_b/(float)n1-sum2_b/(float)n2)) ;

		    float dist_sq = (dist_r*dist_r+dist_g*dist_g+dist_b*dist_b) ;

		    float  t_r = 1. ;
		    float  t_g = 1. ;
		    float  t_b = 1. ;

		    t_r = means_t_test(n1, sum1_r,  sum1sq_r, n2, sum2_r, sum2sq_r, 0) ;
		    t_g = means_t_test(n1, sum1_g,  sum1sq_g, n2, sum2_g, sum2sq_g, 0) ;
		    t_b = means_t_test(n1, sum1_b,  sum1sq_b, n2, sum2_b, sum2sq_b, 0) ;
		    t_r = fabsf(t_r) ;
		    t_g = fabsf(t_g) ;
		    t_b = fabsf(t_b) ;

		    float max_t = t_r ;
		    if(max_t < t_g) max_t = t_g ;
		    if(max_t < t_b) max_t = t_b ;

		    // but if the euclidean distance is small, it doesn't matter so set max_t to a low value
		    if(dist_sq < .02*.02) max_t =sqrtf(dist_sq) ;
		    max_t=log(1.+max_t)*0.1 ;
		    if(max_t > 1.0) max_t=1.0 ;
		    max_t = 1.-max_t ;

		    t_r = max_t ;
		    t_g = max_t ;
		    t_b = max_t ;

		    t_r *= MAX16f ;
		    t_g *= MAX16f ;
		    t_b *= MAX16f ;

		    tvals_r[y] = t_r ;
		    tvals_g[y] = t_g ;
		    tvals_b[y] = t_b ;

		}

		for(int y = delta ; y < ymax ; y++) {
		    tif_set4c(image,x,y, tvals_r[y], tvals_g[y], tvals_b[y], MAX16);
		}

	    }

	    free(tvals_r) ;
	    free(tvals_g) ;
	    free(tvals_b) ;
	}

	fprintf(stderr, "EOS test end processing\n") ;

	return 1 ;

    } else {
	fprintf(stderr, "Did not find eos_test.prms, default end of sky processing\n") ;
	return 0 ;
    }
}
