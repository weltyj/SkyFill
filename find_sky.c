#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>

#include "skyfill_tif.h"
#include "colorspace_conversions.h"
#include "pixel_tests.h"
#include "find_sky.h"

// some shared variables in this module
static float *sky_hue, *sky_val, *sky_sat  ;
static uint8_t *is_clear_sky ;
static int n_sky ;

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


void get_sky_mean_var(tdata_t *image, int x, int y0, int n, float mean[], float var[])
{
    float sum[3] = {0., 0., 0.} ;  // for r,g,b
    float sumsq[3] = {0., 0., 0.} ;  // for r,g,b
    int y ;

    for(y = y0 ; y < y0+n ; y++) {
	uint16_t r,g,b ;
	tif_get3c(image,x,y,r,g,b) ;

	sum[0] += r ;
	sum[1] += g ;
	sum[2] += b ;
	sumsq[0] += r*r ;
	sumsq[1] += g*g ;
	sumsq[2] += b*b ;
    }

    mean[0] = sum[0]/(float)n ;
    mean[1] = sum[1]/(float)n ;
    mean[2] = sum[2]/(float)n ;

    var[0] = sumsq[0] - sum[0]*mean[0] ; // equiv to sumsq - sum*sum/n ;
    var[1] = sumsq[1] - sum[1]*mean[1] ; // equiv to sumsq - sum*sum/n ;
    var[2] = sumsq[2] - sum[2]*mean[2] ; // equiv to sumsq - sum*sum/n ;
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
	m_hue1 += sky_hue[iy]*is_clear_sky[iy] ;
	m_sat1 += sky_sat[iy]*is_clear_sky[iy] ;
	m_val1 += sky_val[iy]*is_clear_sky[iy] ;
	sum_wgts += is_clear_sky[iy] ;
    }

    m_hue1 /= sum_wgts ;
    m_sat1 /= sum_wgts ;
    m_val1 /= sum_wgts ;

    float huediff = fabs(m_hue0 - m_hue1) ;
    float satdiff = fabs(m_sat0 - m_sat1) ;
    float valdiff = fabs(m_val0 - m_val1) ;

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

    return fabs(Gx_r) + fabs(Gx_g) + fabs(Gx_b) + (fabs(Gy_r) + fabs(Gy_g) + fabs(Gy_b))/2. ;
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

int get_sobel_eos(tdata_t *image, int x, SKYFILL_DATA_t *pData)
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

	    float dr = fabs(r1-r0) ;
	    float dg = fabs(g1-g0) ;
	    float db = fabs(b1-b0) ;

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



void find_end_of_sky(int16_t *end_of_sky,int16_t *start_of_sky,int w,int h,tdata_t *image,int fix_edges,int fix_slivers, SKYFILL_DATA_t *pData)
{
    int16_t x, y ;

    is_clear_sky = (uint8_t *)calloc(IMAGE_HEIGHT, sizeof(uint8_t)) ;

    sky_hue = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_val = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_sat = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;

    // look for end of sky now

    for(x = 0 ; x < w ; x++) {

	if(pData->column_mask[x] == 1) {
	    pData->end_of_sky[x] = 0 ;
	    continue ;
	}

	float hdiff, sdiff, vdiff ;

	// start looking at maximum start of sky in local area
	int first_y = pData->raw_start_of_sky[x] ;
	int x0 = x-20 ;
	if(x0 < 0) x0 = 0 ;
	int x1 = x+20 ;
	if(x1 > IMAGE_WIDTH-1) x1 = IMAGE_WIDTH-1 ;

	for(int dx=x0 ; dx <= x1 ; dx++) {
	    if(pData->column_mask[dx] == 1) continue ;

	    if(first_y < pData->raw_start_of_sky[dx])
		first_y = pData->raw_start_of_sky[dx]  ;
	}

	first_y += 10 ;


	//first_y = end_of_sky[x]-6 ;  // temporarily see how this looks 
	float top_mean[3], top_var[3] ;
	float bot_mean[3], bot_var[3] ;

	//fprintf(stderr, "\nNEW EOS x:%5d eos:%5d\n", x, end_of_sky[x]) ;

	for(y = first_y ; y < IMAGE_HEIGHT-10 ; y++) {
	    get_sky_mean_var(image, x, y-10, 10, top_mean, top_var) ;
	    get_sky_mean_var(image, x, y+5, 5, bot_mean, bot_var) ;
	    float rr = fabs(top_mean[0]/bot_mean[0]-1.) ;
	    float rg = fabs(top_mean[1]/bot_mean[1]-1.) ;
	    float rb = fabs(top_mean[2]/bot_mean[2]-1.) ;

	    // do we have the same color above and below the test area?
	    if(rr > .1 || rg > .1 || rb > .1) {
		for(; y < IMAGE_HEIGHT-10 ; y++) {
		    is_clear_sky[y]=1 ;
		}
		break ;
	    }

	    is_clear_sky[y]=1 ;

	    // if transparent or black, must be at bottom of image
	    uint16_t r,g,b,a ;
	    tif_get4c(image,x,y+10,r,g,b,a) ;

	    if(a < HALF16) break ;
	    if( r < BLK16 && g < BLK16 && b < BLK16 ) break ;

	    int last_y = y ;


	    for(int y0 = y ; y0 < y+5 ; y0++) {
		// try to skip small highlights (jet trails) in sky
		tif_get3c(image,x,y,r,g,b) ;

		float rtop=(float)r/top_mean[0] ;
		float gtop=(float)g/top_mean[1] ;
		float btop=(float)b/top_mean[2] ;
		float rbot=(float)r/bot_mean[0] ;
		float gbot=(float)g/bot_mean[1] ;
		float bbot=(float)b/bot_mean[2] ;

		if( rtop > 1. && gtop > 1. && btop > 1. && rbot > 1. && gbot > 1. && bbot > 1.) {
		    is_clear_sky[y0] = 0 ;
		    last_y = y0 ;
		} else {
		    is_clear_sky[y0] = 1 ;
		}
	    }

	    y = last_y ;

	}


	// load up arrays

	for(y = first_y ; y < first_y+10 ; y++) {
	    uint16_t r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;
	    float l = (.299*(float)r+.587*(float)g+.114*(float)b)/MAX16f ;
	    sky_hue[y] = h ;
	    sky_sat[y] = s ;
	    sky_val[y] = l ;
	}

	for(y = first_y ; y < h-10 ; y++) {
	    end_of_sky[x] = y+4 ;

	    if(is_end_of_sky(y+4,5,&hdiff,&sdiff,&vdiff,pData)) {
		break ;
	    }

	    uint16_t r,g,b ;
	    tif_get3c(image,x,y+10,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;
	    float l = (.299*(float)r+.587*(float)g+.114*(float)b)/MAX16f ;
	    sky_hue[y+10] = h ;
	    sky_sat[y+10] = s ;
	    sky_val[y+10] = l ;
	}

	// fine tune the detection -- increase y as long as vdiff becomes larger
	float vdiff_prev=vdiff ;

	for(y = pData->end_of_sky[x]+1 ; y < h-10 ; y++) {
	    uint16_t r,g,b ;
	    tif_get3c(image,x,y+5,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    float l = (.299*(float)r+.587*(float)g+.114*(float)b)/MAX16f ;
	    sky_hue[y+5] = h ;
	    sky_sat[y+5] = s ;
	    sky_val[y+5] = l ;

	    tif_get3c(image,x,y+6,r,g,b) ;

	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    l = (.299*(float)r+.587*(float)g+.114*(float)b)/MAX16f ;
	    sky_hue[y+6] = h ;
	    sky_sat[y+6] = s ;
	    sky_val[y+6] = l ;

	    is_end_of_sky(y,5,&hdiff,&sdiff,&vdiff,pData) ;

	    //if(x == 232)
		//printf("fty:sdiff:vdiff %d:%f:%f\n", y, sdiff, vdiff) ;

	    if(vdiff <= vdiff_prev)
		break ;

	    vdiff_prev = vdiff ;

	    pData->end_of_sky[x]++ ;
	}

	pData->end_of_sky[x]-=2 ;
	pData->end_of_sky[x] = get_sobel_eos(image,x,pData) ;
    }

    //look for and fix slivers
    fix_slivers=0 ;
    while(fix_slivers >0) {
	fprintf(stderr, "Fixing slivers\n") ;

	for(x = 1 ; x < w-3 ; x++) {
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

    if(! pData->full_sky_replacement && fix_edges == 1) {
	fprintf(stderr, "Fixing end of sky edges\n") ;

	// eliminate long vertical edges in the sky end, which can create
	// visible edges in the feathering result
	for(x = 0 ; x < w-1 ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    if(pData->column_mask[x+1] == 1) continue ;

	    if( (pData->end_of_sky[x+1] - pData->end_of_sky[x]) > 2) {
		//fprintf(stderr, "set EOS %d to %d\n", x+1, pData->end_of_sky[x]+2) ;
		pData->end_of_sky[x+1] = pData->end_of_sky[x]+2 ;
	    }
	}

	for(x = w-1 ; x > 0 ; x--) {
	    if(pData->column_mask[x] == 1) continue ;
	    if(pData->column_mask[x-1] == 1) continue ;

	    if( (pData->end_of_sky[x-1] - pData->end_of_sky[x]) > 2) {
		pData->end_of_sky[x-1] = pData->end_of_sky[x]+2 ;
	    }
	}
    }

/*      ZZZZ  */
    if(0 && pData->full_sky_replacement && fix_edges == 1) {
	fprintf(stderr, "Fixing DOWNWARD sky edges\n") ;

	// eliminate long vertical edges in the sky end, which can create
	// visible edges in the feathering result
	// but only correct "downward" slivers
	for(x = 0 ; x < w-1 ; x++) {
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
    for(x = 0 ; x < w-1 ; x++) {
	if(pData->column_mask[x] == 1) continue ;
	float py = 1.-(float)pData->end_of_sky[x]/(float)(h-1) ;

	if(pData->lowest_sky_py_found > py) {
	    pData->lowest_sky_py_found = py ;
	}
    }

/*      print_sobel(image,w/2) ;  */

    free(is_clear_sky) ;
    free(sky_hue) ;
    free(sky_sat) ;
    free(sky_val) ;
}
