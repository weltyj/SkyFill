#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>

#include "skyfill_tif.h"
#include "feather_factor.h"

// how much of the new estimate to use at pixel x, y
float raw_compute_feather_at_y(int x, int y, int feather_length, uint16_t b, float feather_factor, SKYFILL_DATA_t *pData)
{
    //if(x == 0) {
	//fprintf(stderr, "INFO: CFAY 1, x:%d y:%d sos:%d eos:%d fl:%d ff:%f\n", x, y, raw_start_of_sky[x], end_of_sky[x], feather_length, feather_factor) ;
    //} 

    // compute feathering amount
    float fp0 = 1. ; // proportion of prediction to use

    if(y <= pData->raw_start_of_sky[x]) return 1. ;
    if(y > pData->raw_start_of_sky[x]+feather_length) return 0. ;

    if(y >= pData->raw_start_of_sky[x])
	fp0 = 1.-((float)y-(float)pData->raw_start_of_sky[x])/(float)feather_length ;

    if(b < pData->floatBLK*MAX16) {
	// pixel is black, use 100% of predicted value
	fp0 = 1.0 ;
    }

    float fp1 = 1.-fp0 ;
    fp1 *= feather_factor ; // adjust amount of feathering by feather_factor (useful for debugging)
    fp0 = 1.-fp1 ;

    if(fp0 < 0.0) {
	fprintf(stderr, "FATAL: raw CFAY 1, fp0:%f x:%d y:%d sos:%d fl:%d ff:%f\n", fp0, x, y, pData->raw_start_of_sky[x], feather_length, feather_factor) ;
	exit(1) ;
    }

    if(y == pData->end_of_sky[x]) {
	// make the end of sky pixel half estimated
/*  	fp0 = .5 ;  */
/*  	fp1 = .5 ;  */
    }

    if(fp0 < 0.0) {
	fprintf(stderr, "FATAL: raw CFAY 2, fp0:%f x:%d y:%d sos:%d fl:%d ff:%f\n", fp0, x, y, pData->raw_start_of_sky[x], feather_length, feather_factor) ;
	exit(1) ;
    }

    fp0 = powf(fp0,pData->nonlinear_feather) ;
    return fp0 ;
}

int raw_compute_feather_length(int x, int *pFeather_end_y, float depth_of_fill, float feather_factor, int extra, SKYFILL_DATA_t *pData)
{
    float sky_height ;

    if(pData->depth_of_fill_is_absolute) {
	sky_height = pData->depth_of_fill_absolute_y - pData->raw_start_of_sky[x] ;
    } else {
	sky_height = pData->end_of_sky[x] - pData->start_of_sky[x] ;
    }

    *pFeather_end_y = pData->raw_start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) + extra ;

    if(*pFeather_end_y > pData->end_of_sky[x]) *pFeather_end_y = pData->end_of_sky[x] ;

    if(pData->estimate_only) {
	*pFeather_end_y = IMAGE_HEIGHT-1 ;
	feather_factor=0. ;
    }

    return *pFeather_end_y - pData->raw_start_of_sky[x] ;
}

// how much of the new estimate to use at pixel x, y
float compute_feather_at_y(int x, int y, int feather_length, uint16_t pixel_blue_value, float feather_factor, SKYFILL_DATA_t *pData)
{
    //if(x == 0) {
	//fprintf(stderr, "INFO: CFAY 1, x:%d y:%d sos:%d eos:%d fl:%d ff:%f\n", x, y, start_of_sky[x], end_of_sky[x], feather_length, feather_factor) ;
    //} 

    // compute feathering amount
    float fp0 = 1. ; // proportion of prediction to use

    if(y >= pData->start_of_sky[x])
	fp0 = 1.-((float)y-(float)pData->start_of_sky[x])/(float)feather_length ;

    if(pixel_blue_value < pData->floatBLK*MAX16) {
	// pixel is black, use 100% of predicted value
	fp0 = 1.0 ;
    }

    float fp1 = 1.-fp0 ;
    fp1 *= feather_factor ; // adjust amount of feathering by feather_factor (useful for debugging)
    fp0 = 1.-fp1 ;

    if(fp0 < 0.0) {
	fprintf(stderr, "FATAL: CFAY 1, fp0:%f x:%d y:%d sos:%d fl:%d ff:%f\n", fp0, x, y, pData->start_of_sky[x], feather_length, feather_factor) ;
	exit(1) ;
    }

    if(y == pData->end_of_sky[x]) {
	// make the end of sky pixel half estimated
/*  	fp0 = .5 ;  */
/*  	fp1 = .5 ;  */
    }

    if(fp0 < 0.0) {
	fprintf(stderr, "FATAL: CFAY 2, fp0:%f x:%d y:%d sos:%d fl:%d ff:%f\n", fp0, x, y, pData->start_of_sky[x], feather_length, feather_factor) ;
	exit(1) ;
    }

    fp0 = powf(fp0,pData->nonlinear_feather) ;
    return fp0 ;
}

int compute_feather_length_with_eos(int x, int *pFeather_end_y, float depth_of_fill, float feather_factor, int extra, SKYFILL_DATA_t *pData, int end_of_sky)
{
    float sky_height ;

    if(pData->depth_of_fill_is_absolute) {
	sky_height = pData->depth_of_fill_absolute_y - pData->start_of_sky[x] ;
    } else {
	sky_height = end_of_sky - pData->start_of_sky[x] ;
    }

    *pFeather_end_y = pData->start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) + extra ;

    if(*pFeather_end_y > end_of_sky) *pFeather_end_y = end_of_sky ;

    if(pData->estimate_only) {
	*pFeather_end_y = IMAGE_HEIGHT-1 ;
	feather_factor=0. ;
    }

    return *pFeather_end_y - pData->start_of_sky[x] ;
}

int compute_feather_length(int x, int *pFeather_end_y, float depth_of_fill, float feather_factor, int extra, SKYFILL_DATA_t *pData)
{
    return compute_feather_length_with_eos(x, pFeather_end_y, depth_of_fill, feather_factor, extra, pData, pData->end_of_sky[x]) ;
}
