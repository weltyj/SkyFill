#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>

#include "skyfill_tif.h"
#include "colorspace_conversions.h"
#include "pixel_tests.h"

int xy_is_opaque_pixel(tdata_t *image, int32_t x, int32_t y)
{
    /* if alpha channel is nonzero then hugin has placed image data here */
	if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+3] < HALF16)
	    return 0 ;

    return 1 ;
}

int xy_has_nonblack_pixel(tdata_t *image, int32_t x, int32_t y)
{
    /* if alpha channel is nonzero then hugin has placed image data here */
    if(IMAGE_HAS_ALPHA) {
	if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+3] < HALF16)
	    return 0 ;
    }

    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+0] > BLK16)
	return 1 ;
    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+1] > BLK16)
	return 1 ;
    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+2] > BLK16)
	return 1 ;

    return 0 ;
}

int wrong_hue_or_black(tdata_t *image, int32_t x, int32_t y, float hue_sky, float sat_sky, float val_sky)
{
    uint16_t r,g,b,a ;
    tif_get4c(image,x,y,r,g,b,a) ;

    /* if alpha channel is nonzero then hugin has placed image data here */
    if(a < HALF16) return 0 ;

    float h,s,v ;
    rgb2hsv16(r,g,b,&h,&s,&v) ;
    float e_h = fabs((h-hue_sky)/hue_sky) ;
    float e_s = fabs((s-sat_sky)/sat_sky) ;
    float e_v = fabs((v-val_sky)/val_sky) ;

    if(e_h > .05) return 1 ;
    if(e_s > .20) return 1 ;
    if(e_v > .20) return 1 ;

    return 0 ;
}

int xy_is_transparent(tdata_t *image, int32_t x, int32_t y)
{
    /* if alpha channel is nonzero then hugin has placed image data here */
    if(IMAGE_HAS_ALPHA) {
	if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+3] < 256*64)
	    return 1 ;
    }

    return 0 ;
}


int xy_is_black_pixel(tdata_t *image, int32_t x, int32_t y)
{
    /* if alpha channel is nonzero then hugin has placed image data here */
    if(IMAGE_HAS_ALPHA) {
	if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+3] < 256*64)
	    return 0 ;
    }

    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+0] > 256*15)
	return 0 ;
    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+1] > 256*15)
	return 0 ;
    if(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+2] > 256*15)
	return 0 ;

    return 1 ;
}
