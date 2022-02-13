#ifndef REPAIR_SKY_H
#define REPAIR_SKY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "tiffio.h"

#include "skyfill_tif.h"
#include "mstat.h"
#include "colorspace_conversions.h"

#define MAX_PHASES 10
void repair_sky_hue(tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky, SKYFILL_DATA_t *pData) ;
void repair_alpha(tdata_t *image, SKYFILL_DATA_t *pData) ;
void set_minmax_sky_values(SKYFILL_DATA_t *pData) ;
int get_mean_rgb(tdata_t *image, int xc,int yc,double rcoef[],double gcoef[], double bcoef[],int search_width, int clip_eos_flag, int repair_bad_points_flag, int recursion_level, SKYFILL_DATA_t *pData) ;
void repair_top_of_sky(tdata_t *image, int end_of_sky_is_known_flag, int phase, SKYFILL_DATA_t *pData) ;
#endif
