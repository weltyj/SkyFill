#ifndef FIND_SKY_H
#define FIND_SKY_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "tiffio.h"

#include "skyfill_tif.h"

void simple_find_start_of_sky(int16_t *start_of_sky,uint16_t w,uint16_t h,tdata_t *image, SKYFILL_DATA_t *pData) ;
void find_start_of_sky(int16_t *start_of_sky,uint16_t w,uint16_t h,tdata_t *image,int fix_SOS_edges, SKYFILL_DATA_t *pData) ;
void get_sky_mean_var(tdata_t *image, int x, int y0, int n, float mean[], float var[]) ;
int is_end_of_sky(uint16_t y, int n_samples, float *hdiff, float *sdiff, float *vdiff, SKYFILL_DATA_t *pData) ;
float sobel(tdata_t *image, int xc, int yc) ;
void print_sobel(tdata_t *image, int x, SKYFILL_DATA_t *pData) ;
int get_sobel_eos(tdata_t *image, int x, SKYFILL_DATA_t *pData) ;
void find_end_of_sky(int16_t *end_of_sky,int16_t *start_of_sky,int w,int h,tdata_t *image,int fix_edges,int fix_slivers, SKYFILL_DATA_t *pData) ;
#endif
