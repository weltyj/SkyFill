#ifndef ESTIMATE_SKY_H
#define ESTIMATE_SKY_H
void estimate_sky(int x0,int x1,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky,int16_t *final_end_of_sky,
		    float depth_of_sample, int h, int w, float feather_factor, int debug, float depth_of_fill, int extra_feather_length, SKYFILL_DATA_t *pData) ;
#endif
