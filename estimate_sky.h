#ifndef ESTIMATE_SKY_H
#define ESTIMATE_SKY_H

struct sky_pixel_data {
    float hsv_hat[3] ;
    uint16_t rgb_hat[3] ;
    float rgb_err[3] ;
    float hsv_of_e[3] ;
    float hsv[3] ;
    float max_err ;
    float p_hat ;
    float hsv_err ;
    float prob_sky ;
    uint16_t rgb[3], a ;
} ;
float get_prob_sky(tdata_t *image, int x,int y, SKYFILL_DATA_t *pData, struct sky_pixel_data *p) ;

void estimate_sky(int x0,int x1,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky,int16_t *final_end_of_sky,
		    float depth_of_sample, int h, int w, float feather_factor, int debug, float depth_of_fill, int extra_feather_length, SKYFILL_DATA_t *pData) ;
#endif
