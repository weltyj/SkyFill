#ifndef SKYFILL_TIF_H
#define SKYFILL_TIF_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "tiffio.h"

// in case Windows MSVC, or 
#define _USE_MATH_DEFINES
#include <math.h>

// A few truly global variables
extern int32_t IMAGE_HEIGHT ;
extern int32_t IMAGE_WIDTH ;
extern uint16_t IMAGE_NSAMPLES ; // must be 16 bit, needed by call to tif library
extern int IMAGE_HAS_ALPHA ;
extern float p_half_image_width ; // scaled value of half the image width (use as X origin) ;

#define SKY_PIXEL_TYPE 0x01
#define ABOVE_SKY_PIXEL_TYPE 0x02
#define NON_SKY_PIXEL_TYPE 0x04
#define SKY2_PIXEL_TYPE 0x08  // kmeans classifier, next likely sky pixel cluster

// is the sky CIE model also used?
extern int uses_CIE_model ;

// normalized 3D vector
struct V3D {
    float A,B,C ;
    } ;

typedef struct rect_extent { uint16_t l,r,t,b; } RECT_EXTENT_t ; // rectangular extent, left, right, top, bottom



extern struct V3D V_zenith;
extern struct V3D V_sun;

/* coordinate system:
	TIF files are read.  data is stored in image[] array.  Top left is (0,0), bottom right is (width-1,height-1) ;
*/

#define tif_b(image,x,y) (((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2])

// macros to get or set rgb and alpha 
#define tif_get3c(image,x,y,r,g,b) {\
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		}
#define tif_get3cv(image,x,y,array) {\
		array[0] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		array[1] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		array[2] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		}
#define tif_get4cv(image,x,y,array) {\
		array[0] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		array[1] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		array[2] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		array[3] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] ; \
		}
#define dtif_get3c(image,x,y,r,g,b) {\
		fprintf(stderr, "get image data at %d, %d\n", x, y) ; \
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		}
#define dtif_get4c(image,x,y,r,g,b,a) {\
		fprintf(stderr, "get image data at %d, %d\n", x, y) ; \
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		a = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] ; \
		}

#define tif_get4c(image,x,y,r,g,b,a) {\
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		a = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] ; \
		}

#define tif_get_alpha(image,x,y,a) {\
		a = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] ; \
		}

#define tif_set3c(image,x,y,r,g,b) {\
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] = (uint16_t)r; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] = (uint16_t)g; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] = (uint16_t)b; \
		}

#define tif_set4c(image,x,y,r,g,b,a) {\
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] = (uint16_t)r; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] = (uint16_t)g; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] = (uint16_t)b; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] = (uint16_t)a; \
		}

#define MAX16 65535
#define HALF16 32767
#define MAX16f 65535.f
#define FMAX 1.e30f
// value which is black
#define BLK16 8621

#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )

#define IMAGE_PIXEL_X_TO_RELATIVE(x) ((float)( ((float)(x)/(float)(IMAGE_WIDTH-1)) -p_half_image_width))
#define IMAGE_PIXEL_Y_TO_RELATIVE(y) ((float)(1.f-(float)(y)/(float)(IMAGE_HEIGHT-1)))

#define IMAGE_RELATIVE_TO_PIXEL_X(px) (int)(((px)+p_half_image_width)*(float)(IMAGE_WIDTH-1)+0.5)
#define IMAGE_RELATIVE_TO_PIXEL_Y(py) (int)((1.-(py))*(float)(IMAGE_HEIGHT-1)+0.5)

// PX_WGT is for experiments only!!, not actually useful for model
// heavy weight to samples in center of image
//#define PX_WGT(px) (1./(.001+px*px)) 
// heavy weight to samples at left edge of image
//#define PX_WGT(px) (pow((.5-px),8))

// normal use -- even weights for all samples
#define PX_WGT(px) (1.f)

void predict_sky_hsv(float px, float py, float *hsv_hat) ;
void predict_sky_huesat_from_val(float vhat, float *hsv_hat, float px, float py) ;
void predict_sky_modelled_huesat_from_val(float vhat_CIE_sun, float *hsv_hat, float px, float py) ;
void predict_sky_color(float px, float py, uint16_t *rgb_hat) ;
float py_adjusted_for_horizon_curvature(float px, float py) ;
void dump_parameters_to_stdout(void) ;
struct V3D image_relative_pixel_to_V3(float px, float py) ;
struct V3D image_relative_pixel_to_V3_no_horizon(float px, float py) ;
struct V3D image_relative_pixel_to_V3_noclip(float px, float py) ;
float sun_x_angle2px(float sun_x_angle) ;

float find_maximum_CIE_vhat(void) ;
void set_FOV_factor() ;
float F_CIE2003(float px, float py, float *pGamma, float *pTheta, float *pCos_gamma) ;
float F_CIE2003_sun_only(float px, float py, float *pGamma, float *pTheta, float *pCos_gamma) ;

void predict_sky_h(float px, float py_hc, float *pH) ;
void predict_sky_s(float px, float py_hc, float *pS) ;
void predict_sky_v(float px, float py_hc, float *pV) ;
struct V3D P3toV3(float x, float y, float z) ;

// half space coefficients for 2D
struct HS {
    float A, B, C ;
} ;


// construct half space equation orthoganal to two points
// x2,y2 is on the line, x1,y1 is a postive distance from the line
struct HS half_space_from_pts(float x1, float y1, float x2, float y2) ;

// more than 100 and this is probably not a clear enough sky for this app to work
#define MAX_TEST_MASKS 100

typedef struct skyfill_data 
    {
    RECT_EXTENT_t test_extent_mask[MAX_TEST_MASKS] ; // skip over these extents for searching end of sky, or repairing bad pixels in the repair_top_of_sky steps
    int n_test_masks ;

    float floatBLK ; // value at which we assume hugin has masked out a portion of the image and placed black there instead of a alpha of 0 ;
    int verbose ;
    int fix_SOS_edges ;  // by default, don't fill missing data on edges around the start of the sky
    int fill_top_of_sky ;  // by default, don't make start of the sky values equal in first stage of repairing the top of sky
    int full_sky_replacement ; // set to 1, will attempt to replace everything down to end of sky, needs df set to 1
    float full_sky_replacement_thresh ; // threhold at which probability of sky pixel causes replacement of pixels
    float full_sky_replacement_ramp ; // speed at which probability of sky pixel changes
    float final_hsv_factor[3] ;
    float reduce_spline_value_range ; // in splined mode, reduce value range in spine curve. 0 is original curve, 1 is completely flat


    int horizon_was_set ;
    int sat_prediction_method ; // flag to change saturation prediction method for sky
    int full_hsv_model_level ;
    int val_model_full ; // set to 1, will use quadratic term in py for predicting sky value 
    int noclip_end_of_sky ; // set to 1, in full sky replacement mode, will replace pixels all the way to end of sky
    int repair_sky_slope_correct ;  // normally 1, set to 0 turnss of slope correction at the top of raw sky

    // global arrays, will indicate what was detected for each column (x) for these values
    int16_t *start_of_sky ;
    int16_t *end_of_sky ;
    int16_t *final_end_of_sky ;
    int16_t *manual_end_of_sky ; // loaded from file "inputfile".eos, if found
    int16_t *raw_start_of_sky ;  // raw start of sky is as detected before any repairs on image
    uint8_t **pixel_class ; // a 2d array with bit masks that classify the pixels

    // set to 1 to grid search for starting sun and horizon position
    int estimate_only ;
    int show_raw_prediction ; // old model, shows the sky predicted values above the unchanged original image
    int show_sky_prediction ; // shows probability (greyscale) the pixel is a sky pixel
    int show_raw_error ;
    int replace_only_transparent ;  // on final pass, only alter completely transparent (alpha==0), pixels ;
    int CIE_sky_index ; // if 3, will only use CIE coefs to predict given CIE index grid search
    int allowed_sky_type ;
    uint32_t n_estimate_sky_clipped ; // number of pixels at MAX16 during output of estimated pixel values
    int final_box_filter ; // run a box filter over sky area just before dithering

    int fix_sky_hue ;  // try to repair blown out sky areas
    int min_sky_hue_mask ;
    int max_sky_hue_mask ;
    uint32_t n_sky_hue_fixed_clipped ; // number of pixels at MAX16 during output of repaired sky hue pixel values
    int depth_of_fill_is_absolute ;
    int depth_of_fill_absolute_y ;
    float horizon_curvature ; // a proportion, up or down, at the edge of the image relative to the center of the image
    float horizon_slope ; // a factor applied to px to adjust py

    int end_of_sky_method ;


    // defaults for detecting a change from sky pixels to end of sky pixels
    float sky_hue_tolerance ;
    float sky_sat_tolerance ;
    float sky_val_tolerance ;

    // Jan 30 2022, new sky detection logic
    float sky_min_err_squared[3] ; // R,G,B
    float sky_ratio_threshold ;  // if individual pixel predicted (err**2/mean_err**2) > sky_ratio_threshold, it is classified as non sky


    float min_sky_end_p ; // the minimum allowed for end of sky, 0 is bottom of image, 1 is top
    float lowest_sky_py_found ;  // after sky is detected, this is set
    float repair_tos_thresh ;  // how far off a r,g, or b pixel is off from mean before it is not used in local
				    // sky color calculation -- the get_mean_rgb() function.

    int model_is_being_fit ;
    int val_model_is_linear ; // flag for which model is used to predict value
    float sample_based_v_correction ; // in the CIE_MODEL_FULL, factor applied to predicted value.

    float exposure_factor ; // exposure factor

    int have_pto_fov ;

    // horizon_py is at x=0 ;
    float sun_x, sun_py, horizon_py ;
    float hue_horizon, sat_depth, sun_lum ;

    // as determined from the sample points
    float hue_sky ;
    float val_sky ;
    float sat_sky ;

    float horizon_lift_angle_radians ;

    // not currently used, but leave in as comment for now
    /*  float horizon_slope=.05 ;  */

    float FOV_horizontal ;
    float FOV_vertical ;
    float proportion_to_radian_factor_x ; // how many radians per proportion of full image in X dimension ?
    float proportion_to_radian_factor_y ; // how many radians per proportion of full image in Y dimension ?
    float vertical_asin_angle_factor ; // = sin(FOV_vert_rad/2.)/0.5, for computing vertical angle relative to image coordianates
    float maximum_CIE_vhat ; // maximum vhat in sky dome
    float minimum_CIE_vhat ; // minimum vhat in sky dome

    float angle_factor ;  // in HSV sky model, this is 1. to 2. => changes the angular rate on the horizontal
    float valsat_phase_shift ;

    int max_end_of_sky ;
    int min_end_of_sky ;
    int max_start_of_sky ;
    int min_start_of_sky ;

    // 2021, try nonlinear feathering
    float nonlinear_feather ;

    // if > 0., will cause all colors to eventually be a constant at y==0 (aka the zenith) during final estimate sky coloring
    float zenith_blend_depth ;
    float zenith_blend_end_factor ;

    // These are the parameters for the CIE sky model
    // for clear blue sky with some haze this will produce
    // a result with high luminosity around the sun corona
    // as well as some lightening at the horizon
    // these will all be reset anyway by the time they are used
    float perez_A ; // horizon-zenith gradient, -5 to 5
    float perez_B ; // gradient intensity, -10 to 0
    float perez_C ; // circumsolar intensity, 0 to 25
    float perez_D ; // circumsolar radius, -10 to 0
    float perez_E ; // backscattering effect, -1 to 5
    float perez_F ; // width of solar disc in reduced model (Jeff Welty)
    float perez_G ; // width of solar disc in reduced model (Jeff Welty)

    unsigned char *column_mask ; // user requested masked columns for any alteration at all
    unsigned char *column_repair_mask ; // user requested masked columns for repair of bad pixels
    unsigned char *column_sample_mask ; // user requested masked columns for sampling of sky pixels
    unsigned char *column_nonsky_repair_mask ; // user requested masked columns for repair of pixels identified as nonsky

    int output_sample_data ;

} SKYFILL_DATA_t ;

// model lens flare
int lens_flare(char *filename, tdata_t *image, int sun_x_pixel, int sun_y_pixel)  ;

// one more global
extern SKYFILL_DATA_t *pData_fit ;

static inline float px_to_AOY(float px) {
    return -px*pData_fit->proportion_to_radian_factor_x ;
}

static inline float py_to_AOX(float py) {
    float l = (py-0.5)*pData_fit->vertical_asin_angle_factor ;
    if(l > 1.) {
	l=2.-l ; // py went up to the zenith, and over to other side of sky, mirror the angle
	return M_PI - asinf(l) ;
    } else {
	return asinf(l) ;
    }
    //return py*pData_fit->proportion_to_radian_factor_y ;
}

static inline uint16_t is_in_test_mask(int x, int y)
{

    for(int i=0 ; i< pData_fit->n_test_masks ; i++) {

	if( x >= pData_fit->test_extent_mask[i].l && x <= pData_fit->test_extent_mask[i].r
	 && y >= pData_fit->test_extent_mask[i].t && y <= pData_fit->test_extent_mask[i].b) {
	    return 1 ;
	}
    }

    return 0 ;
}
static inline uint16_t max_test_mask(int x, int y)
{

    if( pData_fit->column_mask[x] == 1)
	return IMAGE_HEIGHT ;

    int ymax=y ;
    for(int i=0 ; i< pData_fit->n_test_masks ; i++) {

	if( x >= pData_fit->test_extent_mask[i].l && x <= pData_fit->test_extent_mask[i].r
	 && y >= pData_fit->test_extent_mask[i].t && y <= pData_fit->test_extent_mask[i].b) {

	    // x,y is in the extent, is ymax less than the bottom of the extent?
	    if(ymax < pData_fit->test_extent_mask[i].b)
		ymax = pData_fit->test_extent_mask[i].b ;
	}
    }

    // returns y if pixel is not in extent, bottom of lowest bounding extent if in an extent
    return ymax ;
}

// model types for sky value
#define CIE_NO_MODEL 0x00
#define CIE_SUN_MODEL 0x01
#define CIE_FULL_MODEL 0x02


#define MAX_OPT_PARMS 18

// restrict CIE sky types used in optimization
#define ALL_CIE_SKIES 0x01
#define UNIFORM_CIE_SKIES 0x02
#define NONUNIFORM_CIE_SKIES 0x03

// these are the "standard" CIE parameters for various types of skies.
#define N_CIE_STD_PARMS 16

struct OPT_PARM {
    float *variable_address ;
    char abreviation[8] ; // parameter name abreviation
    char name[80] ; // parameter name
    float default_value ; // parameter value
    float lo ; // parameter minimum allowed value
    float hi ; // parameter maximum allowed value
    int optimize_flag ;  // if 1, include in optimization, if 0 do not include ;
    int grid_optimize_flag ;  // if 1, include in grid search optimization, if 0 do not include ;
    char *units ;
    int used ;
    int is_CIE ;
} ;

struct CIE_STD_PARMS {
    int type ;
    int gradation ;
    float A ;
    float B ;
    float C ;
    float D ;
    float E ;
    char description[120] ;
    } ;

#endif
