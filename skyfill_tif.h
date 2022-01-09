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

// four macros to get or set rgb and alpha 
#define tif_get3c(image,x,y,r,g,b) {\
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		}

#define tif_get4c(image,x,y,r,g,b,a) {\
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
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

#define IMAGE_RELATIVE_TO_PIXEL_X(px) (int)(((px)+p_half_image_width)*(float)(IMAGE_WIDTH-1)+0.5) ;
#define IMAGE_RELATIVE_TO_PIXEL_Y(py) (int)((1.-(py))*(float)(IMAGE_HEIGHT-1)+0.5) ;

// PX_WGT is for experiments only!!, not actually useful for model
// heavy weight to samples in center of image
//#define PX_WGT(px) (1./(.001+px*px)) 
// heavy weight to samples at left edge of image
//#define PX_WGT(px) (pow((.5-px),8))

// normal use -- even weights for all samples
#define PX_WGT(px) (1.f)

void predict_sky_hsv(float px, float py, float *pH, float *pS, float *pV) ;
void predict_sky_huesat_from_val(float vhat, float *pH, float *pS, float *pV, float px, float py) ;
void predict_sky_color(float px, float py, uint16_t *rhat, uint16_t *ghat, uint16_t*bhat) ;
float py_adjusted_for_horizon_curvature(float px, float py) ;
void dump_parameters_to_stdout(void) ;
struct V3D image_relative_pixel_to_V3(float px, float py) ;
struct V3D image_relative_pixel_to_V3_noclip(float px, float py) ;
float sun_x_angle2px(float sun_x_angle) ;
float find_maximum_CIE_vhat(void) ;
void set_FOV_factor() ;

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
    float final_saturation_factor ;
    int horizon_was_set ;
    int sat_prediction_method ; // flag to change saturation prediction method for sky
    int full_hsv_model_level ;
    int val_model_full ; // set to 1, will use quadratic term in py for predicting sky value 
    int noclip_end_of_sky ; // set to 1, in full sky replacement mode, will replace pixels all the way to end of sky

    // global arrays, will indicate what was detected for each column (x) for these values
    int16_t *start_of_sky ;
    int16_t *end_of_sky ;
    int16_t *final_end_of_sky ;
    int16_t *manual_end_of_sky ; // loaded from file "inputfile".eos, if found
    int16_t *raw_start_of_sky ;  // raw start of sky is as detected before any repairs on image

    // set to 1 to grid search for starting sun and horizon position
    int estimate_only ;
    int show_raw_prediction ; // old model, shows the sky predicted values above the unchanged original image
    int show_sky_prediction ; // shows probability (greyscale) the pixel is a sky pixel
    int show_raw_error ;
    int CIE_sky_index ; // if 3, will only use CIE coefs to predict given CIE index grid search
    int allowed_sky_type ;

    int fix_sky_hue ;  // try to repair blown out sky areas
    int min_sky_hue_mask ;
    int max_sky_hue_mask ;
    int depth_of_fill_is_absolute ;
    int depth_of_fill_absolute_y ;
    float horizon_curvature ; // a proportion, up or down, at the edge of the image relative to the center of the image


    // defaults for detecting a change from sky pixels to end of sky pixels
    float sky_hue_tolerance ;
    float sky_sat_tolerance ;
    float sky_val_tolerance ;
    float min_sky_end_p ; // the minimum allowed for end of sky, 0 is bottom of image, 1 is top
    float lowest_sky_py_found ;  // after sky is detected, this is set
    float repair_tos_thresh ;  // how far off a r,g, or b pixel is off from mean before it is not used in local
				    // sky color calculation -- the get_mean_rgb() function.

    int model_is_being_fit ;

    float exposure_factor ; // exposure factor

    int have_pto_fov ;

    // horizon_py is at x=0 ;
    float sun_x, sun_py, horizon_py ;
    float hue_horizon, sat_depth, sky_lum ;

    // as determined from the sample points
    float hue_sky ;
    float val_sky ;
    float sat_sky ;

    float horizon_lift_angle_radians ;

    // not currently used, but leave in as comment for now
    /*  float horizon_slope=.05 ;  */

    float FOV_horizontal ;
    float proportion_to_radian_factor_x ; // how many radians per proportion of full image in X dimension ?
    float proportion_to_radian_factor_y ; // how many radians per proportion of full image in Y dimension ?
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

} SKYFILL_DATA_t ;

// one more global
extern SKYFILL_DATA_t *pData_fit ;

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


#define MAX_OPT_PARMS 17

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
