#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "tiffio.h"


// in case Windows MSVC, or 
#define _USE_MATH_DEFINES
#include <math.h>

#include "skyfill_tif.h"
#include "amoeba_06.h"
#include "string.h"
#include "mstat.h"
#include "find_sky.h"
#include "repair_sky.h"
#include "colorspace_conversions.h"
#include "sample_and_fit_sky_model.h"
#include "feather_factor.h"
#include "estimate_sky.h"
#include "optimize.h"

// A few truly global variables
int32_t IMAGE_HEIGHT ;
int32_t IMAGE_WIDTH ;
uint16_t IMAGE_NSAMPLES ; // must be 16 bit, needed by call to tif library
int IMAGE_HAS_ALPHA ;
float p_half_image_width ; // scaled value of half the image width (use as X origin) ;

// is the sky CIE model also used?
int uses_CIE_model ;

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

// a globally available structure with data for the sky model.  For now -- would be good to eliminate this in the future
SKYFILL_DATA_t *pData_fit ;

void initialize_skyfill_data_struct(SKYFILL_DATA_t *pData)
{

    pData->n_test_masks=0 ;
    pData->floatBLK=.04 ; // value at which we assume hugin has masked out a portion of the image and placed black there instead of a alpha of 0 ;
    pData->verbose=0 ;
    pData->fix_SOS_edges=0 ;  // by default, don't fill missing data on edges around the start of the sky
    pData->fill_top_of_sky=0 ;  // by default, don't make start of the sky values equal in first stage of repairing the top of sky
    pData->full_sky_replacement=0 ; // set to 1, will attempt to replace everything down to end of sky, needs df set to 1
    pData->full_sky_replacement_thresh=.9 ; // threhold at which probability of sky pixel causes replacement of pixels
    pData->full_sky_replacement_ramp=200. ;
    pData->final_saturation_factor=1.0 ;
    pData->horizon_was_set=0 ;
    pData->sat_prediction_method=1 ;
    pData->val_model_full = 1 ; // set to one, will use quadratic term in py for predicting sky value 

    pData->estimate_only=0 ;
    pData->show_raw_prediction=0 ; // old model
    pData->show_raw_error=0 ;
    pData->show_sky_prediction=0 ;
    pData->CIE_sky_index=-1 ; // if 3, will only use CIE coefs to predict given CIE index grid search
    pData->allowed_sky_type = ALL_CIE_SKIES ;

    pData->fix_sky_hue =0 ;  // try to repair blown out sky areas
    pData->min_sky_hue_mask = 0 ;
    pData->max_sky_hue_mask = -1 ;
    pData->depth_of_fill_is_absolute = 0 ;
    pData->horizon_curvature=0. ; // a proportion, up or down, at the edge of the image relative to the center of the image


    // defaults for detecting a change from sky pixels to end of sky pixels
    pData->sky_hue_tolerance=360. ;
    pData->sky_sat_tolerance=.3 ;
    pData->sky_val_tolerance=.015 ;
    pData->min_sky_end_p=1.0 ; // the minimum allowed for end of sky, 0 is bottom of image, 1 is top
    pData->lowest_sky_py_found=1.0 ;  // after sky is detected, this is set
    pData->repair_tos_thresh=0.02 ;  // how far off a r,g, or b pixel is off from mean before it is not used in local
    // sky color calculation -- the get_mean_rgb() function.

    pData->model_is_being_fit=1 ;
    pData->full_hsv_model_level=1 ;

    pData->exposure_factor=1. ; // exposure factor

    pData->have_pto_fov=0 ;

    pData->sun_x = 0. ;
    pData->sun_py = 1.5 ;
    pData->horizon_py = .5 ;
    pData->hue_horizon = 200 ; // not used
    pData->sat_depth = 1.0 ;
    pData->sky_lum=.9 ;

    // as determined from the sample points
    pData->hue_sky = 212 ;
    pData->val_sky=1. ;
    pData->sat_sky=.7 ;

    pData->horizon_lift_angle_radians=0. ;

    // not currently used, but leave in as comment for now
    /*  float horizon_slope=.05 ;  */

    pData->FOV_horizontal=120. ;
    pData->maximum_CIE_vhat = 1. ; // maximum vhat in sky dome
    pData->minimum_CIE_vhat = 0. ; // minimum vhat in sky dome

    pData->angle_factor = 1. ;  // in HSV sky model, this is 1. or 2. => changes the angular rate on the horizontal
    pData->valsat_phase_shift = -1000. ; // the value and saturations models will apply a phase shift.  This is an invalid value which will cause it to be automatically detected

    pData->max_end_of_sky = -1 ;
    pData->min_end_of_sky = -1 ;
    pData->max_start_of_sky = -1 ;
    pData->min_start_of_sky = -1 ;

    // 2021, try nonlinear feathering
    pData->nonlinear_feather=1.0 ;

    // These are the parameters for the CIE sky model
    // for clear blue sky with some haze this will produce
    // a result with high luminosity around the sun corona
    // as well as some lightening at the horizon
    // these will all be reset anyway by the time they are used
    pData->perez_A=-1.5 ; // horizon-zenith gradient, -5 to 5
    pData->perez_B=-0.80 ; // gradient intensity, -10 to 0
    pData->perez_C=2. ; // circumsolar intensity, 0 to 25
    pData->perez_D=-2.4 ; // circumsolar radius, -10 to 0
    pData->perez_E=-0.15 ; // backscattering effect, -1 to 5
    pData->perez_F=0.5 ; // blending amount of sun model into sky hsv model -- Jeff Welty
    pData->perez_G=10. ; // falloff factor of blending amount -- Jeff Welty
}


struct OPT_PARM opt_parms[MAX_OPT_PARMS] ;

void fill_opt_parms(SKYFILL_DATA_t *p)
{
    struct OPT_PARM tmp_parms[MAX_OPT_PARMS] =  {
	{&p->sun_x, "sx", "sun_x", 0., -180., 180., 1, 1, "degrees: 0 is middle of image", 0, 1}, // 0 is straight ahead, -180 or 180 is straight behind
	{&p->sun_py, "sy", "sun_py", .99, 0.0, 3., 1, 1, "proportion: 0 is bottom, 1 is top", 0, 1}, // 0 is bottom of image, 1 is top of image
	{&p->horizon_py, "hy", "horizon_py", .5, -.2, 1.5, 1, 1, "proportion: 0 is bottom, 1 is top"}, // 0 is bottom of image, 1 is top of image
	{&p->horizon_curvature, "hc", "horizon curvature", 0.0, -.2,   .2, 1, 1, "proportion", 0, 0},
	{&p->hue_sky, "hs", "hue sky", 212., 180., 250., 1, 0, "Hue from HSV model", 0, 0},
	{&p->hue_horizon, "hh", "hue_horizon", 212., 180., 250., 1, 0, "Horizon hue, not currenly used", 0, 0},
	{&p->sat_depth, "sd", "sat_depth", .8, 0.01, 1.1, 1, 1, "Amount of saturation to apply, 0 to 1", 0, 0},
	{&p->sky_lum, "sl", "sky_lum", 1., 0.01, 20.1, 1, 1, "Factor applied to predicted sky intensity", 0, 0},
	{&p->FOV_horizontal, "fov", "FOV_horizontal", 90., 1., 360., 1, 0, "The horizontal field of view of the image, degrees", 0, 0},
	{&p->perez_A, "A", "perez_A", .001, -5., 5., 1, 0, "The A parameter of the perez CIE sky model", 0, 1},
	{&p->perez_B, "B", "perez_B", -1., -10., -0.001, 1, 0, "The B parameter of the perez CIE sky model", 0, 1},
	{&p->perez_C, "C", "perez_C", 2., 0., 25., 1, 0, "The C parameter of the perez CIE sky model", 0, 1},
	{&p->perez_D, "D", "perez_D", -1.5, 0.0001, 0.1, 1, 0, "D in sun model, aprox degrees/100.", 0, 1},
	{&p->perez_E, "E", "perez_E", 0.15, -1., 5., 1, 0, "The E parameter of the perez CIE sky model", 0, 1},
	{&p->perez_F, "F", "perez_F", 1.0, .01, 1., 1, 0, "The F parameter of the reduced perez CIE sky model", 0, 1},
	{&p->perez_G, "G", "perez_G", 1.0, .1, 100., 1, 0, "The G parameter of the reduced perez CIE sky model", 0, 1},
	{&p->exposure_factor, "ef", "exposure factor",    1.0,  0.1, 10.0, 1, 1, "Final exposure appled to sky intensity", 0, 0},
    } ;

    for(int i = 0 ; i < MAX_OPT_PARMS ; i++) {
	opt_parms[i] = tmp_parms[i] ;
    }

}

struct V3D V_zenith;
struct V3D V_sun;

void dump_parameters_to_stdout()
{
    int i ;

    for(i=0 ; i < MAX_OPT_PARMS ; i++) {
	if(!strcmp("A",opt_parms[i].abreviation)) continue ; // A no longer used
	if(!strcmp("B",opt_parms[i].abreviation)) continue ; // B no longer used
	if(!strcmp("E",opt_parms[i].abreviation)) continue ; // E no longer used
	printf("%d:%12s %f\n", opt_parms[i].optimize_flag, opt_parms[i].name, *opt_parms[i].variable_address) ;
    }

    printf("Maximum CIE vhat=%f\n", pData_fit->maximum_CIE_vhat) ;
    printf("Minimum CIE vhat=%f\n", pData_fit->minimum_CIE_vhat) ;

    printf("command line settings:\n  ") ;

    for(i=0 ; i < MAX_OPT_PARMS ; i++) {
	if(!strcmp("A",opt_parms[i].abreviation)) continue ; // A no longer used
	if(!strcmp("B",opt_parms[i].abreviation)) continue ; // B no longer used
	if(!strcmp("E",opt_parms[i].abreviation)) continue ; // E no longer used
	printf(" -%s %f", opt_parms[i].abreviation, *opt_parms[i].variable_address) ;
    }
    printf("\n") ;
}

// given sun_x (as an angle from the center of the image), convert
// to px, where px=0 is left side of image, px=1. is right side of image
float sun_x_angle2px(float sun_x_angle)
{
    return sun_x_angle/pData_fit->FOV_horizontal ;
}


void read_tif_image(tdata_t *image,TIFF *tif,int h,uint16_t spp, uint16_t tif_config)
{
    int32_t y ;
    /* read a row */
    for(y = 0 ; y < h ; y++) {
	if(tif_config == PLANARCONFIG_CONTIG) {
	    TIFFReadScanline(tif,image[y],y,0) ;
	} else if(tif_config == PLANARCONFIG_SEPARATE) {
	    uint16_t s ;

	    for(s = 0 ; s < spp ; s++)
		TIFFReadScanline(tif,image[y],y,s) ;

	} else {
	    fprintf(stderr, "Cannot read planar configuration of tiff image, config=%d!\n", (int)tif_config) ;
	    exit(1) ;
	}
    }
}

void copy_pixel(tdata_t *image, int32_t xsrc, int32_t ysrc, int32_t xdest, int32_t ydest)
{
    ((uint16_t *)(image[ydest]))[IMAGE_NSAMPLES*xdest+0] = ((uint16_t *)(image[ysrc]))[IMAGE_NSAMPLES*xsrc+0] ;
    ((uint16_t *)(image[ydest]))[IMAGE_NSAMPLES*xdest+1] = ((uint16_t *)(image[ysrc]))[IMAGE_NSAMPLES*xsrc+1] ;
    ((uint16_t *)(image[ydest]))[IMAGE_NSAMPLES*xdest+2] = ((uint16_t *)(image[ysrc]))[IMAGE_NSAMPLES*xsrc+2] ;

    if(IMAGE_HAS_ALPHA) {
	((uint16_t *)(image[ydest]))[IMAGE_NSAMPLES*xdest+3] = ((uint16_t *)(image[ysrc]))[IMAGE_NSAMPLES*xsrc+3] ;
    }
}


float py_adjusted_for_horizon_curvature(float px, float py)
{

#define OLD_HC
#ifdef OLD_HC
    float tmp = 1. - fabs(pData_fit->horizon_curvature) ;
    float f = 1. - tmp*tmp ;
    // *4 -> px is only -0.5 to 0.5, we want curvature for py as if px runs -1 to 1
    // but the user (or algorithm) considers the horizon curvature at px=0.5 
    tmp = 1. - f*px*px*4. ; 
    if(tmp < 0.) {
/*  	fprintf(stderr, "FATAL: pxpy_toV3D, negative adjusted px\n") ;  */
/*  	fprintf(stderr, "FATAL: px:%f\n", px) ;  */
/*  	fprintf(stderr, "FATAL: py:%f\n", py) ;  */
/*  	fprintf(stderr, "FATAL: horizon_curvature:%f\n", horizon_curvature) ;  */
/*  	fprintf(stderr, "FATAL: f:%f\n", f) ;  */
/*  	fprintf(stderr, "FATAL: tmp:%f\n", tmp) ;  */
	tmp = 0. ;
    }

    float hc_py = 1. - sqrt(tmp) ;
    if(pData_fit->horizon_curvature < 0.) hc_py = -hc_py ;

    // hc_py is where we expect the horizon to actually be at, in image coordinates
    // so if we find the distance between the image given py and hc_py, this will
    // correspond to a 'py' with curvature 0.  (with some obvious issues, but should
    // be good enough for very low curvatures

    return py-hc_py ;
#else

    // bug in this somewhere, causing NAN's in F_CIE2003()

    if(fabs(pData_fit->horizon_curvature) < .00001) {
	return py ;
    } else {
	float R = (pData_fit->horizon_curvature*pData_fit->horizon_curvature+.25)/(2.*pData_fit->horizon_curvature) ;
	R = fabs(R) ;
	float hc_py = R - sqrt(R*R-px*px) ;

	if(pData_fit->horizon_curvature > 0.) {
	    return py-hc_py ;
	} else {
	    return py+hc_py ;
	}
    }
#endif
}


// convert point to normalized vector from origin (0.,0.,0.) to point
struct V3D P3toV3(float x, float y, float z)
{
    struct V3D v ;
    float sum = sqrt(x*x + y*y + z*z) ;
    v.A = x/sum ;
    v.B = y/sum ;
    v.C = z/sum ;

    return v ;
}

void set_FOV_factor()
{
    float FOV_horz_rad = M_PI/180.0*pData_fit->FOV_horizontal ;
    float FOV_vert_rad = FOV_horz_rad * (float)IMAGE_HEIGHT/(float)IMAGE_WIDTH ;

    pData_fit->proportion_to_radian_factor_x = FOV_horz_rad ;
    pData_fit->proportion_to_radian_factor_y = FOV_vert_rad ;
}

float angle_between(struct V3D *v1, struct V3D *v2, float *dotproduct)
{
    *dotproduct = v1->A*v2->A + v1->B*v2->B + v1->C*v2->C ;

    // normally divide by magnitude of each vector
    // but these vectors have been already normalized so their magnitude is 1.0
    // dotproduct /= mag(v1)*mag(v2) ;
    if(*dotproduct > 1.00001f) {
	fprintf(stderr, "AB FATAL, dp(%f) > 1.  (%f,%f,%f) and (%f,%f,%f)\n",
	    *dotproduct,
	    v1->A,
	    v1->B,
	    v1->C,
	    v2->A,
	    v2->B,
	    v2->C) ;
	    exit(1) ;
    }
    if(*dotproduct < -1.00001f) {
	fprintf(stderr, "AB FATAL, dp(%f) > 1.  (%f,%f,%f) and (%f,%f,%f)\n",
	    *dotproduct,
	    v1->A,
	    v1->B,
	    v1->C,
	    v2->A,
	    v2->B,
	    v2->C) ;
	    exit(1) ;
    }
    if(*dotproduct > 1.f) *dotproduct=1.f ; // catch rounding error, acos(> 1.) is NAN
    if(*dotproduct < -1.f) *dotproduct=-1.f ; // catch rounding error, acos(> 1.) is NAN

    float angle = acos(*dotproduct) ;
    return angle ;
}

struct V3D pxpy_toV3D(float px, float py) {
    float py_above_horizon=py_adjusted_for_horizon_curvature(px,py)  ;

    if(0) {
	// for debugging
	pData_fit->horizon_curvature = .1 ;
	py=.1 ;
	for(px = -0.5 ; px < 0.51 ; px += .05) {
	    py_above_horizon=py_adjusted_for_horizon_curvature(px,py)  ;
	    fprintf(stderr, "px:%5.2f py:%5.2f py_above:%5.2f\n", px, py, py_above_horizon) ;
	}
	exit(1) ;
    }

    float angle_about_Y = -px*pData_fit->proportion_to_radian_factor_x ;
    float angle_about_X = py_above_horizon*pData_fit->proportion_to_radian_factor_y ;


    float X = 0. ;
    float Y = 0. ;
    float Z = -1. ;

    float C = cos(angle_about_X) ;
    float S = sin(angle_about_X) ;
    float X1 = X ;
    float Y1 = Y*C-Z*S ;
    float Z1 = Y*S+Z*C ;

    C = cos(angle_about_Y) ;
    S = sin(angle_about_Y) ;
    Y = Y1 ;
    X = X1*C+Z1*S ;
    Z = -X1*S+Z1*C ;

    if(pData_fit->verbose) {
	struct V3D view = P3toV3(X,Y,Z) ;
	float dp_sun, dp_zenith ;
	fprintf(stderr, "pxpy: %5.2f %5.2f, AOX:%7.2f AOY:%7.2f ", px, py, angle_about_X*180./M_PI, angle_about_Y*180./M_PI) ;

	float Gamma = angle_between(&view,&V_sun,&dp_sun) ;
	float Theta = angle_between(&view,&V_zenith,&dp_zenith) ;
	fprintf(stderr, "Asun:%7.2f Azenith:%7.2f ", 180./M_PI*Gamma, 180./M_PI*Theta) ;

	return view ;
    }

    return P3toV3(X,Y,Z) ;
}

// give a image relative pixel (where the bottom of the image is py=0, top is py=1, left is -p_half_image_width, right is p_half_image_width)
// create a 3d vector
struct V3D image_relative_pixel_to_V3(float px, float py)
{
    py -= pData_fit->horizon_py ;
    if(py < 0.) py = 0. ;

    return pxpy_toV3D(px,py) ;
}

struct V3D image_relative_pixel_to_V3_noclip(float px, float py)
{
    py -= pData_fit->horizon_py ;

    return pxpy_toV3D(px,py) ;
}

float F_CIE2003(float px, float py, float *pGamma, float *pTheta, float *pCos_gamma)
{
    /* assume observer is at origin */
    struct V3D view= image_relative_pixel_to_V3(px, py) ;
    float dp_zenith, dp_sun ; // dot product of vectors to view

    *pGamma = angle_between(&view,&V_sun,&dp_sun) ;

    *pCos_gamma = dp_sun ;

    // check if we are within the sun diameter
    if(*pGamma < M_PI/180.f*0.533f/2.) {
	*pGamma = 0.f ;
	*pCos_gamma = 1.f ;
    }

    *pTheta = angle_between(&view,&V_zenith,&dp_zenith) ;

    if(0) {
	// JJW reduce angle between zenith and horizon,
	// which will raise the level of horizon glow (via A & B coefs and Perez formula)
	float horizon_angle_factor = (M_PI/2.)/(M_PI/2.-pData_fit->horizon_lift_angle_radians) ;

	*pTheta *= horizon_angle_factor ;
    }

    if(*pTheta > M_PI/2.) *pTheta = M_PI/2. ;

    // clamp angle between view and zenith at PI/2 so prediction
    // values don't fall below horizon which is not good at all
    if(*pTheta > M_PI/2.-0.01) *pTheta=M_PI/2.-.01 ;

    if(isnan(*pGamma)) {
	fprintf(stderr, "NAN in F_CIE2003, sun_x=%f sun_py=%f px=%f py=%f horizon_py=%f theta=%f gamma=%f\n",
	    pData_fit->sun_x, pData_fit->sun_py, px, py, pData_fit->horizon_py, *pTheta, *pGamma) ;
	fprintf(stderr, "View: %f %f %f\n", view.A, view.B, view.C) ;
	fprintf(stderr, "Sun: %f %f %f\n", V_sun.A, V_sun.B, V_sun.C) ;
	float dotproduct = view.A*V_sun.A + view.B*V_sun.B + view.C*V_sun.C ;
	fprintf(stderr, "dotproduct %f\n", dotproduct) ;
	fprintf(stderr, "acos %f\n", acos(dotproduct)) ;
	exit(1) ;
    }

/*      float L = (1.+perez_A*exp(perez_B/cos(*pTheta)))*(1.+perez_C*(exp(perez_D*(*pGamma))-exp(perez_D*M_PI/2.)) + perez_E*(*pCos_gamma)*(*pCos_gamma)) ;  */

    // reduced model, only adds direct glow from sun
/*      float x = powf(*pGamma*pData_fit->perez_F,pData_fit->perez_G) ;  */
/*      float L = (1.+pData_fit->perez_C*(exp(pData_fit->perez_D*x)-exp(pData_fit->perez_D*M_PI/2.))) ;  */

    // note by using the negative inverse of D, the sun diameter as rendered is linear with D
    float L = (1.+pData_fit->perez_C*(exp((-1./pData_fit->perez_D)*(*pGamma))-exp((-1./pData_fit->perez_D)*M_PI/2.))) ;

    if(isnan(L)) {
	fprintf(stderr, "NAN in F_CIE2003, sun_x=%f sun_py=%f px=%f py=%f horizon_py=%f theta=%f gamma=%f\n",
	    pData_fit->sun_x, pData_fit->sun_py, px, py, pData_fit->horizon_py, *pTheta, *pGamma) ;
	fprintf(stderr, "View: %f %f %f\n", view.A, view.B, view.C) ;
	fprintf(stderr, "Zenith: %f %f %f\n", V_zenith.A, V_zenith.B, V_zenith.C) ;
	float dotproduct = view.A*V_zenith.A + view.B*V_zenith.B + view.C*V_zenith.C ;
	fprintf(stderr, "dotproduct %f\n", dotproduct) ;
	fprintf(stderr, "acos %f\n", acos(dotproduct)) ;
	exit(1) ;
    }

    if(pData_fit->verbose) {
	fprintf(stderr, "sky L:%5.2f\n", L) ;
    }

    return L ;
}

float find_maximum_CIE_vhat(void)
{
    float gamma,theta,cos_gamma ;
    float px, py ;

    if(uses_CIE_model) {
	float sun_px = sun_x_angle2px(pData_fit->sun_x) ;
	float max_vhat = 0. ;
	pData_fit->minimum_CIE_vhat = 1.e30 ;

	if(sun_px >= -0.5 && sun_px <= 0.5 && pData_fit->sun_py >= 0.0 && pData_fit->sun_py <= 1.0) {
	    // if the sun is in the image, get vhat directly at the sun
	    max_vhat = F_CIE2003(sun_x_angle2px(pData_fit->sun_x), pData_fit->sun_py, &gamma, &theta, &cos_gamma) ;
	}


	for(px = -0.5f ; px < 0.501f ; px += .05f) {
	    for(py = 0.f ; py < 1.001f ; py += .05f) {
		float vhat_here = F_CIE2003(px, py, &gamma, &theta, &cos_gamma) ;
		if(vhat_here > max_vhat) max_vhat = vhat_here ;
		if(vhat_here < pData_fit->minimum_CIE_vhat) pData_fit->minimum_CIE_vhat = vhat_here ;
	    }
	}
	return max_vhat ;
    } else {
	pData_fit->minimum_CIE_vhat = 0. ;
	return 1. ;
    }

}

// factor to make mean v from predicted "hsv" match mean v from actual
/*  float sample_based_v_correction = 1.f ;  */

void predict_sky_huesat_from_val(float vhat_CIE_sun, float *pH, float *pS, float *pV, float px, float py)
{

    float py_hc=py_adjusted_for_horizon_curvature(px,py)  ;
    float vhat_final ;

    float hhat     = h_from_pxpy(px, py_hc, pData_fit->full_hsv_model_level) ;
    float shat     = s_from_pxpy(px, py_hc, pData_fit->full_hsv_model_level, pData_fit->horizon_py) ;
    float vhat_raw = v_from_pxpy(px, py_hc, pData_fit->full_hsv_model_level) ;

    if(uses_CIE_model) {
	if(vhat_raw > vhat_CIE_sun)
	    vhat_final = vhat_raw ;
	else
	    vhat_final = vhat_CIE_sun ;

	vhat_CIE_sun*=pData_fit->maximum_CIE_vhat ;
	if(pData_fit->sky_lum > 0.) {
	    // vhat_CIE_sun now runs between max and min, want to scale it to 0.5 to 0.0
	    float vhat_add = pData_fit->sky_lum*(vhat_CIE_sun - pData_fit->minimum_CIE_vhat)/(pData_fit->maximum_CIE_vhat-pData_fit->minimum_CIE_vhat) ;

	    if(vhat_add < 0.0) vhat_add = 0.0 ;

	    vhat_final = vhat_add + vhat_raw ;
	    shat -= vhat_add ;

	} else {
	    // vhat_CIE_sun now runs between max and min, want to scale it to 0.5 to 0.0
	    // since sky_lum is < 0, use -sky_lum and only display CIE vhat
	    float vhat_add = -pData_fit->sky_lum*(vhat_CIE_sun - pData_fit->minimum_CIE_vhat)/(pData_fit->maximum_CIE_vhat-pData_fit->minimum_CIE_vhat) ;

	    if(vhat_add < 0.0) vhat_add = 0.0 ;
	    vhat_final = vhat_add ;
	}

    } else {
	vhat_final = vhat_raw ;
    }

    if(vhat_final > 1.f) vhat_final = 1.f ;
    if(vhat_final < 0.f) vhat_final = 0.f ;


    if(shat > 1.0) shat = 1.0 ;
    if(shat < 0.0) shat = 0.0 ;

    shat *= pData_fit->sat_depth ;

    *pH = hhat ;
    *pS = shat ;
    *pV = vhat_final ;
}

void predict_sky_hsv(float px, float py, float *pH, float *pS, float *pV)
{
    float vhat_sun = 0. ;

    float gamma ; // angle between view and sun vector
    float theta ; // angle between view and zenith vector
    float cos_gamma ;

    if(uses_CIE_model) {
	vhat_sun = F_CIE2003(px, py, &gamma, &theta, &cos_gamma)/pData_fit->maximum_CIE_vhat ;

	if(vhat_sun < 0.) {
/*  	    fprintf(stderr, "WARNING: F_CIE2003 returns negative vhat, fixing!\n") ;  */
	    vhat_sun = 0. ;
	}

	if(vhat_sun > 1.f) vhat_sun = 1.f ;
	if(vhat_sun < 0.f) vhat_sun = 0.f ;
    }

    predict_sky_huesat_from_val(vhat_sun, pH, pS, pV, px, py) ;

    *pV *= pData_fit->exposure_factor ;
}

void predict_sky_color(float px, float py, uint16_t *rhat, uint16_t *ghat, uint16_t*bhat)
{
    float hhat, shat, vhat ;

    predict_sky_hsv(px,py,&hhat,&shat,&vhat) ;

    hsv2rgb16(hhat,shat,vhat,rhat,ghat,bhat) ;
}


// a function to reduce horizontal variation of sky r,g,b values
void make_sky_uniform(tdata_t *image, int flag, int x0, int x1, SKYFILL_DATA_t *pData)
{
    // flag=0 => run smoother from top of image(y0) to lowest end of sky(y1)
    // flag=1 => run smoother from top of image(y0) to lowest start of sky(y1)
    // flag=2 => run smoother from highest start of sky(y0) to lowest end of sky(y1)

    // if flag is 0 or 1, this smoother is run after the entire sky is predicted
    // if flag is 2, this smoother is run before the optimization and sky prediction

/*      int min_esky_y=IMAGE_HEIGHT-1 ;  */
    int min_ssky_y=IMAGE_HEIGHT-1 ;
    int max_esky_y=0 ;
    int max_ssky_y=0 ;

    for(int x = x0 ; x <= x1 ; x++) {
/*  	if(min_esky_y > end_of_sky[x]) min_esky_y = end_of_sky[x] ;  */
	if(min_ssky_y > pData->start_of_sky[x]) min_ssky_y = pData->start_of_sky[x] ;
	if(max_esky_y < pData->end_of_sky[x]) max_esky_y = pData->end_of_sky[x] ;
	if(max_ssky_y < pData->start_of_sky[x]) max_ssky_y = pData->start_of_sky[x] ;
    }

    int y0 = 0 ;
    if(flag == 2) y0 = min_ssky_y ;

    int y1=max_esky_y ;
    if(flag == 1) y1 = max_ssky_y ;

    fprintf(stderr, "make sky uniform, y %d to %d\n", y0, y1) ;

    for(int y = y0 ; y <= y1 ; y++) {
	float sum[3] = {0.,0.,0.} ;
	int n=0 ;
	for(int x = x0 ; x <= x1 ; x++) {
	    if(pData->end_of_sky[x] >= y) {
		if(flag == 2 && y < pData->start_of_sky[x]) continue ;

		uint16_t r,g,b ;
		tif_get3c(image,x,y,r,g,b) ;

		sum[0] += r ;
		sum[1] += g ;
		sum[2] += b ;
		n++ ;
	    }
	}

	if(n > 0) {
	    sum[0] /= (float)n ;
	    sum[1] /= (float)n ;
	    sum[2] /= (float)n ;

    // flag == 0, blend is 0 at eos[x], to 1 at sos[x]
    // flag == 1, blend is 0 y1, to 1 at sos[x]
    // flag == 2, blend is 0 y1, to 1 at (eos[x]+sos[x])/2

	    for(int x = x0 ; x <= x1 ; x++) {
		if(y <= pData->start_of_sky[x]) {
		    tif_set4c(image,x,y, sum[0], sum[1], sum[2], MAX16);
		    continue ;
		}

		if(flag == 0) {
		    int blend_y0 = (int)((pData->start_of_sky[x]+pData->end_of_sky[x])/2.+0.5) ;
		    int blend_y1 = pData->end_of_sky[x] ;
		    if(y > pData->end_of_sky[x]) continue ;

		    if(y <= blend_y0) {
			tif_set4c(image,x,y, sum[0], sum[1], sum[2], MAX16);
			continue ;
		    }

		    float pmean = 1.-(float)(y-blend_y0)/(float)(blend_y1-blend_y0) ;
		    float r = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+0] ;
		    float g = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+1] ;
		    float b = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+2] ;
		    r = pmean*sum[0] + (1.-pmean)*r ;
		    g = pmean*sum[1] + (1.-pmean)*g ;
		    b = pmean*sum[2] + (1.-pmean)*b ;
		    tif_set4c(image,x,y,r,g,b, MAX16) ;
		    continue ;
		}

		float blend_y0 = pData->start_of_sky[x] ;
		float blend_y1 = y1 ;

		if(flag == 2) blend_y0 = (float)(pData->end_of_sky[x]+pData->start_of_sky[x])/2. ;
		if(flag == 1 && blend_y1 > (float)pData->end_of_sky[x]) blend_y1 = (float)pData->end_of_sky[x] ;

		if(blend_y1 > pData->end_of_sky[x]) blend_y1 = (float)pData->end_of_sky[x] ;
		if(flag == 1) blend_y1 = (float)(pData->end_of_sky[x]+pData->start_of_sky[x])/2. ;
		if(flag == 1) blend_y1 = pData->end_of_sky[x] ;

		float pmean = 1.-((float)y-blend_y0)/(float)(blend_y1-blend_y0) ; //at y == y0, use 100% of mean value, linear change to 0% of mean at y == y1

		if(pmean > 0.) {
		    float r = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+0] ;
		    float g = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+1] ;
		    float b = ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+2] ;
		    r = pmean*sum[0] + (1.-pmean)*r ;
		    g = pmean*sum[1] + (1.-pmean)*g ;
		    b = pmean*sum[2] + (1.-pmean)*b ;
		    tif_set4c(image,x,y,r,g,b, MAX16) ;
		}
	    }

	}
    }

    //should we move start of sky up?
    if(flag == 2) {
	// since this happened *before* sky sampling/optimization
	for(int x = x0 ; x <= x1 ; x++) {
	    pData->start_of_sky[x] = y0 ;
	}
    }
}

void usage(char *msg, char *msg2)
{
    fprintf(stderr, "skyfill_tif: %s %s\n", msg, msg2) ;

    fprintf(stderr, "usage: skyfill_tif <inputfile> <outputfile> [options]\n") ;
    fprintf(stderr, "options are:\n") ;
    fprintf(stderr, "    -ds <DOS> -- sample sky values from sky start to (sky end-sky start)*DOS+sky_start in a column, DOS is a proportion, default 0.5\n") ;
    fprintf(stderr, "    -cdf <DOF> -- fill sky values from O to (min sky end- max sky start)*DOF+sky_start in all columns, DOF is a proportion, default 0.5\n") ;
    fprintf(stderr, "              useful when there is a lot of small peaks at the skyline which are visible in the sky merge.  Try if -df doesn't look good\n") ;
    fprintf(stderr, "    -df <DOF> -- fill sky values from O to (sky end-sky_start)*DOF+sky_start in a column, DOF is a proportion, default 0.5\n") ;
    fprintf(stderr, "    -ds <depth_of_sample> -- proportion, control how far in a column to find samples to compare to predicted sky\n") ;
    fprintf(stderr, "             -- affects final blending of predicted sky to original sky at that column\n") ;
    fprintf(stderr, "    -ff <FF> -- feather factor.  Default 1.0 Normally replaced sky is completely estimated at sky start, and\n") ;
    fprintf(stderr, "              uses 100%% of the orginal value at the final row (determined by DOF above\n") ;
    fprintf(stderr, "              setting FF to 0.0 would cause 100%% of the estimated value to be use on every\n") ;
    fprintf(stderr, "              pixel replaced\n") ;
    fprintf(stderr, "    -fse <0|1>  -- turn off|on algorithm to repair vertical edges with missing pixels near the start of the sky\n") ;
    fprintf(stderr, "    -fsg <gamma>  -- a typical image gamma function applied to predicted value of sky (v = v**(1./gamma))\n") ;
    fprintf(stderr, "    -fsh <left> <right> -- attempt to repair sky hue in a portion of the image, left and right are proportions\n") ;
    fprintf(stderr, "                       fsh 0.5 1.0 would attempt the repair on the right half of the image\n") ;
    fprintf(stderr, "    -ftos <0|1> -- brings top of sky up to a constant value before modeling sky values if set to 1\n") ;
    fprintf(stderr, "    -efl <number_pixels>  -- extend computed feather length of all columns by this many pixels \n") ;
    fprintf(stderr, "    -tol <th> <ts> <tv> -- tolerance for determining edge of skyline at horizon, default 360 .06 0.015\n") ;
    fprintf(stderr, "    -im -- interpolate sky in masked columns\n") ;
    fprintf(stderr, "    -hf -- run horizontal filter on sky area to smooth artifacts\n") ;


    fprintf(stderr, "    -fm <l:q><l:q><l:q> -- fit model type for hsv, l is linear, q is quadratic, default is qqq (quadratic fit for h,s,v components\n") ;
    fprintf(stderr, "    -m <l> <r> -- mask out columns from l-r from sky sampling or filling (up to 20 masks may be specified)\n") ;
    fprintf(stderr, "    -d -- debug, output image will have start and end of sky marked\n") ;
    fprintf(stderr, "    -G -- run the grid optimizer, using fixed parameters supplied on command line\n") ;
    fprintf(stderr, "    -rm <l> <r> -- mask out columns from l-r from having pixels replaced i.e. repaird because they appear to be outliers\n") ;
    fprintf(stderr, "    -r_tos_thresh <thresh> -- sets the tolerance for finding \"outlier\" pixels in a small section of the sky, if the pixel\n") ;
    fprintf(stderr, "                       is more than this fraction away from the mean, it is classified as an outlier.  Default thresh is 0.02\n") ;
    fprintf(stderr, "    -spm <0|1> -- saturation prediction method 0 will just use a mean value, 1 will predict saturation as a function of value\n") ;
    fprintf(stderr, "    -verbose -- adds a little information at the end of the run, really only useful for debuggin\n") ;
    fprintf(stderr, "    -t -- quick test most, scales the output image by 50%% to speed up processing\n") ;
    fprintf(stderr, "    -EO -- show estimated sky only, will replace the entire image with the optimized estimate of the sky\n") ;
    fprintf(stderr, "    -SR -- show estimated sky only for the sky area to be filled -- it will not have final blending applied\n") ;



    fprintf(stderr, "\n") ;
    fprintf(stderr, "    -O <one or more parms> optimize given parameters of sky dome, possible parameters are:\n") ;

    int i ;

    for(i = 0 ; i < MAX_OPT_PARMS ; i++) {
	fprintf(stderr, "               %s (%s) range is %lf to %lf --  %s\n",
	    opt_parms[i].abreviation,
	    opt_parms[i].name,
	    opt_parms[i].lo,
	    opt_parms[i].hi,
	    opt_parms[i].units) ;
    }
    fprintf(stderr, "\n               More than one -O may be specified, causing multiple optimizations to be performed\n") ;
    fprintf(stderr, "               sequentially with individual control over which parameters are optimized\n\n") ;

    fprintf(stderr, "    -any parameter may be specified by prepending with a '-', i.e. \"-sx 90\" sets the sun x position 90 degrees to the right\n") ;

    exit(1) ;
}

#ifdef TEST_HSI_MODEL
// temporary, see if implementation of HSI model conversion is good
void test_I_model(void)
{
    uint16_t r0,g0,b0 ;
    uint16_t r1,g1,b1 ;
    float H, S, I ;
    uint16_t values[11] ;

    for(int i=0 ; i <= 10 ; i++) {
	values[i] = (uint16_t)((float)i/10.*MAX16f+0.5) ;
    }
    int i_r,i_g,i_b ;
    for(i_r=0 ; i_r <= 10 ; i_r++) {
    for(i_g=0 ; i_g <= 10 ; i_g++) {
    for(i_b=0 ; i_b <= 10 ; i_b++) {
	r0 = values[i_r] ;
	g0 = values[i_g] ;
	b0 = values[i_b] ;
	rgb2hsv16(r0,g0,b0,&H,&S,&I) ;
	hsv2rgb16(H,S,I,&r1,&g1,&b1) ;
	int good=1 ;
	if( abs(r0-r1) > 1.) good=0 ;
	if( abs(g0-g1) > 1.) good=0 ;
	if( abs(b0-b1) > 1.) good=0 ;
	if(! good) {
	    float m = I*(1.-S)*MAX16f ;
	    fprintf(stderr, "HSI fails, rgb in (%d,%d,%d), HSI (%f,%f,%f), rgb out(%d,%d,%d), m:%f\n", r0,g0,b0,H,S,I,r1,g1,b1,m) ;
	}
    }
    }
    }
    exit(1) ;
}
#endif

// setup an array (0 to IMAGE_WIDTH-1, given a number of specified regions to mask
void setup_mask_array(unsigned char **array, int *mask_l, int *mask_r, int n_masks, int input_W, int quick_test)
{
    // apply column masks to array

    *array = (unsigned char *)calloc(input_W, sizeof(unsigned char)) ;
    for(int x = 0 ; x < input_W ; x++) {
	(*array)[x/(quick_test+1)] = 0 ;
    }
    for(int m = 0 ; m < n_masks ; m++) {
	for(int x = mask_l[m] ; x <= mask_r[m] ; x++) {
	    if(x < 0) continue ;
	    if(x > input_W-1) continue ;
	    (*array)[x/(quick_test+1)] = 1 ;
	}
    }
}

int main(int argc, char* argv[])
{
    float depth_of_sample = 0.5 ; // how far into the sky to sample for hue,val and sat
    float depth_of_fill = 0.5 ; // how far into existing sky to fill
    float feather_factor = 1.0 ; // how much feathering to do
    int extra_feather_length = 0 ; // extend feathering length by this many pixels
    int run_hf=0 ; // run horizontal filter
    int debug=0 ;
#ifdef TEST_HSI_MODEL
    test_I_model() ;
#endif

    // data structure to hold data & variables necessary for program
    SKYFILL_DATA_t skyfill_data ;
    SKYFILL_DATA_t *pData = &skyfill_data ;

    // set the global pointer
    pData_fit = pData ;

    initialize_skyfill_data_struct(pData) ;
    samples = NULL ; // make sure the samples array pointer is NULL
    fill_opt_parms(pData) ;

    int n_masks=0 ;
    int mask_l[20] ;
    int mask_r[20] ;

    int n_repair_masks=0 ;
    int repair_mask_l[20] ;
    int repair_mask_r[20] ;
    int n_nonsky_repair_masks=0 ;
    int nonsky_repair_mask_l[20] ;
    int nonsky_repair_mask_r[20] ;

    int seam_line_smooth_size = -1 ;  // if greater than 0, run a smoother along the start of sky seam line

    int n_sample_masks=0 ;
    int sample_mask_l[20] ;
    int sample_mask_r[20] ;

    char infilename[512] ;
    char outfilename[512] ;
    char ptofilename[512] ;
    char eosfilename[512] ;
    int interpolate_masked_columns=0 ;
    int verbose_flag=0 ;
    int use_exiftool=1 ;
    int make_sky_uniform_flag=-1 ;
    int quick_test=0 ;
    int df_was_set=0 ;
    float fmin_sky_hue_mask = 2. ;
    float fmax_sky_hue_mask = 2. ;


    if(argc < 3)
	usage("Input or output file not specified", "") ;

    // strip final \n off the last input argument
    {
	int l = strlen(argv[argc-1]) ;

	if(argv[argc-1][l-1] == '\n')
	    argv[argc-1][l-1] = '\0' ;
    }

    // make sure input tif is not the same as output tif
    if(! strcmp(argv[1], argv[2]) ) {
	usage("Input file name must not be same as output filename", "") ;
    }


    strcpy(infilename, argv[1]) ;
    strcpy(outfilename, argv[2]) ;
    strcpy(ptofilename, argv[1]) ;
    strcpy(eosfilename, argv[1]) ;
    int pto_fov=-1 ;

    // can we open a corresponding hugin pto file, and get the horizontal FOV from it?
    {
	int i ;

	for(i=0 ; i < strlen(ptofilename) ; i++) {
	    if(ptofilename[i] == '.')
		break ;
	}

	ptofilename[i+1] = 'p' ;
	ptofilename[i+2] = 't' ;
	ptofilename[i+3] = 'o' ;
	ptofilename[i+4] = '\0' ;
	eosfilename[i+1] = 'e' ;
	eosfilename[i+2] = 'o' ;
	eosfilename[i+3] = 's' ;
	eosfilename[i+4] = '\0' ;
	printf("ptofilename is %s\n", ptofilename) ;
	printf("eosfilename is %s\n", eosfilename) ;
    }

    FILE *fp = fopen(ptofilename,"r") ;

    if(fp != NULL) {
	char s[1000] ;
	while(fgets(s, 999, fp)) {
	    int j ;

	    if(s[0] != 'p')
		continue ;

	    for(j=0 ; j < strlen(s) ; j++) {
		if(s[j] == 'v')
		    break ;
	    }

	    if(s[j] == 'v' && s[j-1] == ' ') {
		sscanf(&s[j+1], "%d", &pto_fov) ;
		printf("PTO fov is %d\n", pto_fov) ;
	    }
	}

	fclose(fp) ;
    }


    TIFF* tif = TIFFOpen(infilename, "r");
    argc -= 3 ;
    argv += 3 ;
    int n_custom_optimizations=0 ;
    int force_grid_optimization = 0 ;

    int optimization_vars[10][MAX_OPT_PARMS+1] ;

    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &IMAGE_NSAMPLES);

    if(IMAGE_NSAMPLES < 4) {
	fprintf(stderr, "FATAL:  skyfill only works with 4 channel (rgba) images produced by hugin\n") ;
	exit(1) ;
    }

    uint32_t input_W, input_H ;
    uint16_t BPS, SPP, TIF_CONFIG ;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &input_W);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &input_H);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &BPS);
    TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &TIF_CONFIG);

    if(BPS != 16) {
	fprintf(stderr, "FATAL:  skyfill only works with 16 bit tiff images\n") ;
	exit(1) ;
    }

    uses_CIE_model=0 ;


    while(argc > 0) {
	int found_set_parm=0 ;
	int i ;

	for(i = 0 ; i < MAX_OPT_PARMS ; i++) {
	    char set_string[20]  ;
	    snprintf(set_string, 19, "-%s", opt_parms[i].abreviation) ;

	    if(!strcmp(argv[0], set_string)) {

		if(opt_parms[i].is_CIE)
		    uses_CIE_model=1 ;

		if(!strcmp(opt_parms[i].abreviation, "D")) {
		    // D parameter specifies Degrees of view for sun
		    // to make it approximately correct in the formula need to divide by 100
		    *opt_parms[i].variable_address = atof(argv[1])/100. ;
		} else {
		    *opt_parms[i].variable_address = atof(argv[1]) ;
		}
		opt_parms[i].grid_optimize_flag = 0 ;
		opt_parms[i].optimize_flag = 0 ;
		printf("setting %s to %f\n", opt_parms[i].name, *opt_parms[i].variable_address) ;
		found_set_parm=1 ;

		if(!strcmp(opt_parms[i].abreviation, "hy"))
		    pData->horizon_was_set=1 ;

		argv++ ;
		argc-- ;
		argv++ ;
		argc-- ;
		break ;
	    }
	}

	if(found_set_parm) continue ;


	if(!strcmp(argv[0], "-O")) {
	    // optimization requested, which variables?
	    argc-- ;
	    argv++ ;
	    int n_vars=0 ;

	    while(argc > 0) {
		if(argv[0][0] != '-') {
		    for(i = 0 ; i < MAX_OPT_PARMS ; i++) {
			if(!strcmp(argv[0], opt_parms[i].abreviation)) {

			    if(opt_parms[i].is_CIE)
				uses_CIE_model=1 ;

			    optimization_vars[n_custom_optimizations][n_vars] = i ;
			    n_vars++ ;
			    optimization_vars[n_custom_optimizations][n_vars] = -1 ;
			    printf("optimization %d, optimize %s\n", n_custom_optimizations, opt_parms[i].name) ;
			    argc-- ;
			    argv++ ;
			    break ;
			}
		    }
		} else {
		    break ;
		}
	    }
	    n_custom_optimizations++ ;
	    continue ;
	}

	if(!strcmp(argv[0], "-CIE")) {
	    pData->CIE_sky_index = atoi(argv[1]) ;
	    set_cie_parms(pData->CIE_sky_index-1) ;
	    uses_CIE_model=1 ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-CIE_types")) {

	    if(argv[1][0] == 'a')
		pData->allowed_sky_type = ALL_CIE_SKIES ;
	    else if(argv[1][0] == 'n')
		pData->allowed_sky_type = NONUNIFORM_CIE_SKIES ;
	    else if(argv[1][0] == 'u')
		pData->allowed_sky_type = UNIFORM_CIE_SKIES ;
	    else {
		usage("unknown sky type for -CIE_types, allowed are 'a', 'n' and 'u'", "") ;
	    }

	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-d")) {
	    if(debug < 1) debug=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-t")) {
	    quick_test=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}

	if(strlen(argv[0]) == 4 || strlen(argv[0]) == 3) {
	    if(argv[0][0] == '-' && argv[0][1] == 'd') {
		if(argv[0][2] >= '0' && argv[0][2] <= '9') {
		    debug = atoi(&(argv[0][2])) ;
		    fprintf(stderr, "debug set to %d\n", debug) ;
		    argc -= 1 ;
		    argv += 1 ;
		    continue ;
		}
	    }
	}

	if(!strcmp(argv[0], "-cdf")) {
	    depth_of_fill = atof(argv[1]) ;
	    pData->depth_of_fill_is_absolute = 1 ;
	    df_was_set=1 ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-df")) {
	    depth_of_fill = atof(argv[1]) ;
	    df_was_set=1 ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-ds")) {
	    depth_of_sample = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-efl")) {
	    extra_feather_length = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fhml")) {
	    pData->full_hsv_model_level= atoi(argv[1]) ;
	    fprintf(stderr, "full_hsv_model_level set to %d\n", pData->full_hsv_model_level) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-ff")) {
	    feather_factor = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fs")) {
	    pData->final_saturation_factor = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fse")) {
	    pData->fix_SOS_edges = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fsr")) {
	    pData->full_sky_replacement = 1 ;
	    if(! df_was_set) {
		depth_of_fill= 0.5 ;
		pData->nonlinear_feather = 1. ;
	    }

	    pData->noclip_end_of_sky = 0 ;

	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fsrt")) {
	    pData->full_sky_replacement_thresh = atof(argv[1]) ;
	    pData->full_sky_replacement_ramp = atof(argv[2]) ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fsh")) {
	    pData->fix_sky_hue=1 ;
	    fmin_sky_hue_mask = atof(argv[1]) ;
	    fmax_sky_hue_mask = atof(argv[2]) ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-ftos")) {
	    pData->fill_top_of_sky = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-hf")) {
	    run_hf=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-im")) {
	    interpolate_masked_columns = 1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-m")) {
	    mask_l[n_masks] = atoi(argv[1]) ;
	    mask_r[n_masks] = atoi(argv[2]) ;
	    n_masks++ ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-tm")) {
	    if(pData->n_test_masks == MAX_TEST_MASKS) {
		fprintf(stderr, "FATAL: too many test masks specified, maximum is %d\n", MAX_TEST_MASKS) ;
		exit(1) ;
	    }
	    pData->test_extent_mask[pData->n_test_masks].l = atoi(argv[1]) ;
	    pData->test_extent_mask[pData->n_test_masks].r = atoi(argv[2]) ;
	    pData->test_extent_mask[pData->n_test_masks].t = atoi(argv[3]) ;
	    pData->test_extent_mask[pData->n_test_masks].b = atoi(argv[4]) ;
	    pData->n_test_masks++ ;
	    argc -= 5 ;
	    argv += 5 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-mse")) {
	    pData->min_sky_end_p = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-msu")) {
	    make_sky_uniform_flag = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-nceos")) {
	    pData->noclip_end_of_sky = 1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-nf")) {
	    pData->nonlinear_feather = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-nsrm")) {
	    nonsky_repair_mask_l[n_nonsky_repair_masks] = atoi(argv[1]) ;
	    nonsky_repair_mask_r[n_nonsky_repair_masks] = atoi(argv[2]) ;
	    n_nonsky_repair_masks++ ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-r_tos_thresh")) {
	    pData->repair_tos_thresh = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-rm")) {
	    repair_mask_l[n_repair_masks] = atoi(argv[1]) ;
	    repair_mask_r[n_repair_masks] = atoi(argv[2]) ;
	    n_repair_masks++ ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-sls")) {
	    seam_line_smooth_size = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-sm")) {
	    sample_mask_l[n_repair_masks] = atoi(argv[1]) ;
	    sample_mask_r[n_repair_masks] = atoi(argv[2]) ;
	    n_sample_masks++ ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-spm")) {
	    pData->sat_prediction_method = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-tol")) {
	    pData->sky_hue_tolerance = atof(argv[1]) ;
	    pData->sky_sat_tolerance = atof(argv[2]) ;
	    pData->sky_val_tolerance = atof(argv[3]) ;
	    argc -= 4 ;
	    argv += 4 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-verbose")) {
	    verbose_flag = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-EO")) {
	    pData->estimate_only=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-G")) {
	    force_grid_optimization=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-SR")) {
	    pData->show_raw_prediction=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-SRE")) {
	    pData->show_raw_error=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-SSP")) {
	    pData->show_sky_prediction=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-vm")) {
	    pData->val_model_full = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}
	usage("Unknown option", argv[0]) ;
    }

    if(pto_fov > 1) {
	pData->FOV_horizontal = (float)pto_fov ;
	pData->have_pto_fov=1 ;
    }
    if (!tif) {
	fprintf(stderr, "failed to open %s\n", infilename) ;
	exit(1) ;
    }


    TIFF* tifout = TIFFOpen(outfilename, "w");

    tsize_t ROWSIZE ;
    int32_t x, y, w, h ;

    // prepare output image to have same tags as input image
    uint16_t alltags[] = {
			TIFFTAG_IMAGEWIDTH,
			TIFFTAG_IMAGELENGTH,
			TIFFTAG_BITSPERSAMPLE,
			TIFFTAG_SAMPLESPERPIXEL,
			TIFFTAG_ROWSPERSTRIP,
			TIFFTAG_ORIENTATION,
			TIFFTAG_PLANARCONFIG,
			TIFFTAG_PHOTOMETRIC,
			TIFFTAG_SAMPLEFORMAT,
			TIFFTAG_COMPRESSION,
			} ;
    int i ;

    for(i = 0 ; i < 10 ; i++) {
	uint16_t tag = alltags[i] ;
	uint32_t tmp ;

	if(tag == TIFFTAG_IMAGEWIDTH && quick_test) {
	    if(TIFFGetField(tif, tag, &tmp)) {
		tmp /= 2 ;
		TIFFSetField(tifout, tag, tmp);
	    }
	} else if(tag == TIFFTAG_IMAGELENGTH && quick_test) {
	    if(TIFFGetField(tif, tag, &tmp)) {
		tmp /= 2 ;
		TIFFSetField(tifout, tag, tmp);
	    }
	} else {
	    if(TIFFGetField(tif, tag, &tmp)) {
		fprintf(stderr, "setting tag number %d in output file\n", tag) ;
		TIFFSetField(tifout, tag, tmp);
	    }
	}
    }

    if(!use_exiftool) {
	{
	    uint16_t tag = TIFFTAG_ICCPROFILE ;

	    uint32_t tmp[4] ;
	    unsigned char buf[4][2000] ;
	    if(TIFFGetField(tif, tag, tmp, buf)) {
		fprintf(stderr, "setting ICCPROFILE tag in output file, tmp is %d,%d,%d\n", (int)tmp[0],(int)tmp[1],(int)tmp[2]) ;
		TIFFSetField(tifout, tag, tmp[0], buf);
	    }
	    //exit(1) ;
	}
    }


    tdata_t *image ;


    ROWSIZE = TIFFScanlineSize(tif) ;
    SPP = IMAGE_NSAMPLES ;

    IMAGE_HAS_ALPHA = IMAGE_NSAMPLES > 3 ? 1 : 0 ;


    // make signed int values for width and height
    h = (int32_t)input_H ;
    w = (int32_t)input_W ;

    if(quick_test) {
	h /= 2 ;
	w /= 2 ;
    }

    pData->min_sky_hue_mask = fmin_sky_hue_mask*w ;
    pData->max_sky_hue_mask = fmax_sky_hue_mask*w ;

    IMAGE_HEIGHT = (int32_t)h ; // set global variable
    IMAGE_WIDTH = (int32_t)w ; // set global variable

    p_half_image_width = 0.5 ;
    V_zenith= P3toV3(0., 1., 0.) ;
    V_sun= P3toV3(0., 1., 0.) ; // default with sun directly overhead ;
    pData->maximum_CIE_vhat = find_maximum_CIE_vhat() ;

    pData->manual_end_of_sky = (int16_t *)calloc(IMAGE_WIDTH, sizeof(int16_t *)) ;
    for(int i = 0 ; i < IMAGE_WIDTH ; i++)
	pData->manual_end_of_sky[i] = -1 ;

    fp = NULL ;
    fp = fopen(eosfilename, "r") ;

    if(fp != NULL) {
	fscanf(fp, "%*s %*s") ; // skip "X Y"
	int x0=-1 ;
	int y0=-1 ;
	int x1, y1 ;

	while(fscanf(fp, "%d %d", &x1, &y1) == 2) {
	    if(quick_test) {
		x1 /= 2 ;
		y1 /= 2 ;
	    }

	    if(x0 != -1 && x1 != -1) {
		float dy = (float)(y1-y0)/(float)(x1-x0) ;
		float fy = y0 ;
		for(int x=x0 ; x<= x1 ; x++) {
		    int y = (int)(fy+0.5) ;
		    pData->manual_end_of_sky[x] = y ;
/*  		    printf("MEOS x:%d y:%d\n", x, y) ;  */
		    fy += dy ;
		}
	    }

	    x0 = x1 ;
	    y0 = y1 ;
	}

	fclose(fp) ;
    }

    fprintf(stderr, "input image is %d x %d, %d bits, %d channels, %d = %d\n", (int)input_W, (int)input_H, (int)BPS, (int)SPP, (int)(input_W*BPS/8*SPP), (int)ROWSIZE) ;

    // apply column masks to column_mask array
    setup_mask_array(&pData->column_mask, mask_l, mask_r, n_masks, input_W, quick_test) ;
    setup_mask_array(&pData->column_repair_mask, repair_mask_l, repair_mask_r, n_repair_masks, input_W, quick_test) ;
    setup_mask_array(&pData->column_nonsky_repair_mask, nonsky_repair_mask_l, nonsky_repair_mask_r, n_nonsky_repair_masks, input_W, quick_test) ;
    setup_mask_array(&pData->column_sample_mask, sample_mask_l, sample_mask_r, n_sample_masks, input_W, quick_test) ;


    int x0 ;

    for(x0 = 0 ; x0 < w ; x0++) {
	if(pData->column_mask[x0] == 0) break ;
    }

    if(x0 >= w) {
	fprintf(stderr, "All columns are masked!, no output generated\n") ;
	exit(1) ;
    }



    // apply column repair masks to column_repair_mask array

/*      pData->column_repair_mask = (unsigned char *)calloc(w, sizeof(unsigned char)) ;  */
/*      for(x = 0 ; x < w ; x++) {  */
/*  	pData->column_repair_mask[x/(quick_test+1)] = 0 ;  */
/*      }  */
/*      for(int m = 0 ; m < n_repair_masks ; m++) {  */
/*  	for(x = repair_mask_l[m] ; x <= repair_mask_r[m] ; x++) {  */
/*  	    if(x < 0) continue ;  */
/*  	    if(x > input_W-1) continue ;  */
/*  	    pData->column_repair_mask[x/(quick_test+1)] = 1 ;  */
/*  	}  */
/*      }  */

    // Read the input tif image

    if((int)ROWSIZE != (int)(input_W*BPS/8*SPP)) {
	fprintf(stderr, "ROWSIZE is WRONG!, fixing...\n") ;
	ROWSIZE = input_W*BPS/8*SPP ;
	fprintf(stderr, "ROWSIZE is now %d\n", (int)ROWSIZE) ;
    }

    /* allocate memory for image */
    image = (tdata_t *)calloc(input_H, sizeof(tdata_t *)) ;

    for(y = 0 ; y < input_H ; y++) {
	image[y] = (tdata_t) _TIFFmalloc(ROWSIZE);
	if(image[y] == NULL) {
	    fprintf(stderr, "Failed to alloc memory for image at row %d\n", y) ;
	    exit(1) ;
	}
    }

    read_tif_image(image,tif,input_H,SPP, TIF_CONFIG) ;


    if(quick_test) {
	//scale image by 1/2
	fprintf(stderr, "Scale image (%d,%d) by half using nearest pixel method to (%d,%d)\n", input_W,input_H,IMAGE_WIDTH,IMAGE_HEIGHT) ;

	for(x = 0 ; x < IMAGE_WIDTH  ; x++) {
	    int x2 = x*2 ;
	    if(x2 > input_W) {
		fprintf(stderr, "uh-oh, x:%d * 2 (%d) > input width %d\n", x, x*2, input_W) ;
		exit(1) ;
	    }


	    for(y = 0 ; y < IMAGE_HEIGHT ; y++) {
		int y2 = y*2 ;
		if(y2 > input_H) {
		    fprintf(stderr, "uh-oh, y:%d * 2 (%d) > input height %d\n", y, y*2, input_H) ;
		    exit(1) ;
		}
/*  		fprintf(stderr, "y:%d * 2 (%d) > input height %d\n", y, y*2, input_H) ;  */
		uint16_t r,g,b,a ;
		tif_get4c(image,x2,y2,r,g,b,a) ;
		tif_set4c(image,x,y,r,g,b,a) ;
	    }
	}

    }

    // input image is read, and possibly scaled by 1/2, ready for image processing

    fprintf(stderr, "repair alpha\n") ;
    repair_alpha(image, pData) ;

    /* find start of real image data */
    pData->raw_start_of_sky = (int16_t *)calloc(w, sizeof(int16_t *)) ;
    pData->start_of_sky = (int16_t *)calloc(w, sizeof(int16_t *)) ;
    pData->end_of_sky = (int16_t *)calloc(w, sizeof(int16_t *)) ;
    pData->final_end_of_sky = (int16_t *)calloc(w, sizeof(int16_t *)) ;


    simple_find_start_of_sky(pData->start_of_sky,w,h,image,pData) ;

    for(x = 0 ; x < w ; x++) {
	pData->raw_start_of_sky[x] = pData->start_of_sky[x] ;
    }

    find_start_of_sky(pData->start_of_sky,w,h,image,pData->fix_SOS_edges,pData) ;

    int fix_slivers=1 ;
    int fix_edges=1 ;
    find_end_of_sky(pData->end_of_sky,pData->start_of_sky,w,h,image,fix_edges,fix_slivers,pData) ;

    repair_top_of_sky(image, 1, 1, pData) ;  // 1-> end of sky is known, 1 means phase 1
    if(debug == 1) goto writeout ;

    find_start_of_sky(pData->start_of_sky,w,h,image,pData->fix_SOS_edges,pData) ;
    find_end_of_sky(pData->end_of_sky,pData->start_of_sky,w,h,image,fix_edges,fix_slivers,pData) ;
    if(debug == 2) goto writeout ;

    set_minmax_sky_values(pData) ;

    if(pData->depth_of_fill_is_absolute) {
	float length = pData->min_end_of_sky - pData->max_start_of_sky ;
	pData->depth_of_fill_absolute_y = (int)((float)pData->max_start_of_sky + depth_of_fill*length+0.5) ;
    }

    set_FOV_factor() ;
    set_xy_constraints(opt_parms) ;

    int n_samples = sample_sky_points(10, 200,image,pData->start_of_sky,pData->end_of_sky, pData) ;
    fprintf(stderr, "sampled %d points in image\n", n_samples) ;

    if(pData->fix_sky_hue) {
	// mainly corrects for blown out regions where hue is messed up,
	// ** Note, the starting hue_sky is used as estimate of actual hue
	repair_sky_hue(image,pData->start_of_sky,pData->end_of_sky,pData) ;

	// resample points with corrected sky hue
	// Not going to do this -- the sample points are already constrained to
	// not include points in the area of the sky identified as needing hue repaired
/*  	fix_sky_hue=0 ;  */
/*  	n_samples = sample_sky_points(10, 200,image,start_of_sky,end_of_sky) ;  */
    }

    for(int phase = 3 ; phase <= MAX_PHASES ; phase++) {
	repair_top_of_sky(image, 1, phase, pData) ;  // 1 means end of sky is known
	if(debug == phase) goto writeout ;
    }

    if(make_sky_uniform_flag == 2) make_sky_uniform(image,make_sky_uniform_flag,0,IMAGE_WIDTH-1, pData) ;

    if(debug == MAX_PHASES+1) goto writeout ;

    fit_sky_modelv2(force_grid_optimization, n_custom_optimizations, optimization_vars, n_samples, pData,opt_parms) ;

estimate_sky:
    pData->maximum_CIE_vhat = find_maximum_CIE_vhat() ;
    pData->model_is_being_fit=0 ;

    estimate_sky(0,IMAGE_WIDTH-1,image,pData->start_of_sky,pData->end_of_sky,pData->final_end_of_sky,depth_of_sample,h,w,feather_factor,debug,depth_of_fill,extra_feather_length,pData) ;

    fprintf(stderr, "back in main\n") ;

    if(pData->show_raw_error || pData->show_raw_prediction || pData->estimate_only || pData->show_sky_prediction)
	goto writeout ;

    int min_sky=h ;

    if(interpolate_masked_columns || run_hf) {
	// find minimum end of sky for all x
	for(x = 0 ; x < w ; x++) {
	    if(pData->column_mask[x] == 1) continue ;
	    if(min_sky > pData->end_of_sky[x]) min_sky = pData->end_of_sky[x] ;
	}
	fprintf(stderr, "min sky is %d\n", min_sky) ;
    }

    if(interpolate_masked_columns) {
	fprintf(stderr, "main:interpolate masked columns\n") ;
	for(int m = 0 ; m < n_masks ; m++) {
	    if(mask_l[m] > 0 && mask_r[m] < w-1) {
		int32_t x0 = mask_l[m]-1 ;
		int32_t x1 = mask_r[m]+1 ;

		int max_y = pData->start_of_sky[x0] ;

		for(x = x0+1 ; x < x1 ; x++) {
		    if(max_y < pData->start_of_sky[x]) max_y = pData->start_of_sky[x] ;
		}

		max_y += 10 ;

		for(y = 0 ; y < max_y ; y++) {
		    /* interpolate value from start to end */
		    float start_r, start_g, start_b ;
		    float   end_r, end_g, end_b ;

		    tif_get3c(image,x0,y,start_r,start_g,start_b) ;
		    tif_get3c(image,x1,y,end_r,end_g,end_b) ;

		    for(x = x0+1 ; x < x1 ; x++) {
			float pStart = 1.f - ((float)x-((float)x0))/((float)x1-(float)x0+1.f) ;
			float pEnd = 1.f - pStart ;
			tif_set4c(image,x,y, (uint16_t)(pStart*start_r + pEnd*end_r), (uint16_t)(pStart*start_g + pEnd*end_g), (uint16_t)(pStart*start_b + pEnd*end_b), MAX16) ;
		    }
		}

	    }
	}
    }

    if(make_sky_uniform_flag == 0 || make_sky_uniform_flag == 1) make_sky_uniform(image,make_sky_uniform_flag, 0,IMAGE_WIDTH-1, pData) ;

    if(run_hf) {
	// run a Horizontal filter on the sky to smooth edges
	fprintf(stderr, "main:horizontal filter\n") ;


	// alloc a buffer to hold output data
#define HBW 50
	float *buf = (float *) calloc(w, sizeof(float)) ;

	for(int ch=0 ; ch < 3 ; ch++) {
	    for(y = 0 ; y < min_sky ; y++) {
		// read values into buf
		for(x = 0 ; x < w ; x++) {
		    buf[x] = (float)((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+ch] ;
		}

		for(x=HBW ; x < w-HBW ; x++) {
		    float sum=0. ;
		    float n =0.;
		    int dx ;

		    for(dx=x-HBW ; dx < x+HBW ; dx++) {
			if(pData->column_mask[dx] == 1) continue ;
			sum += buf[dx] ;
			n += 1. ;
		    }
		    ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+ch] = (uint16_t)(sum/n) ;

		    if(x == HBW) {
			for(dx=0 ; dx < HBW ; dx++) {
			    if(pData->column_mask[dx] == 1) continue ;
			    ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*dx+ch] = (uint16_t)(sum/(float)(HBW*2+1)) ;
			}
		    } else if(x == w-HBW-1) {
			for(dx=x ; dx < x+HBW ; dx++) {
			    if(pData->column_mask[dx] == 1) continue ;
			    ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*dx+ch] = (uint16_t)(sum/(float)(HBW*2+1)) ;
			}
		    } else {
			if(pData->column_mask[x] == 1) continue ;
			((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+ch] = (uint16_t)(sum/(float)(HBW*2+1)) ;
		    }

		}

	    }
	}

	free(buf) ;
    }


    if(seam_line_smooth_size > 1) {
	// this needs to be in {}'s because the "goto writeout" isn't happy when
	// other variables are declared in the scope between the goto and the label
	uint16_t **output_buf = (uint16_t **) calloc(IMAGE_HEIGHT, sizeof(uint16_t *)) ;

	// alloc a buffer to hold output data
	for(i = 0 ; i < IMAGE_HEIGHT ; i++) {
	    output_buf[i] = (uint16_t *) calloc(IMAGE_WIDTH, sizeof(uint16_t)) ;
	}

	for(int ch=0 ; ch < 3 ; ch++) {
	    int ox, oy ;

	    for(ox = 0 ; ox < IMAGE_WIDTH ; ox++) {
		if(pData->column_mask[ox] == 1) continue ;

		for(oy = pData->start_of_sky[ox]-seam_line_smooth_size ; oy < pData->start_of_sky[x]+seam_line_smooth_size ; oy++) {
		    if(oy > pData->end_of_sky[ox]) break ;

		    float sum=0. ;
		    float n =0.;

		    for(x=ox-seam_line_smooth_size ; x <= ox+seam_line_smooth_size ; x++) {
			if(pData->column_mask[x] == 1) continue ;
			if(x < 0) continue ;
			if(x > IMAGE_WIDTH-1) continue ;

			for(y=oy-seam_line_smooth_size ; y <= oy+seam_line_smooth_size ; y++) {

			    if(y < 0) continue ;
			    if(y > pData->end_of_sky[x]) break ;

			    sum += ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+ch] ;
			    n += 1. ;
			}
		    }

		    output_buf[oy][ox] = (uint16_t)(sum/n+0.5) ;
		}
	    }

	    for(ox = 0 ; ox < IMAGE_WIDTH ; ox++) {
		if(pData->column_mask[ox] == 1) continue ;

		for(oy = pData->start_of_sky[ox]-seam_line_smooth_size ; oy < pData->start_of_sky[x]+seam_line_smooth_size ; oy++) {
		    if(oy > pData->end_of_sky[ox]) break ;

		    ((uint16_t *)(image[oy]))[IMAGE_NSAMPLES*ox+ch] = output_buf[oy][ox] ;
		}
	    }
	}


	for(i = 0 ; i <= pData->max_end_of_sky ; i++) {
	    free(output_buf[i]) ;
	}
	free(output_buf) ;
    }

#define DBW 4
#define DBH 6

    if(pData->full_sky_replacement==0) {
	int filter_passes=2 ;

	// this needs to be in {}'s because the "goto writeout" isn't happy when
	// other variables are declared in the scope between the goto and the label
	uint16_t **output_buf = (uint16_t **) calloc(IMAGE_HEIGHT, sizeof(uint16_t *)) ;
	float save_depth_of_fill=depth_of_fill ;
	depth_of_fill = 1. ;
	//extra_feather_length=5 ;

	// run a Box filter on the sky to smooth edges
	if(filter_passes > 0) {

	    fprintf(stderr, "main:box filter size %d X %d\n", DBW, DBH) ;
	    int i ;

	    // alloc a buffer to hold output data
	    for(i = 0 ; i < IMAGE_HEIGHT ; i++) {
		output_buf[i] = (uint16_t *) calloc(IMAGE_WIDTH, sizeof(uint16_t)) ;
	    }
	}

	while(filter_passes > 0) {
	    int i ;
	    filter_passes-- ;

	    fprintf(stderr, "main:box filter\n") ;

	    for(int ch=0 ; ch < 3 ; ch++) {
		int ox, oy ;

		for(ox = 0 ; ox < IMAGE_WIDTH ; ox++) {
		    if(pData->column_mask[ox] == 1) continue ;
		    int feather_end_y ;
		    int feather_length = raw_compute_feather_length(ox, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length, pData) ;
		    int maxy = pData->start_of_sky[ox]+feather_length ;

		    if(maxy < pData->raw_start_of_sky[ox]) maxy = pData->raw_start_of_sky[ox] ;
		    if(maxy > IMAGE_HEIGHT-1) maxy = IMAGE_HEIGHT-1 ;

		    if(ox == 0 && ch == 0) fprintf(stderr, "main:box filter, start y loop\n") ;

		    for(oy = 0 ; oy < maxy ; oy++) {

			float sum=0. ;
			float n =0.;

			for(x=ox-DBW ; x <= ox+DBW ; x++) {
			    if(pData->column_mask[x] == 1) continue ;
			    if(x < 0) continue ;
			    if(x > IMAGE_WIDTH-1) continue ;

			    for(y=oy-DBH ; y <= oy+DBH ; y++) {
				if(y < 0) continue ;
				if(y > maxy) break ;
				if(y > IMAGE_HEIGHT-1) break ;
				sum += ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+ch] ;
				n += 1. ;
			    }
			}

			output_buf[oy][ox] = (uint16_t)(sum/n+0.5) ;
		    }
		}

		for(ox = 0 ; ox < IMAGE_WIDTH ; ox++) {
		    if(pData->column_mask[ox] == 1) continue ;
		    int feather_end_y ;
		    int feather_length = raw_compute_feather_length(ox, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length, pData) ;
		    int maxy = pData->start_of_sky[ox]+feather_length ;

		    if(maxy < pData->raw_start_of_sky[ox]) maxy = pData->raw_start_of_sky[ox] ;
		    if(maxy > IMAGE_HEIGHT-1) maxy = IMAGE_HEIGHT-1 ;
		    if(ox == 0 && ch == 0) fprintf(stderr, "main:box filter, start second y loop\n") ;

		    for(oy = 0 ; oy < maxy ; oy++) {

			uint16_t b = ((uint16_t *)(image[oy]))[IMAGE_NSAMPLES*ox+2] ;
			float fp0 = raw_compute_feather_at_y(ox, oy, feather_length, b, feather_factor, pData) ;
			uint16_t new = output_buf[oy][ox] ;
			uint16_t old = ((uint16_t *)(image[oy]))[IMAGE_NSAMPLES*ox+ch] ;
			((uint16_t *)(image[oy]))[IMAGE_NSAMPLES*ox+ch] = (uint16_t)(fp0*(float)new + (1.-fp0)*(float)old + 0.5) ;

/*  			uint16_t old = ((uint16_t *)(image[oy]))[IMAGE_NSAMPLES*ox+ch] = output_buf[oy][ox] ;  */
		    }
		}
	    }


	    if(filter_passes == 0) {
		for(i = 0 ; i <= pData->max_end_of_sky ; i++) {
		    free(output_buf[i]) ;
		}
		free(output_buf) ;
	    }

	}

	depth_of_fill = save_depth_of_fill ;
    }


    int dither_sky=1 ;
    if(dither_sky) {
	fprintf(stderr, "main:dithering the sky\n") ;
	srand(0) ;
	/* dither sky */
	for(x = 0 ; x < IMAGE_WIDTH  ; x++) {
	    if(pData->column_mask[x] == 1) continue ;

	    float sky_height = pData->end_of_sky[x] - pData->start_of_sky[x] ;
	    int feather_end_y = pData->start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) ;

	    for(y = 0 ; y < feather_end_y ; y++) {
		for(int ch=0 ; ch < 3 ; ch++) {
		    float f16 = (float)(((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+ch]) ;
		    float dither = (float)(((rand()%128) - 64)*8) ;
		    f16 += dither ;
		    if(f16 < 0.0) f16 = 0.0 ;
		    if(f16 > MAX16f-1.) f16 = MAX16f-1. ;
		    ((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+ch] = (uint16_t)f16 ;
		}

#ifdef OLD
		// red
		int16_t dither = ((rand()%128) - 64)*8 ;
		((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+0] += dither ;

		// green
		dither = ((rand()%128) - 64)*8 ;
		((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+1] += dither ;

		// blue
		dither = ((rand()%128) - 64)*8 ;
		((uint16_t *)(image[y]))[IMAGE_NSAMPLES*x+2] += dither ;
#endif
	    }
	}
    }

writeout:
    if(debug) {
	for(x = 0 ; x < IMAGE_WIDTH  ; x++) {

	    // indicate manual end of sky as red
	    y = pData->manual_end_of_sky[x] ;
	    if(y >= 0 && y <h) {
		int yset = y+2 ; // set the color of the pixel two below EOS to green, so you can see the transition pixels above

		if(yset > h-1) yset = h-1 ;

		tif_set4c(image,x,yset,MAX16,0,0,MAX16) ;
	    }

	    if(pData->column_mask[x] == 1) continue ;

	    float sky_height = pData->end_of_sky[x] - pData->start_of_sky[x] ;
	    int feather_end_y = pData->start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) ;

	    // indicate start of sky as red line
	    for(y = pData->start_of_sky[x] ; y < pData->start_of_sky[x]+1 ; y++) {
		if(y < 0 || y >=h) {
		    fprintf(stderr, "debug: sos y=%d\n", y) ;
		    exit(1) ;
		}
		if(y > 0)
		    tif_set4c(image,x,y-1,MAX16,0,0,MAX16) ;
	    }

	    // indicate raw_start of sky as yellow line
	    for(y = pData->raw_start_of_sky[x] ; y < pData->raw_start_of_sky[x]+1 ; y++) {
		if(y < 0 || y >=h) {
		    fprintf(stderr, "debug: sos y=%d\n", y) ;
		    exit(1) ;
		}

		if(y > 0)
		    tif_set4c(image,x,y-1,HALF16,HALF16,0,MAX16) ;
	    }

	    if(debug > 3) {
		// indicate end of fill as grey
		for(y = feather_end_y ; y < feather_end_y+1 ; y++) {
		    if(y < 0 || y >=h) {
			fprintf(stderr, "debug: feather y=%d\n", y) ;
			exit(1) ;
		    }

		    tif_set4c(image,x,y,HALF16,HALF16,HALF16,MAX16) ;
		}
	    }

	    // indicate end of detected sky as green
	    for(y = pData->end_of_sky[x] ; y < pData->end_of_sky[x]+1 ; y++) {
		if(y < 0 || y >=h) {
		    fprintf(stderr, "debug: eos x=%d, y=%d\n", x, y) ;
		    exit(1) ;
		}

		int yset = y+2 ; // set the color of the pixel two below EOS to green, so you can see the transition pixels above

		if(yset > h-1) yset = h-1 ;

		tif_set4c(image,x,yset,0,MAX16,0,MAX16) ;
	    }
	}

	if(samples != NULL) {
	    // mark location of sky samples with black pixels
	    for(int i=0 ; i < n_samples ; i++) {
		int xc = samples[i].x ;
		int yc = samples[i].y ;
/*  		fprintf(stderr, "sample %d, x,y:%d,%d\n", i,xc,yc) ;  */
		float fact = samples[i].abs_v_error[0] * 10. ;
		if(fact > 1.) fact = 1. ;
		uint16_t f = (uint16_t)(fact * (float)MAX16+0.5) ;
		tif_set4c(image,xc,yc,f,f,f,MAX16) ;
	    }
	}
    }



    fprintf(stderr, "main:save\n") ;
    fprintf(stderr, "sky detect tolerances (hsv) are %lf, %lf, %lf\n", pData->sky_hue_tolerance, pData->sky_sat_tolerance, pData->sky_val_tolerance) ;
/*      fprintf(stderr, "sample_based_v_correction is %lf\n", sample_based_v_correction) ;  */
    printf("Maximum CIE vhat=%f\n", pData->maximum_CIE_vhat) ;
    printf("Minimum CIE vhat=%f\n", pData->minimum_CIE_vhat) ;

    pData->verbose=verbose_flag ;

    if(pData->verbose) {
	//some useful debugging information
	{
	    float sun_px = sun_x_angle2px(pData->sun_x) ;
	    float gamma,theta,cos_gamma ;
	    printf("horizon_py:%f, sun_x_angle:%f sun(px,py):(%6.2f,%6.2f)\n\n", pData->horizon_py, pData->sun_x, sun_px, pData->sun_py) ;
	    printf("======== image_relative pixel_to_V3 for sun ============\n") ;
	    F_CIE2003(sun_px, pData->sun_py, &gamma, &theta, &cos_gamma) ;
	    printf("======== image_relative pixel_to_V3 for sun ============\n\n") ;
	}

	float px, py ;
	py = pData->sun_py ;
	for(px = -0.5 ; px < 0.51 ; px += .5) {
	    float gamma,theta,cos_gamma ;
	    float vhat = F_CIE2003(px, py, &gamma, &theta, &cos_gamma) ;
	    if(0) {
		printf("vhat:%f\n", vhat) ;
	    }
	}
    }

    srand(0) ;
    /* save image data */
    for(y = 0 ; y < h ; y++) {
/*  	if(! (y%200)) {  */
/*  	    fprintf(stderr, "main:save row:%d\n", (int)y) ;  */
/*  	}  */
	if(TIF_CONFIG == PLANARCONFIG_CONTIG) {
	    TIFFWriteScanline(tifout,image[y],y,0) ;
	} else if(TIF_CONFIG == PLANARCONFIG_SEPARATE) {
	    uint16_t s ;
	    for(s = 0 ; s < SPP ; s++)
		TIFFWriteScanline(tifout,image[y],y,s) ;
	} else {
	    fprintf(stderr, "Cannot read planar configuration of tiff image, config=%d!\n", (int)TIF_CONFIG) ;
	    exit(1) ;
	}
    }

    TIFFClose(tif);
    TIFFClose(tifout);

    fprintf(stderr, "free-ing image memory\n") ;

    /* free memory for image */

    for(y = 0 ; y < h ; y++) {
	_TIFFfree(image[y]) ;
    }

    fprintf(stderr, "free-ing image\n") ;
    free(image) ;
    free(pData->raw_start_of_sky) ;
    free(pData->start_of_sky) ;
    free(pData->end_of_sky) ;
    free(pData->manual_end_of_sky) ;
    free(pData->final_end_of_sky) ;
    free(pData->column_mask) ;
    free(pData->column_repair_mask) ;
    free(pData->column_nonsky_repair_mask) ;
    free(pData->column_sample_mask) ;
    if(samples != NULL)
	free(samples) ;

    // finally, copy exif tags from original to output image
    // The ICC profile seems to only be recognized in exif tags, not tiff tags, by at least
    // Geegie, and the linux image viewer "xviewer"
    if(use_exiftool) {
	char cmd[1000] ;
	if(0) {
	    //this doesn't work, maybe a bug in exiftool, it seems to create a malformed
	    //tag
/*  	    snprintf(cmd,999,"exiftool -tagsFromFile %s %s\n", infilename, outfilename) ;  */
/*  	    printf("%s", cmd) ;  */
/*  	    system(cmd) ;  */
	} else {
	    // extract the ICC profile to external file
	    snprintf(cmd,999,"exiftool -icc_profile -b %s > skyfill_tmp.icc\n", infilename) ;
	    printf("%s", cmd) ;
	    int retval = system(cmd) ;

	    // add the ICC profile from the external file to the new tiff image
	    snprintf(cmd,999,"exiftool \"-icc_profile<=skyfill_tmp.icc\" %s\n", outfilename) ;
	    printf("%s", cmd) ;
	    retval = system(cmd) ;

	    // system command to delete the file
	    remove("skyfill_tmp.icc") ;
	}
    }

    exit(0) ;
}
