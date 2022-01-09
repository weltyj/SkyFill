#ifndef SAMPLE_AND_FIT_SKY_MODEL_H
#define SAMPLE_AND_FIT_SKY_MODEL_H
struct sample_point {
    float px, py ;
    float h,s,v ;
    float v_hat ;
    float abs_error[2] ; // absolute error from model of r+g+b, by location index
    float abs_h_error[2] ; // absolute error from model of hue, by location index
    float abs_s_error[2] ; // absolute error from model of sat, by location index
    float abs_v_error[2] ; // absolute error from model of val, by location index
    uint16_t r,g,b,x,y ;
} ;

// these have to be global
// because the nonlinear fitting function doesn't pass data, only
// parameters to the function to be minimized
extern struct sample_point *samples ;
extern int n_samples_to_optimize ;

#define FULL_RAW_HSV_MODEL

// coefs to model sky HSV from input file sky area
// model is H|S|V = B0 + B1*px + B2*py + B3*px*py ;
struct HSV_model_coefs {
    float h_coefs[9];
    float s_coefs[9];
    float v_coefs[9];
    } ;

extern struct HSV_model_coefs raw_hsv_coefs[2] ;

// weighting for all points when fitting model for single function for entire sky
static inline float px_wgt_single(float px, float abs_error)
{
    return 1./(1.+abs_error)  ;
}

#define LOWER_PEAK
// weighting for all points when fitting model for left half of sky
static inline float px_wgt_left(float px, float abs_error)
{
#ifdef LOWER_PEAK
    px *= 1.175 ;
    px -= .25 ;
    if(px < 0.) px = 0. ;
    return (1.-3.*px*px+2.*px*px*px) / (1.+abs_error) ;
#else
    px = px*1.5 -.25;

    if(px < 0.) px = 0. ;
    if(px > 1.) px = 1. ;
    float lr_wgt = (1.-3.*px*px+2.*px*px*px) ;
    return lr_wgt*lr_wgt*lr_wgt*lr_wgt / (1.+abs_error) ;
#endif
}

// weighting for all points when fitting model for right half of sky
static inline float px_wgt_right(float px, float abs_error)
{
#ifdef LOWER_PEAK
    px = (1.-px)*1.175 -.25;
    if(px < 0.) px = 0. ;
    return (1.-3.*px*px+2.*px*px*px) / (1.+abs_error) ;
#else
    px = (1.-px)*1.5 -.25;

    if(px < 0.) px = 0. ;
    if(px > 1.) px = 1. ;
    float lr_wgt = (1.-3.*px*px+2.*px*px*px) ;
    return lr_wgt*lr_wgt*lr_wgt*lr_wgt / (1.+abs_error) ;
#endif
}


static inline float local_h_from_pxpy(float px, float py_hc, float *coefs)
{
    float half=coefs[0] ;
    float spread=coefs[1] ;
    float speed=coefs[2] ;
    float center=coefs[3] ;
    return spread*(py_hc-center)*speed/(1.+fabs(py_hc-center)*speed) + half ;

    //return p->h_coefs[0] + p->h_coefs[1]*px + p->h_coefs[2]*log(0.+py_hc) ;
/*      return p->h_coefs[0]*exp(p->h_coefs[1]/(0.001+py_hc)) ;  */
}

/*  static inline float local_ty_from_pxpy(float px, float py_hc, struct HSV_model_coefs *p, float horizon_py)  */
/*  {  */
/*      py_hc -= horizon_py ;  */
/*    */
/*      return 1./(1.+exp(-py_hc*p->S_inv_coef)) ; // really a simple logistic centered on horizon_py  */
/*  }  */

static inline float local_s_from_pxpy(float px, float py_hc, struct HSV_model_coefs *p)
{
    // nonlinear model now, fit with amoeba (Nelder Meade)
    float angle_factor = p->s_coefs[5] ;
    float x_mapped = ((pData_fit->FOV_horizontal*px+pData_fit->valsat_phase_shift)*angle_factor-180.)
          *M_PI/180. ;
    float x_adjust= (p->s_coefs[6]-1./(cos(x_mapped)+ 1.+1./p->s_coefs[6]))/10. ;

    return p->s_coefs[0]/(1.+exp(-(py_hc-p->s_coefs[1])*p->s_coefs[2]))
	    +p->s_coefs[3] +p->s_coefs[4]*(py_hc-p->s_coefs[1]) + x_adjust ;
}

static inline float local_v_from_pxpy(float px, float py_hc, float *coefs)
{
#define VAL_INV_ADD 1.5
    float angle_factor = coefs[0]+coefs[5]*py_hc ;
    float x_mapped = ((pData_fit->FOV_horizontal*px+pData_fit->valsat_phase_shift)*angle_factor-180.)
                      * M_PI/180. ;
/*      return coefs[1]/(cos(x_mapped)*coefs[2]+VAL_INV_ADD)+coefs[3]+coefs[4]*py_hc + coefs[5]*py_hc*py_hc ;  */
    return coefs[1]/(cos(x_mapped)*coefs[2]+VAL_INV_ADD)
	    +coefs[4]*py_hc
	    +coefs[3] ;
}

#define hsv_blended(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,&raw_hsv_coefs[0]) ;\
    float rvalue = f(px,py,&raw_hsv_coefs[1]) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }

#define v_blended(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,raw_hsv_coefs[0].v_coefs) ;\
    float rvalue = f(px,py,raw_hsv_coefs[1].v_coefs) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }

#define h_blended(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,raw_hsv_coefs[0].h_coefs) ;\
    float rvalue = f(px,py,raw_hsv_coefs[1].h_coefs) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }

// saturation needs to pass horizon_py as an additional argument
#define hsv_blended_hpy(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,&raw_hsv_coefs[0]) ;\
    float rvalue = f(px,py,&raw_hsv_coefs[1]) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }


static inline float h_from_pxpy(float px, float py_hc, int full_hsv_model_level)
{
    if(full_hsv_model_level==1)
	return local_h_from_pxpy(px,py_hc,raw_hsv_coefs[0].h_coefs) ;

    float h ;
    h_blended(local_h_from_pxpy, px, py_hc, h) ;
    return h ;
}

static inline float s_from_pxpy(float px, float py_hc, int full_hsv_model_level, float horizon_py)
{
    if(full_hsv_model_level==1)
	return local_s_from_pxpy(px,py_hc,&raw_hsv_coefs[0]) ;

    float s ;
    hsv_blended_hpy(local_s_from_pxpy, px, py_hc, s) ;
    return s ;
}

static inline float v_from_pxpy(float px, float py_hc, int full_hsv_model_level)
{
    if(full_hsv_model_level==1)
	return local_v_from_pxpy(px,py_hc,raw_hsv_coefs[0].v_coefs) ;

    float v ;
    v_blended(local_v_from_pxpy, px, py_hc, v) ;
    return v ;
}

float fit_full_HSV_model(int n_samples, int location_index, int calc_error, SKYFILL_DATA_t *pData) ;

void reset_sample_errors(int n_samples) ;

int sample_sky_points(int n_per_column, int n_columns,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky, SKYFILL_DATA_t *pData) ;

#endif
