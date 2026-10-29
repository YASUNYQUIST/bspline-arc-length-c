#ifndef BSPLINE_H
#define BSPLINE_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define BS_MAX_DEGREE 16
#define BS_MAX_POINTS 1000000
typedef enum { BS_OK=0, BS_INVALID, BS_MEMORY, BS_LIMIT, BS_NOT_CONVERGED,
               BS_UNREACHABLE } BsStatus;
typedef struct { double x,y; } BsPoint;
typedef struct {
    BsPoint a,o,b,c,ab;
    double radius,direction,theta0,s_ab,s_circle,s_total;
} BsGeometry;
/* Owning structs must be zero-initialized; free before reusing as outputs. */
typedef struct { size_t count; int degree; BsPoint *points; double *knots; } BsCurve;
typedef struct {
    double d,actual_step,max_error;
    size_t evaluations;
    BsGeometry geometry;
    BsCurve curve;
} BsFit;
typedef struct {
    BsCurve derivative;
    double *breaks,*prefix;
    size_t spans;
    double tolerance,estimated_error;
    int n_start,n_max;
} BsArc;
typedef struct {
    double u,length,abs_error;
    BsPoint point;
    int iterations,bisection_steps;
} BsInverse;
const char *bs_status_string(BsStatus status);
BsStatus bs_geometry(BsPoint a,BsPoint o,double radius,int branch,BsGeometry *out);
BsStatus bs_sample(const BsGeometry *geo,double d,int degree,BsCurve *out,
                   double *actual_step);
BsStatus bs_knots(BsCurve *curve);
void bs_curve_free(BsCurve *curve);
BsStatus bs_eval(const BsCurve *curve,double u,BsPoint *out);
BsStatus bs_derivative(const BsCurve *curve,BsCurve *out);
double bs_theoretical_error(BsPoint point,const BsGeometry *geo);
BsStatus bs_max_error(const BsCurve *curve,const BsGeometry *geo,
                      size_t samples,double *out);
/* Compatibility bisection. Failure retains the last valid candidate in out. */
BsStatus bs_fit(const BsGeometry *geo,double err_min,double err_max,double d_init,
                int degree,size_t samples,int max_iter,double tol_d,BsFit *out);
void bs_fit_free(BsFit *fit);
/* Gauss n/2n differences are error estimates, not rigorous upper bounds. */
BsStatus bs_arc_init(const BsCurve *curve,double tolerance,int n_start,int n_max,
                     BsArc *out);
void bs_arc_free(BsArc *arc);
BsStatus bs_arc_cumulative(const BsArc *arc,double u,double *out);
BsStatus bs_arc_length(const BsArc *arc,double a,double b,double *out);
BsStatus bs_inverse(const BsCurve *curve,const BsArc *arc,double u0,double length,
                    double tolerance,int max_iter,BsInverse *out);
#ifdef __cplusplus
}
#endif
#endif
