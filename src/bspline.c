#include "bspline.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.141592653589793238462643383279502884
static BsPoint add(BsPoint a,BsPoint b) { return (BsPoint){a.x+b.x,a.y+b.y}; }
static BsPoint sub(BsPoint a,BsPoint b) { return (BsPoint){a.x-b.x,a.y-b.y}; }
static BsPoint mul(BsPoint a,double t) { return (BsPoint){a.x*t,a.y*t}; }
static double dot(BsPoint a,BsPoint b) { return a.x*b.x+a.y*b.y; }
static double norm(BsPoint a) { return hypot(a.x,a.y); }
static int finite_point(BsPoint a) { return isfinite(a.x)&&isfinite(a.y); }
const char *bs_status_string(BsStatus s) {
    const char *names[]={"ok","invalid input","allocation failure","resource limit",
                         "not converged","unreachable target"};
    return s>=BS_OK && s<=BS_UNREACHABLE ? names[s] : "unknown status";
}
void bs_curve_free(BsCurve *c) {
    if(c) { free(c->points); free(c->knots); memset(c,0,sizeof(*c)); }
}
void bs_fit_free(BsFit *f) { if(f) { bs_curve_free(&f->curve); memset(f,0,sizeof(*f)); } }
BsStatus bs_geometry(BsPoint a,BsPoint o,double r,int branch,BsGeometry *g) {
    if(!g||!finite_point(a)||!finite_point(o)||!isfinite(r)||r<=0) return BS_INVALID;
    BsPoint v=sub(a,o); double d=norm(v);
    if(!isfinite(d)||d<=r+1e-12) return BS_INVALID;
    BsPoint e=mul(v,1/d), n={-e.y,e.x};
    double ct=r/d, st=sqrt(fmax(0,1-ct*ct));
    BsPoint radial=mul(add(mul(e,ct),mul(n,branch>=0?st:-st)),r);
    memset(g,0,sizeof(*g));
    g->a=a; g->o=o; g->radius=r; g->b=add(o,radial); g->ab=sub(g->b,a);
    BsPoint ccw={-radial.y,radial.x};
    g->direction=dot(ccw,g->ab)<0 ? 1:-1;
    g->c=add(g->b,g->ab); g->theta0=atan2(radial.y,radial.x);
    g->s_ab=norm(g->ab); g->s_circle=2*PI*r; g->s_total=2*g->s_ab+g->s_circle;
    return isfinite(g->s_total)&&finite_point(g->c) ? BS_OK:BS_INVALID;
}
BsStatus bs_knots(BsCurve *c) {
    if(!c||!c->points||!c->knots||c->degree<1||c->degree>BS_MAX_DEGREE||
       c->count<(size_t)c->degree+1||c->count>BS_MAX_POINTS) return BS_INVALID;
    size_t count=c->count; int p=c->degree;
    double *t=calloc(count,sizeof(double));
    if(!t) return BS_MEMORY;
    for(size_t i=0;i<count;i++) {
        if(!finite_point(c->points[i])) { free(t); return BS_INVALID; }
        if(i) t[i]=t[i-1]+norm(sub(c->points[i],c->points[i-1]));
    }
    double total=t[count-1];
    if(!isfinite(total)) { free(t); return BS_INVALID; }
    for(size_t i=0;i<count;i++) t[i]=total<1e-12 ? (double)i/(count-1):t[i]/total;
    memset(c->knots,0,(count+p+1)*sizeof(double));
    for(size_t j=count;j<count+p+1;j++) c->knots[j]=1;
    for(size_t j=p+1;j<count;j++) {
        double sum=0; for(size_t i=j-p;i<j;i++) sum+=t[i];
        c->knots[j]=fmax(c->knots[j-1],sum/p);
    }
    free(t); return BS_OK;
}
BsStatus bs_sample(const BsGeometry *g,double d,int p,BsCurve *out,double *step) {
    if(!g||!out||!isfinite(d)||d<=0||p<1||p>BS_MAX_DEGREE||
       !isfinite(g->s_total)||g->s_total<=0||g->s_ab<=0||g->radius<=0) return BS_INVALID;
    double segments=fmax(1,ceil(g->s_total/d-1e-12));
    if(segments>BS_MAX_POINTS-1) return BS_LIMIT;
    size_t n=(size_t)segments;
    if(n<(size_t)p) return BS_INVALID;
    BsCurve c={n+1,p,NULL,NULL};
    c.points=calloc(c.count,sizeof(BsPoint)); c.knots=calloc(c.count+p+1,sizeof(double));
    if(!c.points||!c.knots) { bs_curve_free(&c); return BS_MEMORY; }
    double h=g->s_total/n;
    for(size_t i=0;i<=n;i++) {
        double s=fmin(i*h,g->s_total);
        if(s<=g->s_ab) c.points[i]=add(g->a,mul(g->ab,s/g->s_ab));
        else if(s<=g->s_ab+g->s_circle) {
            double angle=g->theta0+g->direction*(s-g->s_ab)/g->radius;
            c.points[i]=add(g->o,(BsPoint){g->radius*cos(angle),g->radius*sin(angle)});
        } else c.points[i]=add(g->b,mul(g->ab,(s-g->s_ab-g->s_circle)/g->s_ab));
    }
    BsStatus status=bs_knots(&c);
    if(status!=BS_OK) { bs_curve_free(&c); return status; }
    *out=c; if(step) *step=h; return BS_OK;
}
/* Public curves must have nondecreasing finite knots and a valid domain. */
static int valid_curve(const BsCurve *c) {
    return c&&c->points&&c->knots&&c->degree>=0&&c->degree<=BS_MAX_DEGREE&&
        c->count>(size_t)c->degree&&c->count<=BS_MAX_POINTS;
}
BsStatus bs_eval(const BsCurve *c,double u,BsPoint *out) {
    if(!out||!valid_curve(c)||!isfinite(u)) return BS_INVALID;
    int p=c->degree; const double *U=c->knots; size_t n=c->count-1,k;
    if(u<U[p]||u>U[n+1]) return BS_INVALID;
    if(u>=U[n+1]) k=n;
    else if(u<=U[p]) k=(size_t)p;
    else {
        size_t lo=p,hi=n+1; k=(lo+hi)/2;
        while(u<U[k]||u>=U[k+1]) {
            if(u<U[k]) hi=k; else lo=k;
            k=(lo+hi)/2;
        }
    }
    BsPoint work[BS_MAX_DEGREE+1];
    for(int j=0;j<=p;j++) work[j]=c->points[k-p+j];
    for(int r=1;r<=p;r++) for(int j=p;j>=r;j--) {
        size_t idx=k-p+j; double denom=U[idx+p-r+1]-U[idx];
        double alpha=denom>1e-12 ? (u-U[idx])/denom:0;
        work[j]=add(mul(work[j-1],1-alpha),mul(work[j],alpha));
    }
    *out=work[p]; return finite_point(*out)?BS_OK:BS_INVALID;
}
BsStatus bs_derivative(const BsCurve *c,BsCurve *out) {
    if(!out||!valid_curve(c)||c->degree==0) return BS_INVALID;
    BsCurve q={c->count-1,c->degree-1,NULL,NULL};
    q.points=calloc(q.count,sizeof(BsPoint)); q.knots=calloc(q.count+q.degree+1,sizeof(double));
    if(!q.points||!q.knots) { bs_curve_free(&q); return BS_MEMORY; }
    for(size_t i=0;i<q.count;i++) {
        double denom=c->knots[i+c->degree+1]-c->knots[i+1];
        if(denom>1e-12) q.points[i]=mul(sub(c->points[i+1],c->points[i]),c->degree/denom);
    }
    memcpy(q.knots,c->knots+1,(q.count+q.degree+1)*sizeof(double));
    *out=q; return BS_OK;
}
double bs_theoretical_error(BsPoint point,const BsGeometry *g) {
    BsPoint v=sub(g->c,g->a),w=sub(point,g->a); double denom=dot(v,v);
    double line=denom<1e-15?norm(w):norm(sub(w,mul(v,dot(w,v)/denom)));
    return fmin(line,fabs(norm(sub(point,g->o))-g->radius));
}
BsStatus bs_max_error(const BsCurve *c,const BsGeometry *g,size_t samples,double *out) {
    if(!valid_curve(c)||!g||!out||samples<2) return BS_INVALID;
    double result=0,lo=c->knots[c->degree],hi=c->knots[c->count];
    for(size_t i=0;i<samples;i++) {
        BsPoint pt; BsStatus s=bs_eval(c,i==samples-1?hi:lo+(hi-lo)*i/(samples-1),&pt);
        if(s!=BS_OK) return s;
        result=fmax(result,bs_theoretical_error(pt,g));
    }
    *out=result; return BS_OK;
}
static BsStatus candidate(const BsGeometry *g,double d,int p,size_t samples,BsFit *out) {
    BsFit f={0}; f.geometry=*g; f.d=d;
    BsStatus s=bs_sample(g,d,p,&f.curve,&f.actual_step);
    if(s==BS_OK) s=bs_max_error(&f.curve,g,samples,&f.max_error);
    if(s!=BS_OK) { bs_fit_free(&f); return s; }
    size_t evaluations=out->evaluations+1;
    bs_fit_free(out); *out=f; out->evaluations=evaluations; return BS_OK;
}
BsStatus bs_fit(const BsGeometry *g,double emin,double emax,double d,int p,
                size_t samples,int max_iter,double tol,BsFit *out) {
    if(!out||!g||!isfinite(emin)||!isfinite(emax)||emin<=0||emax<emin||
       !isfinite(tol)||tol<=0||max_iter<=0||samples<2) return BS_INVALID;
    BsStatus s=candidate(g,d,p,samples,out); if(s!=BS_OK) return s;
    if(out->max_error>=emin&&out->max_error<=emax) return BS_OK;
    double lo,hi; int bracketed=0;
    int shrinking=out->max_error>emax;
    if(shrinking) { hi=d; lo=d/2; } else { lo=d; hi=d*2; }
    for(int i=0;i<64;i++) {
        s=candidate(g,shrinking?lo:hi,p,samples,out);
        if(s!=BS_OK) return s==BS_INVALID?BS_NOT_CONVERGED:s;
        if((shrinking&&out->max_error<=emax)||(!shrinking&&out->max_error>=emin)) {
            bracketed=1; break;
        }
        if(shrinking) { hi=lo; lo/=2; } else { lo=hi; hi*=2; }
    }
    if(!bracketed) return BS_NOT_CONVERGED;
    for(int i=0;i<max_iter;i++) {
        double mid=0.5*(lo+hi);
        s=candidate(g,mid,p,samples,out);
        if(s!=BS_OK) return s;
        if(out->max_error>=emin&&out->max_error<=emax) return BS_OK;
        if(out->max_error>emax) hi=mid; else lo=mid;
        if(hi-lo<tol) break;
    }
    return BS_NOT_CONVERGED;
}
/* Roots of Legendre P_n by Newton; pairwise symmetry reduces work. */
static BsStatus gauss(const BsCurve *derivative,double a,double b,int n,double *out) {
    double sum=0,half=(b-a)/2,center=(a+b)/2;
    for(int i=0;i<(n+1)/2;i++) {
        double z=cos(PI*(i+0.75)/(n+0.5)),dp=0;
        int converged=0;
        for(int iter=0;iter<64;iter++) {
            double p1=1,p2=0;
            for(int j=1;j<=n;j++) { double p3=p2; p2=p1; p1=((2*j-1)*z*p2-(j-1)*p3)/j; }
            dp=n*(z*p1-p2)/(z*z-1);
            double next=z-p1/dp;
            if(fabs(next-z)<=4*DBL_EPSILON) { z=next; converged=1; break; }
            z=next;
        }
        if(!converged) return BS_NOT_CONVERGED;
        /* Recompute derivative at the accepted root for accurate weights. */
        double p1=1,p2=0;
        for(int j=1;j<=n;j++) { double p3=p2; p2=p1; p1=((2*j-1)*z*p2-(j-1)*p3)/j; }
        dp=n*(z*p1-p2)/(z*z-1);
        double weight=2/((1-z*z)*dp*dp);
        BsPoint left,right;
        BsStatus s=bs_eval(derivative,center-half*z,&left); if(s!=BS_OK) return s;
        double value=norm(left);
        if(i!=n-1-i) {
            s=bs_eval(derivative,center+half*z,&right); if(s!=BS_OK) return s;
            value+=norm(right);
        }
        sum+=weight*value;
    }
    *out=half*sum; return isfinite(*out)?BS_OK:BS_INVALID;
}
static BsStatus integrate_recursive(const BsArc *arc,double a,double b,double tolerance,
                                    double *out,double *estimate,int depth) {
    if(b<=a) { *out=0; *estimate=0; return BS_OK; }
    int n=arc->n_start; double prev=0,curr=0;
    BsStatus s=gauss(&arc->derivative,a,b,n,&prev); if(s!=BS_OK) return s;
    while(n<arc->n_max) {
        int next=n>arc->n_max/2?arc->n_max:n*2;
        s=gauss(&arc->derivative,a,b,next,&curr); if(s!=BS_OK) return s;
        *estimate=fabs(curr-prev); *out=curr;
        if(*estimate<tolerance) return BS_OK;
        prev=curr; n=next;
    }
    /* Reversals can put a zero of speed inside a span. Increasing GL order
       alone then converges slowly; subdivide while retaining the error budget. */
    double mid=a+(b-a)/2;
    if(depth>=24||mid==a||mid==b) return BS_NOT_CONVERGED;
    double left,right,le,re;
    s=integrate_recursive(arc,a,mid,tolerance/2,&left,&le,depth+1);
    if(s!=BS_OK) return s;
    s=integrate_recursive(arc,mid,b,tolerance/2,&right,&re,depth+1);
    if(s!=BS_OK) return s;
    *out=left+right; *estimate=le+re; return BS_OK;
}
static BsStatus integrate(const BsArc *arc,double a,double b,double tolerance,
                           double *out,double *estimate) {
    return integrate_recursive(arc,a,b,tolerance,out,estimate,0);
}
void bs_arc_free(BsArc *a) {
    if(a) { bs_curve_free(&a->derivative); free(a->breaks); free(a->prefix); memset(a,0,sizeof(*a)); }
}
BsStatus bs_arc_init(const BsCurve *c,double tol,int n_start,int n_max,BsArc *out) {
    if(!out||!valid_curve(c)||c->degree<1||!isfinite(tol)||tol<=0||
       n_start<1||n_max<=n_start||n_max>256) return BS_INVALID;
    /* Validate once; bs_eval is deliberately O(log N + p^2). */
    for(size_t i=0;i<c->count+c->degree+1;i++)
        if(!isfinite(c->knots[i])||(i&&c->knots[i]<c->knots[i-1])) return BS_INVALID;
    if(c->knots[c->count]<=c->knots[c->degree]) return BS_INVALID;
    for(size_t i=0;i<c->count;i++) if(!finite_point(c->points[i])) return BS_INVALID;
    BsArc a={0}; a.tolerance=tol; a.n_start=n_start; a.n_max=n_max;
    BsStatus s=bs_derivative(c,&a.derivative); if(s!=BS_OK) return s;
    a.breaks=calloc(c->count+1,sizeof(double)); a.prefix=calloc(c->count+1,sizeof(double));
    if(!a.breaks||!a.prefix) { bs_arc_free(&a); return BS_MEMORY; }
    a.breaks[0]=c->knots[c->degree];
    for(size_t i=c->degree+1;i<=c->count;i++)
        if(c->knots[i]>a.breaks[a.spans]) a.breaks[++a.spans]=c->knots[i];
    for(size_t i=0;i<a.spans;i++) {
        double value=0,error=0;
        /* Reserve half the budget for partial-span queries. */
        s=integrate(&a,a.breaks[i],a.breaks[i+1],tol/(2*a.spans),&value,&error);
        if(s!=BS_OK) { bs_arc_free(&a); return s; }
        a.prefix[i+1]=a.prefix[i]+value; a.estimated_error+=error;
    }
    *out=a; return BS_OK;
}
BsStatus bs_arc_cumulative(const BsArc *a,double u,double *out) {
    if(!a||!a->breaks||!a->prefix||!a->spans||!out||!isfinite(u)||
       u<a->breaks[0]||u>a->breaks[a->spans]) return BS_INVALID;
    if(u==a->breaks[a->spans]) { *out=a->prefix[a->spans]; return BS_OK; }
    size_t lo=0,hi=a->spans;
    while(hi-lo>1) { size_t mid=(lo+hi)/2; if(a->breaks[mid]<=u) lo=mid; else hi=mid; }
    double partial=0,error=0;
    BsStatus s=integrate(a,a->breaks[lo],u,a->tolerance/2,&partial,&error);
    if(s==BS_OK) *out=a->prefix[lo]+partial;
    return s;
}
BsStatus bs_arc_length(const BsArc *a,double u,double v,double *out) {
    if(!out) return BS_INVALID;
    double first,last; BsStatus s=bs_arc_cumulative(a,u,&first);
    if(s==BS_OK) s=bs_arc_cumulative(a,v,&last);
    if(s==BS_OK) *out=fabs(last-first);
    return s;
}
BsStatus bs_inverse(const BsCurve *c,const BsArc *a,double u0,double target,
                    double tol,int max_iter,BsInverse *out) {
    if(!out||!valid_curve(c)||!a||!isfinite(target)||target<0||!isfinite(tol)||tol<=0||
       max_iter<1) return BS_INVALID;
    memset(out,0,sizeof(*out)); out->u=u0;
    double base; BsStatus s=bs_arc_cumulative(a,u0,&base); if(s!=BS_OK) return s;
    s=bs_eval(c,u0,&out->point); if(s!=BS_OK) return s;
    if(target==0) return BS_OK;
    double end=a->breaks[a->spans],remaining=a->prefix[a->spans]-base;
    if(target>remaining+tol) return BS_UNREACHABLE;
    if(fabs(target-remaining)<=tol) {
        out->u=end; out->length=remaining; out->abs_error=fabs(target-remaining);
        return bs_eval(c,end,&out->point);
    }
    if(remaining<=0) return BS_UNREACHABLE;
    double lo=u0,hi=end,u=u0+(end-u0)*(target/remaining);
    for(int i=0;i<max_iter;i++) {
        double cumulative;
        s=bs_arc_cumulative(a,u,&cumulative); if(s!=BS_OK) return s;
        double length=cumulative-base,f=length-target;
        out->u=u; out->length=length; out->abs_error=fabs(f); out->iterations=i+1;
        s=bs_eval(c,u,&out->point); if(s!=BS_OK) return s;
        if(fabs(f)<=tol) return BS_OK;
        if(f<0) lo=u; else hi=u;
        BsPoint derivative; s=bs_eval(&a->derivative,u,&derivative); if(s!=BS_OK) return s;
        double speed=norm(derivative),next=NAN;
        if(speed>=1e-12) next=u-f/speed;
        /* Safeguarded Newton: retain a bracket at every iteration. */
        if(!isfinite(next)||next<=lo||next>=hi||fabs(next-u)>(hi-lo)*0.75) {
            next=lo+(hi-lo)/2; out->bisection_steps++;
        }
        if(next==u) return BS_NOT_CONVERGED;
        u=next;
    }
    return BS_NOT_CONVERGED;
}
