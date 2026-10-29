#include "bspline.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int number(const char *text,double *out) {
    char *end; errno=0; *out=strtod(text,&end);
    return end!=text&&!*end&&!errno&&isfinite(*out);
}
static void point(BsPoint p) { printf("[%.17g,%.17g]",p.x,p.y); }
int main(int argc,char **argv) {
    double args[]={500,-150,0,0,150,-1,20,0.04,0.05};
    int json=0,continuous=0,pos=0,result=1;
    BsFit fit={0}; BsArc arc={0}; BsCurve derivative={0};
    for(int i=1;i<argc;i++) {
        if(!strcmp(argv[i],"--json")) json=1;
        else if(!strcmp(argv[i],"--continuous")) continuous=1;
        else if(!strcmp(argv[i],"--help")) {
            puts("bspline_demo [--json] [--continuous] [Ax Ay Ox Oy R branch d_init err_min err_max]");
            return 0;
        } else if(pos>=9||!number(argv[i],&args[pos++])) { fputs("Invalid argument\n",stderr); return 2; }
    }
    if(pos!=0&&pos!=9) { fputs("Expected all 9 numeric arguments\n",stderr); return 2; }
    BsGeometry geo;
    BsStatus status=bs_geometry((BsPoint){args[0],args[1]},(BsPoint){args[2],args[3]},
                               args[4],args[5]>=0?1:-1,continuous,&geo);
    if(status!=BS_OK) goto cleanup;
    status=bs_fit(&geo,args[7],args[8],args[6],4,1000,60,1e-5,&fit);
    if(status!=BS_OK) goto cleanup;
    /* Tighter than the original 1e-6 full-length request; reused for inversions. */
    status=bs_arc_init(&fit.curve,1e-10,4,256,&arc); if(status!=BS_OK) goto cleanup;
    status=bs_derivative(&fit.curve,&derivative); if(status!=BS_OK) goto cleanup;
    BsInverse inv;
    status=bs_inverse(&fit.curve,&arc,0.5,1.2,1e-6,100,&inv); if(status!=BS_OK) goto cleanup;
    if(!json) {
        printf("mode: %s\nd: %.12g\nactual step: %.12g\ncontrol points: %zu\n"
               "sampled max error: %.12g\narc length: %.12g\nGL difference estimate: %.3e\n"
               "u0: 0.5\nu1: %.15g\npoint: (%.12g, %.12g)\n"
               "requested length: 1.2\ncomputed length: %.12g\nresidual: %.3e\n",
               continuous?"tangent-continuous":"original-compatible",fit.d,fit.actual_step,
               fit.curve.count,fit.max_error,arc.prefix[arc.spans],arc.estimated_error,
               inv.u,inv.point.x,inv.point.y,inv.length,inv.abs_error);
    } else {
        printf("{\"d\":%.17g,\"actual_step\":%.17g,\"max_error\":%.17g,"
               "\"length\":%.17g,\"estimated_error\":%.17g,\"direction\":%.17g,\"controls\":[",
               fit.d,fit.actual_step,fit.max_error,arc.prefix[arc.spans],arc.estimated_error,geo.direction);
        for(size_t i=0;i<fit.curve.count;i++) { if(i) putchar(','); point(fit.curve.points[i]); }
        printf("],\"knots\":[");
        for(size_t i=0;i<fit.curve.count+5;i++) { if(i) putchar(','); printf("%.17g",fit.curve.knots[i]); }
        printf("],\"samples\":[");
        for(int i=0;i<=256;i++) {
            double u=i/256.0; BsPoint value,dv; double cumulative;
            status=bs_eval(&fit.curve,u,&value); if(status!=BS_OK) goto cleanup;
            status=bs_eval(&derivative,u,&dv); if(status!=BS_OK) goto cleanup;
            status=bs_arc_cumulative(&arc,u,&cumulative); if(status!=BS_OK) goto cleanup;
            if(i) putchar(',');
            printf("[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]",u,value.x,value.y,dv.x,dv.y,cumulative);
        }
        printf("],\"inverse\":{\"u0\":0.5,\"target\":1.2,\"u\":%.17g,\"point\":",inv.u);
        point(inv.point);
        printf(",\"length\":%.17g,\"residual\":%.17g,\"iterations\":%d,\"bisection_steps\":%d},\"inverses\":[",
               inv.length,inv.abs_error,inv.iterations,inv.bisection_steps);
        int emitted=0;
        for(int i=0;i<4;i++) for(int j=0;j<3;j++) {
            double u0=i*0.23,remaining;
            status=bs_arc_length(&arc,u0,1,&remaining); if(status!=BS_OK) goto cleanup;
            double target=remaining*(j+1)/4;
            status=bs_inverse(&fit.curve,&arc,u0,target,1e-8,100,&inv); if(status!=BS_OK) goto cleanup;
            if(emitted++) putchar(',');
            printf("[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]",u0,target,inv.u,
                   inv.point.x,inv.point.y,inv.length);
        }
        puts("]}");
    }
    result=0;
cleanup:
    if(status!=BS_OK) fprintf(stderr,"%s (last d=%.17g, sampled error=%.17g)\n",
                             bs_status_string(status),fit.d,fit.max_error);
    bs_curve_free(&derivative); bs_arc_free(&arc); bs_fit_free(&fit);
    return result;
}
