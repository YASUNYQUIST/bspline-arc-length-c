#include "bspline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(expr) do { if(!(expr)) { fprintf(stderr,"line %d: %s\n",__LINE__,#expr); return 1; } } while(0)
int main(void) {
    /* Quartic Bezier with equally spaced collinear controls: C(u)=(3u,4u). */
    BsPoint points[]={{0,0},{.75,1},{1.5,2},{2.25,3},{3,4}};
    double knots[]={0,0,0,0,0,1,1,1,1,1};
    BsCurve line={5,4,points,knots}; BsArc arc={0}; BsPoint pt; BsInverse inv;
    CHECK(bs_arc_init(&line,1e-10,4,256,&arc)==BS_OK);
    CHECK(fabs(arc.prefix[arc.spans]-5)<1e-12);
    CHECK(bs_eval(&line,1,&pt)==BS_OK&&pt.x==3&&pt.y==4);
    CHECK(bs_eval(&line,-.1,&pt)==BS_INVALID);
    CHECK(bs_eval(&line,NAN,&pt)==BS_INVALID);
    for(int i=0;i<100;i++) {
        double u=((i*37)%100)/100.0,s;
        CHECK(bs_arc_cumulative(&arc,u,&s)==BS_OK&&fabs(s-5*u)<1e-12);
    }
    CHECK(bs_inverse(&line,&arc,.2,2,1e-9,100,&inv)==BS_OK);
    CHECK(fabs(inv.u-.6)<1e-12&&fabs(inv.point.x-1.8)<1e-12);
    CHECK(bs_inverse(&line,&arc,0,5,1e-9,100,&inv)==BS_OK&&inv.u==1);
    CHECK(bs_inverse(&line,&arc,1,0,1e-9,100,&inv)==BS_OK&&inv.u==1);
    CHECK(bs_inverse(&line,&arc,0,6,1e-9,100,&inv)==BS_UNREACHABLE);
    CHECK(bs_inverse(&line,&arc,0,-1,1e-9,100,&inv)==BS_INVALID);
    double len;
    CHECK(bs_arc_length(&arc,.8,.2,&len)==BS_OK&&fabs(len-3)<1e-12);
    bs_arc_free(&arc);
    /* A degree-one curve on a non-normalized domain [2,5]. */
    BsPoint segment_points[]={{0,0},{3,4}};
    double segment_knots[]={2,2,5,5};
    BsCurve segment={2,1,segment_points,segment_knots};
    CHECK(bs_arc_init(&segment,1e-10,4,256,&arc)==BS_OK);
    CHECK(bs_inverse(&segment,&arc,2,2.5,1e-9,100,&inv)==BS_OK);
    CHECK(fabs(inv.u-3.5)<1e-12);
    bs_arc_free(&arc);
    /* A stationary first span followed by a moving span. */
    BsPoint plateau_points[]={{0,0},{0,0},{3,4}};
    double plateau_knots[]={0,0,.5,1,1};
    BsCurve plateau={3,1,plateau_points,plateau_knots};
    CHECK(bs_arc_init(&plateau,1e-10,4,256,&arc)==BS_OK);
    CHECK(bs_inverse(&plateau,&arc,0,2.5,1e-9,100,&inv)==BS_OK);
    CHECK(fabs(inv.u-.75)<1e-12);
    bs_arc_free(&arc);
    /* Zero speed at u=0: C(u)=(u^4,0), analytic length=u^4. */
    BsPoint slow[]={{0,0},{0,0},{0,0},{0,0},{1,0}};
    BsCurve power={5,4,slow,knots};
    CHECK(bs_arc_init(&power,1e-11,4,256,&arc)==BS_OK);
    CHECK(bs_inverse(&power,&arc,0,.0001,1e-12,100,&inv)==BS_OK);
    CHECK(fabs(inv.u-.1)<1e-8);
    CHECK(bs_inverse(&power,&arc,0,.0001,1e-14,1,&inv)==BS_NOT_CONVERGED);
    bs_arc_free(&arc);
    BsPoint constant[]={{2,3},{2,3},{2,3},{2,3},{2,3}};
    BsCurve zero={5,4,constant,knots};
    CHECK(bs_arc_init(&zero,1e-10,4,256,&arc)==BS_OK);
    CHECK(arc.prefix[arc.spans]==0);
    CHECK(bs_inverse(&zero,&arc,.5,1,1e-9,100,&inv)==BS_UNREACHABLE);
    bs_arc_free(&arc);
    BsGeometry geo; BsFit fit={0}; BsCurve c={0}; double step;
    CHECK(bs_geometry((BsPoint){1,0},(BsPoint){0,0},1,1,0,&geo)==BS_INVALID);
    CHECK(bs_geometry((BsPoint){500,-150},(BsPoint){0,0},150,-1,0,&geo)==BS_OK);
    CHECK(bs_sample(&geo,0,4,&c,&step)==BS_INVALID);
    CHECK(bs_sample(&geo,1e-20,4,&c,&step)==BS_LIMIT);
    CHECK(bs_sample(&geo,geo.s_total,4,&c,&step)==BS_INVALID);
    CHECK(bs_fit(&geo,.04,.05,20,4,1000,60,1e-5,&fit)==BS_OK);
    CHECK(fit.max_error>=.04&&fit.max_error<=.05);
    CHECK(bs_arc_init(&fit.curve,1e-10,4,256,&arc)==BS_OK);
    CHECK(bs_inverse(&fit.curve,&arc,.5,1.2,1e-6,100,&inv)==BS_OK);
    bs_arc_free(&arc); bs_fit_free(&fit);
    CHECK(bs_fit(&geo,.0400000000001,.0400000000002,20,4,1000,60,1e-5,&fit)==BS_NOT_CONVERGED);
    bs_fit_free(&fit);
    puts("All C checks passed (analytic curves, endpoints, stationary speed, failure paths, fit).");
    return 0;
}
