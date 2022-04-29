//Adapted from fdom.c
//These methods deal with basic plane geometry. The main purpose is to implement Mobius maps acting on circles/lines.

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "apol.h"
#endif


//STATIC METHOD DECLARATIONS

//BASIC LINE, CIRCLE, AND POINT OPERATIONS
//static GEN arc_init(GEN c, GEN p1, GEN p2, int dir, long prec);
static GEN arc_midpoint(GEN c, GEN p1, GEN p2, GEN tol, long prec);
static GEN midpoint(GEN p1, GEN p2);
static GEN mobius(GEN M, GEN c, GEN tol, long prec);
static GEN mobius_arcseg(GEN M, GEN c, int isarc, GEN tol, long prec);
static GEN mobius_circle(GEN M, GEN c, GEN tol, long prec);
static GEN mobius_line(GEN M, GEN l, GEN tol, long prec);
static GEN perpbis(GEN p1, GEN p2, GEN tol, long prec);
static GEN radialangle(GEN c, GEN p, GEN tol, long prec);
static GEN slope(GEN p1, GEN p2, GEN tol, long prec);

//INTERSECTION OF LINES/CIRCLES
static GEN arc_int(GEN c1, GEN c2, GEN tol, long prec);
static GEN arcseg_int(GEN c, GEN l, GEN tol, long prec);
static GEN circle_int(GEN c1, GEN c2, GEN tol, long prec);
static GEN circleline_int(GEN c, GEN l, GEN tol, long prec);
static GEN line_int(GEN l1, GEN l2, GEN tol, long prec);
static int onarc(GEN c, GEN p, GEN tol, long prec);
static int onseg(GEN l, GEN p, GEN tol, long prec);
static GEN seg_int(GEN l1, GEN l2, GEN tol, long prec);

//GEOMETRIC HELPER METHODS
static GEN anglediff(GEN ang, GEN bot, GEN tol, long prec);
static int geom_check(GEN c);
static GEN shiftangle(GEN ang, GEN bot, GEN tol, long prec);
static int tolcmp(GEN x, GEN y, GEN tol, long prec);
static int toleq(GEN x, GEN y, GEN tol, long prec);
static int toleq0(GEN x, GEN tol, long prec);

//Circle is stored as [centre, radius, curvature]. radius>0, but curvature can be <0 to mean the outside
//Circle arc is [centre, radius, curvature, start pt, end pt, start angle, end angle, dir]. It is the arc counterclockwise from startpt to endpt, and dir=1 means oriented counterclockwise, and dir=-1 means oriented clockwise. This can also be left uninitialized if arc is not oriented. The final 0 represents that we have an arc and not a segment.
//Line is stored as [slope, intercept], where the line is y=slope*x+intercept unless slope=oo, where it is x=intercept instead.
//Line segment is stored as [slope, intercept, startpt, endpt, ooendptor, dir]. The final 1 is to signal a segment. dir=1 means we go from startpt to endpt in the upper half plane, and -1 means we go through infinity (only when neither endpoint is infinite). If one endpoint is oo, ooendptor stores which way we get to it. If ooendptor=1, this means the segment travels either vertically up or right, and -1 means the arc is vertically down or left.

//GEN tol -> The tolerance.



//BASIC LINE, CIRCLE, AND POINT OPERATIONS


//Creates the arc on circle c going from p1 to p2 counterclockwise. if dir=1, the arc is oriented counterclockwise, else it is clockwise (i.e. clockwise from p2 to p1). If dir=0, we take it to be unoriented
GEN arc_init(GEN c, GEN p1, GEN p2, int dir, long prec){
  pari_sp top=avma;
  GEN ang2=radialangle(c, p2, gen_0, prec);//No tolerance need
  GEN arc=cgetg(9, t_VEC);
  gel(arc, 1)=gcopy(gel(c, 1));//Centre
  gel(arc, 2)=gcopy(gel(c, 2));///Radius
  gel(arc, 3)=gcopy(gel(c, 3));//Curvature
  gel(arc, 4)=gcopy(p1);//Start point
  gel(arc, 5)=gcopy(p2);//End point
  gel(arc, 6)=radialangle(c, p1, gen_0, prec);//Start angle
  gel(arc, 7)=shiftangle(ang2, gel(arc, 5), gen_0, prec);//End angle; no need for tolerance.
  if(dir==1) gel(arc, 8)=gen_1;
  else if(dir==0) gel(arc, 8)=gen_0;
  else gel(arc, 8)=gen_m1;
  return gerepileupto(top, arc);
}

//Returns the midpoint of the arc between p1 and p2 (counterclockwise) on c.
static GEN arc_midpoint(GEN c, GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  GEN pts=circleline_int(c, perpbis(p1, p2, tol, prec), tol, prec);
  GEN a1=radialangle(c, p1, gen_0, prec);
  GEN a2=shiftangle(radialangle(c, p2, gen_0, prec), a1, gen_0, prec);//No tolerance concerns
  GEN angint1=shiftangle(radialangle(c, gel(pts, 1), gen_0, prec), a1, gen_0, prec);//the angle formed by pts[1] to c with base a1. Again, no tolerance need.
  if(gcmp(a2, angint1)>0) return gerepilecopy(top, gel(pts, 1));//No need for tolerance, as if this is an issue our points p1 and p2 would be equal up to tolerance.
  return gerepilecopy(top, gel(pts, 2));
}

//Circle with centre cent passing through p
GEN circle_fromcp(GEN cent, GEN p, long prec){
  pari_sp top=avma;
  GEN abspmcent=gabs(gsub(p, cent), prec);
  return gerepilecopy(top, mkvec3(cent, abspmcent, divoo(gen_1, abspmcent)));
}

//Circle through 3 points (with one allowed to be oo, making a line instead)
GEN circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol, long prec){
  if(typ(p1)==t_INFINITY) return line_frompp(p2, p3, tol, prec);
  if(typ(p2)==t_INFINITY) return line_frompp(p1, p3, tol, prec);
  if(typ(p3)==t_INFINITY) return line_frompp(p1, p2, tol, prec);//Lines
  pari_sp top=avma;
  GEN l1=perpbis(p1, p2, tol, prec), l2=perpbis(p1, p3, tol, prec);
  GEN centre=line_int(l1, l2, tol, prec);//centre is intersection of perp bisectors.
  if(typ(centre)==t_INFINITY) return gerepileupto(top, line_frompp(p1, p2, tol, prec));//p1, p2, p3 collinear
  return gerepileupto(top, circle_fromcp(centre, p1, prec));//The circle!
}

//The line through p with slope s
GEN line_fromsp(GEN s, GEN p){
  if(typ(s)==t_INFINITY) retmkvec2(mkoo(), greal(p));//oo slope
  pari_sp top=avma;
  GEN strealp=gmul(s, real_i(p));
  return gerepilecopy(top, mkvec2(s, gsub(imag_i(p), strealp)));
}

//Line through two points
GEN line_frompp(GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  return gerepileupto(top, line_fromsp(slope(p1, p2, tol, prec), p1));
}

//Mx, where M is a 2x2 matrix and x is complex or infinite.
GEN mat_eval(GEN M, GEN x){
  pari_sp top = avma;
  if(typ(x)==t_INFINITY) return gerepileupto(top, divoo(gcoeff(M, 1, 1), gcoeff(M, 2, 1)));
  return gerepileupto(top, divoo(gadd(gmul(gcoeff(M, 1, 1), x), gcoeff(M, 1, 2)), gadd(gmul(gcoeff(M, 2, 1), x), gcoeff(M, 2, 2))));
}

//Midpoint of p1 and p2
static GEN midpoint(GEN p1, GEN p2){
  pari_sp top=avma;
  return gerepileupto(top, gdivgs(gadd(p1, p2), 2));
}

//Mobius transform accessible in GP
GEN mobius_gp(GEN M, GEN c, long prec){
  pari_sp top=avma;
  if(typ(M)!=t_MAT || lg(M)!=3 || lg(gel(M, 1))!=3) pari_err_TYPE("M needs to be a 2x2 matrix", M);
  GEN tol=deftol(prec);
  return gerepileupto(top, mobius(M, c, tol, prec));
}

//Returns M(c), for c a circle/line/arc/segment
static GEN mobius(GEN M, GEN c, GEN tol, long prec){
  int i=geom_check(c);
  switch(i){
    case 0: return mobius_circle(M, c, tol, prec);
    case 1: return mobius_line(M, c, tol, prec);
    case 2: return mobius_arcseg(M, c, 1, tol, prec);
    case 3: return mobius_arcseg(M, c, 0, tol, prec);
  }
  pari_err_TYPE("Please input a circle/line/arc/segment", c);
  return gen_0;//Never reach this far
}

//Mobius map acting on an arc or segment
static GEN mobius_arcseg(GEN M, GEN c, int isarc, GEN tol, long prec){
  pari_sp top=avma;//We start by finding the new circle/line(ignoring the arc/segment bits)
  GEN endpt1, endpt2, extpt;//Endpoints and an extra point (used to have 3 points to translate)
  endpt1=mat_eval(M, gel(c, 4));
  endpt2=mat_eval(M, gel(c, 5));
  //Now we must find an extra point extpt
  if(isarc==1) extpt=mat_eval(M, arc_midpoint(c, gel(c, 4), gel(c, 5), tol, prec));//arc
  else{//segment
    if(typ(gel(c, 3))==t_INFINITY){//Start point infinity
      GEN u;
      if(gequal1(gel(c, 5))) u=gen_m2;//segment goes vertically up or right
      else u=gen_2;//segment goes vertically down or left
      if(typ(gel(c, 1))==t_INFINITY) extpt=mat_eval(M, gadd(gel(c, 4), mulcxI(u)));//Vertical line, M(c[4]+U*I)
      else extpt=mat_eval(M, gadd(gel(c, 4), gmul(u, gaddsg(1, mulcxI(gel(c, 1))))));//non-vertical line, M(c[4]+u+u*c[1]*I)
    }
    else if(typ(gel(c, 4))==t_INFINITY){//End point infinity
      GEN u;
      if(gequal1(gel(c, 5))) u=gen_2;//segment goes vertically up or right
      else u=gen_m2;//segment goes vertically down or left
      if(typ(gel(c, 1))==t_INFINITY) extpt=mat_eval(M, gadd(gel(c, 3), mulcxI(u)));//Vertical line, M(c[3]+U*I)
      else extpt=mat_eval(M, gadd(gel(c, 3), gmul(u, gaddsg(1, mulcxI(gel(c, 1))))));//non-vertical line, M(c[3]+u+u*c[1]*I)
    }
    else{//Start/end points in the plane
      if(gequal1(gel(c, 6))) extpt=mat_eval(M, midpoint(gel(c, 3), gel(c, 4)));//Does not go through oo, can take midpoint
      else extpt=mat_eval(M, mkoo());//Use oo, since the line goes through there
    }
  }
  //Now we have the 3 new points used to define the new arc/segment. Let's finish the process.
  GEN ret;//The returned arc/segment
  GEN newcirc=circle_fromppp(endpt1, endpt2, extpt, tol, prec);//The new circle/line
  if(lg(newcirc)==4){//Circle
    ret=vec_lengthen(newcirc, 8);
    gel(ret, 4)=endpt1;//Start point
    gel(ret, 5)=endpt2;//end point. These may be in the wrong order, so will fix later if so.
    gel(ret, 6)=radialangle(newcirc, endpt1, tol, prec);//angle 1
    gel(ret, 7)=shiftangle(radialangle(newcirc, endpt2, tol, prec), gel(ret, 6), tol, prec);//angle 2
    if(isarc) gel(ret, 8)=gel(c, 8);//Temporary; match the old to the new direction
    else gel(ret, 8)=gen_1;//Temporary; since for segments we have a bona fide start and end, we start with the direction being 1.
    if(!onarc(ret, extpt, tol, prec)){//Must swap start/endpoints, angles and the direction
      gel(ret, 4)=endpt2;
      gel(ret, 5)=endpt1;
      GEN tempang=gel(ret, 6);//angle to endpt1
      gel(ret, 6)=shiftangle(gel(ret, 7), gen_0, tol, prec);//The angle to endpt2 shifted to [0, 2*Pi)
      gel(ret, 7)=shiftangle(tempang, gel(ret, 6), tol, prec);//Angle to endpt1 shifted with a base of the angle to endpt2
      gel(ret, 8)=gneg(gel(ret, 8));//We now run backwards
    }
  }
  else{//Line
    ret=vec_lengthen(newcirc, 6);
    gel(ret, 3)=endpt1;//Start point
    gel(ret, 4)=endpt2;//end point. These may be in the wrong order, so will fix later if so.
    if(isarc==1 && gequalm1(gel(c, 8))){//We need to reverse the order of the points, because the arc ran backwards.
      gel(ret, 3)=endpt2;
      gel(ret, 4)=endpt1;
    }
    if(typ(endpt1)==t_INFINITY || typ(endpt2)==t_INFINITY){//oo endpoint
      gel(ret, 5)=gen_1;//Temporary, assume we go up/right
      gel(ret, 6)=gen_0;//Not used
      if(!onseg(ret, extpt, tol, prec)) gel(ret, 5)=gen_m1;//We were wrong, and go down/left.
    }
    else{//Both endpoints finite
      gel(ret, 5)=gen_0;
      gel(ret, 6)=gen_1;//Temporary, assume we go through the plane only
      if(!onseg(ret, extpt, tol, prec)) gel(ret, 6)=gen_m1;//We were wrong, and go through oo.
    }
  }
  return gerepilecopy(top, ret);
}

//Mobius map acting on circle
static GEN mobius_circle(GEN M, GEN c, GEN tol, long prec){
  pari_sp top=avma;
  GEN p1=mat_eval(M, gadd(gel(c, 1), gel(c, 2)));//M(c[1]+c[2])
  GEN p2=mat_eval(M, gadd(gel(c, 1), mulcxI(gel(c, 2))));//M(c[1]+c[2]*I)
  GEN p3=mat_eval(M, gsub(gel(c, 1), gel(c, 2)));//M(c[1]-c[2])
  return gerepileupto(top, circle_fromppp(p1, p2, p3, tol, prec));
}

//Mobius map acting on line
static GEN mobius_line(GEN M, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  GEN p1, p2, p3;
  if(typ(gel(l, 1))==t_INFINITY){//Vertical line
    p1=mat_eval(M, gel(l, 2));//M(x-intercept)
    p2=mat_eval(M, mkcomplex(gel(l, 2), gen_1));//M(x-intercept+I)
    p3=mat_eval(M, mkcomplex(gel(l, 2), gen_m1));//M(x-intercept-I)
  }
  else{//Non-vertical line
    GEN slopeIp1=mkcomplex(gen_1, gel(l, 1));//gaddgs(mulcxI(gel(l, 1)), 1);//1+Slope*I
    GEN p1base=mulcxI(gel(l, 2));//y-intercept
    GEN p2base=gadd(p1base, slopeIp1);//y-intercept+1+slope*I
    GEN p3base=gadd(p2base, slopeIp1);//y-intercept+2+2*slope*I
    p1=mat_eval(M, p1base);p2=mat_eval(M, p2base);p3=mat_eval(M, p3base);
  }
  return gerepileupto(top, circle_fromppp(p1, p2, p3, tol, prec));
}

//Perpendicular bisector of distinct points
static GEN perpbis(GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  return gerepileupto(top, line_fromsp(divoo(gen_m1, slope(p1, p2, tol, prec)), midpoint(p1, p2)));
}

//Angle between p and the centre of c, in the range [0, 2*Pi)
static GEN radialangle(GEN c, GEN p, GEN tol, long prec){
  pari_sp top=avma;
  return gerepileupto(top, shiftangle(garg(gsub(p, gel(c, 1)), prec), gen_0, tol, prec));
}

//The slope of the line through p1, p2
static GEN slope(GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  GEN p2mp1=gsub(p2, p1);
  GEN ftop=imag_i(p2mp1);
  GEN fbot=real_i(p2mp1);
  if(toleq0(fbot, tol, prec)) fbot=gen_0;
  if(toleq0(ftop, tol, prec)) ftop=gen_0;
  return gerepileupto(top, divoo(ftop, fbot));
}



//INTERSECTION OF LINES/CIRCLES


//Returns the intersection points of two arcs
static GEN arc_int(GEN c1, GEN c2, GEN tol, long prec){
  pari_sp top=avma;
  GEN ipts=circle_int(c1, c2, tol, prec);
  if(lg(ipts)==1) return gc_0vec(top);//No intersection
  if(lg(ipts)==2){//One intersection point (tangent circles)
    if(!onarc(c1, gel(ipts, 1), tol, prec)) return gc_0vec(top);//Not on arc 1
    if(!onarc(c2, gel(ipts, 1), tol, prec)) return gc_0vec(top);//Not on arc 2
    return gerepilecopy(top, ipts);//On arc
  }
  //Two intersections
  int i1=onarc(c1, gel(ipts, 1), tol, prec);
  if(i1==1) i1=onarc(c2, gel(ipts, 1), tol, prec);//Now i1==1 iff the ipts[1] is on both c1 and c2
  int i2=onarc(c1, gel(ipts, 2), tol, prec);
  if(i2==1) i2=onarc(c2, gel(ipts, 2), tol, prec);//Now i2==1 iff the ipts[2] is on both c1 and c2
  if(i1==1){
    if(i2==1) return gerepilecopy(top, ipts);//Both pts on the arcs
    return gerepilecopy(top, mkvec(gel(ipts, 1)));//Just point 1
  }
  //Now i1=0
  if(i2==0) return gc_0vec(top);//Not on either arc
  return gerepilecopy(top, mkvec(gel(ipts, 2)));//Just point 2
}

//Returns the intersection points of an arc and a segment
static GEN arcseg_int(GEN c, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  GEN ipts=circleline_int(c, l, tol, prec);
  if(lg(ipts)==1) return gc_0vec(top);//No intersection
  if(lg(ipts)==2){//One intersection point (tangent circles)
    if(!onarc(c, gel(ipts, 1), tol, prec)) return gc_0vec(top);//Not on arc
    if(!onseg(l, gel(ipts, 1), tol, prec)) return gc_0vec(top);//Not on segment
    return gerepilecopy(top, ipts);//On both
  }
  //Two intersections
  int i1=onarc(c, gel(ipts, 1), tol, prec);
  if(i1==1) i1=onseg(l, gel(ipts, 1), tol, prec);//Now i1==1 iff the ipts[1] is on both c and l
  int i2=onarc(c, gel(ipts, 2), tol, prec);
  if(i2==1) i2=onseg(l, gel(ipts, 2), tol, prec);//Now i2==1 iff the ipts[2] is on both c and l
  if(i1==1){
    if(i2==1) return gerepilecopy(top, ipts);//Both pts on both
    return gerepilecopy(top, mkvec(gel(ipts, 1)));//Just point 1
  }
  //Now i1=0
  if(i2==0) return gc_0vec(top);//Not on either
  return gerepilecopy(top, mkvec(gel(ipts, 2)));//Just point 2
}

//Returns the set of points in the intersection of circles c1, c2
static GEN circle_int(GEN c1, GEN c2, GEN tol, long prec){
  pari_sp top=avma;
  GEN a1=real_i(gel(c1, 1)), b1=imag_i(gel(c1, 1)), r1=gel(c1, 2);//x, y coords and radius of c1
  GEN a2=real_i(gel(c2, 1)), b2=imag_i(gel(c2, 1)), r2=gel(c2, 2);//x, y coords and radius of c2
  GEN a1ma2=gsub(a1, a2), b1mb2=gsub(b1, b2), x1, x2, y1, y2;
  int oneint=0;
  if(gcmp(gabs(a1ma2, prec), gabs(b1mb2, prec))>=0){//We want to divide by the larger of the two quantities to maximize precision and avoid errors when the centres are on the same line.
    if(toleq0(a1ma2, tol, prec)==1) return gc_0vec(top);//Same centre, cannot intersect.
    //u=(r1^2-r2^2+b2^2-b1^2+a2^2-a1^2)/(2*a2-2*a1)-a1;
    GEN u=gsub(gdiv(gsub(gadd(gsub(gadd(gsub(gsqr(r1), gsqr(r2)), gsqr(b2)), gsqr(b1)), gsqr(a2)), gsqr(a1)), gmulgs(a1ma2, -2)), a1);
    GEN v=gneg(gdiv(b1mb2, a1ma2));//v=(b1-b2)/(a2-a1), and x=a1+u+vy
    GEN uvmb1=gsub(gmul(u, v), b1);//uv-b1
    GEN vsqrp1=gaddgs(gsqr(v), 1);//v^2+1
    GEN rtpart=gsub(gsqr(uvmb1), gmul(vsqrp1, gadd(gsqr(b1), gsub(gsqr(u), gsqr(r1)))));//(u*v-b1)^2-(v^2+1)*(b1^2+u^2-r1^2)
    oneint=tolcmp(rtpart, gen_0, tol, prec);//Comparing rtpart to 0
    if(oneint==-1) return gc_0vec(top);//rtpart must be square rooted, so if it's negative the circles do not intersect
    if(oneint==0){//One intersection, so we take rtpart=0. This is CRUCIAL, as taking the square root kills our precision if we don't do this here.
      y1=gdiv(gneg(uvmb1), vsqrp1);//y1=(b1-u*v)/(1*v^2+1)
      x1=gadd(gadd(a1, u), gmul(v, y1));//x1=a1+u+v*y1
    }
    else{
      y1=gdiv(gadd(gneg(uvmb1), gsqrt(rtpart, prec)), vsqrp1);//y1=(b1-u*v+sqrt((u*v-b1)^2-*(v^2+1)*(b1^2+u^2-r1^2)))/(1*v^2+1)
      y2=gadd(gneg(y1), gdiv(gmulgs(uvmb1, -2), vsqrp1));//y2=-y1+(2*b1-2*u*v)/(v^2+1)
      GEN a1pu=gadd(a1, u);
      x1=gadd(a1pu, gmul(v, y1));//x1=a1+u+v*y1
      x2=gadd(a1pu, gmul(v, y2));//x1=a1+u+v*y2
    }
  }
  else{
    if(toleq0(b1mb2, tol, prec)==1) return gc_0vec(top);//Same centre, cannot intersect.
    //u=(r1^2-r2^2+b2^2-b1^2+a2^2-a1^2)/(2*b2-2*b1)-b1;
    GEN u=gsub(gdiv(gsub(gadd(gsub(gadd(gsub(gsqr(r1), gsqr(r2)), gsqr(b2)), gsqr(b1)), gsqr(a2)), gsqr(a1)), gmulgs(b1mb2, -2)), b1);
    GEN v=gneg(gdiv(a1ma2, b1mb2));//v=(a1-a2)/(b2-b1), and y=b1+u+vx
    GEN uvma1=gsub(gmul(u, v), a1);//uv-a1
    GEN vsqrp1=gaddgs(gsqr(v), 1);//v^2+1
    GEN rtpart=gsub(gsqr(uvma1), gmul(vsqrp1, gadd(gsqr(a1), gsub(gsqr(u), gsqr(r1)))));//(u*v-a1)^2-(v^2+1)*(a1^2+u^2-r1^2))
    oneint=tolcmp(rtpart, gen_0, tol, prec);//Comparing rtpart to 0
    if(oneint==-1) return gc_0vec(top);//rtpart must be square rooted, so if it's negative the circles do not intersect
    if(oneint==0){//One intersection, so we take rtpart=0. This is CRUCIAL, as taking the square root kills our precision if we don't do this here.
      x1=gdiv(gneg(uvma1), vsqrp1);//x1=(a1-u*v)/(v^2+1);
      y1=gadd(gadd(b1, u), gmul(v, x1));//y1=b1+u+v*x1;
    }
    else{
      x1=gdiv(gadd(gneg(uvma1), gsqrt(rtpart, prec)), vsqrp1);//x1=(a1-u*v+sqrt((u*v-a1)^2-(v^2+1)*(a1^2+u^2-r1^2)))/(v^2+1);
      x2=gadd(gneg(x1), gdiv(gmulgs(uvma1, -2), vsqrp1));//x2=-x1+(2*a1-2*u*v)/(v^2+1)
      GEN b1pu=gadd(b1, u);
      y1=gadd(b1pu, gmul(v, x1));//y1=b1+u+v*x1;
      y2=gadd(b1pu, gmul(v, x2));//y2=b1+u+v*x2;
    }
  }
  GEN P1=mkcomplex(x1, y1);
  if(oneint==0) return gerepilecopy(top, mkvec(P1));//One point of intersection (0 pts of intersection was already dealt with and returned)
  return gerepilecopy(top, mkvec2(P1, mkcomplex(x2, y2)));
}

//Returns the intersection points of c and l
static GEN circleline_int(GEN c, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  if(typ(gel(l, 1))==t_INFINITY){
    GEN x1=gel(l, 2);
    GEN rtpart=gsub(gsqr(gel(c, 2)), gsqr(gsub(x1, real_i(gel(c, 1)))));//c[2]^2-(x1-real(c[1]))^2
    if(gsigne(rtpart)==-1) return gc_0vec(top);//No intersections.
    GEN y1=gadd(imag_i(gel(c, 1)), gsqrt(rtpart, prec));//y1=imag(c[1])+sqrt(c[2]^2-(x1-real(c[1]))^2)
    if(toleq0(rtpart, tol, prec)) return gerepilecopy(top, mkvec(mkcomplex(x1, y1)));
    //Two intersection points
    GEN y1py2=gmulgs(imag_i(gel(c, 1)), 2);//2*imag(c[1])
    return gerepilecopy(top, mkvec2(mkcomplex(x1, y1), mkcomplex(x1, gsub(y1py2, y1))));
  }
  //Now y=mx+b with m finite
  GEN A=gaddgs(gsqr(gel(l, 1)), 1);//l[1]^2+1
  GEN l2mic1=gsub(gel(l, 2), imag_i(gel(c, 1)));//l[2]-imag(c[1])
  GEN B=gadd(gmulgs(real_i(gel(c, 1)), -2), gmulsg(2, gmul(gel(l, 1), l2mic1)));//-2*real(c[1])+2*l[1]*(l[2]-imag(c[1]))
  GEN C=gadd(gsqr(real_i(gel(c, 1))), gsub(gsqr(l2mic1), gsqr(gel(c, 2))));//real(c[1])^2+(l[2]-imag(c[1]))^2-c[2]^2
  GEN rtpart=gsub(gsqr(B), gmulsg(4, gmul(A, C)));
  int rtpartsig=tolcmp(rtpart, gen_0, tol, prec);
  if(rtpartsig==-1) return gc_0vec(top);//No intersection
  if(rtpartsig==0){//One root, and rtpart=0
    GEN x1=gdiv(B, gmulgs(A, -2));//-B/(2A)
    GEN y1part=gmul(gel(l, 1), x1);//l[1]*x1
    return gerepilecopy(top, mkvec(mkcomplex(x1, gadd(y1part, gel(l, 2)))));//y1=l[1]*x1+l[2];
  }
  //Two roots
  GEN x1=gdiv(gsub(gsqrt(rtpart, prec), B), gmulgs(A, 2));//x1=(-B+sqrt(B^2-4*A*C))/(2*A);
  GEN x2=gsub(gneg(gdiv(B, A)), x1);//-B/A-x1
  GEN y1part=gmul(gel(l, 1), x1);//l[1]*x1
  GEN y2part=gmul(gel(l, 1), x2);//l[1]*x2
  GEN y1=gadd(y1part, gel(l, 2));//l[1]*x1+l[2];
  GEN y2=gadd(y2part, gel(l, 2));//l[1]*x2+l[2];
  return gerepilecopy(top, mkvec2(mkcomplex(x1, y1), mkcomplex(x2, y2)));
}

//Returns the intersection of c1 and c2, where c1 and c2 are circles/lines/arcs/segments
GEN geom_int(GEN c1, GEN c2, long prec){
  pari_sp top=avma;
  GEN tol=deftol(prec);
  int i1=geom_check(c1), i2=geom_check(c2);
  switch(i1){
    case 0://c1=circle
      switch(i2){
        case 0: return gerepileupto(top, circle_int(c1, c2, tol, prec));
        case 1: return gerepileupto(top, circleline_int(c1, c2, tol, prec));
        case 2: return gerepileupto(top, arc_int(c1, c2, tol, prec));
        case 3: return gerepileupto(top, circleline_int(c1, c2, tol, prec));
      }
    case 1://c1=line
      switch(i2){
        case 0: return gerepileupto(top, circleline_int(c2, c1, tol, prec));
        case 1: return gerepileupto(top, line_int(c1, c2, tol, prec));
        case 2: return gerepileupto(top, arcseg_int(c2, c1, tol, prec));
        case 3:;
          GEN ret=seg_int(c1, c2, tol, prec);
          if(ret) return gerepilecopy(top, mkvec(ret));
          return gc_0vec(top);//No intersection on the segments.
      }
    case 2://c1=arc
      switch(i2%2){
        case 0: return gerepileupto(top, arc_int(c1, c2, tol, prec));
        case 1: return gerepileupto(top, arcseg_int(c1, c2, tol, prec));
      }
    case 3://c1=segment
      switch(i2%2){
        case 0: return gerepileupto(top, arcseg_int(c2, c1, tol, prec));
        case 1:;
          GEN ret=seg_int(c1, c2, tol, prec);
          if(ret) return gerepilecopy(top, mkvec(ret));
          return gc_0vec(top);//No intersection on the segments.
    }
  }
  return gc_0vec(top);//Bad input.
}

//The intersection of two lines
static GEN line_int(GEN l1, GEN l2, GEN tol, long prec){
  GEN s1=gel(l1, 1), s2=gel(l2, 1);//Slopes
  if(toleq(s1, s2, tol, prec)) return mkoo();//Parallel or equal
  pari_sp top=avma;
  if(typ(s1)==t_INFINITY){//l1 vertical
    GEN ypart=gmul(s2, gel(l1, 2));//s2*l1[2]
    return gerepilecopy(top, mkcomplex(gel(l1, 2), gadd(ypart, gel(l2, 2))));//y=s2*l1[2]+l2[2]
  }
  if(typ(s2)==t_INFINITY){//l2 vertical
    GEN ypart=gmul(s1, gel(l2, 2));//s1*l2[2]
    return gerepilecopy(top, mkcomplex(gel(l2, 2), gadd(ypart, gel(l1, 2))));//y=s1*l2[2]+l1[2]
  }
  GEN x=gdiv(gsub(gel(l2, 2), gel(l1, 2)), gsub(s1, s2));//(l2[2]-l1[2])/(s1-s2)
  GEN ypart=gmul(s1, x);//s1*x
  return gerepilecopy(top, mkcomplex(x, gadd(ypart, gel(l1, 2))));//y=s1*x+l1[2]
}

//p is assumed to be on the circle defined by c; this checks if it is actually on the arc (running counterclockwise from c[3] to c[4]).
static int onarc(GEN c, GEN p, GEN tol, long prec){
  if(lg(c)==4) return 1;//Allow input of just a circle, so the return is trivially 1
  if(toleq(gel(c, 4), p, tol, prec)) return 1;//p=the start point. We have this done seperately in case rounding errors take the angle to <c[6], as this will cause issues with the shifting angle.
  pari_sp top=avma;
  GEN angle=shiftangle(radialangle(c, p, tol, prec), gel(c, 6), tol, prec);//Getting the angle in the range [c[6], c[6]+2*Pi)
  if(tolcmp(angle, gel(c, 7), tol, prec)<=0) return gc_int(top, 1);//On the arc
  return gc_int(top, 0);//Beyond the arc.
}

//p is assumed to be on the line defined by l; this checks if it is actually on the segment l
static int onseg(GEN l, GEN p, GEN tol, long prec){
  if(lg(l)==3) return 1;//Allow input of a line, so return is trivially 1
  if(typ(p)==t_INFINITY){//p is the point at oo
    if(!gequal0(gel(l, 5)) || gequalm1(gel(l, 6))) return 1;//oo is an endpoint OR the seg passes through oo
    return 0;//If not, does not pass through oo
  }
  pari_sp top=avma;
  //Okay, now p is not oo and l is a line segment
  if(typ(gel(l, 1))==t_INFINITY){//Vertical line
    if(typ(gel(l, 3))==t_INFINITY){//Start point in oo
      if(equali1(gel(l, 5))){//p must lie BELOW l[4]
          if(tolcmp(imag_i(p), imag_i(gel(l, 4)), tol, prec)<=0) return gc_int(top, 1);//Lies below l[4]
          return gc_int(top, 0);//lies above l[4]
      }
      //p must lie ABOVE l[4]
      if(tolcmp(imag_i(p), imag_i(gel(l, 4)), tol, prec)>=0) return gc_int(top, 1);//Lies above l[4]
      return gc_int(top, 0);//lies below l[4]
    }
    if(typ(gel(l, 4))==t_INFINITY){//End point is oo
      if(equali1(gel(l, 5))){//p must lie ABOVE l[3]
          if(tolcmp(imag_i(p), imag_i(gel(l, 3)), tol, prec)>=0) return gc_int(top, 1);//Lies above l[3]
          return gc_int(top, 0);//lies below l[3]
      }
      //p must lie BELOW l[3]
      if(tolcmp(imag_i(p), imag_i(gel(l, 3)), tol, prec)<=0) return gc_int(top, 1);//Lies below l[3]
      return gc_int(top, 0);//lies above l[3]
    }
    //Start and end points are finite
    int i1=tolcmp(imag_i(gsub(p, gel(l, 3))), gen_0, tol, prec);//sign of imag(p)-imag(l[3])
    int i2=tolcmp(imag_i(gsub(p, gel(l, 4))), gen_0, tol, prec);//sign of imag(p)-imag(l[4])
    set_avma(top);
    if(i1==0 || i2==0) return 1;//endpoint
    if(i1==i2){//p on the same side of l[3] and l[4], so return 1 iff l passes through oo
        if(gequal1(gel(l, 6))) return 0;//Not through oo
        return 1;//through oo
    }
    //p is between l[3] and l[4], so return 1 iff l does not pass through oo
    if(gequal1(gel(l, 6))) return 1;//not through oo
    return 0;//through oo
  }
  //Non-vertical line
  if(typ(gel(l, 3))==t_INFINITY){//Start point in oo
    if(equali1(gel(l, 5))){//p must lie LEFT OF l[4]
      if(tolcmp(real_i(p), real_i(gel(l, 4)), tol, prec)<=0) return gc_int(top, 1);//Lies left of l[4]
      return gc_int(top, 0);//lies right of  l[4]
    }
    //p must lie RIGHT OF l[4]
    if(tolcmp(real_i(p), real_i(gel(l, 4)), tol, prec)>=0) return gc_int(top, 1);//Lies right of l[4]
    return gc_int(top, 0);//lies left of l[4]
  }
  if(typ(gel(l, 4))==t_INFINITY){//End point is oo
    if(equali1(gel(l, 5))){//p must lie RIGHT OF l[3]
      if(tolcmp(real_i(p), real_i(gel(l, 3)), tol, prec)>=0) return gc_int(top, 1);//Lies right of l[3]
      return gc_int(top, 0);//lies below l[3]
    }
    //p must lie LEFT OF l[3]
    if(tolcmp(real_i(p), real_i(gel(l, 3)), tol, prec)<=0) return gc_int(top, 1);//Lies left of l[3]
    return gc_int(top, 0);//lies above l[3]
  }
  //Start and end points are finite
  int i1=tolcmp(real_i(gsub(p, gel(l, 3))), gen_0, tol, prec);//sign of real(p)-real(l[3])
  int i2=tolcmp(real_i(gsub(p, gel(l, 4))), gen_0, tol, prec);//sign of real(p)-real(l[4])
  set_avma(top);
  if(i1==0 || i2==0) return 1;//endpoint
  if(i1==i2){//p on the same side of l[3] and l[4], so return 1 iff l passes through oo
    if(gequal1(gel(l, 6))) return 0;//Not through oo
    return 1;//through oo
  }
  //p is between l[3] and l[4], so return 1 iff l does not pass through oo
  if(gequal1(gel(l, 6))) return 1;//not through oo
  return 0;//through oo
}

static GEN seg_int(GEN l1, GEN l2, GEN tol, long prec){
  pari_sp top=avma;
  GEN p=line_int(l1, l2, tol, prec);
  if(!onseg(l1, p, tol, prec)) return gc_NULL(top);
  if(!onseg(l2, p, tol, prec)) return gc_NULL(top);
  return gerepilecopy(top, mkvec(p));
}



//GEOMETRIC HELPER METHODS


//Returns the angle ang-bot in the range [0, 2*Pi)
static GEN anglediff(GEN ang, GEN bot, GEN tol, long prec){
  pari_sp top=avma;
  GEN twopi=Pi2n(1, prec);
  GEN angdiff=gmod(gsub(ang, bot), twopi);
  if(toleq(angdiff, twopi, tol, prec) || toleq0(angdiff, tol, prec)){set_avma(top);return gen_0;}
  return gerepileupto(top, angdiff);
}

//Returns 0 if c is a circle, 1 if c is a line, 2 if c is a circle arc, 3 if c is line segment, and -1 if none of the above
static int geom_check(GEN c){
  if(typ(c)!=t_VEC) return -1;
  switch(lg(c)){
    case 3: return 1;//Line
    case 4: return 0;//Circle
    case 7: return 3;//Line segment
    case 9: return 2;//Arc
  }
  return -1;
}

//Returns the default tolerance given the precision. This may need to be updated to check for 32/64 bit systems? Not sure.
GEN deftol(long prec){
  pari_sp top=avma;
  return gerepileupto(top, gtofp(powis(gen_2, 32*(2-prec)), prec));
}

//Shifts the given angle ang by multiples of 2*Pi into the range [bot, bot+2*Pi).
static GEN shiftangle(GEN ang, GEN bot, GEN tol, long prec){
  pari_sp top=avma;
  GEN diff=anglediff(ang, bot, tol, prec);
  if(gequal0(diff)){set_avma(top);return gcopy(bot);}
  return gerepileupto(top, gadd(bot, diff));
}

//Returns -1 if x<y, 0 if x==y, 1 if x>y (x, y are t_REAL). Accounts for the tolerance, so will deem x==y if they are equal up to tol AND at least one is inexact
static int tolcmp(GEN x, GEN y, GEN tol, long prec){
  if(typ(x)==t_INFINITY || typ(y)==t_INFINITY) return gcmp(x, y);//No precision concerns
  pari_sp top=avma;
  GEN d=gsub(x, y);
  if(precision(d)==0) return gc_int(top, gsigne(d));//Exact objects
  if(gcmp(gabs(d, prec), tol)<0) return gc_int(top, 0);//Within tolerance
  return gc_int(top, gsigne(d));
}

//Returns 1 if x==y, and 0 if x!=y. If x or y is not a precise objects (e.g. t_REAL), will return 1 if they are equal up to the tolerance tol.
static int toleq(GEN x, GEN y, GEN tol, long prec){
  if(typ(x)==t_INFINITY || typ(y)==t_INFINITY) return gequal(x, y);//Do oo case separately.
  pari_sp top=avma;
  GEN d=gsub(x, y);
  return gc_int(top, toleq0(d, tol, prec));//Just compare d with 0.
}

//Returns 1 if x==0, and 0 if x!=y. If x or y is not a precise objects (e.g. t_REAL or t_COMPLEX), will return 1 if they are equal up to the tolerance tol.
static int toleq0(GEN x, GEN tol, long prec){
  if(gequal0(x)) return 1;//Deemed equal already
  if(precision(x)==0) return 0;//Exact object, and already checked if it's 0
  pari_sp top=avma;
  if(gcmp(gabs(x, prec), tol)<0) return gc_int(top, 1);//Within tolerance
  return gc_int(top, 0);//Not within tolerance
}

