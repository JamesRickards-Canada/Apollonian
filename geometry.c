/*These methods implement Mobius maps acting on circles/lines. There are no restrictions on the types of numbers that can appear; real, rational, and integral are all OK.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"

/*STATIC DECLARATIONS*/

/*SECTION 1: BASIC GEOMETRIC OPERATIONS*/
static GEN circle_fromcp(GEN centre, GEN p, long prec);
static GEN circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol);
static GEN line_fromsp(GEN s, GEN p);
static GEN line_frompp(GEN p1, GEN p2, GEN tol);
static GEN line_int(GEN l1, GEN l2, GEN tol);
static GEN midpoint(GEN p1, GEN p2);
static GEN perpbis(GEN p1, GEN p2, GEN tol);
static GEN slope(GEN p1, GEN p2, GEN tol);

/*SECTION 2: MOBIUS*/
static int geom_check(GEN x);
static GEN mobius_i(GEN M, GEN x, GEN tol);
static GEN mobius_circle(GEN M, GEN c, GEN tol);
static GEN mobius_line(GEN M, GEN l, GEN tol);
static GEN mobius_point(GEN M, GEN x, GEN tol);

/*SECTION 3: TOLERANCE*/
static GEN toldiv(GEN a, GEN b, GEN tol);
static int toleq0(GEN x, GEN tol);


/*MAIN BODY*/

/* FORMATTING
LINE:
    [slope, intercept]
    y=slope*x+intercept unless slope=oo, where it is x=intercept instead.
CIRCLE:
    [centre, radius, curvature]
    Radius > 0, but curvature can be < 0 to mean the "outside" is the inside. This does NOT get preserved under mobius.
TOLERANCE:
    We deem two points as equal if they are equal up to tolerance. The default value is half the precision, i.e. if we work with 38-digit precision, then they are equal if 19 digits agree.
*/


/*SECTION 1: BASIC GEOMETRIC OPERATIONS*/

/*Circle with centre centre passing through p. This must be a circle, so p!=cent, and neither can be infinite. Not stack clean or gerepileupto suitable.*/
static GEN
circle_fromcp(GEN centre, GEN p, long prec)
{
  GEN abspmcent = gabs(gsub(p, centre), prec);
  return mkvec3(centre, abspmcent, ginv(abspmcent));
}

/*Circle through 3 points (with one allowed to be oo or all 3 collinear, making a line instead). Not stack clean or gerepileupto suitable.*/
static GEN
circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol)
{
  if (typ(p1) == t_INFINITY) return line_frompp(p2, p3, tol);/*Lines*/
  if (typ(p2) == t_INFINITY) return line_frompp(p1, p3, tol);
  if (typ(p3) == t_INFINITY) return line_frompp(p1, p2, tol);
  GEN l1 = perpbis(p1, p2, tol), l2 = perpbis(p1, p3, tol);
  GEN centre = line_int(l1, l2, tol);/*The centre is the intersection of the perpendicular bisectors.*/
  if (typ(centre) == t_INFINITY) return line_frompp(p1, p2, tol);/*p1, p2, p3 collinear*/
  return circle_fromcp(centre, p1, lg(tol));
}

/*The line through p with slope s, which can be oo. Not stack clean or gerepileupto suitable.*/
static GEN
line_fromsp(GEN s, GEN p)
{
  if (typ(s) == t_INFINITY) retmkvec2(mkoo(), greal(p));
  GEN strealp = gmul(s, real_i(p));
  return mkvec2(s, gsub(imag_i(p), strealp));
}

/*Line through two points. Not stack clean or gerepileupto suitable.*/
static GEN
line_frompp(GEN p1, GEN p2, GEN tol)
{
  return line_fromsp(slope(p1, p2, tol), p1);
}

/*Intersection point of two lines, infinite if parallel or equal. Not stack clean or gerepileupto suitable.*/
static GEN
line_int(GEN l1, GEN l2, GEN tol)
{
  GEN s1 = gel(l1, 1), m1 = gel(l1, 2), s2 = gel(l2, 1), m2 = gel(l2, 2);
  if (toleq(s1, s2, tol)) return mkoo();/*Parallel or equal*/
  if (typ(s1) == t_INFINITY) {/*l1 vertical, l2 is not.*/
    GEN ypart = gmul(s2, m1);
    return mkcomplex(m1, gadd(ypart, m2));
  }
  if (typ(s2) == t_INFINITY) {/*l2 vertical, l1 is not.*/
    GEN ypart = gmul(s1, m2);
    return mkcomplex(m2, gadd(ypart, m1));
  }
  GEN x = gdiv(gsub(m2, m1), gsub(s1, s2));
  GEN ypart = gmul(s1, x);
  return mkcomplex(x, gadd(ypart, m1));
}

/*Midpoint of p1 and p2. Not stack clean.*/
static GEN
midpoint(GEN p1, GEN p2)
{
  return gdivgs(gadd(p1, p2), 2);
}

/*Perpendicular bisector of distinct points. Not stack clean or gerepileupto suitable.*/
static GEN
perpbis(GEN p1, GEN p2, GEN tol)
{
  return line_fromsp(toldiv(gen_m1, slope(p1, p2, tol), tol), midpoint(p1, p2));
}

/*The slope of the line through p1, p2, which are not oo. Not stack clean*/
static GEN
slope(GEN p1, GEN p2, GEN tol)
{
  GEN p2mp1 = gsub(p2, p1);
  return toldiv(imag_i(p2mp1), real_i(p2mp1), tol);
}


/*SECTION 2: MOBIUS*/

/*Returns 0 if x is a point, 1 if x is a line, 2 if x is a circle, and -1 if none of the above*/
static int
geom_check(GEN x)
{
  switch (typ(x)) {
    case t_VEC: break;
    case t_INT: case t_FRAC: case t_REAL: case t_COMPLEX: case t_INFINITY: return 0;
    default: return -1;
  }
  switch (lg(x)) {
    case 3: return 1;/*Line*/
    case 4: return 2;/*Circle*/
  }
  return -1;
}

/*Mobius transformation with default tolerance.*/
GEN
mobius(GEN M, GEN x, long prec)
{
  pari_sp av = avma;
  GEN tol = deftol(prec);
  return gerepilecopy(av, mobius_i(M, x, tol));
}

/*Returns Mx, where x is a point, circle, or line. Not stack clean or gerepileupto suitable.*/
static GEN
mobius_i(GEN M, GEN x, GEN tol)
{
  int i = geom_check(x);
  switch (i) {
    case 0: return mobius_point(M, x, tol);
    case 1: return mobius_line(M, x, tol);
    case 2: return mobius_circle(M, x, tol);
  }
  pari_err_TYPE("Please input a point/circle/line", x);
  return gen_0;/*So the compiler does not complain.*/
}

/*Mobius map acting on a circle. Not stack clean or gerepileupto suitable.*/
static GEN
mobius_circle(GEN M, GEN c, GEN tol)
{
  GEN ctr = gel(c, 1), r = gel(c, 2);
  GEN p1 = mobius_point(M, gadd(ctr, r), tol);
  GEN p2 = mobius_point(M, gadd(ctr, mulcxI(r)), tol);
  GEN p3 = mobius_point(M, gsub(ctr, r), tol);
  return circle_fromppp(p1, p2, p3, tol);
}

/*Mobius map acting on a line. Not stack clean or gerepileupto suitable.*/
static GEN
mobius_line(GEN M, GEN l, GEN tol)
{
  GEN p1, p2, p3;
  GEN m = gel(l, 1), b = gel(l, 2);
  if (typ(m) == t_INFINITY) {/*Vertical line*/
    p1 = mobius_point(M, b, tol);
    p2 = mobius_point(M, mkcomplex(b, gen_1), tol);
    p3 = mobius_point(M, mkcomplex(b, gen_m1), tol);
  }
  else {/*Non-vertical line*/
    GEN slopeIp1 = mkcomplex(gen_1, m);
    GEN p1base = mulcxI(b);
    GEN p2base = gadd(p1base, slopeIp1);
    GEN p3base = gadd(p2base, slopeIp1);
    p1 = mobius_point(M, p1base, tol);
    p2 = mobius_point(M, p2base, tol);
    p3 = mobius_point(M, p3base, tol);
  }
  return circle_fromppp(p1, p2, p3, tol);
}

/*Mobius map acting on a point. Not stack clean or gerepileupto suitable.*/
static GEN
mobius_point(GEN M, GEN x, GEN tol)
{
  if (typ(x) == t_INFINITY) return toldiv(gcoeff(M, 1, 1), gcoeff(M, 2, 1), tol);
  GEN numer = gadd(gmul(gcoeff(M, 1, 1), x), gcoeff(M, 1, 2));
  GEN denom = gadd(gmul(gcoeff(M, 2, 1), x), gcoeff(M, 2, 2));
  return toldiv(numer, denom, tol);
}


/*SECTION 3: TOLERANCE*/

/*Returns the default tolerance given the precision, which is saying that x==y if they are equal up to half of the precision.*/
GEN
deftol(long prec)
{
  return real2n((BITS_IN_LONG >> 1)*(2 - prec), prec);
}

/*Divides a and b, and allows for oo and division by 0, with this checked up to tolerance. Returns oo for x/0 and oo/x (any x, even 0 or oo), 0 for x/oo with x!=oo, and treats all infinities as +oo.*/
static GEN
toldiv(GEN a, GEN b, GEN tol)
{
  if (toleq0(b, tol)) return mkoo();
  if (typ(a) == t_INFINITY) return mkoo();
  if (typ(b) == t_INFINITY) return gen_0;
  if (toleq0(a, tol)) return gen_0;
  return gdiv(a, b);
}

/*If x is of type t_INT or t_FRAC, returns 1 iff x==0. Otherwise, x must be of type t_REAL or t_COMPLEX, and returns 1 iff x=0 up to tolerance tol.*/
static int
toleq0(GEN x, GEN tol)
{
  switch (typ(x)) {
    case t_REAL:
      if (abscmprr(x, tol) < 0) return 1;/*|x|<tol*/
      return 0;
    case t_COMPLEX:;
      long i;
      for (i = 1; i <= 2; i++) {
        switch (typ(gel(x, i))) {
          case t_REAL:
            if (abscmprr(gel(x, i), tol) >= 0) return 0;/*Too large*/
            break;
          case t_INT:
            if (signe(gel(x, i))) return 0;
            break;
          case t_FRAC:/*Fraction component, cannot be 0*/
            return 0;
          default:/*Illegal input*/
            pari_err_TYPE("Tolerance equality only valid for type t_INT, t_FRAC, t_REAL, t_COMPLEX", x);
        }
      }
      return 1;/*We passed*/
    case t_INT:/*Given exactly*/
      return !signe(x);
    case t_FRAC: case t_INFINITY:/*t_FRAC and t_INFINITY cannot be 0*/
      return 0;
  }
  pari_err_TYPE("Tolerance equality only valid for type t_INT, t_FRAC, t_REAL, t_COMPLEX", x);
  return 0;/*So that there is no warning*/
}

/*Returns 1 if x==y up to tolerance tol, where two infinities are declared as equal (ignores sign). If x and y are both t_INT/t_FRAC, will only return 1 if they are exactly equal.*/
int
toleq(GEN x, GEN y, GEN tol)
{
  pari_sp av = avma;
  if (typ(x) == t_INFINITY) {
    if (typ(y) == t_INFINITY) return 1;
    return 0;
  }
  if (typ(y) == t_INFINITY) return 0;
  GEN d = gsub(x, y);
  return gc_int(av, toleq0(d, tol));/*Just compare d with 0.*/
}

