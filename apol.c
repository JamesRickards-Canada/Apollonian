/*Basic methods to deal with Apollonian circle packings. This part of the package is not designed to be as efficient as possible: it should be good, but we use pari GENS which are slower than C longs. Efficient methods should be written in C and placed in apol_fast. We allow integral and real packings, but NOT fractions (and cannot mix real and integral).*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"


/*STATIC DECLARATIONS*/

/*SECTION 1: BASIC METHODS*/
static int apol_check_primitive(GEN v);
static GEN apol_fix(GEN v, int *isint, long prec);
static void apol_move_1i(GEN v, int ind);
static void apol_move_1r(GEN v, int ind);
static void apol_move_batchi(GEN v, GEN bat);
static void apol_move_batchr(GEN v, GEN bat);
static GEN apol_red0(GEN v, int seq, void (*move1)(GEN, int), int (*cmp)(GEN, GEN));
static GEN apol_red_partial0(GEN v, long maxsteps, void (*move1)(GEN, int), int (*cmp)(GEN, GEN));

/*SECTION 2: CREATION OF ACPS*/
static GEN apol_make_n(GEN q, GEN n, int red);


static int ismp(GEN x);



static GEN apol_search_bound(GEN v, GEN bound, int countsymm, void *info, GEN (*getdata)(GEN, int, GEN, void*, int), GEN (*nextquad)(GEN, int, void*), GEN (*retquad)(GEN), int overridestrip);
static GEN apol_search_depth(GEN v, int depth, GEN bound, void *info, GEN (*getdata)(GEN, int, GEN, void*, int), GEN (*nextquad)(GEN, int, void*), GEN (*retquad)(GEN));
static GEN apol_circles_getdata(GEN vdat, int ind, GEN reps, void *nul, int state);
static int circequaltol(GEN c1, GEN c2, GEN tol, long prec);
static GEN apol_thirdtangent_line1(GEN circ1, GEN circ2, GEN c3, GEN c4, int right);
static GEN apol_circles_nextquad(GEN vdat, int ind, void* nul);
static GEN apol_circles_retquad(GEN vdat);
static GEN apol_generic_nextquad(GEN vdat, int ind, void* nul);
static GEN apol_generic_retquad(GEN vdat);
static GEN apol_curvatures_getdata(GEN vdat, int ind, GEN reps, void *nul, int state);
static GEN apol_find_getdata(GEN vdat, int ind, GEN reps, void *N, int state);
static GEN apol_layer_getdata(GEN vdat, int ind, GEN reps, void *nul, int state);
static GEN apol_layer_nextquad(GEN vdat, int ind, void *maxlayers);
static GEN apol_primes_getdata(GEN vdat, int ind, GEN reps, void *nul, int state);
static GEN apol_primes_layer_getdata(GEN vdat, int ind, GEN reps, void *nul, int state);

static GEN apol_stairs_nextquad(GEN vdat, int ind, void *info);
static GEN apol_stairs_getdata(GEN vdat, int ind, GEN reps, void *info, int state);


static GEN ZV_copy(GEN v);

/*MAIN BODY*/


/*SECTION 1: BASIC METHODS*/

/*Returns the possible obstructions modulo 24 of a primitive ACP, sorted lexicographically.*/
GEN
apol_admissiblesets()
{
  pari_sp av = avma;
  GEN ret = cgetg(7, t_VEC);
  gel(ret, 1) = mkvecsmalln(6, 0L, 1L, 4L, 9L, 12L, 16L);
  gel(ret, 2) = mkvecsmalln(6, 0L, 5L, 8L, 12L, 20L, 21L);
  gel(ret, 3) = mkvecsmalln(6, 0L, 4L, 12L, 13L, 16L, 21L);
  gel(ret, 4) = mkvecsmalln(6, 0L, 8L, 9L, 12L, 17L, 20L);
  gel(ret, 5) = mkvecsmalln(8, 3L, 6L, 7L, 10L, 15L, 18L, 19L, 22L);
  gel(ret, 6) = mkvecsmalln(8, 2L, 3L, 6L, 11L, 14L, 15L, 18L, 23L);
  long i;
  for (i = 1; i <= 6; i++) gel(ret, i) = zv_to_ZV(gel(ret, i));/*Lazy solution.*/
  return gerepilecopy(av, ret);
}

/*Checks if v gives 4 circles that generate an Apollonian circle packing, up to tolerance. Returns 1 if it is, 0 if not.*/
int
apol_check(GEN v, long prec)
{
  pari_sp av = avma;
  if (typ(v) != t_VEC || lg(v) != 5) pari_err_TYPE("Must be a length 4 vector", v);
  long i;
  for (i = 1; i <= 4; i++) {
	if (!ismp(gel(v, i))) pari_err_TYPE("Each entry must be integral or real.", v);
  }
  GEN L = gen_0, R = gen_0;
  for (i = 1; i <= 4; i++) L = mpadd(L, mpsqr(gel(v, i)));
  L = mpshift(L, 1);
  for (i = 1; i <= 4; i++) R = mpadd(R, gel(v, i));
  R = mpsqr(R);
  return gc_int(av, toleq(L, R, deftol(prec)));
}

/*Checks if v is a primitive integral packing. Returns 1 if it is, 0 if not.*/
static int
apol_check_primitive(GEN v)
{
  pari_sp av = avma;
  int isint;
  GEN w = apol_fix(v, &isint, 3);/*Don't care about the precision.*/
  if (!w || !isint) return gc_int(av, 0);
  GEN g = ZV_content(w);
  return gc_int(av, equali1(g));
}

/*Returns chi of the ACP. We start by finding a pair of coprime tangent circles.*/
long
apol_chi(GEN v)
{
  pari_sp av = avma;
  if (!apol_check_primitive(v)) pari_err_TYPE("must be a primitive integral packing", v);
  GEN vcop = shallowcopy(v);
  GEN a = gel(vcop, 1), b;
  if (gequal0(a)) return gc_long(av, 1);/*Strip packing is [0, 0, 1, 1].*/
  if (equali1(gcdii(a, gel(vcop, 2)))) b = gel(vcop, 2);
  else if (equali1(gcdii(a, gel(vcop, 3)))) b = gel(vcop, 3);
  else if (equali1(gcdii(a, gel(vcop, 4)))) b = gel(vcop, 4);
  else {/*Do S3 and S4 until we get a coprime one.*/
	for (;;) {
	  apol_move_1i(vcop, 3);
	  if (equali1(gcdii(a, gel(vcop, 3)))) { b = gel(vcop, 3); break; }
	  apol_move_1i(vcop, 4);
	  if (equali1(gcdii(a, gel(vcop, 4)))) { b = gel(vcop, 4); break; }
	}
  }
  long m4 = Mod4(a);/*a can be negative so need Mod4 not mod4.*/
  if (m4 <= 1) return gc_long(av, kronecker(b, a));
  else if (m4 == 2) return gc_long(av, kronecker(negi(b), shifti(a, -1)));
  return gc_long(av, kronecker(shifti(b, 1), a));
}

/*Given three curvatures, finds the Descartes quadruple containing them. We pick the smaller of the two possible curvatures, and sort the output.*/
GEN
apol_complete(GEN a, GEN b, GEN c, long prec)
{
  pari_sp av = avma;
  long t = typ(a);
  if (t == t_VEC || t == t_COL) {
	c = gel(a, 3); b = gel(a, 2); a = gel(a, 1);
  }
  int ta = ismp(a);
  if (!ta) pari_err_TYPE("Each input must be integral or real.", a);
  int tb = ismp(b);
  if (!tb) pari_err_TYPE("Each input must be integral or real.", b);
  if (ta != tb) {
	if (ta == 1) a = gtofp(a, prec);/*Fix to both be real.*/
	else b = gtofp(b, prec);
  }
  int tc = ismp(c);
  if (!tc) pari_err_TYPE("Each input must be integral or real.", c);
  if (tb != tc) {
	if (tb == 1) { a = gtofp(a, prec); b = gtofp(b, prec); }
	else c = gtofp(c, prec);
  }
  GEN bpc = mpadd(b, c);
  GEN bc = mpmul(b, c);
  GEN tort = mpadd(mpmul(a, bpc), bc);/*ab+ac+bc*/
  GEN tol = deftol(prec), rt;
  if (toleq0(tort, tol)) rt = gen_0;/*To account for rounding to be slightly negative.*/
  else {
	if (signe(tort) < 0) pari_err_TYPE("There does not exist a Descartes quadruple with these curvatures. Are two of them negative, or is the negative curvature too small?", mkvec3(a, b, c));
	if (typ(tort) == t_INT) {
	  GEN r;
	  rt = sqrtremi(tort, &r);
	  if (!isintzero(r)) {/*Square root not integral.*/
		rt = gsqrt(tort, prec);
		a = gtofp(a, prec);/*Must update them all.*/
		b = gtofp(b, prec);
		c = gtofp(c, prec);
	  }
	}
	else rt = sqrtr(tort);
  }
  GEN term2 = mpshift(rt, 1);/*2*sqrt(ab+ac+bc)*/
  GEN d = mpsub(mpadd(a, bpc), term2);/*The smaller of the roots.*/
  return gerepileupto(av, sort(mkvec4(a, b, c, d)));
}

/*Returns the external depth of v, i.e. the minimal number of swaps required to reach a quadruple with non-positive curvature.*/
long
apol_extdepth(GEN v, long prec)
{
  pari_sp av = avma;
  int isint;
  v = apol_fix(v, &isint, 3);
  if (!v) pari_err_TYPE("Not a Descartes quadruple", v);
  v = shallowcopy(v);
  long ind = vecindexmin(v), step = 0;
  void (*move1)(GEN, int);
  if (isint) move1 = &apol_move_1i;
  else move1 = &apol_move_1r;
  for (;;) {
    if (signe(gel(v, ind)) != 1) return gc_long(av, step);/*An index is <=0.*/
    step++;
    ind = vecindexmax(v);
    move1(v, ind);
  }
}

/*Returns NULL if not a Descartes quadruple. Otherwise, makes sure all entries are t_INT or t_REAL, and returns the updated vector. If no updates required, it is shallow and returns v without copying. Updates isint to 1 if it is integral, 0 if not (if not passed as NULL).*/
static GEN
apol_fix(GEN v, int *isint, long prec)
{
  pari_sp av = avma;
  if (typ(v) != t_VEC || lg(v) != 5) return NULL;
  int t = ismp(gel(v, 1));
  if (!t) pari_err_TYPE("Each entry must be integral or real.", v);
  int onlyint = t, mixing = 0, i;
  for (i = 2; i <= 4; i++) {
	t = ismp(gel(v, i));
	if (!t) pari_err_TYPE("Each entry must be integral or real.", v);
	if (t != onlyint) mixing = 1;/*Go back and fix*/
  }
  GEN L = gen_0, R = gen_0;
  for (i = 1; i <= 4; i++) L = mpadd(L, mpsqr(gel(v, i)));
  L = mpshift(L, 1);
  for (i = 1; i <= 4; i++) R = mpadd(R, gel(v, i));
  R = mpsqr(R);
  if (!toleq(L, R, deftol(prec))) return gc_NULL(av);/*Descartes equation not satisfied.*/
  if (!mixing) {/*Uniform type.*/
    if (isint) *isint = 2 - onlyint;
	set_avma(av);
	return v;
  }
  if (isint) *isint = 0;/*Must be real, as we mixed the two.*/
  GEN w = cgetg(5, t_VEC);
  for (i = 1; i <= 4; i++) gel(w, i) = gtofp(gel(v, i), prec);
  return gerepileupto(av, w);
}

/*Returns [S1, S2, S3, S4, K], where Si generate the Apollonian group, and K*[n,A,B,C]~=theta([A, B, C]) (see Staircase paper for theta description)*/
GEN
apol_matrices()
{
  pari_sp av = avma;
  GEN S1 = mkmat4(mkcol4s(-1, 0, 0, 0), mkcol4s(2, 1, 0, 0), mkcol4s(2, 0, 1, 0), mkcol4s(2, 0, 0, 1));
  GEN S2 = mkmat4(mkcol4s(1, 2, 0, 0), mkcol4s(0, -1, 0, 0), mkcol4s(0, 2, 1 ,0), mkcol4s(0, 2, 0, 1));
  GEN S3 = mkmat4(mkcol4s(1, 0, 2, 0), mkcol4s(0, 1, 2, 0), mkcol4s(0, 0, -1, 0), mkcol4s(0, 0, 2, 1));
  GEN S4 = mkmat4(mkcol4s(1, 0, 0, 2), mkcol4s(0, 1, 0, 2), mkcol4s(0, 0, 1, 2), mkcol4s(0, 0, 0, -1));
  GEN K = mkmat4(mkcol4s(1, -1, -1, -1), mkcol4s(0, 1, 0, 1), mkcol4s(0, 0, 0, -1), mkcol4s(0, 0, 1, 1));
  return gerepilecopy(av, mkvec5(S1, S2, S3, S4, K));
}

/*Returns the set of admissible residues modulo 24 of a primitive packing. There are 6 possible primitive sets: 
[0, 1, 4, 9, 12, 16]; primes are 1 mod 24                 1  1  1
[0, 4, 12, 13, 16, 21]; primes are 13 mod 24              1 -1  1
[0, 5, 8, 12, 20, 21]; primes are 5 mod 24                1 -1 -1
[0, 8, 9, 12, 17, 20]; primes are 17 mod 24               1  1 -1
[2, 3, 6, 11, 14, 15, 18, 23]; primes are 11, 23 mod 24  -1 -1  1   -1  1  1
[3, 6, 7, 10, 15, 18, 19, 22]; primes are 7, 19 mod 24   -1  1 -1   -1 -1 -1
You ONLY need to go to depth 3 to find which class we are in (proven by brute force check).*/
GEN
apol_mod24(GEN v)
{
  pari_sp av = avma;
  if (!apol_check_primitive(v)) pari_err_TYPE("must be a primitive integral packing", v);
  GEN tw4 = stoi(24);
  GEN orb = apol_curvatures_depth(v, 3, gen_0);/*Only need depth 3*/
  long lo = lg(orb), i;
  for (i = 1; i < lo; i++) gel(orb, i) = Fp_red(gel(orb, i), tw4);/*Reduce modulo 24.*/
  return gerepileupto(av, ZV_sort_uniq(orb));
}

/*Calls apol_move_(1/batch)(i/r), returning a clean result without affecting v. We also first check that v is a Descartes quadruple.*/
GEN
apol_move(GEN v, GEN command, long prec)
{
  pari_sp av = avma;
  int isint;
  v = apol_fix(v, &isint, prec);
  if (!v) pari_err_TYPE("v is not a Descartes quadruple.", v);
  GEN vcop = shallowcopy(v);
  long t = typ(command);
  switch (t) {
    case t_INT:
	  if (isint) apol_move_1i(vcop, itos(command));
	  else apol_move_1r(vcop, itos(command));
	  return gerepilecopy(av, vcop);
	case t_VEC:
	case t_COL:
	  command = gtovecsmall(command);
    case t_VECSMALL:
	  if (isint) apol_move_batchi(vcop, command);
	  else apol_move_batchr(vcop, command);
	  return gerepilecopy(av, vcop);
  }
  pari_err_TYPE("Input does not represent a valid move (or series of moves)", command);
  return NULL;
}

/*Replace the index ind in the integral ACP. Not gerepileupto safe, shallow, leaves garbage.*/
static void
apol_move_1i(GEN v, int ind)
{
  GEN a = gel(v, 1), b = gel(v, 2), c = gel(v, 3), d = gel(v, 4);
  switch (ind) {
	case 1:
	  gel(v, 1) = subii(shifti(addii(addii(b, c), d), 1), a);
	  return;
	case 2:
	  gel(v, 2) = subii(shifti(addii(addii(a, c), d), 1), b);
	  return;
	case 3:
	  gel(v, 3) = subii(shifti(addii(addii(a, b), d), 1), c);
	  return;
	case 4:
	  gel(v, 4) = subii(shifti(addii(addii(a, b), c), 1), d);
  }
}

/*Replace the index ind in the real ACP. Not gerepileupto safe, shallow, leaves garbage.*/
static void
apol_move_1r(GEN v, int ind)
{
  GEN a = gel(v, 1), b = gel(v, 2), c = gel(v, 3), d = gel(v, 4);
  switch (ind) {
	case 1:
	  gel(v, 1) = subrr(shiftr(addrr(addrr(b, c), d), 1), a);
	  return;
	case 2:
	  gel(v, 2) = subrr(shiftr(addrr(addrr(a, c), d), 1), b);
	  return;
	case 3:
	  gel(v, 3) = subrr(shiftr(addrr(addrr(a, b), d), 1), c);
	  return;
	case 4:
	  gel(v, 4) = subrr(shiftr(addrr(addrr(a, b), c), 1), d);
  }
}

/*Does apol_move_1i for v and bat[1], bat[2], ... The input bat needs to be a Vecsmall. Not gerepileupto safe, shallow, leaves garbage.*/
static void
apol_move_batchi(GEN v, GEN bat)
{
  long lbat = lg(bat), i;
  for (i = 1; i < lbat; i++) apol_move_1i(v, bat[i]);
}

/*Does apol_move_1r for v and bat[1], bat[2], ... The input bat needs to be a Vecsmall. Not gerepileupto safe, shallow, leaves garbage.*/
static void
apol_move_batchr(GEN v, GEN bat)
{
  long lbat = lg(bat), i;
  for (i = 1; i < lbat; i++) apol_move_1r(v, bat[i]);
}

/*Returns the intrgral quadratic form whose primitive values are v[ind]+curvatures touching the ind^th circle. The formula is [a+b,a+b+c-d, a+c] if v=[a,b,c,d] and ind=1.*/
GEN
apol_qf(GEN v, int ind)
{
  pari_sp av = avma;
  int isint;
  v = apol_fix(v, &isint, 3);
  if (!v) pari_err_TYPE("v is not a Descartes quadruple.", v);
  if (!isint) pari_err_TYPE("v is not integral", v);
  GEN is = cgetg(5, t_VECSMALL);/*Used to shift around things since ind might not be 1.*/
  is[1] = ind;
  int i;
  for (i = 2; i <= 4; i++) is[i] = is[i - 1]%4 + 1;
  GEN a = gel(v, is[1]), apb = addii(a, gel(v, is[2]));/*a+b*/
  GEN c = gel(v, is[3]), apbpc = addii(apb, c);/*a+b+c*/
  GEN D = shifti(negi(sqri(a)), 2);/*-4a^2*/
  return gerepilecopy(av, mkqfb(apb, subii(apbpc, gel(v, is[4])), addii(a, c), D));
}

/*Returns the reduction of v. If seq=1, also returns a VECSMALL of the sequence of indices swapped to reduce.*/
GEN
apol_red(GEN v, int seq, long prec)
{
  pari_sp av = avma;
  int isint;
  v = apol_fix(v, &isint, prec);
  if (!v) pari_err_TYPE("v is not a Descartes quadruple.", v);
  if (isint) return gerepileupto(av, apol_red0(v, seq, &apol_move_1i, &cmpii));
  return gerepileupto(av, apol_red0(v, seq, &apol_move_1r, &cmprr));
}

/*Reduce the ACP given the move and comparison type.*/
static GEN
apol_red0(GEN v, int seq, void (*move1)(GEN, int), int (*cmp)(GEN, GEN))
{
  pari_sp av = avma;
  v = shallowcopy(v);
  long ind;
  GEN dold;
  if (!seq) {
    do {
      ind = vecindexmax(v);
      dold = gel(v, ind);
      move1(v, ind);
    } while(cmp(gel(v, ind), dold) < 0);
	move1(v, ind);/*Must go back one!*/
    return gerepilecopy(av, v);
  }
  long i = 0, len = 40;
  GEN S = cgetg(len + 1, t_VECSMALL);
  do {
    ind = vecindexmax(v);
    dold = gel(v, ind);
    move1(v, ind);
	i++;
	if (i > len) {
	  len <<= 1;
	  S = vecsmall_lengthen(S, len);
	}
	S[i] = ind;
  } while(cmp(gel(v, ind), dold) < 0);
  S = vecsmall_shorten(S, i - 1);/*Remove last move.*/
  move1(v, ind);
  return gerepilecopy(av, mkvec2(v, S));
}

/*Reduces v, where we go ONLY at most maxsteps steps.*/
GEN
apol_red_partial(GEN v, long maxsteps, long prec)
{
  pari_sp av = avma;
  int isint;
  v = apol_fix(v, &isint, prec);
  if (!v) pari_err_TYPE("v is not a Descartes quadruple.", v);
  if (isint) return gerepileupto(av, apol_red_partial0(v, maxsteps, &apol_move_1i, &cmpii));
  return gerepileupto(av, apol_red_partial0(v, maxsteps, &apol_move_1r, &cmprr));
}

/*Reduces v, where we go ONLY at most maxsteps steps.*/
static GEN
apol_red_partial0(GEN v, long maxsteps, void (*move1)(GEN, int), int (*cmp)(GEN, GEN))
{
  pari_sp av = avma;
  if (!maxsteps) return ZV_copy(v);
  long ind, step;
  GEN dold;
  for (step = 1; step <= maxsteps; step++) {
    ind = vecindexmax(v);
    dold = gel(v, ind);
	move1(v, ind);
    if (cmp(gel(v, ind), dold) >= 0) { move1(v, ind); break; }/*Must go back one.*/
  }
  return gerepilecopy(av, v);
}

/*Returns the "type" of v, i.e. the number of resiudes modulo 24 and the smallest residue coprime to 6, which uniquely identifies it. If chi = 1, also includes the chi value.*/
GEN
apol_type(GEN v, int chi)
{
  pari_sp av = avma;
  GEN m24 = apol_mod24(v);
  long second = itos(gel(m24, 2));/*Uniquely identified by the second element*/
  set_avma(av);
  if (chi) {
	switch (second) {
      case 1: return mkvec3s(6, 1, apol_chi(v));
	  case 3: return mkvec3s(8, 11, apol_chi(v));
	  case 4: return mkvec3s(6, 13, apol_chi(v));
	  case 5: return mkvec3s(6, 5, apol_chi(v));
	  case 6: return mkvec3s(8, 7, apol_chi(v));
	  case 8: return mkvec3s(6, 17, apol_chi(v));
    }
  }
  switch (second) {
    case 1: return mkvec2s(6, 1);
	case 3: return mkvec2s(8, 11);
	case 4: return mkvec2s(6, 13);
	case 5: return mkvec2s(6, 5);
	case 6: return mkvec2s(8, 7);
	case 8: return mkvec2s(6, 17);
  }
  pari_err(e_MISC, "We didn't find one of the 6 possible admissible sets, are you sure you inputted a primitive Apollonian circle packing?");
  return gen_0;
}


/*SECTION 2: CREATION OF ACPS*/

/*Given a qfb q of discriminant -4n^2, this gives the corresponding Descartes quadruple. if pos=-1, we give the quadruple with -n, and if pos=1, we give the quadruple with +n. If red=1 we reduce, else we don't.*/
GEN
apol_make(GEN q, int pos, int red)
{
  pari_sp av = avma;
  if (typ(q) != t_QFB) pari_err_TYPE("q must be an integral binary quadratic form of discriminant -4n^2 for some integer n.", q);
  GEN D = gel(q, 4);
  if (signe(D) != -1) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  long rem;
  GEN nsqr = divis_rem(D, -4, &rem);
  if (rem) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  GEN sqrtrem;
  GEN n = sqrtremi(nsqr, &sqrtrem);
  if (!gequal0(sqrtrem)) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  n = pos ? n : negi(n);/*+/-n*/
  return gerepileupto(av, apol_make_n(q, n, red));
}

/*apol_make, but we are given the value of n, where disc(q)=-4n^2.*/
static GEN
apol_make_n(GEN q, GEN n, int red)
{
  pari_sp av = avma;/*q=[A, B, C]->[n, A-n, C-n, A+C-n-B]*/
  GEN Amn = subii(gel(q, 1), n);
  GEN v = mkvec4(n, Amn, subii(gel(q, 3), n), addii(Amn, subii(gel(q, 3), gel(q, 2))));/*The ACP*/
  if (red) v = apol_red0(v, 0, &apol_move_1i, &cmpii);
  return gerepileupto(av, ZV_sort(v));
}

/*Computes qfbnarrowlex(-4n^2), and output the ACP's created from these quadratic forms with apol_make. We only count ambiguous forms ONCE, and we take pos=sign(n).*/
GEN
apol_makeall(GEN n, int red, long prec)
{
  pari_sp av = avma;
  GEN D = shifti(negi(sqri(n)), 2);
  GEN forms = gel(qfbnarrowlex(D, prec), 3);
  long lf = lg(forms), i;
  GEN quads = vectrunc_init(lf);
  for (i = 1; i < lf; i++) {/*If we have [A, B, C] with B<0 we do not count it.*/
    GEN q = gel(forms, i);
    if (signe(gel(q, 2)) < 0) continue;
    vectrunc_append(quads, apol_make_n(q, n, red));
  }
  return gerepileupto(av, quads);
}



//SEARCHING FOR CURVATURES


/*
The outline of the searching methods is:
    We use depth first search on the Apollonian graph.
    The current partial path is stored in the variable W. Each element of W is normally a Descartes quadruple, but there may be extra information.
    We have functions to:
        Retrieve the data we want from a new element (the curvature, entire quadruple, etc.), or returning NULL if we don't want to count it (e.g. only counting quadruples with certain properties)
        Move 1 layer deeper in the tree, and return the new element to be stored in W (again, normally a Descartes quadruple)
        Retrieve the Descartes quadruple from the current element of W
    These functions are (respectively):
        static GEN getdata(GEN vdat, int ind, GEN reps, void *info, int state)
        static GEN nextquad(GEN vdat, int ind, void *info)
        static GEN retquad(GEN vdat)
        
    getdata:
        vdat is the newest node, and keeps track of the Descartes quadruple and whatever extra info we need
        1<=ind<=4 is the index of the circle we last swapped
        reps passes in the current set of data found. This is used if W stores indices of elements in reps, and is also used at the end to clean up the data (we may want to sort the final answer, or do something else, etc.)
        info tracks extra information you may want, and may not be modified. It is not touched in apol_search_(type). If not needed, pass it as NULL.
        state tracks what we are doing with this method.
            state=1 means we are adding a new piece of information potentially. Your method should compute the data, return it if you want it added, and return NULL if we don't want to count it.
            state=2 is for the initial 4 circles, and you can do something special there if you wish. We always run the method on these first 4 circles in order at the start.
            state=0 happens at the end, and is used as a wrap-up. You should return a modified version of reps that is gerepileupto compatible, e.g. sort it, just copy it, etc.
    
    nextquad:
        essentially does apol_move(vdat, ind), but formats it correctly in case vdat is not just the Descartes quadruple. Pass in ind=0 for the initial setup. Return NULL to say "too far, go backwards".
  
    retquad:
        Returns the Descartes quadruple from vdat. Normally this is just {return vdat;}, but if we keep track of more it is not.
*/


//Starting at the integral Descartes quadruple v, we search through the circle packing, finding all circles with curvature <=bound. For each such circle we compute some data (from getdata), and return the final set. If countsymm=0, we do not double count when there is symmetry, otherwise we do. For the strip packing, we force countsymm=0, since otherwise it would be infinite. It can be shown that all symmetries are realized for the reduced form (i.e. [a, b, c, d] -> [a, b, c, d] or [a, b, c, c], if this happens in a packing, it happens for the reduced form, so we can just check it there). There are situtations where we want to override the strip packing not counting symmetries (when we sort this out in nextquad), and we can pass this in as an argument.
static GEN apol_search_bound(GEN v, GEN bound, int countsymm, void *info, GEN (*getdata)(GEN, int, GEN, void*, int), GEN (*nextquad)(GEN, int, void*), GEN (*retquad)(GEN), int overridestrip){
  pari_sp top=avma, mid;
  v=apol_red(v, 0, 0);//Start by reducing v, since we want curvatures up to a given bound.
  ZV_sort_inplace(v);//May as well make the minimal curvature first.
  int ind=1;//We reuse ind to track which depth we are going towards.
  int depth=10;//Initially we go up to depth 10, but this may change if we need to go deeper.
  GEN W=zerovec(depth);//Tracks the sequence of Descartes quadruples; W[ind] is at ind-1 depth
  GEN vdat=nextquad(v, 0, info);//The first one. Stores the quadruple + data
  if(!vdat) return gerepilecopy(top, cgetg(1, t_VEC));//Our initialization failed.
  gel(W, 1)=vdat;
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of indices of the replacements
  int forward=1;//Tracks if we are going forward or not.
  long Nreps=400;//Initial size of the return.
  GEN reps=vectrunc_init(Nreps);//We may need to extend the length of reps later on.
  int isstrip=gequal0(gel(v, 1));//If we have the strip packing or not, which causes some subtle issues.
  if(isstrip){
	if(overridestrip) isstrip=0;//Override
	else countsymm=0;//Strip packing, must not count symmetry
  }
  if(countsymm){
    for(int i=1;i<=4;i++){
      if(cmpii(gel(v, i), bound)<=0){//First 4 reps
        GEN dat=getdata(vdat, i, reps, info, 2);
        if(dat) vectrunc_append(reps, dat);//We have data to store.
      }
    }
  }
  else{//Must ignore symmetry. v is already sorted.
    for(int i=1;i<=4;i++){
      if(cmpii(gel(v, i), bound)<=0 && (i==1 || !equalii(gel(v, i-1), gel(v, i)))){//First 4 reps
        GEN dat=getdata(vdat, i, reps, info, 2);
        if(dat) vectrunc_append(reps, dat);//We have data to store.
      }
    }
  }
  do{//1<=ind<=depth is assumed.
    if(gc_needed(top, 1)){//Garbage day!
      GEN oldreps=reps;
      mid=avma;
      W=gcopy(W);
      I=gcopy(I);
      reps=vectrunc_init(Nreps);
      for(long i=1;i<lg(oldreps);i++) vectrunc_append(reps, gcopy(gel(oldreps, i)));//Copying reps; can't use vectrunc_append_batch.
      gerepileallsp(top, mid, 3, &W, &I, &reps);
    }
    if(lg(reps)==Nreps){//We don't have enough space! Double the possible length of reps.
      Nreps=2*Nreps-1;
      GEN newreps=vectrunc_init(Nreps);
      vectrunc_append_batch(newreps, reps);//Append the old list
      reps=newreps;
    }
    //Now lets try to move along by updating I.
    I[ind]=forward? 1:I[ind]+1;//The new index
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat consecutively
    if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
    //At this point, we can go on with valid and new inputs.
    v=retquad(gel(W, ind));//Current Descartes quadruple
    if(isstrip || (ind==1 && !countsymm)){//Make sure we haven't seen this element before. If this happens at some point, then it ALSO happens at the reduced form!
      int goback=0;
      for(int j=1;j<I[ind];j++){
        if(equalii(gel(v, j), gel(v, I[ind]))){goback=1;break;}
      }
      if(goback){forward=0;continue;}//We've already done this element in a previous move.
    }
    GEN newvdat=nextquad(gel(W, ind), I[ind], info);//The data for the move.
    if(!newvdat){forward=0;continue;}//We must go backwards.
    GEN newcurv=gel(retquad(newvdat), I[ind]);
    if(cmpii(newcurv, bound)>0){forward=0;continue;}//Must go back, elt too big
    if(ind==1 && !countsymm){//The replacement being the same can ONLY happen for the reduced form, even for the strip packing.
      if(equalii(gel(v, I[ind]), newcurv)){forward=0;continue;}//Must go back, same thing (e.g. the strip packing)
    }
    GEN dat=getdata(newvdat, I[ind], reps, info, 1);//Retrieve the data
    if(dat) vectrunc_append(reps, dat);//Add the new piece of data if we are allowed.
    ind++;
    if(ind>depth){//Reached maximum depth, must extend.
      int newdepth=2*depth;
      GEN newW=zerovec(newdepth), newI=vecsmall_ei(newdepth, 1);
      for(long i=1;i<=depth;i++){gel(newW, i)=gel(W, i);newI[i]=I[i];}//Copy them over.
      W=newW;
      I=newI;
      depth=newdepth;
    }
    gel(W, ind)=newvdat;
    forward=1;
  }while(ind>0);
  return gerepileupto(top, getdata(NULL, 0, reps, NULL, 0));
}

//Starting at the integral Descartes quadruple v, we search through the circle packing, finding all circles of depth at most depth (maximum number of circle replacements), optionally also keeping those with curvature <=curvature.
static GEN apol_search_depth(GEN v, int depth, GEN bound, void *info, GEN (*getdata)(GEN, int, GEN, void*, int), GEN (*nextquad)(GEN, int, void*), GEN (*retquad)(GEN)){
  pari_sp top=avma, mid;
  int ind=1;//We reuse ind to track which depth we are going towards.
  GEN W=zerovec(depth);//Tracks the sequence of Descartes quadruples; W[ind] is at ind-1 depth
  GEN vdat=nextquad(v, 0, info);//The first one. Stores the quadruple + data
  gel(W, 1)=vdat;
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1, usebound=1-gequal0(bound);;//Tracks if we are going forward or not, and whether or not to use the bound.
  long Nreps=400;
  GEN reps=vectrunc_init(Nreps);//We may need to extend the length of reps later on.
  for(int i=1;i<=4;i++){
    if(!usebound || gcmp(gel(v, i), bound)<=0){//First 4 reps
      GEN dat=getdata(vdat, i, reps, info, 2);
      if(dat) vectrunc_append(reps, dat);//We have data to store.
    }
  }  
  do{//1<=ind<=depth is assumed.
    if(gc_needed(top, 2)){//Garbage day!
      GEN oldreps=reps;
      mid=avma;
      W=gcopy(W);
      I=gcopy(I);
      reps=vectrunc_init(Nreps);
      for(long i=1;i<lg(oldreps);i++) vectrunc_append(reps, gcopy(gel(oldreps, i)));//Copying reps; can't use vectrunc_append_batch.
      gerepileallsp(top, mid, 3, &W, &I, &reps);
    }
    if(lg(reps)==Nreps){//We don't have enough space! Double the possible length of reps.
      Nreps=2*Nreps-1;
      GEN newreps=vectrunc_init(Nreps);
      vectrunc_append_batch(newreps, reps);//Append the old list
      reps=newreps;
    }
    I[ind]=forward? 1:I[ind]+1;//The new index
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
    if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
    //At this point, we can go on with valid and new inputs.
    v=retquad(gel(W, ind));//Current Descartes quadruple
    GEN newvdat=nextquad(gel(W, ind), I[ind], info);//Make the move
    if(!newvdat){forward=0;continue;}//We must go backwards.
    GEN newcurv=gel(retquad(newvdat), I[ind]);//The new element
    if(usebound && gcmp(newcurv, bound)>0){forward=0;continue;}//Must go back, new curvature too big
    GEN dat=getdata(newvdat, I[ind], reps, info, 1);//Retrieve the data
    if(dat) vectrunc_append(reps, dat);//Add the new piece of data if we are allowed.
    ind++;
    if(ind>depth){forward=0;continue;}//Reached maximum depth, must go back
    gel(W, ind)=newvdat;
    forward=1;
  }while(ind>0);
  return gerepileupto(top, getdata(NULL, 0, reps, NULL, 0));
}


//vdat=[v, [indices in reps of the circles]]. We want to return the circle, i.e. [centre, radius, curvature]. radius>=0 always, but curvature may be negative.
static GEN apol_circles_getdata(GEN vdat, int ind, GEN reps, void *nul, int state){
  if(state==0) return gcopy(reps);
  GEN v=gel(vdat, 1);//The new quadruple
  if(state==2){//Initial 4 circles
    switch(ind){
      case 1:;
        GEN c1=gel(v, 1);//First curvature
		if(gequal0(c1)) return mkvec2(gen_0, gen_0);//Line!
		if(gsigne(c1)==-1) return mkvec3(gen_0, gdivsg(-1, c1), c1);//Outer circle. First circle has negative curvature necessarily, hence why r1=-1/c1
		return mkvec3(gen_0, gdivsg(1, c1), c1);//First circle has positive curvature, may be a full plane packing!
      case 2:;
        GEN c2=gel(v, 2);
		if(gequal0(c2)){//Line!
		  if(lg(gel(reps, 1))==3){//First circle was also a line!
			return mkvec2(gen_0, gdivsg(2, gel(v, 3)));
		  }
		  //First circle was a circle. Place line below.
		  return mkvec2(gen_0, gsub(imag_i(gmael(reps, 1, 1)), gmael(reps, 1, 2)));
		}
		//c2 is a circle.
        GEN r2=gdivsg(1, c2);
		if(lg(gel(reps, 1))==3){//First circle was a line! It must be the x-axis
		  return mkvec3(mkcomplex(gen_0, r2), r2, c2);
		}
		if(gsigne(gmael(reps, 1, 3))==-1){//First circle was outer one.
          return mkvec3(mkcomplex(gen_0, gsub(gmael(reps, 1, 2), r2)), r2, c2);//first inner circle, placed vertically at the top. Note gmael(reps, 1, 2)=r1>0.
		}
		//First circle was not the outer one.
		return mkvec3(mkcomplex(gen_0, gadd(gmael(reps, 1, 2), r2)), r2, c2);
      case 3:
        return apol_thirdtangent(gel(reps, 1), gel(reps, 2), gel(v, 3), gel(v, 4), 0);//Third circle goes left.
      case 4:
        return apol_thirdtangent(gel(reps, 2), gel(reps, 3), gel(v, 4), gel(v, 1), 0);//Fourth circle is left of circ2 ->circ3.
    }
  }
  //Now we are at a normal place, i.e state=1.
  GEN prevind=gel(vdat, 2);
  int is[3]={0, 0, 0};
  int is_ind=0;
  for(int i=1;i<=4;i++){
    if(ind==i) continue;
    is[is_ind]=i;
    is_ind++;
  }//is are the three non-ind indices in {1, 2, 3, 4}.
  GEN oldcirc1=gel(reps, prevind[is[0]]);//One of the old circles
  GEN oldcirc2=gel(reps, prevind[is[1]]);//Another of the old circles
  GEN newc=gel(v, ind);
  GEN newcirc=apol_thirdtangent(oldcirc1, oldcirc2, newc, gel(v, is[2]), 1);//The new circle, if it is to the right of oldcirc1 ->oldcirc2.
  
  GEN prevcirc=gel(reps, prevind[ind]);//The circle we are "replacing"
  GEN tol=deftol(3);//Default tolerance
  if(circequaltol(newcirc, prevcirc, tol, 3)) newcirc=apol_thirdtangent(oldcirc1, oldcirc2, newc, gel(v, is[2]), 0);//If the two curvatures were the same, this could trigger.
  else{
    GEN oldcirc3=gel(reps, prevind[is[2]]);//The unused old circle. Our newcirc must be tangent to it.
    GEN rsums=gsqr(gadd(gel(oldcirc3, 2), gel(newcirc, 2)));//(r1+r2)^2
    GEN dcentres=gnorm(gsub(gel(oldcirc3, 1), gel(newcirc, 1)));//dist(centres)^2
    if(!toleq(rsums, dcentres, tol)) newcirc=apol_thirdtangent(oldcirc1, oldcirc2, newc, gel(v, is[2]), 0);//Must be the other side.
  }
  prevind[ind]=lg(reps);//Updating the location of the circle.
  return newcirc;
}

static int circequaltol(GEN c1, GEN c2, GEN tol, long prec){//Returns 1 if the circles are equal, 0 else.
  if(lg(c1)==3){
	if(lg(c2)!=3) return 0;//1 line 1 circle
	if(toleq(gel(c1, 2), gel(c2, 2), tol)) return 1;//Equal up to tolerance!
	return 0;//Not equal up to tolerance
  }
  if(lg(c2)==3) return 0;//1 circle 1 line
  if(!toleq(gel(c1, 1), gel(c2, 1), tol)) return 0;
  if(!toleq(gel(c1, 3), gel(c2, 3), tol)) return 0;
  return 1;//No need to check radii if curvatures are okay.
}

//Store a circle as [centre, radius, curvature]. Given two tangent circles and a third curvature, this finds this third circle that is tangent to the first two. For internal tangency, we need a positive radius and negative curvature. There are always 2 places to put the circle: left or right of the line from circ1 to circ2. If right=1, we put it right, else we put it left. c4 is one of the curvatures to complete an Apollonian quadruple (supplying it allows us to always work with exact numbers in the case of integral ACPs). In the case of a horizontal line, we ASSUME that it has slope 0.
GEN apol_thirdtangent(GEN circ1, GEN circ2, GEN c3, GEN c4, int right){
  if(lg(circ1)==3) return apol_thirdtangent_line1(circ1, circ2, c3, c4, right);//circ1 is a line!
  if(lg(circ2)==3) return apol_thirdtangent_line1(circ2, circ1, c3, c4, 1-right);//circ2 is a line!
 pari_sp top=avma;//Now neither circ1 or circ2 is a line.
  if(gequal0(c3)){//c3 is a line.
	if(right==1){//Below
	  return gerepilecopy(top, mkvec2(gen_0, gsub(imag_i(gel(circ1, 1)), gel(circ1, 2))));
	}
	return gerepilecopy(top, mkvec2(gen_0, gadd(imag_i(gel(circ1, 1)), gel(circ1, 2))));//Above
  }
  //Now circ1, circ2, and c3 are all circles.
  long prec=3;//Does not matter, things here are normally exact, and when not prec=3 is fine.
  //The centres form a triangle with sides r1+r2, r1+r3, r2+r3, or -r1-r2, -r1-r3, r2+r3 (if internal tangency, where r1<0). Let theta be the angle at the centre of c1.
  GEN c1=gel(circ1, 3), c2=gel(circ2, 3);//Curvatures
  GEN c1pc2=gadd(c1, c2), c1pc3=gadd(c1, c3);
  GEN denom=gmul(c1pc2, c1pc3);
  GEN costheta=gsubsg(1, gdiv(gmulsg(2, gsqr(c1)), denom));//1-2c1^2/((c1+c2)(c1+c3))
  GEN sintheta=gdiv(gabs(gmul(c1, gadd(c1pc2, gsub(c3, c4))), prec), denom);//|c1(c1+c2+c3-c4)|/((c1+c2)(c1+c3))
  GEN r1=gel(circ1, 2), r2=gel(circ2, 2), r3=gdivsg(1, c3);//The radii
  if(gsigne(c1)<0) r1=gneg(r1);
  if(gsigne(c2)<0) r2=gneg(r2);//Correcting r1 and r2 to be negative.
  GEN x1=real_i(gel(circ1, 1)), y1=imag_i(gel(circ1, 1)), x2=real_i(gel(circ2, 1)), y2=imag_i(gel(circ2, 1));
  //We need alpha, the angle to the centre of c2 from the centre of c1.
  GEN r1pr2=gadd(r1, r2);
  GEN cosalpha=gdiv(gsub(x2, x1), r1pr2);
  GEN sinalpha=gdiv(gsub(y2, y1), r1pr2);//cos and sin of alpha.
  //If right=1, we have angle alpha-theta, and if right=0, we have angle alpha+theta. We derive the new coordinates from the (co)sine addition/subtraction formulae.
  GEN relcos, relsin;
  if(right==1){//cos(alpha-theta), sin(alpha-theta)
    relcos=gadd(gmul(cosalpha, costheta), gmul(sinalpha, sintheta));
    relsin=gsub(gmul(sinalpha, costheta), gmul(cosalpha, sintheta));
  }
  else{//cos(alpha+theta), sin(alpha+theta)
    relcos=gsub(gmul(cosalpha, costheta), gmul(sinalpha, sintheta));
    relsin=gadd(gmul(sinalpha, costheta), gmul(cosalpha, sintheta));
  }
  GEN r1pr3=gadd(r1, r3);
  GEN x=gadd(x1, gmul(r1pr3, relcos));
  GEN y=gadd(y1, gmul(r1pr3, relsin));
  return gerepilecopy(top, mkvec3(mkcomplex(x, y), r3, c3));
}

//apol_thirdtangent when circ1 is actually a line. We assume that it is horizontal.
static GEN apol_thirdtangent_line1(GEN circ1, GEN circ2, GEN c3, GEN c4, int right){
  pari_sp top=avma;
  if(lg(circ2)==3){//circ2 is also a line. We will place the third circle on the y-axis. This may cause issues if you input something funny, but if we start at [0, 0, 1, 1] then it should be OK I think.
    GEN r=gdivgs(gsub(gel(circ1, 2), gel(circ2, 2)), 2);//The radius
	if(gsigne(r)==-1){//circ2 above circ1
	  r=gneg(r);
	  return gerepilecopy(top, mkvec3(mkcomplex(gen_0, gadd(gel(circ1, 2), r)), r, gdivsg(1, r)));
	}
	return gerepilecopy(top, mkvec3(mkcomplex(gen_0, gadd(gel(circ2, 2), r)), r, gdivsg(1, r)));
  }
  //Now circ2 is a circle.
  if(gequal0(c3)){//We are making the other tangent line.
    GEN y=imag_i(gel(circ2, 1));//height
	if(gcmp(gel(circ1, 2), y)>0) return gerepilecopy(top, mkvec2(gen_0, gsub(y, gel(circ2, 2))));//tangent line circ1 above
	return gerepilecopy(top, mkvec2(gen_0, gadd(y, gel(circ2, 2))));//tangent line circ1 below
  }
  //Now circ3 is a circle too.
  GEN r2=gel(circ2, 2);
  GEN r3=gdivsg(1, c3), x;
  if(typ(gel(circ2, 3))==t_INT && typ(c3)==t_INT && typ(c4)==t_INT){//Must be the strip packing, so sqrt(r2r3) is rational
    GEN c2c3=mulii(gel(circ2, 3), c3);
	x=Qdivii(gen_2, sqrti(c2c3));//x=2sqrt(r2r3)=2/sqrt(c2c3)
  }
  else x=gmulsg(2, gsqrt(gmul(r2, r3), 3));//2sqrt(r2r3) is the x-distance we move. Uses lowest precision.
  if(gcmp(gel(circ1, 2), imag_i(gel(circ2, 1)))>0){//Tangent line circ1 above
	if(right==1){//-x
	  return gerepilecopy(top, mkvec3(mkcomplex(gsub(real_i(gel(circ2, 1)), x), gsub(gel(circ1, 2), r3)), r3, c3));
	}
	else{//+x
	  return gerepilecopy(top, mkvec3(mkcomplex(gadd(real_i(gel(circ2, 1)), x), gsub(gel(circ1, 2), r3)), r3, c3));
	}
  }
  //Tangent line circ1 below
  if(right==1){//+x
	return gerepilecopy(top, mkvec3(mkcomplex(gadd(real_i(gel(circ2, 1)), x), gadd(gel(circ1, 2), r3)), r3, c3));
  }
  //-x
  return gerepilecopy(top, mkvec3(mkcomplex(gsub(real_i(gel(circ2, 1)), x), gadd(gel(circ1, 2), r3)), r3, c3));
}

//vdat=[v, [indices in reps of previous circles]], unless ind=0, when vdat=v
static GEN apol_circles_nextquad(GEN vdat, int ind, void *nul){
  if(ind==0) return mkvec2(vdat, mkvecsmall4(1, 2, 3, 4));//The initial vdat.
  GEN v=gel(vdat, 1);//The actual quadruple
  GEN newv=apol_move(v, stoi(ind), 3);//The new v
  return mkvec2(newv, zv_copy(gel(vdat, 2)));//We do NOT update the index in ind, which will update to lg(reps). This will be done in the getdata method, since we still need the old one for now.
}

//vdat=[v, [indices in reps of previous circles]]. This returns the first element of vdat
static GEN apol_circles_retquad(GEN vdat){return gel(vdat, 1);}

//Given a bounded integral Descartes quadruple, this returns equations for the circles with curvatures <=maxcurv at depth<=depth. An equation takes the form [curvature, radius, x, y]. The outer circle has centre (0, 0).
GEN apol_circles(GEN v, GEN maxcurv){
  return apol_search_bound(v, maxcurv, 1, NULL, &apol_circles_getdata, &apol_circles_nextquad, &apol_circles_retquad, 0);
}

//apol_circles, but goes by depth instead.
GEN apol_circles_depth(GEN v, int depth, GEN maxcurv){
  return apol_search_depth(v, depth, maxcurv, NULL, &apol_circles_getdata, &apol_circles_nextquad, &apol_circles_retquad);
}

//vdat=v=Descartes quadruple. This returns the next one
static GEN apol_generic_nextquad(GEN vdat, int ind, void *nul){return apol_move(vdat, stoi(ind), 3);}

//vdat=v=Descartes quadruple. This returns it
static GEN apol_generic_retquad(GEN vdat){return vdat;}

//Helper method for apol_curvatures, to feed into apol_search_bound. vdat=v=Descartes quadruple
static GEN apol_curvatures_getdata(GEN vdat, int ind, GEN reps, void *nul, int state){
  if(state==0) return ZV_sort(reps);//Sort it and return.
  return gel(vdat, ind);//state=1/2, we want the new curvature
}

//Returns the curvatures in the packing v up to bound.
GEN apol_curvatures(GEN v, GEN bound, int countsymm){
  return apol_search_bound(v, bound, countsymm, NULL, &apol_curvatures_getdata, &apol_generic_nextquad, &apol_generic_retquad, 0);
}

//Returns the curvatures up to depth depth from v. If bound>0, we only count those at most bound.
GEN apol_curvatures_depth(GEN v, int depth, GEN bound){
  if(depth<=0) depth=1;//To avoid errors from bad input.
  return apol_search_depth(v, depth, bound, NULL, &apol_curvatures_getdata, &apol_generic_nextquad, &apol_generic_retquad);
}

//Finds all curvatures in layer at most maxlayers with curvature at most bound. The layer of a circle is how many tangencies it is away from the outer circle.
GEN apol_curvatures_layer(GEN v, int maxlayers, GEN bound, int countsymm){
  return apol_search_bound(v, bound, countsymm, &maxlayers, &apol_layer_getdata, &apol_layer_nextquad, &apol_circles_retquad, 0);
}

//Helper method for apol_find, to feed into apol_search_bound.
static GEN apol_find_getdata(GEN vdat, int ind, GEN reps, void *N, int state){
  if(state==0) return gcopy(reps);//Nothing to do.
  if(equalii(gel(vdat, ind), *(GEN *)N)) return vdat;//We have found N!
  return NULL;//This was not N, do not return anything.
}

//Searches for all circles with curvature N, and returns the corresponding quadruples. If countsymm=1, we may have repeats coming from the symmetry.
GEN apol_find(GEN v, GEN N, int countsymm){
  return apol_search_bound(v, N, countsymm, &N, &apol_find_getdata, &apol_generic_nextquad, &apol_generic_retquad, 0);
}

//Sort the final list or return the new curvature.
static GEN apol_layer_getdata(GEN vdat, int ind, GEN reps, void *nul, int state){
  if(state==0) return ZV_sort(reps);//Sort it and return.
  return gmael(vdat, 1, ind);//state=1/2, we want the new curvature
}

//An element of vdat is [v, current layer], unless ind=0, when vdat=v. The current layer is [L(v[1]), L(v[2]), L(v[3]), L(v[4]), min layer, freq of min layer], where L(v[i]) is the layer of the ith circle in v.
static GEN apol_layer_nextquad(GEN vdat, int ind, void *maxlayers){
  pari_sp top=avma;
  if(ind==0){
    if(gequal0(gel(vdat, 1))) return mkvec2(vdat, mkvecsmalln(6, 0L, 0L, 1L, 1L, 0L, 2L));//Strip packing; layers 0, 0, 1, 1; min 0 1 time
    return mkvec2(vdat, mkvecsmalln(6, 0L, 1L, 1L, 1L, 0L, 1L));//First circles are in layers 0, 1, 1, 1; min 0 1 time
  }
  GEN curlayer=gel(vdat, 2), nextlayer=vecsmall_copy(curlayer);//Current layer, next layer
  if(curlayer[ind]==curlayer[5]){//Currently min. If max, it MUST be repeated, and does not change anything.
    if(curlayer[6]==1){nextlayer[ind]=nextlayer[ind]+2;nextlayer[5]++;nextlayer[6]=3;}//go from x, x+1, x+1, x+1 to x+2, x+1, x+1, x+1.
    else{nextlayer[ind]++;nextlayer[6]--;}//Repeated minimum just moves up 1.
  }
  if(nextlayer[ind]>*(int *)maxlayers) return gc_NULL(top);//Too many layers deep.
  GEN newv=apol_move(gel(vdat, 1), stoi(ind), 3);//Make the move
  return mkvec2(newv, nextlayer);//The new element.
}

//apol_curvatures, but only returns primes
GEN apol_primes(GEN v, GEN bound, int countsymm){
  return apol_search_bound(v, bound, countsymm, NULL, &apol_primes_getdata, &apol_generic_nextquad, &apol_generic_retquad, 0);
}

//Checks if the new curvature is prime.
static GEN apol_primes_getdata(GEN vdat, int ind, GEN reps, void *nul, int state){
  if(state==0) return ZV_sort(reps);//Sort it and return.
  GEN p=gel(vdat, ind);//State 1/2, we want the new curvature
  if(isprime(p)) return p;
  return NULL;
}

//Does apol_primes, but searches up to a maximum number of layers in.
GEN apol_primes_layer(GEN v, int maxlayers, GEN bound, int countsymm){
  return apol_search_bound(v, bound, countsymm, &maxlayers, &apol_primes_layer_getdata, &apol_layer_nextquad, &apol_circles_retquad, 0);
}

//Checks if the new curvature is prime.
static GEN apol_primes_layer_getdata(GEN vdat, int ind, GEN reps, void *nul, int state){
  if(state==0) return ZV_sort(reps);//Sort it and return.
  GEN p=gmael(vdat, 1, ind);//State 1/2, we want the new curvature
  if(isprime(p)) return p;
  return NULL;
}



/*Runs C code to find the missing curvatures up to the given bound, then returns them in a vector. If family=1, removes the known families.*/
GEN
apol_missing(GEN v, GEN B, int family, int load)
{
  pari_sp av = avma;
  GEN w = ZV_sort(apol_red(v, 0, 0));
  GEN modres = apol_mod24(w);
  char *torun = pari_sprintf("./missing_curvatures %d %Pd %Pd %Pd %Pd %Pd", family, B, gel(w, 1), gel(w, 2), gel(w, 3), gel(w, 4));
  long lmod = lg(modres), i;
  for (i = 1; i < lmod; i++) torun = pari_sprintf("%s %Pd", torun, gel(modres, i));/*Making the command to run.*/
  int s = system(torun);
  if (s == -1) {
	pari_err(e_MISC, "problem finding the missing curvatures.");
	return gc_const(av, gen_m1);
  }
  if (!load) return gc_const(av, gen_0);/*Do not load.*/
  set_avma(av);
  return apol_missing_load(v, B, family);
}

/*Loads the saved file of curvatures.*/
GEN
apol_missing_load(GEN v, GEN B, int family)
{
  pari_sp av = avma;
  v = ZV_sort(apol_red(v, 0, 0));
  char *fname;
  if (signe(gel(v, 1)) < 0) fname = pari_sprintf("m%Pd", negi(gel(v, 1)));
  else fname = pari_sprintf("%Pd", gel(v, 1));
  long i;
  for (i = 2; i <= 4; i++) fname = pari_sprintf("%s_%Pd", fname, gel(v, i));
  if (family) fname = pari_sprintf("%s_%Pd_fam.dat", fname, B);
  else fname = pari_sprintf("%s_%Pd.dat", fname, B);
  if (pari_is_dir("missing")) {
    fname = pari_sprintf("missing/%s", fname);
  }
  set_avma(av);
  return gp_readvec_file(fname);/*Load them up!*/
}



//STRIP PACKING METHODS

//Returns [centre, radius, curvature]/[slope, intercept], where the depth pairing corresponding to L is given by corresponding circle/line. If L is an integer, this corresponds to Id_L. If L is a vecsmall/vector, this corresponds to (S_L[1]*...*S_L[n], L[1]).
GEN apol_depthelt_circle(GEN L){
  pari_sp top=avma;
  int t=typ(L);
  switch(t){
    case t_INT:;//Required to stop error with int sL
      int sL=itos(L);
      switch(sL){
        case 1:
          return gerepilecopy(top, mkvec2(gen_0, gen_0));//y=0
        case 2:
          return gerepilecopy(top, mkvec2(gen_0, gen_1));//y=1
        case 3:
          return gerepilecopy(top, mkvec3(mkcomplex(gen_0, ghalf), ghalf, gen_2));//Not sure I need to copy, but will do so to be safe.
        default://i.e. case 4
          return gerepilecopy(top, mkvec3(mkcomplex(gen_m1, ghalf), ghalf, gen_2));
      }
    case t_VEC:
      L=gtovecsmall(L);//Make it a vecsmall
    case t_VECSMALL:
      break;
    default:
      pari_err_TYPE("L needs to be an integer from 1 to 4, or a vecsmall/vector of such integers", L);
  }
  GEN M=apol_matrices();
  GEN W=matid(4);
  for(long i=1;i<lg(L);i++) W=ZM_mul(W, gel(M, L[i]));
  W=ZM_mul(W, gel(M, 5));//Times k at the end
  GEN mtuvw=row(W, L[1]);//[-t, u, v, w]
  GEN twow=shifti(gel(mtuvw, 4), 1);//2w=curvature
  GEN a=Qdivii(gel(mtuvw, 3), gel(mtuvw, 4));
  GEN b=Qdivii(negi(gel(mtuvw, 1)), twow);
  GEN r=Qdivii(gen_1, twow);
  return gerepilecopy(top, mkvec3(mkcomplex(a, b), r, twow));
}

//Returns the set of apol_farey_qf(p, q) for -q/2<p<0 and gcd(p^3-p, q)=1
GEN apol_farey_allqf(GEN q){
  pari_sp top=avma;
  GEN maxp=shifti(subis(q, 1), -1);
  GEN v=vectrunc_init(itos(q)+1);
  for(GEN p=gen_1;cmpii(p, maxp)<=0;p=addis(p, 1)){
    if(!equali1(gcdii(subis(sqri(p), 1), q))) continue;
    GEN qf=apol_farey_qf(p, q);
    if(!gequal0(qf)) vectrunc_append(v, qf);
  }
  return gerepilecopy(top, v);
}

//Returns the quadratic form corresponding to the PSL(2, Z) orbit of the upside down Farey circle at (p, q), where p, q are coprime and q>0. This might not be primitive (it is not if gcd(p^2-1, q)>1 or q is even).
GEN apol_farey_qf(GEN p, GEN q){
  pari_sp top=avma;
  if(!equali1(gcdii(p, q))) return gen_0;
  GEN qsqr=sqri(q);
  return gerepilecopy(top, mkvec3(qsqr, shifti(mulii(p, q), 1), subis(addii(sqri(p), qsqr), 1)));//[q^2, 2pq, p^2+q^2-1]
}

//Returns the data associated to the stair corresponding to the depth element L. If format=1, returns [t, a_W], where t is a positive integer and a_W\in{2, 6, 12} is as in my paper. Returns 0 if the depth element does not intersect the fundamental domain, or if t=2. If format=0, we return [cutoff, height], which has the formula [t-sqrt{t^2-1}, a_W/sqrt(t^2-1)], and works for all depth elements (if we don't intersect the fundamental domain, returns [0, 0]).
GEN apol_stair(GEN L, int format, long prec){
  pari_sp top=avma;
  int ty=typ(L);
  switch(ty){
    case t_INT:;//Required to stop error with int sL
      int sL=itos(L);
      switch(sL){
        case 2:
          return format? gen_0:gerepilecopy(top, mkvec2(gen_1, gdivsg(3, mppi(prec))));
        default://i.e. case 1/3/4
		  return format? gen_0:mkvec2(gen_0, gen_0);
      }
    case t_VEC:
      L=gtovecsmall(L);//Make it a vecsmall
    case t_VECSMALL:
      break;
    default:
      pari_err_TYPE("L needs to be an integer from 1 to 4, or a vecsmall/vector of such integers", L);
  }
  //Now we have a non-identity element.
  long lenL=lg(L)-1;
  int curval=1, isbot=1;
  if(L[lenL]!=1) return format? gc_const(top, gen_0):gerepileupto(top, mkvec2(gen_0, gen_0));//length>=1 means we end with S_1
  lenL--;
  if(lenL>0 && L[lenL]!=4) return format? gc_const(top, gen_0):gerepileupto(top, mkvec2(gen_0, gen_0));//length>=2 means we end with S_4S_1
  lenL--;
  for(long i=lenL;i>0;i--){//Go backwards, we must have the backwards word start with S_1S_4S_1... S_1 or S_4 then S_3 then anything
    if(L[i]==curval){curval=5-curval;continue;}//We are marching along the bottom.
	if(L[i]==3){isbot=0;break;}//We have entered the area, and are not along the bottom. No need to go further here.
	return format? gc_const(top, gen_0):gerepileupto(top, mkvec2(gen_0, gen_0));//We must have L[i]=2, and have left the fundamental domain.
  }
  //At this point we are guaranteed an intersection with the fundamental domain.
  if(lg(L)==2){//S_1
	if(format) return gerepilecopy(top, mkvec2(stoi(7), gen_2));//We do this case explicitly.
	GEN rt12=gsqrt(stoi(12), prec);
	return gerepilecopy(top, mkvec2(gsubsg(7, gmulgs(rt12, 2)), gdivsg(1, rt12)));//S_1
  }
  //Now we must find t. Code copied from apol_depthelt_circle
  GEN M=apol_matrices();
  GEN W=matid(4);
  for(long i=1;i<lg(L);i++) W=ZM_mul(W, gel(M, L[i]));
  W=ZM_mul(W, gel(M, 5));//Times k at the end
  GEN t=negi(gcoeff(W, L[1], 1));//t
  long aw;
  if(isbot==1) aw=6;//(S_4S_1)^k or S_1(S_4S_1)^k
  else aw=12;
  if(format) return gerepilecopy(top, mkvec2(t, stoi(aw)));
  GEN rt=gsqrt(gsubgs(gsqr(t), 1), prec);
  return gerepilecopy(top, mkvec2(gsub(t, rt), gdivsg(aw, rt)));
}

//vdat=[v, matrix, onbot], info=[apol_matrices(), -tmax];
static GEN apol_stairs_nextquad(GEN vdat, int ind, void *info){
  pari_sp top=avma;
  if(ind==0) return mkvec3(vdat, gmael((GEN)info, 1, 5), gen_1);//vdat=[0, 0, 1, 1]. We must initialize vdat=[v, mtx, onbot].
  GEN onbot=gel(vdat, 3);
  if(equali1(onbot)){//Are we still on the bottom? If we aren't, we have nothing to worry about.
    if(ind==2) return gc_NULL(top);//We must swap between S_1 and S_4 until we eventually get to S_3.
	else if(gequal0(gmael(vdat, 1, 1))){//[0, 0, 1, 1]. Must force S_1
	  if(ind!=1) return gc_NULL(top);//Didn't start with S_1.
	}
	else if(equalis(gmael(vdat, 1, 4), 1)){//[4, 0, 1, 1]. Must force S_4
	  if(ind!=4) return gc_NULL(top);//Didn't start with S_4S_1
	}//At this point we are OK, we have an S_4S_1 at the end, and are either continuing or are doing S_3.
	else if(ind==3) onbot=gen_0;//Now we are no longer on the bottom.
  }
  GEN mnew=ZM_mul(gmael((GEN)info, 1, ind), gel(vdat, 2));//Multiply on the left by the correct matrix.
  GEN vnew=apol_move(gel(vdat, 1), stoi(ind), 3);//Make the move
  if(cmpii(gcoeff(mnew, ind, 1), gel((GEN)info, 2))<0) return gc_NULL(top);//t is too large! (we store -tmax and get -t, so we test if -t<-tmax)
  return mkvec3(vnew, mnew, onbot);//Return the new set of data.
}

//Returns the new pair [t, aW]
static GEN apol_stairs_getdata(GEN vdat, int ind, GEN reps, void *nul, int state){
  if(state==0) return lexsort(reps);//Sort it and return.
  if(state==2) return NULL;//Nothing to add at the start.
  GEN t=negi(gcoeff(gel(vdat, 2), ind, 1));//Retrieve t.
  if(gequal0(gel(vdat, 3))) return mkvec2(t, stoi(12));//in the middle, so get [t, 12]
  if(equalis(t, 7)) return mkvec2(t, stoi(2));//S_1, the starting one.
  return mkvec2(t, stoi(6));//On the boundary and not S_1.
}

//Returns the stairs in the packing up to tmax. We use the format of apol_stair with format=1, i.e. [t, a_W], hence we skip the identity element. We also don't combine stairs of the same height.
GEN apol_stairs(GEN tmax){
  pari_sp top=avma;
  GEN M=apol_matrices();
  GEN info=mkvec2(M, negi(tmax));//The info we want to pass to the methods.
  GEN v=mkvec4s(0, 0, 1, 1);//Starting quadruple.
  long rm;
  GEN maxcurv=divis_rem(mulis(tmax, 13), 24, &rm);//13/24*tmax is the maximal curvature, since t=2*curvature*y-coordinate (the 2 since we are using [0, 0, 1, 1], not [0, 0, 2, 2] which is what is actually going on), and y-coord>=12/13
  return gerepileupto(top, apol_search_bound(v, maxcurv, 0, info, &apol_stairs_getdata, &apol_stairs_nextquad, &apol_circles_retquad, 1));
}

/*Returns the quadratic form corresponding to the circle in the stip packing designated by L. If L is an integer, this corresponds to Id_L. If L is a vecsmall/vector, this corresponds to S_L[1]*...*S_L[n]. We can't have L=1 or 2, this doesn't give a circle.*/
GEN
apol_strip_qf(GEN L, int red)
{
  pari_sp av = avma;
  GEN c = apol_depthelt_circle(L);/*The corresponding circle=[centre, r, v].*/
  if (lg(c) == 3) pari_err_TYPE("Please give a circle instead, i.e. don't input L=1 or 2", L);
  GEN C = shifti(gel(c, 3), -1);/*v/2, the curvature of the corresponding circle in [0,0,1,1]*/
  GEN B = gmul(gel(c, 3), real_i(gel(c, 1)));/*B=curvature*Real(centre)*/
  GEN A = shifti(gsub(gmul(gnorm(gel(c, 1)), gel(c, 3)), gel(c, 2)), -1);/*A=cocurvature/2=(|centre|^2*curvature-radius)/2*/
  GEN q = Qfb0(A, B, C);
  if(red) return gerepileupto(av, qfbred(q));
  return gerepilecopy(av, q);
}



//VISUALIZATION


//Given a list of circles, this prints them to the screen in a format suitable for Desmos.
void printcircles_desmos(GEN c){
  for(long i=1;i<lg(c);i++){
    if(lg(gel(c, i))==3) pari_printf("y=%Ps\n", gmael(c, i, 2));//Horizontal line
    else pari_printf("(x-%Ps)^2+(y-%Ps)^2=1/(%Ps)^2\n", real_i(gmael(c, i, 1)), imag_i(gmael(c, i, 1)), gmael(c, i, 3));
  }
}

//Given a list of circles, this prints them to the tex file images/build/imagename_build.tex using tikz. If compile=1, we compile and move the output up to images/imagename.pdf. If open=1, we also open the file, assuming we are working with WSL. Assume the largest circle occurs first. This can take horizontal lines as well. Can supply the scalingfactor, or have the program auto-set it if NULL.
GEN printcircles_tex(GEN c, char *imagename, int addnumbers, int modcolours, int compile, int open, long prec){
  pari_sp top=avma;
  if(!pari_is_dir("images/build")){
    int s=system("mkdir -p images/build");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY images/build");
  }
  char *autofilestr=stack_sprintf("images/build/%s_build.tex", imagename);
  FILE *f=fopen(autofilestr, "w");//Now we have created the output file f.
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{anyfontsize, pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\tikzset{external/force remake}\n  \\pgfplotsset{compat=1.16}\n");
  if(modcolours>0){//Define some colours!
    GEN col=tex_makecolours(modcolours);
    if(typ(gel(col, 1))==t_STR){
      for(long i=1;i<=modcolours;i++){
        pari_fprintf(f, "\\definecolor{col%d}{HTML}{%Ps}\n", i-1, gel(col, i));//Print the custom colours.
      }
    }
    else{//Too many colours, use RGB.
      for(long i=1;i<=modcolours;i++){
        pari_fprintf(f, "\\definecolor{col%d}{rgb}{%P.4f,%P.4f,%P.4f}\n", i-1, gmael(col, i, 1), gmael(col, i, 2), gmael(col, i, 3));//Print the custom colours.
      }
    }
  }
  pari_fprintf(f, "\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n", imagename);
  
  //Now we treat the circles:
  long lc;
  GEN cscale=cgetg_copy(c, &lc);//Scale it so the first circle has radius 3in and centre at (0, 0). The first circle is supposed to be the biggest, having negative curvature, and centre 0, 0. If it is not the largest, we have some work to do.
  GEN largestcirc=stoi(3);//Radius of largest circle
  if(lg(gel(c, 1))==4 && gcmpgs(gmael(c, 1, 3), 0)<0){//First circle neg curvature, assume all circles
    GEN scalingfactor=gdiv(largestcirc, gmael(c, 1, 2));//Scaling factor.
    gel(cscale, 1)=mkvec3(largestcirc, gen_0, gen_0);
    for(long i=2;i<lc;i++){
      gel(cscale, i)=mkvec3(gmul(gmael(c, i, 2), scalingfactor), gmul(real_i(gmael(c, i, 1)), scalingfactor), gmul(imag_i(gmael(c, i, 1)), scalingfactor));//r, x, y
    }//Circles have been scaled!
	pari_printf("Scale:%P.20f\nHorizontal shift: 0\n", scalingfactor);
  }
  else{//Strip packing, OR the largest curvature does not come first. The width should be at least as much as the height.
    GEN minx=mkoo(), maxx=mkmoo();
    for(long i=1;i<lc;i++){
      if(lg(gel(c, i))==3) continue;//Horizontal line.
      GEN circxmin=gsub(real_i(gmael(c, i, 1)), gmael(c, i, 2));
      GEN circxmax=gadd(real_i(gmael(c, i, 1)), gmael(c, i, 2));
      minx=gmin_shallow(minx, circxmin);
      maxx=gmax_shallow(maxx, circxmax);
    }
    if(typ(minx)==t_INFINITY || typ(maxx)==t_INFINITY) pari_err_TYPE("You don't have any circles", c);
    GEN scalingfactor=gmulgs(gdiv(largestcirc, gsub(maxx, minx)), 2);
    GEN horizshift=gdivgs(gadd(maxx, minx), 2);
    for(long i=1;i<lc;i++){
      if(lg(gel(c, i))==3){//Horizontal line
        gel(cscale, i)=mkvec3(mkoo(), largestcirc, gmul(gmael(c, i, 2), scalingfactor));//oo, "radius", y-intercept ("radius"=r means going from [-r, y] to [r, y])
        continue;
      }
      gel(cscale, i)=mkvec3(gmul(gmael(c, i, 2), scalingfactor), gmul(gsub(real_i(gmael(c, i, 1)), horizshift), scalingfactor), gmul(imag_i(gmael(c, i, 1)), scalingfactor));//r, x, y
    }
	pari_printf("Scale:%P.20f\nHorizontal shift: %P.20f\n", scalingfactor, gneg(horizshift));
  }
  
  //Time to draw the circles
  if(modcolours==0){
    char *drawoptions="[ultra thin]";
    for(long i=1;i<lc;i++){
      if(typ(gmael(cscale, i, 1))==t_INFINITY) pari_fprintf(f, "  \\draw%s (%P.10fin, %P.10fin) -- (%P.10fin, %P.10fin);\n", drawoptions, gneg(gmael(cscale, i, 2)), gmael(cscale, i, 3), gmael(cscale, i, 2), gmael(cscale, i, 3));
      else pari_fprintf(f, "  \\draw%s (%P.10fin, %P.10fin) circle (%P.10fin);\n", drawoptions, gmael(cscale, i, 2), gmael(cscale, i, 3), gmael(cscale, i, 1));
    }
  }
  else{//Must add colouring. Not sure how this works for the strip packing.
    for(long i=1;i<lc;i++){
      if(typ(gmael(cscale, i, 1))==t_INFINITY) pari_fprintf(f, "  \\draw[ultra thin, fill=col%d] (%P.10fin, %P.10fin) -- (%P.10fin, %P.10fin);\n", smodis(gmael(c, i, 3), modcolours), gneg(gmael(cscale, i, 2)), gmael(cscale, i, 3), gmael(cscale, i, 2), gmael(cscale, i, 3));
      else{
        if(signe(gmael(c, i, 3))==-1) pari_fprintf(f, "  \\draw[ultra thin] (%P.10fin, %P.10fin) circle (%P.10fin);\n", gmael(cscale, i, 2), gmael(cscale, i, 3), gmael(cscale, i, 1));
        else pari_fprintf(f, "  \\draw[ultra thin, fill=col%d] (%P.10fin, %P.10fin) circle (%P.10fin);\n", smodis(gmael(c, i, 3), modcolours), gmael(cscale, i, 2), gmael(cscale, i, 3), gmael(cscale, i, 1));
      }
    }
  }
  GEN ten=stoi(10);
  if(addnumbers){
      char *options;
      if(modcolours>0) options=", white";
      else options="";
    for(long i=1;i<lc;i++){
      if(lg(gel(c, i))==3) continue;//No numbers for horizontal lines
      GEN curv=gmael(c, i, 3), scaleby;
      if(signe(curv)!=1) continue;//Don't add numbers for curvatures <0.
      long ndigits=logint(curv, ten)+1;
      switch(ndigits){//Scaling the font.
        case 1:
          scaleby=dbltor(1.4);break;
        case 2:
          scaleby=gen_1;break;
        case 3:
          scaleby=dbltor(0.8);break;
        case 4:
          scaleby=dbltor(0.6);break;
        case 5:
          scaleby=dbltor(0.5);break;
        default:
          scaleby=gdivgs(gen_2, ndigits);
      }
      pari_fprintf(f, "  \\node[align=center%s] at (%P.10fin, %P.10fin) {\\fontsize{%P.10fin}{0in}\\selectfont %Pd};\n", options, gmael(cscale, i, 2), gmael(cscale, i, 3), gmul(gmael(cscale, i, 1), scaleby), curv);
    }
  }
  
  //Ending stuff.
  pari_fprintf(f, "\\end{tikzpicture}\n\\end{document}");
  fclose(f);
  tex_compile(imagename, open);
  return gerepilecopy(top, mkvec2(strtoGENstr(imagename), stoi(open)));
}



//SUPPORTING METHODS


//Returns all 4*3^(d-1) reduced words of depth d.
GEN apol_words(int d){
  if(d<=0) return cgetg(1, t_VEC);
  pari_sp top=avma;
  int ind=1;//We reuse ind to track which depth we are going towards.
  GEN I=vecsmall_ei(d, 1);//Tracks the words
  int forward=1;//Tracks if we are going forward or not.
  long Nreps=4*(itos(powuu(3, d-1)))+1;//4*3^(d-1)+1
  GEN reps=vectrunc_init(Nreps);
  do{//1<=ind<=d is assumed.
    I[ind]=forward? 1:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
    if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
    //At this point, we can go on with valid and new inputs
    if(ind==d){
      forward=0;
      vectrunc_append(reps, gcopy(I));//Add in I
      continue;
    }
    //We can keep going forward
    ind++;
    forward=1;
  }while(ind>0);
  return gerepilecopy(top, reps);
}

//Copies an integer vector
static GEN
ZV_copy(GEN v){
  long len=lg(v);
  GEN rvec=cgetg(len, t_VEC);
  for(long i=1;i<len;i++) gel(rvec, i)=icopy(gel(v, i));
  return rvec;
}



//SPECIALIZED METHODS
//These are scripts that are useful to me, but not likely useful in general.


//Computes apol_makeall, and returns the external depths of the corresponding circles of curvature n. This works because we use the reduced quadratic forms: A-n<=C-n<=A+C-B-n. Swapping A+C-B-n gives A+C+B-n, so this cannot decrease it. Thus either we are already reduced (and A-n<=0 necessarily), or it is not the largest term (which thus must be n). Thus n is immediately swapped if depth>0, and this gives the circle depth.
GEN apol_makeall_extdepths(GEN n, long prec){
  pari_sp top=avma;
  GEN forms=apol_makeall(n, 0, prec);//The forms
  long maxdepth=0, lf=lg(forms);
  GEN depths=cgetg(lf, t_VECSMALL);
  for(long i=1;i<lf;i++){//Computing depths
    long d=apol_extdepth(gel(forms, i), prec);
    depths[i]=d;
    if(d>maxdepth) maxdepth=d;
  }
  GEN dcount=zero_zv(maxdepth+1);//Depths 0 through d.
  for(long i=1;i<lf;i++) dcount[depths[i]+1]++;//Counting.
  return gerepileupto(top, dcount);
}

//Computes apol_makeall, and returns the sorted vector of smallest curvature for each example. We do not remove repeats, and output the negative of the curvatures (as they are all negative). If red=0, we actually just take the data as is, and output the smallest curvature in each of the quadruples (no negation).
GEN apol_makeall_small(GEN n, int red, long prec){
  pari_sp top=avma;
  GEN forms=apol_makeall(n, red, prec);
  long lf;
  GEN curves=cgetg_copy(forms, &lf);
  for(long i=1;i<lf;i++){
    gel(curves, i)=gmael(forms, i, 1);
    if(red) togglesign_safe(&gel(curves, i));
  }
  return gerepileupto(top, ZV_sort(curves));
}

//apol_makeall_small, but we only reduce maxsteps steps.
GEN apol_makeall_small_maxsteps(GEN n, long maxsteps, long prec){
  pari_sp top=avma;
  GEN forms=apol_makeall(n, 0, prec);
  long lf;
  GEN curves=cgetg_copy(forms, &lf);
  for(long i=1;i<lf;i++){
    GEN q=apol_red_partial0(gel(forms, i), maxsteps, &apol_move_1i, &cmpii);
    gel(curves, i)=gel(q, vecindexmin(q));
  }
  return gerepileupto(top, ZV_sort(curves));
}


/*Returns 1 if x is t_INT, 2 if t_REAL, 0 else*/
static int
ismp(GEN x)
{
  long t = typ(x);
  if (t == t_INT) return 1;
  if (t == t_REAL) return 2;
  return 0;
}
