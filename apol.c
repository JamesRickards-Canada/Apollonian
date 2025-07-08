/*Basic methods to deal with Apollonian circle packings. This part of the package is not designed to be as efficient as possible: it should be good, but we use pari GENS which are slower than C longs. Efficient methods should be written in C and placed in apol_fast. We allow integral and real packings, but NOT fractions (and cannot mix real and integral).*/

/*INCLUSIONS*/
#include <pari.h>
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

/*SECTION 3: COMPUTING THE CIRCLES*/
static GEN apol_thirdtangent(GEN circ1, GEN circ2, GEN c3, GEN c4, int right, long prec);
static GEN apol_thirdtangent_line1(GEN circ1, GEN circ2, GEN c3, GEN c4, int right, long prec);
static int toleq_circ(GEN c1, GEN c2, GEN tol);

/*SECTION 4: SUPPORTING METHODS*/
static int ismp(GEN x);
static long quarticresidue_int(GEN x, GEN y);
static GEN gaussian_makeodd(GEN x, long *pwr);
static GEN gaussian_makeprimary(GEN x, long *pwr);
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

/*Returns chi_2 of the ACP. We start by finding a pair of coprime tangent circles.*/
long
apol_chi2(GEN v)
{
  pari_sp av = avma;
  if (!apol_check_primitive(v)) pari_err_TYPE("must be a primitive integral packing", v);
  GEN vcop = shallowcopy(v);
  GEN a = gel(vcop, 1), b;
  if (isintzero(a)) return gc_long(av, 1);/*Strip packing has chi_2 = 1.*/
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

/*Returns chi_4 of the ACP of type (6, 1) or (6, 17).*/
GEN
apol_chi4(GEN v)
{
  pari_sp av = avma;
  GEN a = gel(v, 1);
  if (isintzero(a)) return gen_1;/*Strip packing has chi_4 = 1.*/
  GEN q = qfbred(apol_qf(v, 1));
  GEN sos = qfbsos1(q), beta = NULL;
  if (equali1(gcdii(a, gel(q, 1)))) beta = mkcomplex(gcoeff(sos, 1, 1), gcoeff(sos, 2, 1));
  else if (equali1(gcdii(a, gel(q, 3)))) beta = mkcomplex(gcoeff(sos, 1, 2), gcoeff(sos, 2, 2));
  else {
    long x, y;
    for (x = 1; ; x++) {
      for (y = 1; y <= x; y++) {
        GEN term2 = mulis(addii(mulis(gel(q, 2), x), mulis(gel(q, 3), y)), y);/*(Bx+Cy)*y*/
        GEN term1 = mulis(gel(q, 1), x * x);/*Ax^2*/
        GEN term = addii(term1, term2);
        if (equali1(gcdii(a, term))) {
          GEN concoef = addii(mulis(gcoeff(sos, 1, 1), x), mulis(gcoeff(sos, 1, 2), y));
          GEN Icoef = addii(mulis(gcoeff(sos, 2, 1), x), mulis(gcoeff(sos, 2, 2), y));
          beta = mkcomplex(concoef, Icoef);
          break;
        }
      }
      if (beta) break;
    }
  }
  beta = simplify(beta);/*In case we did 4+0*I.*/
  GEN nprime;
  long e = Z_lvalrem(a, 2, &nprime), pw;/*nprime = odd part*/
  if (e) beta = gaussian_makeprimary(beta, &pw);/*Needs to be primary.*/
  GEN qr = quarticresidue(beta, nprime);
  if (!e) return gerepilecopy(av, qr);/*n is odd.*/
  if (e == 2) {
    if (Mod4(nprime) == 1) return gerepilecopy(av, qr);/*n==4(8), n'==1(4).*/
    return gerepileupto(av, gneg(qr));/*n==4(8), n'==3(4).*/
  }
  if (Mod8(mulis(imag_i(beta), e))) return gerepileupto(av, gneg(qr));/*n==0(8)*/
  return gerepilecopy(av, qr);
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
  if (!apol_check(v, prec)) return NULL;
  int t = ismp(gel(v, 1));
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
  long m3, m8, i = 0;
  do { i++; m3 = smodis(gel(v, i), 3); } while (!m3);/*Find non-zero residue modulo 3.*/
  i = 0;
  do { i++; m8 = Mod8(gel(v, i)); } while (!(m8 % 2));/*Find odd residue modulo 8.*/
  GEN m24;
  switch (m8) {
    case 1:
      if (m3 == 1) m24 = mkvecsmalln(6, 0L, 1L, 4L, 9L, 12L, 16L);
      else m24 = mkvecsmalln(6, 0L, 8L, 9L, 12L, 17L, 20L);
      break;
    case 5:
      if (m3 == 1) m24 = mkvecsmalln(6, 0L, 4L, 12L, 13L, 16L, 21L);
      else m24 = mkvecsmalln(6, 0L, 5L, 8L, 12L, 20L, 21L);
      break;
    default:/*(8, k) type*/
      if (m3 == 1) m24 = mkvecsmalln(8, 3L, 6L, 7L, 10L, 15L, 18L, 19L, 22L);
      else m24 = mkvecsmalln(8, 2L, 3L, 6L, 11L, 14L, 15L, 18L, 23L);
  }
  return gerepileupto(av, zv_to_ZV(m24));
}

/*Calls apol_move_(1/batch)(i/r), returning a clean result without affecting v. We also first check that v is a Descartes quadruple.*/
GEN
apol_move(GEN v, GEN command, long prec)
{
  pari_sp av = avma;
  int isint;
  GEN vorig = v;
  v = apol_fix(v, &isint, prec);
  if (!v) pari_err_TYPE("v is not a Descartes quadruple.", vorig);
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
  GEN w = apol_fix(v, &isint, prec);
  if (!w) pari_err_TYPE("v is not a Descartes quadruple.", v);
  if (isint) return gerepileupto(av, apol_red0(w, seq, &apol_move_1i, &cmpii));
  return gerepileupto(av, apol_red0(w, seq, &apol_move_1r, &cmprr));
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

/*Returns the "type" of v, i.e. the number of resiudes modulo 24 and the smallest residue coprime to 6, which uniquely identifies it. If chi = 1, also includes the chi2 value, and if chi=2, includes the chi4 value (0 if not (6, 1) or (6, 17) packings).*/
GEN
apol_type(GEN v, int chi)
{
  pari_sp av = avma;
  GEN m24 = apol_mod24(v);
  long second = itos(gel(m24, 2));/*Uniquely identified by the second element*/
  set_avma(av);
  if (chi == 1) {
    switch (second) {
      case 1: return mkvec3s(6, 1, apol_chi2(v));
      case 3: return mkvec3s(8, 11, apol_chi2(v));
      case 4: return mkvec3s(6, 13, apol_chi2(v));
      case 5: return mkvec3s(6, 5, apol_chi2(v));
      case 6: return mkvec3s(8, 7, apol_chi2(v));
      case 8: return mkvec3s(6, 17, apol_chi2(v));
    }
  }
  if (chi == 2) {
    switch (second) {
      case 1: return gerepilecopy(av, mkvec4(stoi(6), gen_1, stoi(apol_chi2(v)), apol_chi4(v)));
      case 3: return mkvec4s(8, 11, apol_chi2(v), 0);
      case 4: return mkvec4s(6, 13, apol_chi2(v), 0);
      case 5: return mkvec4s(6, 5, apol_chi2(v), 0);
      case 6: return mkvec4s(8, 7, apol_chi2(v), 0);
      case 8: return gerepilecopy(av, mkvec4(stoi(6), stoi(17), stoi(apol_chi2(v)), apol_chi4(v)));
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
  if (isintzero(n)) return gerepilecopy(av, mkvec(mkvec4s(0, 0, 1, 1)));
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


/*SECTION 3: COMPUTING THE CIRCLES*/

/*Finds equations for all circles in the ACP generated by v that are at most the bound B. v can be integral or not. Can also pass in a maximal depth, maximal x value to save at the end, and a minimal curvature.*/
GEN
apol_circles(GEN v, GEN bounds, long depth, GEN maxxval, long prec)
{
  pari_sp av = avma;
  GEN B, Bmin = NULL;
  long t = typ(bounds);
  if (t == t_VEC) {
    Bmin = gel(bounds, 1);
    t = typ(Bmin);
    if (t != t_INT && t != t_FRAC && t != t_REAL) pari_err_TYPE("bounds must be integers, rational numbers, or real numbers.", bounds);
    B = gel(bounds, 2);
    t = typ(B);
  }
  else B = bounds;
  if (t != t_INT && t != t_FRAC && t != t_REAL) pari_err_TYPE("bounds must be integers, rational numbers, or real numbers.", bounds);
  int isint;
  GEN w = apol_fix(v, &isint, prec);
  if (!w) pari_err_TYPE("v is not a Descartes quadruple.", v);
  void (*move1)(GEN, int);
  if (isint) move1 = &apol_move_1i;
  else move1 = &apol_move_1r;
  v = sort(w);/*Sort it, and move back to the variable v.*/
  if (depth < 0) depth = 0;
  if (gequal0(gel(w, 1)) && !depth) pari_err_TYPE("For a strip, half-plane, and full-plane packing, the depth MUST be set.", stoi(depth));
  long maxdepth;
  if (depth) maxdepth = depth + 1;
  else maxdepth = 100;/*Maximal depth, to start.*/
  GEN swaps = const_vecsmall(maxdepth, 0);/*Tracks the sequence of Apollonian moves. We do NOT use swaps[1].*/
  GEN depthseq = cgetg(maxdepth + 1, t_VEC);/*Each entry stores [v, inds]: the current quadruple v, and a Vecsmall of the indices in vfound of the equations for the circles.*/
  gel(depthseq, 1) = mkvec2(v, mkvecsmall4(1, 2, 3, 4));/*Base case is v, with the first 4 entries of vfound as circles.*/
  GEN tol = deftol(prec);
  long maxfound = 3000, foundind = 4;
  GEN vfound = cgetg(maxfound + 1, t_VEC);/*Storing the circle equations found.*/
  
  /*Now, we insert the first four circles.*/
  GEN c1 = gel(v, 1);/*First, or outer circle. It may also be a line.*/
  if (gequal0(c1)) gel(vfound, 1) = mkvec2(gen_0, gen_0);/*Line!*/
  else if (signe(c1) < 0) gel(vfound, 1) = mkvec3(gen_0, gdivsg(-1, c1), c1);/*Outer circle. First circle has negative curvature necessarily, hence why r1=-1/c*/
  else gel(vfound, 1) = mkvec3(gen_0, ginv(c1), c1);/*First circle has positive curvature, may be a full plane packing!*/
  GEN c2 = gel(v, 2);/*Second circle*/
  if (gequal0(c2)) {/*It is a line.*/
    if (gequal0(c1)) gel(vfound, 2) = mkvec2(gen_0, gdivsg(2, gel(v, 3)));/*First circle was also a line*/
    else gel(vfound, 2) = mkvec2(gen_0, gsub(imag_i(gmael(vfound, 1, 1)), gmael(vfound, 1, 2)));/*First circle was a circle. Place line below.*/
  }
  else {/*It is a circle*/
    GEN r2 = ginv(c2);
    if (gequal0(c1)) gel(vfound, 2) = mkvec3(mkcomplex(gen_0, r2), r2, c2);/*First circle was a line! It must be the x-axis*/
    else if (signe(c1) < 0){/*First circle was outer one, negative curvature.*/
      gel(vfound, 2) = mkvec3(mkcomplex(gen_0, gsub(gmael(vfound, 1, 2), r2)), r2, c2);/*first inner circle, placed vertically at the top. Note gmael(vfound, 1, 2) = r1 > 0.*/
    }
    else gel(vfound, 2) = mkvec3(mkcomplex(gen_0, gadd(gmael(vfound, 1, 2), r2)), r2, c2);/*First circle was not the outer one.*/
  }
  gel(vfound, 3) = apol_thirdtangent(gel(vfound, 1), gel(vfound, 2), gel(v, 3), gel(v, 4), 0, prec);/*Third circle goes left*/
  gel(vfound, 4) = apol_thirdtangent(gel(vfound, 2), gel(vfound, 3), gel(v, 4), gel(v, 1), 0, prec);/*Fourth circle is left of circ2 ->circ3*/
  
  /*Depth first search section.*/
  long ind = 2;/*Which depth we are working at. This differs to apol_fast methods as we start at index 0 there and index 1 here.*/
  while (ind > 1) {/*We are coming in trying to swap this circle out.*/
    int cind = ++swaps[ind];/*Increment the swapping index.*/
    if (cind == 5 || (depth && ind == maxdepth)) {/*Overflowed, go back.*/
      swaps[ind] = 0;
      ind--;
      continue;
    }
    long lastind = ind - 1;
    if (cind == swaps[lastind]) continue; /*Same thing twice, so skip it.*/
    GEN newv = shallowcopy(gmael(depthseq, lastind, 1));/*Shallow copy the previous entry.*/
    move1(newv, cind);/*Move it.*/
    GEN newc = gel(newv, cind);
    if (gcmp(newc, B) > 0) continue;/*Too big! go back.*/
    
    GEN newinds = cgetg(5, t_VECSMALL);/*The new indices*/
    GEN prevind = gmael(depthseq, lastind, 2);/*The corresponding circle indices of the previous quadruple.*/
    long i = 1;
    while (i == cind) i++;
    newinds[i] = prevind[i];
    GEN oldcirc1 = gel(vfound, prevind[i]);/*First old circle*/
    i++;
    while (i == cind) i++;
    newinds[i] = prevind[i];
    GEN oldcirc2 = gel(vfound, prevind[i]);/*Second old circle*/
    i++;
    while (i == cind) i++;
    newinds[i] = prevind[i];
    GEN newcirc = apol_thirdtangent(oldcirc1, oldcirc2, newc, gel(newv, i), 1, prec);/*The new circle, if it is to the right of oldcirc1 -> oldcirc2.*/
    GEN prevcirc = gel(vfound, prevind[cind]);/*The circle we are replacing.*/
    /*Now we need to check that we found the new circle on the correct side of the old ones.*/
    if (toleq0(gel(newv, 1), tol) && toleq0(gel(newv, 2), tol)) {/*Strip packing moving along the top.*/
      GEN oldcirc3 = gel(vfound, prevind[i]);/*The unused old circle. Our newcirc must be tangent to it.*/
      GEN ct = gel(oldcirc3, 1);/*Centre of oldcirc3*/
      GEN ctdiff = real_i(gsub(ct, gel(prevcirc, 1)));
      newcirc = mkvec3(gadd(ct, ctdiff), gcopy(gel(oldcirc3, 2)), gcopy(gel(oldcirc3, 3)));
    }
    else if (toleq_circ(newcirc, prevcirc, tol)) newcirc = apol_thirdtangent(oldcirc1, oldcirc2, newc, gel(newv, i), 0, prec);/*If the two curvatures were the same, this could trigger.*/
    else {
      GEN oldcirc3 = gel(vfound, prevind[i]);/*The unused old circle. Our newcirc must be tangent to it.*/
      GEN rsums = gsqr(gadd(gel(oldcirc3, 2), gel(newcirc, 2)));/*(r1+r2)^2*/
      GEN dcentres = gnorm(gsub(gel(oldcirc3, 1), gel(newcirc, 1)));/*dist(centres)^2*/
      if (!toleq(rsums, dcentres, tol)) newcirc = apol_thirdtangent(oldcirc1, oldcirc2, newc, gel(newv, i), 0, prec);/*Must be the other side.*/
    }
    /*Update the data.*/
    foundind++;
    if (foundind > maxfound) {
        maxfound <<= 1;/*Double it*/
        vfound = vec_lengthen(vfound, maxfound);
    }
    gel(vfound, foundind) = newcirc;/*Add the new circle to the list*/
    newinds[cind] = foundind;
    gel(depthseq, ind) = mkvec2(newv, newinds);/*Update the depth sequence.*/
    ind++;
    if (ind > maxdepth) {/*We are going too deep, must double the depth. This cannot trigger if depth=0.*/
      long newdepth = maxdepth << 1;/*Double it.*/
      depthseq = vec_lengthen(depthseq, newdepth);
      swaps = vecsmall_lengthen(swaps, newdepth);
      for (i = maxdepth + 1; i <= newdepth; i++) swaps[i] = 0;
      maxdepth = newdepth;
    }
  }
  if (!maxxval) {
    if (!Bmin) return gerepilecopy(av, vec_shorten(vfound, foundind));
    GEN remaining = vectrunc_init(foundind + 1);
    long i;
    for (i = 1; i <= foundind; i++) {
      GEN cir = gel(vfound, i);
      if (lg(cir) == 4 && gcmp(gel(cir, 3), Bmin) < 0) continue;
      vectrunc_append(remaining, cir);
    }
    return gerepilecopy(av, remaining);
  }
  if (!Bmin) Bmin = mkmoo();
  GEN remaining = vectrunc_init(foundind + 1);
  long i;
  for (i = 1; i <= foundind; i++) {
    GEN cir = gel(vfound, i);
    if (lg(cir) == 4 && gcmp(gel(cir, 3), Bmin) < 0) continue;/*Curvature too small.*/
    GEN xright = gadd(real_i(gel(cir, 1)), gel(cir, 2));
    GEN xleft = gsub(real_i(gel(cir, 1)), gel(cir, 2));
    GEN mx = gmax(gabs(xright, prec), gabs(xleft, prec));
    if (gcmp(mx, maxxval) > 0) continue;/*Went too far.*/
    vectrunc_append(remaining, cir);
  }
  return gerepilecopy(av, remaining);
}

/*Store a circle as [centre, radius, curvature]. Given two tangent circles and a third curvature, this finds this third circle that is tangent to the first two. For internal tangency, we need a positive radius and negative curvature. There are always 2 places to put the circle: left or right of the line from the centre of circ1 to the centre of circ2. If right=1, we put it right, else we put it left. c4 is one of the curvatures to complete an Apollonian quadruple (supplying it allows us to always work with exact numbers in the case of integral ACPs). In the case of a horizontal line, we ASSUME that it has slope 0.*/
static GEN
apol_thirdtangent(GEN circ1, GEN circ2, GEN c3, GEN c4, int right, long prec)
{
  if (lg(circ1) == 3) return apol_thirdtangent_line1(circ1, circ2, c3, c4, right, prec);/*circ1 is a line!*/
  if (lg(circ2) == 3) return apol_thirdtangent_line1(circ2, circ1, c3, c4, 1 - right, prec);/*circ2 is a line!*/
  pari_sp av = avma;/*Now neither circ1 or circ2 is a line.*/
  if (gequal0(c3)) {/*c3 is a line.*/
    if (right) return gerepilecopy(av, mkvec2(gen_0, gsub(imag_i(gel(circ1, 1)), gel(circ1, 2))));/*Below*/
    return gerepilecopy(av, mkvec2(gen_0, gadd(imag_i(gel(circ1, 1)), gel(circ1, 2))));/*Above*/
  }
  /*Now circ1, circ2, and c3 are all circles.*/
  /*The centres form a triangle with sides r1+r2, r1+r3, r2+r3, or -r1-r2, -r1-r3, r2+r3 (if internal tangency, where r1<0). Let theta be the angle at the centre of c1.*/
  GEN c1 = gel(circ1, 3), c2 = gel(circ2, 3);/*Curvatures*/
  GEN c1pc2 = gadd(c1, c2), c1pc3 = gadd(c1, c3);
  GEN denom = gmul(c1pc2, c1pc3);
  GEN costheta = gsubsg(1, gdiv(gmulsg(2, gsqr(c1)), denom));/*1-2c1^2/((c1+c2)(c1+c3))*/
  GEN sintheta=gdiv(gabs(gmul(c1, gadd(c1pc2, gsub(c3, c4))), prec), denom);/*|c1(c1+c2+c3-c4)|/((c1+c2)(c1+c3))*/
  GEN r1 = gel(circ1, 2), r2 = gel(circ2, 2), r3 = ginv(c3);/*The radii*/
  if (signe(c1) < 0) r1 = gneg(r1);
  if (signe(c2) < 0) r2 = gneg(r2);/*Correcting r1 and r2 to be negative.*/
  GEN x1 = real_i(gel(circ1, 1)), y1 = imag_i(gel(circ1, 1)), x2 = real_i(gel(circ2, 1)), y2 = imag_i(gel(circ2, 1));
  /*We need alpha, the angle to the centre of c2 from the centre of c1.*/
  GEN r1pr2 = gadd(r1, r2);
  GEN cosalpha = gdiv(gsub(x2, x1), r1pr2);
  GEN sinalpha = gdiv(gsub(y2, y1), r1pr2);/*cos and sin of alpha.*/
  /*If right=1, we have angle alpha-theta, and if right=0, we have angle alpha+theta. We derive the new coordinates from the (co)sine addition/subtraction formulae.*/
  GEN relcos, relsin;
  if (right) {/*cos(alpha-theta), sin(alpha-theta)*/
    relcos = gadd(gmul(cosalpha, costheta), gmul(sinalpha, sintheta));
    relsin = gsub(gmul(sinalpha, costheta), gmul(cosalpha, sintheta));
  }
  else{/*cos(alpha+theta), sin(alpha+theta)*/
    relcos = gsub(gmul(cosalpha, costheta), gmul(sinalpha, sintheta));
    relsin = gadd(gmul(sinalpha, costheta), gmul(cosalpha, sintheta));
  }
  GEN r1pr3 = gadd(r1, r3);
  GEN x = gadd(x1, gmul(r1pr3, relcos));
  GEN y = gadd(y1, gmul(r1pr3, relsin));
  return gerepilecopy(av, mkvec3(mkcomplex(x, y), r3, c3));
}

/*apol_thirdtangent when circ1 is actually a line. We assume that it is horizontal.*/
static GEN
apol_thirdtangent_line1(GEN circ1, GEN circ2, GEN c3, GEN c4, int right, long prec)
{
  pari_sp av = avma;
  if (lg(circ2) == 3) {/*circ2 is also a line. We will place the third circle on the y-axis. This may cause issues if you input something funny, but if we start at [0, 0, 1, 1] then it should be OK I think.*/
    GEN r = gdivgs(gsub(gel(circ1, 2), gel(circ2, 2)), 2);/*The radius*/
    if (signe(r) < 0){/*circ2 above circ1*/
      r = gneg(r);
      return gerepilecopy(av, mkvec3(mkcomplex(gen_0, gadd(gel(circ1, 2), r)), r, ginv(r)));
    }
    return gerepilecopy(av, mkvec3(mkcomplex(gen_0, gadd(gel(circ2, 2), r)), r, ginv(r)));
  }
  /*Now circ2 is a circle.*/
  if (gequal0(c3)) {/*We are making the other tangent line.*/
    GEN y = imag_i(gel(circ2, 1));/*height*/
    if (gcmp(gel(circ1, 2), y) > 0) return gerepilecopy(av, mkvec2(gen_0, gsub(y, gel(circ2, 2))));/*tangent line circ1 above*/
    return gerepilecopy(av, mkvec2(gen_0, gadd(y, gel(circ2, 2))));/*tangent line circ1 below*/
  }
  /*Now circ3 is a circle too.*/
  GEN r2 = gel(circ2, 2);
  GEN r3 = ginv(c3), x;
  if (typ(gel(circ2, 3)) == t_INT && typ(c3) == t_INT && typ(c4) == t_INT){/*Must be the strip packing, so sqrt(r2r3) is rational*/
    GEN c2c3 = mulii(gel(circ2, 3), c3);
    x = Qdivii(gen_2, sqrti(c2c3));/*x=2sqrt(r2r3)=2/sqrt(c2c3)*/
  }
  else x = gmulsg(2, gsqrt(gmul(r2, r3), prec));/*2sqrt(r2r3) is the x-distance we move. Uses lowest precision.*/
  if (gcmp(gel(circ1, 2), imag_i(gel(circ2, 1))) > 0) {/*Tangent line circ1 above*/
    if (right) return gerepilecopy(av, mkvec3(mkcomplex(gsub(real_i(gel(circ2, 1)), x), gsub(gel(circ1, 2), r3)), r3, c3)); /*-x*/
    return gerepilecopy(av, mkvec3(mkcomplex(gadd(real_i(gel(circ2, 1)), x), gsub(gel(circ1, 2), r3)), r3, c3));/*+x*/
  }
  /*Tangent line circ1 below*/
  if (right) return gerepilecopy(av, mkvec3(mkcomplex(gadd(real_i(gel(circ2, 1)), x), gadd(gel(circ1, 2), r3)), r3, c3));/*+x*/
  return gerepilecopy(av, mkvec3(mkcomplex(gsub(real_i(gel(circ2, 1)), x), gadd(gel(circ1, 2), r3)), r3, c3));/*-x*/
}

/*Returns 1 if the circles/lines are equal, 0 else. Assumes lines are always horizontal.*/
static int
toleq_circ(GEN c1, GEN c2, GEN tol)
{
  if (lg(c1) == 3) {
    if (lg(c2) != 3) return 0;/*line vs circle*/
    return toleq(gel(c1, 2), gel(c2, 2), tol);/*Only need to compare the intercepts, both slopes assumed to be 0.*/
  }
  if (lg(c2) == 3) return 0;/*circle vs line*/
  if (!toleq(gel(c1, 1), gel(c2, 1), tol)) return 0;/*Different centres.*/
  return toleq(gel(c1, 3), gel(c2, 3), tol);/*No need to check radii if curvatures are okay.*/
}

/*Given a list of circles, this prints them to the screen in a format suitable for Desmos.*/
void
printcircles_desmos(GEN c)
{
  long i;
  for (i = 1; i < lg(c); i++) {
    if (lg(gel(c, i)) == 3) pari_printf("y=%P.10f\n", gmael(c, i, 2));/*Horizontal line*/
    else pari_printf("(x-%P.10f)^2+(y-%P.10f)^2=1/(%P.10f)^2\n", real_i(gmael(c, i, 1)), imag_i(gmael(c, i, 1)), gmael(c, i, 3));
  }
}

/*Given a list of circles/lines in an ACP, this prints them to the tex file images/build/imagename_build.tex using tikz. If compile=1, we compile and move the output up to images/imagename.pdf. If open=1, we also open the file, assuming we are working with WSL. We can also supply the radius of the bounding circle (if relevant), defaults to 3in.*/
GEN
printcircles_tex(GEN c, char *imagename, int addnumbers, int modcolours, GEN outerrad, int compile, int open, long prec)
{
  pari_sp av = avma;
  if( !pari_is_dir("images/build")) {
    int s = system("mkdir -p images/build");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY images/build");
  }
  char *autofilestr = stack_sprintf("images/build/%s_build.tex", imagename);
  FILE *f = fopen(autofilestr, "w");/*Now we have created the output file f.*/
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{anyfontsize, pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\tikzset{external/force remake}\n  \\pgfplotsset{compat=1.16}\n");
  long i;
  if (modcolours) {/*Define some colours!*/
    GEN col = tex_makecolours(modcolours);
    if (typ(gel(col, 1)) == t_STR) {
      for (i = 1; i <= modcolours; i++) {
        pari_fprintf(f, "\\definecolor{col%ld}{HTML}{%Ps}\n", i - 1, gel(col, i));/*Print the custom colours.*/
      }
    }
    else {/*Too many colours, use RGB.*/
      for (i = 1; i <= modcolours; i++) {
        pari_fprintf(f, "\\definecolor{col%ld}{rgb}{%P.4f,%P.4f,%P.4f}\n", i - 1, gmael(col, i, 1), gmael(col, i, 2), gmael(col, i, 3));/*Print the custom colours.*/
      }
    }
  }
  pari_fprintf(f, "\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n", imagename);
  
  /*Now we treat the circles:*/
  long lc;
  GEN cscale = cgetg_copy(c, &lc);/*Scale it so the first circle has radius 3in and centre at (0, 0). The first circle is supposed to be the biggest, having negative curvature, and centre 0, 0. If it is not the largest, we have some work to do.*/
  GEN largestcirc;
  if (outerrad) largestcirc = outerrad;
  else largestcirc = stoi(3);/*Radius of largest circle*/
  if (lg(gel(c, 1)) == 4 && gcmpgs(gmael(c, 1, 3), 0) < 0) {/*First circle neg curvature, assume all circles*/
    GEN scalingfactor = gdiv(largestcirc, gmael(c, 1, 2));/*Scaling factor.*/
    gel(cscale, 1) = mkvec3(largestcirc, gen_0, gen_0);
    for (i = 2; i < lc; i++) {
      gel(cscale, i) = mkvec3(gmul(gmael(c, i, 2), scalingfactor), gmul(real_i(gmael(c, i, 1)), scalingfactor), gmul(imag_i(gmael(c, i, 1)), scalingfactor));/*r, x, y*/
    }/*Circles have been scaled!*/
    pari_printf("Scale:%P.20f\nHorizontal shift: 0\n", scalingfactor);
  }
  else {/*Strip packing, OR the largest curvature does not come first. The width should be at least as much as the height.*/
    GEN minx = mkoo(), maxx = mkmoo();
    for (i = 1; i < lc; i++) {
      if (lg(gel(c, i)) == 3) continue;/*Horizontal line.*/
      GEN circxmin = gsub(real_i(gmael(c, i, 1)), gmael(c, i, 2));
      GEN circxmax = gadd(real_i(gmael(c, i, 1)), gmael(c, i, 2));
      minx = gmin_shallow(minx, circxmin);
      maxx = gmax_shallow(maxx, circxmax);
    }
    if (typ(minx) == t_INFINITY || typ(maxx) == t_INFINITY) pari_err_TYPE("You don't have any circles", c);
    GEN scalingfactor = gmulgs(gdiv(largestcirc, gsub(maxx, minx)), 2);
    GEN horizshift = gdivgs(gadd(maxx, minx), 2);
    for (i = 1; i < lc; i++) {
      if (lg(gel(c, i)) == 3) {/*Horizontal line*/
        gel(cscale, i) = mkvec3(mkoo(), largestcirc, gmul(gmael(c, i, 2), scalingfactor));/*oo, "radius", y-intercept ("radius"=r means going from [-r, y] to [r, y])*/
        continue;
      }
      gel(cscale, i) = mkvec3(gmul(gmael(c, i, 2), scalingfactor), gmul(gsub(real_i(gmael(c, i, 1)), horizshift), scalingfactor), gmul(imag_i(gmael(c, i, 1)), scalingfactor));/*r, x, y*/
    }
    pari_printf("Scale:%P.20f\nHorizontal shift: %P.20f\n", scalingfactor, gneg(horizshift));
  }
  
  /*Time to draw the circles*/
  if (!modcolours) {
    char *drawoptions = "[ultra thin]";
    for (i = 1; i < lc; i++) {
      if (typ(gmael(cscale, i, 1)) == t_INFINITY) pari_fprintf(f, "  \\draw%s (%P.10fin, %P.10fin) -- (%P.10fin, %P.10fin);\n", drawoptions, gneg(gmael(cscale, i, 2)), gmael(cscale, i, 3), gmael(cscale, i, 2), gmael(cscale, i, 3));
      else pari_fprintf(f, "  \\draw%s (%P.10fin, %P.10fin) circle (%P.10fin);\n", drawoptions, gmael(cscale, i, 2), gmael(cscale, i, 3), gmael(cscale, i, 1));
    }
  }
  else {/*Must add colouring. For strip, won't colour the horizontal line.*/
    for (i = 1; i < lc; i++) {
      if (typ(gmael(cscale, i, 1)) == t_INFINITY) pari_fprintf(f, "  \\draw[ultra thin] (%P.10fin, %P.10fin) -- (%P.10fin, %P.10fin);\n", gneg(gmael(cscale, i, 2)), gmael(cscale, i, 3), gmael(cscale, i, 2), gmael(cscale, i, 3));
      else {
        if(signe(gmael(c, i, 3)) < 0) pari_fprintf(f, "  \\draw[ultra thin] (%P.10fin, %P.10fin) circle (%P.10fin);\n", gmael(cscale, i, 2), gmael(cscale, i, 3), gmael(cscale, i, 1));
        else pari_fprintf(f, "  \\draw[ultra thin, fill=col%ld] (%P.10fin, %P.10fin) circle (%P.10fin);\n", smodis(gmael(c, i, 3), modcolours), gmael(cscale, i, 2), gmael(cscale, i, 3), gmael(cscale, i, 1));
      }
    }
  }
  GEN ten = stoi(10);
  if (addnumbers) {
    char *options;
    if (modcolours) options = ", white";
    else options = "";
    for (i = 1; i < lc; i++) {
      if (lg(gel(c, i)) == 3) continue;/*No numbers for horizontal lines*/
      GEN curv = gmael(c, i, 3), scaleby;
      if (signe(curv) <= 0) continue;/*Don't add numbers for curvatures <0.*/
      long ndigits;
      if (typ(curv) == t_REAL) ndigits = logint(gfloor(curv), ten) + 3;/*Real numbers to two decimal places.*/
      else ndigits = logint(curv, ten) + 1;
      switch(ndigits){/*Scaling the font.*/
        case 1:
          scaleby = dbltor(1.4);break;
        case 2:
          scaleby = gen_1;break;
        case 3:
          scaleby = dbltor(0.8);break;
        case 4:
          scaleby = dbltor(0.6);break;
        case 5:
          scaleby = dbltor(0.5);break;
        default:
          scaleby = gdivgs(gen_2, ndigits);
      }
      if (typ(curv) == t_REAL) {
        pari_fprintf(f, "  \\node[align=center%s] at (%P.10fin, %P.10fin) {\\fontsize{%P.10fin}{0in}\\selectfont %P.2f};\n", options, gmael(cscale, i, 2), gmael(cscale, i, 3), gmul(gmael(cscale, i, 1), scaleby), curv);
      }
      else {
        pari_fprintf(f, "  \\node[align=center%s] at (%P.10fin, %P.10fin) {\\fontsize{%P.10fin}{0in}\\selectfont %Pd};\n", options, gmael(cscale, i, 2), gmael(cscale, i, 3), gmul(gmael(cscale, i, 1), scaleby), curv);
      }
    }
  }
  /*Ending stuff.*/
  pari_fprintf(f, "\\end{tikzpicture}\n\\end{document}");
  fclose(f);
  tex_compile(imagename, open);
  return gerepilecopy(av, mkvec2(strtoGENstr(imagename), stoi(open)));
}


/*SECTION 4: STRIP PACKING METHODS*/


/*SECTION 5: SUPPORTING METHODS*/

/*Returns all 4*3^(d-1) reduced words of depth d.*/
GEN
apol_words(int d)
{
  if (d <= 0) return cgetg(1, t_VEC);
  pari_sp av = avma;
  int ind = 1;/*We reuse ind to track which depth we are going towards.*/
  GEN I = vecsmall_ei(d, 1);/*Tracks the words*/
  int forward = 1;/*Tracks if we are going forward or not.*/
  long Nreps = ((itos(powuu(3, d - 1))) << 2) + 1;/*4*3^(d-1)+1*/
  GEN reps = vectrunc_init(Nreps);
  do {/*1<=ind<=d is assumed.*/
    I[ind] = forward? 1:I[ind] + 1;
    if (ind>1 && I[ind - 1] == I[ind]) I[ind]++;/*Don't repeat*/
    if (I[ind] > 4) { ind--; continue; }/*Go back. Forward already must =0, so no need to update.*/
    /*At this point, we can go on with valid and new inputs*/
    if (ind == d) {
      forward = 0;
      vectrunc_append(reps, gcopy(I));/*Add in I*/
      continue;
    }
    /*We can keep going forward*/
    ind++;
    forward = 1;
  } while(ind > 0);
  return gerepilecopy(av, reps);
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

/*Returns the quartic residue symbol [x/y]=i^e for coprime Gaussian integers x, y with y odd. Does not verify that x, y are coprime.*/
GEN
quarticresidue(GEN x, GEN y)
{
  long tx = typ(x);
  if (tx == t_COMPLEX) {
    if (typ(gel(x, 1)) != t_INT || typ(gel(x, 2)) != t_INT) pari_err_TYPE("x and y must be Gaussian integers", x);
  }
  else if (tx != t_INT) pari_err_TYPE("x and y must be Gaussian integers", x);
  long ty = typ(y);
  if (ty == t_COMPLEX) {
    if (typ(gel(y, 1)) != t_INT || typ(gel(y, 2)) != t_INT) pari_err_TYPE("x and y must be Gaussian integers", y);
  }
  else if (ty != t_INT) pari_err_TYPE("x and y must be Gaussian integers", y);
  long pow = quarticresidue_int(x, y);
  if (!pow) return gen_1;
  if (pow == 1) return mkcomplex(gen_0, gen_1);
  if (pow == 2) return gen_m1;
  return mkcomplex(gen_0, gen_m1);
}

/*Returns 0<=e<=3 where the quartic residue symbol [x/y]=i^e for coprime Gaussian integers x, y with y odd. Does not check these things.*/
static long
quarticresidue_int(GEN x, GEN y)
{
  pari_sp av = avma;
  GEN sixteen = stoi(16);
  long shift;
  y = gaussian_makeprimary(y, &shift);
  long pwr = 0;/*Tracks the power of i.*/
  for (;;) {
    GEN a = real_i(y), b = imag_i(y);
    GEN z = gdiv(x, y);
    GEN w = mkcomplex(ground(real_i(z)), ground(imag_i(z)));
    GEN delta = gsub(z, w);/*x/y=w+delta with w in Z[i].*/
    if (gequal0(delta)) return gc_long(av, pwr);/*Done, ASSUMING coprime.*/
    x = gmul(delta, y);/*Reduce x mod y to something smaller.*/
    x = gaussian_makeodd(x, &shift);
    if (shift) {/*Update pwr, divided by (1+i)^shift*/
        GEN t = Fp_sub(a, Fp_add(b, Fp_add(Fp_sqr(b, sixteen), gen_1, sixteen), sixteen), sixteen);/*a-b-b^2-1*/
        pwr += (Mod16(t) >> 2) * shift;/*Add the power t/4*/
    }
    x = gaussian_makeprimary(x, &shift);
    if (shift) {/*Update pwr, divided by i^shift*/
      pwr += (Mod8(subsi(1, a)) >> 1) * shift;
    }
    if (Mod4(b) && Mod4(imag_i(x))) pwr += 2;/*x and y both are == (3, 2) mod 4*/
    GEN t = x;
    x = y;
    y = t;
    pwr = pwr % 4;
  }
}

/*Assuming x=a+bi is a non-zero t_COMPLEX in Z[i], divides by (1+i)^k to make it odd, and stores k mod 4 in pwr. Leaves garbage.*/
static GEN
gaussian_makeodd(GEN x, long *pwr)
{
  GEN a = real_i(x), b = imag_i(x);
  long pw = 0;
  while (Mod2(a) == Mod2(b)) {
    GEN t = a;
    a = shifti(addii(a, b), -1);
    b = shifti(subii(b, t), -1);
    pw++;
  }
  *pwr = pw % 4;
  return mkcomplex(a, b);
}

/*Assuming x=a+bi is a t_COMPLEX in Z[i], updates it to its primary associate, and stores the negative of the power of i (mod 4) we multiplied by in pwr. Leaves garbage.*/
static GEN
gaussian_makeprimary(GEN x, long *pwr)
{
  GEN a = real_i(x), b = imag_i(x);
  int a4 = Mod4(a), b4 = Mod4(b);
  switch (a4) {
    case 0:
      if (b4 == 1) { *pwr = 1; return mkcomplex(b, negi(a)); }
      *pwr = 3;
      return mkcomplex(negi(b), a);/*b==3(4)*/
    case 1:
      if (!b4) { *pwr = 0; return x; }
      *pwr = 2;
      return mkcomplex(negi(a), negi(b));/*b==2(4)*/
    case 2:
      if (b4 == 1) { *pwr = 3; return mkcomplex(negi(b), a); }
      *pwr = 1;
      return mkcomplex(b, negi(a));/*b==3(4)*/
  }/*a==3(4)*/
  if (b4) { *pwr = 0; return x; }
  *pwr = 2;
  return mkcomplex(negi(a), negi(b));/*b==0(4)*/
}

/*Copies an integer vector*/
static GEN
ZV_copy(GEN v) {
  long len = lg(v), i;
  GEN rvec = cgetg(len, t_VEC);
  for (i = 1; i < len; i++) gel(rvec, i) = icopy(gel(v, i));
  return rvec;
}

