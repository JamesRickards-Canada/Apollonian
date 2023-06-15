/*Methods to deal with quadratic forms.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"

/*STATIC DECLARATIONS*/

/*SECTION 2: BASIC QUADRATIC FORM METHODS*/
static GEN nf_get_rootD(GEN nf, GEN D);

/*SECTION 3: CLASS GROUP*/
static int qfbequal1(GEN q);
static GEN qfbidentity(GEN D);

/*MAIN BODY*/


/*SECTION 1: DISCRIMINANT METHODS*/

/*Generates a list of discriminants from D1 to D2, can specify if they are fundamental and coprime to a given input.*/
GEN
disclist(GEN D1, GEN D2, int fund, GEN cop)
{
  pari_sp av = avma;
  if (typ(D1) != t_INT) pari_err_TYPE("D1 must be an integer", D1);
  if (typ(D2) != t_INT) pari_err_TYPE("D2 must be an integer", D2);
  if (typ(cop) != t_INT) pari_err_TYPE("cop must be an integer", cop);
  long maxlen = (itos(subii(D2, D1)) >> 1) + 5;/*Add 5 to cover edge cases, just to be safe.*/
  GEN Dlist = vectrunc_init(maxlen);
  GEN D = D1;
  if (!fund) {/*Don't have to be fundamental.*/
    if (gequal0(cop)) {/*No coprimality restrictions.*/
      for (; cmpii(D, D2) <= 0; D = addis(D, 1)) {
        if (isdisc(D)) vectrunc_append(Dlist, D);
      }
    }
    else {
      for (; cmpii(D, D2) <= 0; D = addis(D, 1)) {
        if (equali1(gcdii(cop, D)) && isdisc(D)) vectrunc_append(Dlist, D);
      }
    }
    return gerepilecopy(av, Dlist);
  }
  /*Must be fundamental.*/
  if (gequal0(cop)) {/*Don't have to be fundamental.*/
    for (; cmpii(D, D2) <= 0; D = addis(D, 1)){/*coredisc(0)=0, coredisc(1)=1, but we don't want to count them.*/
      if (equalii(coredisc(D), D) && !gequal0(D) && !equali1(D)) vectrunc_append(Dlist, D);
    }
    return gerepilecopy(av, Dlist);
  }
  for (; cmpii(D, D2) <= 0; D = addis(D, 1)) {
    if (equali1(gcdii(cop, D))) {
      if (equalii(coredisc(D), D) && !gequal0(D) && !equali1(D)) vectrunc_append(Dlist, D);
    }
  }
  return gerepilecopy(av, Dlist);
}

/*Returns 1 if discriminant and 0 if not.*/
int
isdisc(GEN D)
{
  if (typ(D) != t_INT) return 0;
  if (smodis(D, 4) < 2 && !Z_issquare(D)) return 1;
  return 0;
}


/*SECTION 2: BASIC QUADRATIC FORM METHODS*/

/*Returns L^n acting on q.*/
GEN
qfbapplyL(GEN q, GEN n)
{
  pari_sp av = avma;
  GEN D = gel(q, 4);
  GEN An = mulii(gel(q, 1), n);/*An*/
  GEN AnpB = addii(An, gel(q, 2));/*An+B*/
  GEN Bnew = addii(An, AnpB);/*B+2nA*/
  GEN AnsqBn = mulii(AnpB, n);/*An^2+Bn*/
  return gerepilecopy(av, mkqfb(gel(q, 1), Bnew, addii(AnsqBn, gel(q, 3)), D));
}

/*Returns R^n acting on q.*/
GEN
qfbapplyR(GEN q, GEN n)
{
  pari_sp av = avma;
  GEN D = gel(q, 4);
  GEN Cn = mulii(gel(q, 3), n);/*Cn*/
  GEN CnpB = addii(Cn, gel(q, 2));/*Cn+B*/
  GEN Bnew = addii(Cn, CnpB);/*B+2nC*/
  GEN CnsqBn = mulii(CnpB, n);/*Cn^2+Bn*/
  return gerepilecopy(av, mkqfb(addii(CnsqBn, gel(q, 1)), Bnew, gel(q, 3), D));
}

/*Returns S acting on q*/
GEN
qfbapplyS(GEN q)
{
  GEN Sq = cgetg(4, t_QFB);
  gel(Sq, 1) = icopy(gel(q, 3));
  gel(Sq, 2) = negi(gel(q, 2));
  gel(Sq, 3) = icopy(gel(q, 1));
  gel(Sq, 4) = icopy(gel(q, 4));
  return Sq;
}

/*Converts the ideal x in the (necessarily quadratic) number field nf into an integral quadratic form. We are assuming that we are working in the maximal ideal for now. We follow Buell.*/
GEN
idealtoqfb(GEN nf, GEN x)
{
  pari_sp av = avma;
  if (nf_get_degree(nf) != 2) pari_err_TYPE("must be a quadratic number field", nf);
  GEN D = nf_get_disc(nf);
  GEN hnfx = idealhnf(nf, x);
  GEN a1 = gel(hnfx, 1), a2 = gel(hnfx, 2);/*Generators of the ideal, with a1 in Q.*/
  GEN a2tr = nftrace(nf, a2);/*integer or fraction*/
  GEN a2conj = nfsub(nf, a2tr, a2);
  GEN diff = nfsub(nf, a2conj, a2);/*conj(a)-a2*/
  GEN normrt = nfmul(nf, a1, diff);/*a1*conj(a2)-conj(a1)*a2 since a1=conj(a1)*/
  GEN y = nf_get_rootD(nf, D);/*sqrt(D), chosen consistently.*/
  GEN norm = lift(basistoalg(nf, nfdiv(nf, normrt, y)));/*(a1*conj(a2)-conj(a1)*a2)/sqrt(D)=N(a), lives in Q*/
  int swap = 0;
  if (gsigne(norm) < 0) { swap = 1; norm = gneg(norm); }/*We order so that norm>0, so must swap if not.*/
  GEN invnorm = ginv(norm);
  GEN A = gmul(nfnorm(nf, a1), invnorm);
  GEN B = gmul(gmul(lift(basistoalg(nf, a1)), a2tr), invnorm);
  GEN C = gmul(nfnorm(nf, a2), invnorm);
  if (swap) return gerepilecopy(av, mkqfb(C, B, A, D));
  return gerepilecopy(av, mkqfb(A, B, C, D));
}

/*Converts the qfb (positive definite if disc<0) to a fractional ideal in the number field nf. Under SL(2, Z) equivalence, this is inverse to idealtoqfb.*/
GEN
qfbtoideal(GEN nf, GEN q)
{
  pari_sp av = avma;
  GEN D = gel(q, 4), A = gel(q, 1);
  GEN y = nf_get_rootD(nf, D), b, delta;
  if (mod2(D)) {/*D odd*/
    b = shifti(subis(gel(q, 2), 1), -1);/*(B-1)/2*/
    delta = gdivgs(nfsub(nf, gen_1, y), 2);/*1-sqrt(D)/2*/
  }
  else {/*D even*/
    b = shifti(gel(q, 2), -1);/*B/2*/
    delta = gneg(y);/*-sqrt(D)*/
  }
  GEN a1 = mkcol2(A, gen_0);/*The element A in the number field*/
  GEN a2 = nfadd(nf, b, delta);/*b+delta*/
  if (signe(A) < 0) {
    a1 = nfmul(nf, a1, delta);/*A*delta*/
    a2 = nfmul(nf, a2, delta);/*(b+delta)*delta*/
  }
  GEN x = mkmat2(algtobasis(nf, a1), algtobasis(nf, a2));/*a1 and a2 are the generators of the ideal*/
  return gerepilecopy(av, idealhnf(nf, x));
}

/*In a quadratic number field Q(sqrt(D)), returns the algebraic expression for the element which squares to D and has positive coefficient in the number field variable. Not stack clean.*/
static GEN
nf_get_rootD(GEN nf, GEN D)
{
  pari_sp av = avma;
  GEN y = mkcol2(gen_0, gen_1);
  y = nfsub(nf, y, gdivgs(nftrace(nf, y), 2));/*Trace 0, so squares to D*square*/
  GEN ysqr = lift(basistoalg(nf, nfsqr(nf, y)));/*integer or fraction*/
  GEN shift = gdiv(D, ysqr);/*Times the square root of this.*/
  if (typ(shift) == t_INT) shift = sqrti(shift);/*Guaranteed to be a square.*/
  else { gel(shift, 1) = sqrti(gel(shift, 1)); gel(shift, 2) = sqrti(gel(shift, 2)); }/*Must be a square fraction.*/
  y = lift(basistoalg(nf, nfmul(nf, y, shift)));/*Now, y^2=D*/
  GEN firstcoef = polcoef_i(y, 1, nf_get_varn(nf));/*The coefficient of the variable defining nf.*/
  if (gsigne(firstcoef) < 0) y = gneg(y);/*We want y to be the square root of D for which the first coefficient is positive. This is consistent.*/
  return gerepilecopy(av, y);
}


/*SECTION 3: CLASS GROUP*/

/*Given a vector/Vecsmall of positive integers v and an index ind, this returns the Vecsmall [e1, e2, ..., er] which is the index ind output of a lexicographic ordering of [0, v[1]-1] x [0, v[2]-1] x ... x [0, v[r]-1].*/
GEN
lexind(GEN v, long ind)
{
  pari_sp av = avma;
  if (typ(v) == t_VEC) v = ZV_to_zv(v);
  long lord;
  GEN es = cgetg_copy(v, &lord);
  long a = ind - 1, b, n, i;/*we start with ind-1*/
  for (i = lord - 1; i > 0; i--) {/*Go backwards.*/
    n = v[i];
    a = sdivss_rem(a, n, &b);/*a_{i+1}=a_i*n_i+b_i; a_{r+1}=i-1.*/
    es[i] = b;
  }
  return gerepilecopy(av, es);
}

/*This computes the narrow class group*/
GEN
qfbnarrow(GEN D, long prec)
{
  pari_sp av = avma;
  if (!isdisc(D)) pari_err_TYPE("D must be a discriminant, i.e. a non-square integer equivalent to 0 or 1 modulo 4.", D);
  GEN cgp = quadclassunit0(D, 0, NULL, prec);/*Class group*/
  if (signe(D) < 0 || quadunitnorm(D) == -1) {/*Negative disc or norm(fund. unit)=-1*/
    if (equali1(gel(cgp, 1))) {
      /*We modify the second and third entries, as we want [1, Vecsmall([1]), [identity], reg] instead of [1, [], [], reg].*/
      gel(cgp, 2) = mkvecsmall(1);
      gel(cgp, 3) = mkvec(qfbred(qfbidentity(D)));
    }
    else gel(cgp, 2) = ZV_to_zv(gel(cgp, 2));/*Convert to Vecsmall.*/
    return gerepilecopy(av, vec_shorten(cgp, 3));
  }
  GEN minus1;/*-id, the new form to add*/
  if (mod2(D)) minus1 = qfbred(mkqfb(gen_m1, gen_1, shifti(subis(D, 1), -2), D));/*[-1, 1, (D-1)/4]*/
  else minus1 = qfbred(mkqfb(gen_m1, gen_0, shifti(D, -2), D));/*[-1, 0, D/4]*/
  if (equali1(gel(cgp, 1))) return gerepilecopy(av, mkvec3(gen_2, mkvecsmall(2), mkvec(minus1)));/*h^+(D)=2, done separately.*/
  GEN newsize = shifti(gel(cgp, 1), 1);
  GEN orders = ZV_to_zv(gel(cgp, 2));/*Convert to Vecsmall.*/
  GEN gens = gel(cgp, 3);
  long lgens = lg(gens), firstm1 = 0, i;/*Start by checking if any of our current generators power to -1.*/
  for (i = 1; i < lgens; i++) {
    GEN qpow = qfbpows(gel(gens, i), orders[i]);
    if (!qfbequal1(qpow)) { firstm1 = i; break; }/*note: qpow is either -1 or 1*/
  }
  if (!firstm1) {/*Nothing powers to -1, just add in new Z/2Z copy*/
    long lastodd = 0;/*Index of largest odd order*/
    for (i = 1; i < lgens; i++) { if (orders[i] % 2) { lastodd = i; break; } }
    if (!lastodd) {/*All even, add a new Z/2Z at the end.*/
      orders = vecsmall_append(orders, 2);
      gens = vec_append(gens, minus1);
    }
    else {/*Now there is an odd order element. Composing this with minus1 gives the new generator, and is the only change required.*/
      orders[lastodd] <<= 1; /*Double order*/
      gel(gens, lastodd) = qfbcomp_i(gel(gens, lastodd), minus1);/*q*minus1*/
    }
    return gerepilecopy(av, mkvec3(newsize, orders, gens));
  }
  /*Now we have an element powering to -1 and not 1. We must fix the previous and subsequent entries.*/
  GEN g = gel(gens, firstm1);
  long ai = orders[firstm1];/*Initially this is half the order of g, since g^ai=-1*/
  GEN gpow = g;/*gpow stores g^(.5order(g)/ai)*/
  for (i = firstm1 + 1; i < lgens; i++) {/*Fix the later ones by composing with gpow^(.5order(g)/ai) if they power to -1 instead of 1*/
    GEN qpow = qfbpows(gel(gens, i), orders[i]);
    if (qfbequal1(qpow)) continue;/*OK already*/
    gpow = qfbpows(gpow, ai/orders[i]);
    ai = orders[i];/*Updating ai*/
    gel(gens, i) = qfbcomp_i(gel(gens, i), gpow);/*Modifying the generator.*/
  }
  /*We are done fixing the generators past firstm1; they all have the correct order. We shift focus to those before it.*/
  /*In order to keep the divisibility of orders, we must figure out which generator should double in order.*/
  long optplace = 1;/*The place where we should be inserting the 2x in the class group*/
  long v = vals(orders[firstm1]);/*2-adic valuation*/  
  for (i = firstm1 - 1; i > 0; i--) {
    if (vals(orders[i]) > v) { optplace = i + 1; break; }
  }/*If this for method never triggers, we are left with optplace=1=the start, as desired.*/
  /*Now we must modify the elements between firstm1 and optplace to only double the order of optplace and retain the other orders. In fact, we only need to modify these two elements*/
  if (firstm1 != optplace) {/*If equal, there is nothing left to do! If !=, we must shift the elements appropriately*/
    if (!v) {/*The elements between firstm1 and optplace are all ODD powered, so we can adjust firstm1 and optplace by minus1*/
      gel(gens, firstm1) = qfbcomp_i(gel(gens, firstm1), minus1);
      gel(gens, optplace) = qfbcomp_i(gel(gens, optplace), minus1);
    }
    else {/*firstm1 is g1 and has order 2^((v+1)h1), optplace is g2 and has order 2^(vh2) with h1, h2 odd and h1|h2.*/
      long h2 = orders[optplace] >> v;
      /*Replace g1 by g1^2*g2^h2, which has order 2^vh1, and replace g2 by g1*g2, which has order 2^(v+1)h2. These elements generate the same subgroup as well since 2-h2 is coprime to 2h2.*/
      gel(gens, firstm1) = qfbcomp_i(qfbpow(g, gen_2), qfbpows(gel(gens, optplace), h2));/*g1^2*g2^h2*/
      gel(gens, optplace) = qfbcomp_i(g, gel(gens, optplace));/*g1*g2*/
    }
  }
  orders[optplace] <<= 1;/*Doubling the correct place*/
  return gerepilecopy(av, mkvec3(newsize, orders, gens));
}

/*Returns 1 if q is similar to the form [1, D%2, (D%2-D)/4] where D=disc(q), and 0 else.*/
static int
qfbequal1(GEN q)
{
  pari_sp av = avma;
  if (signe(gel(q, 4)) < 0) {/*D<0*/
    GEN qred = qfbred(q);
    return gc_int(av, equali1(gel(qred, 1)));
  }
  GEN sols = qfbsolve(q, gen_1, 0);
  if (lg(sols) > 1) return gc_int(av, 1);
  return gc_int(av, 0);
}

/*Returns the identity element of discriminant D. Not stack clean.*/
static GEN
qfbidentity(GEN D)
{
  if (mod2(D)) return mkqfb(gen_1, gen_1, shifti(subsi(1, D), -2), D);/*[1, 1, (1-D)/4]*/
  return mkqfb(gen_1, gen_0, shifti(negi(D), -2), D);/*[1, 0, -D/4]*/
}

/*qfbnarrow, but the third entry is this list of forms in lexicographic ordering. Can input qfbnarrow(D) instead of D.*/
GEN
qfbnarrowlex(GEN D, long prec)
{
  pari_sp av = avma;
  GEN nar;/*Narrow class group.*/
  if (typ(D) == t_VEC) {/*Fix D/nar*/
    nar = D;
    D = gmael3(nar, 3, 1, 4);
  }
  else nar = qfbnarrow(D, prec);/*Get the narrow class group.*/
  long narclno = itos(gel(nar, 1));/*Narrow class number*/
  GEN ords = gel(nar, 2);
  GEN gens = gel(nar, 3);
  GEN allforms = cgetg(narclno + 1, t_VEC);
  gel(allforms, 1) = qfbred(qfbidentity(D));/*Start with the identity.*/
  long f1, f2 = 1, pow, i, j, k = 2;
  for (i = lg(gel(nar, 2)) - 1; i >= 1; i--) {
    f1 = 1;
    f2 = k - 1;
    for (pow = 1; pow <= ords[i] - 1; pow++) {
      k = f2 + 1;
      for (j = f1; j <= f2; j++) {
        gel(allforms, k) = qfbcomp_i(gel(allforms, j), gel(gens, i));
        k++;
      }
      f1 = f2+1;
      f2 = k - 1;/*update f1, f2*/
    }
  }
  return gerepilecopy(av, mkvec3(gel(nar, 1), ords, allforms));
}

