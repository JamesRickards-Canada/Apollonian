/*Methods to deal with quadratic forms.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"

/*STATIC DECLARATIONS*/

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
	  gel(cgp, 3) = mkvec(qfbidentity(D));
	}
	else gel(cgp, 2) = ZV_to_zv(gel(cgp, 2));/*Convert to Vecsmall.*/
	return gerepilecopy(av, vec_shorten(cgp, 3));
  }
  GEN minus1;/*-id, the new form to add*/
  if (mod2(D)) minus1 = mkqfb(gen_m1, gen_1, shifti(subis(D, 1), -2), D);/*[-1, 1, (D-1)/4]*/
  else minus1 = mkqfb(gen_m1, gen_0, shifti(D, -2), D);/*[-1, 0, D/4]*/
  if (equali1(gel(cgp, 1))) return gerepilecopy(av, mkvec4(gen_2, mkvecsmall(2), mkvec(minus1), gel(cgp, 4)));/*h^+(D)=2, done separately.*/
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
	return gerepilecopy(av, mkvec4(newsize, orders, gens, gel(cgp, 4)));
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




//TUPLE REPS OF PRIMES
static int mod_collapse(GEN L);
static int bqf_tuplevalid_cmp(void *data, GEN x, GEN y);






//Converts the ideal ideal in the (necessarily quadratic) number field numf into an integral quadratic form.
GEN ideal_tobqf(GEN numf, GEN ideal){
  pari_sp top = avma;
  ideal = idealhnf0(numf, ideal, NULL);
  GEN alph1 = gadd(gmul(gcoeff(ideal, 1, 1), gel(member_zk(numf), 1)), gmul(gcoeff(ideal, 2, 1), gel(member_zk(numf), 2)));
  GEN alph2 = gadd(gmul(gcoeff(ideal, 1, 2), gel(member_zk(numf), 1)), gmul(gcoeff(ideal, 2, 2), gel(member_zk(numf), 2)));//alph1,alph2 generate the ideal
  long varno=varn(alph1);//The variable number
  GEN alph2conj=gsubst(alph2, varno, gneg(pol_0(varno)));//Conjugating alph2
  GEN a1 = lift(gmodulo(gmul(alph1, alph2conj), member_pol(numf))); //a1=alph1*conj(alph2)=u*sqrt(D)+b
  GEN A, C;
  if(gcmpgs(polcoef_i(a1, 1, varno), 0) > 0){//if >0 we are not properly ordered (we require (conj(a1)-a1)/sqrt(D)>0), and must swap alph1,alph2; */
    A=nfnorm(numf, alph2);
    C=nfnorm(numf, alph1);
  }
  else{
    A=nfnorm(numf, alph1);
    C=nfnorm(numf, alph2);
  }
  GEN B=gmul(polcoef_i(a1, 0, varno), gen_2);//B/2 may be a half integer, so can't use shifti or mulii sadly
  togglesign_safe(&B);//Negate B in place
  GEN d=gcdii(gcdii(A, B), C);
  GEN Q=cgetg(4, t_VEC);
  if(equali1(d)){
    gel(Q, 1)=icopy(A);
    gel(Q, 2)=icopy(B);
    gel(Q, 3)=icopy(C);
  }
  else{
    gel(Q, 1)=diviiexact(A, d);
    gel(Q, 2)=diviiexact(B, d);
    gel(Q, 3)=diviiexact(C, d);
  }
  return gerepileupto(top, Q);
}





//TUPLE REPS OF PRIMES

//Returns the residue classes modulo D that q could represent.
GEN bqf_primesmod(GEN q){
  pari_sp top=avma, mid;
  GEN D=absi(bqf_disc(q));//WLOG work with a positive number.
  mid=avma;
  GEN divs=divisors(D);
  GEN L=vectrunc_init(lg(divs)*itos(D));//Stores q(x, y) as x loops over divisors of D and 0, and y ranges from 0 to D-1
  if(equali1(gcdii(gel(q, 3), D))) vectrunc_append(L, Fp_red(gel(q, 3), D));//q(0, 1), the only value we care about if x=0.
  for(long i=1;i<lg(divs)-1;i++){//Don't want to include D as a divisor, already treated this case.
    GEN x=gel(divs, i);
	GEN c1=Fp_mul(gel(q, 1), Fp_sqr(x, D), D);//q[1]*x^2
	GEN c2part=Fp_mul(gel(q, 2), x, D);//q[2]*x
	for(GEN y=gen_0;cmpii(y, D)<0;y=addis(y, 1)){
	  if(!equali1(gcdii(x, y))) continue;//gcd(x, y)!=1
	  GEN val=Fp_addmul(c1, Fp_addmul(c2part, gel(q, 3), y, D), y, D);//q(x, y) mod D
	  if(!equali1(gcdii(val, D))) continue;//gcd(D, q(x, y))!=1
	  vectrunc_append(L, val);
	}
  }
  L=gerepileupto(mid, ZV_sort_uniq(L));//All primitive values of q coprime to D are equal to a square times an element of L.
  GEN gens=vectrunc_init(lg(L));//Stores one element per part.
  GEN sq=modsquares(D, 1);//coprime squares mod D
  long lg;
  while(lg(L)>1){
	GEN x=gel(L, 1);//New rep.
	vectrunc_append(gens, x);
	GEN xgen=cgetg_copy(sq, &lg);
	for(long i=1;i<lg;i++) gel(xgen, i)=Fp_mul(gel(sq, i), x, D);
	xgen=ZV_sort(xgen);//No need for uniq, already distinct guaranteed.
	L=setminus(L, xgen);
  }
  long lsq=lg(sq), lgen=lg(gens);
  GEN ret=vectrunc_init((lsq-1)*(lgen-1)+1);
  for(long i=1;i<lsq;i++){
	for(long j=1;j<lgen;j++) vectrunc_append(ret, Fp_mul(gel(sq, i), gel(gens, j), D));//Multiply sq and gens
  }
  return gerepileupto(top, ZV_sort(ret));//Again, no need for uniq
}

//Returns the smallest prime primitively represented by all BQFs in v in the range pmin to pmax
GEN bqf_primetuplereps(GEN v, GEN pmin, GEN pmax){
  pari_sp top=avma;
  long lgv;
  GEN vqfb=cgetg_copy(v, &lgv), Dlist=cgetg_copy(v, &lgv);
  for(long i=1;i<lgv;i++){gel(vqfb, i)=Qfb0(gel(v, i), NULL, NULL);gel(Dlist, i)=negi(bqf_disc(gel(v, i)));}//Convert to Qfbs and get discs
  Dlist=ZV_sort(Dlist);
  int pbigenough=0;
  forprime_t T;
  if(gequal0(pmax)){pmax=pmin;pmin=gen_2;}//In case we only supply a maximum.
  forprime_init(&T, pmin, pmax);
  GEN p;
  GEN fact=mkmat2(mkcol(gen_1), mkcol(gen_1));//Matrix [1 1], will be modified to be [p 1]=factor(p)
  while((p=forprime_next(&T))){
	int success=1;
	if(!pbigenough){//Need to check that p is coprime to all discs.
	  for(long i=1;i<lgv;i++) if(!equali1(gcdii(p, gel(Dlist, i)))){success=0;break;}
	  if(!success) continue;//p not always coprime
	  if(cmpii(gel(Dlist, lgv-1), p)<0) pbigenough=1;//p>all discs, no need to check for gcd=1 as it's guaranteed.
	}
	gcoeff(fact, 1, 1)=p;
	for(long i=1;i<lgv;i++){
	  if(lg(qfbsolve(gel(vqfb, i), fact, 0))==1){success=0;break;}
	}
	if(success) return gerepilecopy(top, p);
  }
  return gc_const(top, gen_0);
}

//Returns 1 if there are residue classes that could contain primes simultaneously represented by all BQF's in v, and 0 if not (some restrictions mod n).
int bqf_tuplevalid(GEN v){
  pari_sp top=avma;
  long lv;
  GEN res=cgetg_copy(v, &lv);
  for(long i=1;i<lv;i++) gel(res, i)=mod_breakdown(bqf_primesmod(gel(v, i)), bqf_disc(gel(v, i)));
  res=shallowconcat1(res);//Put the residues all together
  return gc_int(top, mod_collapse(res));
}

//Given a bunch of residue classes [res, [p, e, p^e]], this returns 1 if there exists a number that obeys 1 congruence from each class, and 0 if not.
static int mod_collapse(GEN L){
  pari_sp top=avma;
  GEN Lsort=gen_sort(L, NULL, &bqf_tuplevalid_cmp);
  long i0=1;
  for(long i=2;i<=lg(L);i++){
	if(i<lg(L) && equalii(gmael3(Lsort, i, 2, 1), gmael3(Lsort, i-1, 2, 1))) continue;//Go until we hit a new prime
	long i1=i-1;//We go from i0 to i1.
	if(i1==i0){i0=i;continue;}//Fine!
	//Go backwards and keep set intersecting
	GEN resbase=gmael(Lsort, i1, 1);//The baseline to search with
	for(long j=i1-1;j>=i0;j--){
	  resbase=FpV_red(resbase, gmael3(Lsort, j, 2, 3));//Mod the lower prime power
	  resbase=ZV_sort_uniq(resbase);//Sort it
	  resbase=setintersect(resbase, gmael(Lsort, j, 1));
	  if(lg(resbase)==1) return gc_int(top, 0);//No can do
	}
	i0=i;
  }
  return gc_int(top, 1);//Made it through, it works.
}

static int bqf_tuplevalid_cmp(void *data, GEN x, GEN y){//[res, [p, e, p^e]]. We sort by p, then by e.
  int a=cmpii(gmael(x, 2, 1), gmael(y, 2, 1));
  if(a!=0) return a;
  return cmpii(gmael(x, 2, 2), gmael(y, 2, 2));
}


