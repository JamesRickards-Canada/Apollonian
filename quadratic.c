/*Methods to deal with quadratic forms.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"

/*STATIC DECLARATIONS*/

/*SECTION 3: CLASS GROUP*/
static int qfbequal1(GEN q);


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

/*This computes the narrow class group*/
GEN
qfbnarrow(GEN D, long prec)
{
  pari_sp av = avma;
  if (!isdisc(D)) pari_err_TYPE("D must be a discriminant, i.e. a non-square integer equivalent to 0 or 1 modulo 4.", D);
  GEN cgp = quadclassunit0(D, 0, NULL, prec);/*Class group*/
  if (signe(D) < 0 || quadunitnorm(D) == -1) return cgp;/*0 or norm(fund. unit)=-1*/
  GEN minus1;/*-id, the new form to add*/
  if (mod2(D)) minus1 = mkqfb(gen_m1, gen_1, shifti(subis(D, 1), -2), D);/*[-1, 1, (D-1)/4]*/
  else minus1 = mkqfb(gen_m1, gen_0, shifti(D, -2), D);/*[-1, 0, D/4]*/
  if (equali1(gel(cgp, 1))) return gerepilecopy(av, mkvec4(gen_2, mkvec(gen_2), mkvec(minus1), gel(cgp, 4)));/*h^+(D)=2, done separately.*/
  GEN newsize = shifti(gel(cgp, 1), 1);
  GEN orders = gel(cgp, 2);
  GEN gens = gel(cgp, 3);
  long lgens = lg(gens), firstm1 = 0, i;/*Start by checking if any of our current generators power to -1.*/
  for (i = 1; i < lgens; i++) {
	GEN qpow = qfbpow(gel(gens, i), gel(orders, i));
    if (!qfbequal1(qpow)) { firstm1 = i; break; }/*note: qfbpow is either -1 or 1*/
  }
  if (!firstm1) {/*Nothing powers to -1, just add in new Z/2Z copy*/
    long lastodd = 0;/*Index of largest odd order*/
    for (i = 1; i < lgens; i++) { if (mod2(gel(orders, i))) { lastodd = i; break; } }
    if (!lastodd) {/*All even, add a new Z/2Z at the end.*/
	  orders = vec_append(orders, gen_2);
	  gens = vec_append(gens, minus1);
    }
	else {/*Now there is an odd order element. Composing this with minus1 gives the new generator, and is the only change required.*/
	  gel(orders, lastodd) = shifti(gel(orders, lastodd), 1);/*Double order*/
	  gel(gens, lastodd) = qfbcomp(gel(gens, lastodd), minus1);/*q*minus1*/
	}
	return gerepilecopy(av, mkvec4(newsize, orders, gens, gel(cgp, 4)));
  }
  /*Now we have an element powering to -1 and not 1. We must fix the previous and subsequent entries.*/
  GEN g = gel(gens, firstm1);
  GEN ai = gel(orders, firstm1);/*Initially this is half the order of g, since g^ai=-1*/
  GEN gpow = g;/*gpow stores g^(.5order(g)/ai)*/
  for (i = firstm1 + 1; i < lgens; i++) {/*Fix the later ones by composing with gpow^(.5order(g)/ai) if they power to -1 instead of 1*/
    GEN qpow = qfbpow(gel(gens, i), gel(orders, i));
    if (qfbequal1(qpow)) continue;/*OK already*/
    gpow = qfbpow(gpow, diviiexact(ai, gel(orders, i)));
    ai = gel(orders, i);/*Updating ai*/
    gel(gens, i) = qfbcomp(gel(gens, i), gpow);/*Modifying the generator.*/
  }
  /*We are done fixing the generators past firstm1; they all have the correct order. We shift focus to those before it.*/
  /*In order to keep the divisibility of orders, we must figure out which generator should double in order.*/
  long optplace = 1;/*The place where we should be inserting the 2x in the class group*/
  long v = vali(gel(orders, firstm1));/*2-adic valuation*/  
  for (i = firstm1 - 1; i > 0; i--) {
    if (vali(gel(orders, i)) > v) { optplace = i + 1; break; }
  }/*If this for method never triggers, we are left with optplace=1=the start, as desired.*/
  /*Now we must modify the elements between firstm1 and optplace to only double the order of optplace and retain the other orders. In fact, we only need to modify these two elements*/
  if (firstm1 != optplace) {/*If equal, there is nothing left to do! If !=, we must shift the elements appropriately*/
    if (!v) {/*The elements between firstm1 and optplace are all ODD powered, so we can adjust firstm1 and optplace by minus1*/
      gel(gens, firstm1) = qfbcomp(gel(gens, firstm1), minus1);
      gel(gens, optplace) = qfbcomp(gel(gens, optplace), minus1);
    }
    else {/*firstm1 is g1 and has order 2^((v+1)h1), optplace is g2 and has order 2^(vh2) with h1, h2 odd and h1|h2.*/
      GEN h2 = shifti(gel(orders, optplace), -v);
	  /*Replace g1 by g1^2*g2^h2, which has order 2^vh1, and replace g2 by g1*g2, which has order 2^(v+1)h2. These elements generate the same subgroup as well since 2-h2 is coprime to 2h2.*/
      gel(gens, firstm1) = qfbcomp(qfbpow(g, gen_2), qfbpow(gel(gens, optplace), h2));/*g1^2*g2^h2*/
	  gel(gens, optplace) = qfbcomp(g, gel(gens, optplace));/*g1*g2*/
    }
  }
  gel(orders, optplace) = shifti(gel(orders, optplace), 1);//Doubling the correct place
  return gerepilecopy(av, mkvec4(newsize, orders, gens, gel(cgp, 4)));
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





//STATIC DECLARATIONS

//COMPOSITION/CLASS GROUP
static GEN bqf_ncgp_nonfundnarrow(GEN cgp, GEN D, GEN rootD);

//TUPLE REPS OF PRIMES
static int mod_collapse(GEN L);
static int bqf_tuplevalid_cmp(void *data, GEN x, GEN y);






//Composes q1, q2. Code adapted from "static void qfb_comp(GEN z, GEN x, GEN y)" in Qfb.c
GEN bqf_comp(GEN q1, GEN q2){
  pari_sp top=avma;
  GEN n, c, d, y1, v1, v2, c3, m, p1, r;
  if(ZV_equal(q1, q2)) return bqf_square(q1);
  n=shifti(subii(gel(q2, 2), gel(q1, 2)), -1);
  v1=gel(q1, 1);
  v2=gel(q2, 1);
  c=gel(q2, 3);
  d=bezout(v2, v1, &y1, NULL);
  if (equali1(d)) m=mulii(y1, n);
  else{
    GEN s = subii(gel(q2, 2), n);
    GEN x2, y2, d1=bezout(s, d, &x2, &y2); /* x2 s + y2 (x1 v1 + y1 v2) = d1 */
    if(!equali1(d1)){
      v1=diviiexact(v1, d1);
      v2=diviiexact(v2, d1); // gcd = 1 iff q1 or q2 primitive
      v1=mulii(v1, gcdii(c, gcdii(gel(q1, 3), gcdii(d1, n))));
      c=mulii(c, d1);
    }
    m=addii(mulii(mulii(y1, y2), n), mulii(gel(q2, 3), x2));
  }
  togglesign(m);
  r=modii(m, v1);
  p1=mulii(r, v2);
  c3=addii(c, mulii(r, addii(gel(q2, 2), p1)));
  p1=shifti(p1, 1);
  GEN z=cgetg(4, t_VEC);
  gel(z,1)=mulii(v1, v2);
  gel(z,2)=addii(gel(q2, 2), p1);
  gel(z,3)=diviiexact(c3, v1); 
  return gerepileupto(top, z);
}

//bqf_comp with reduction. If D<0, can pass rootD as NULL
GEN bqf_comp_red(GEN q1, GEN q2, GEN rootD, int Dsign){
  pari_sp top=avma;
  GEN qcomp=bqf_comp(q1, q2);
  return gerepileupto(top, bqf_red(qcomp, rootD, Dsign, 0));
}

//bqf_comp_red with typecheck
GEN bqf_comp_tc(GEN q1, GEN q2, int tored, long prec){
  pari_sp top=avma;
  if(tored==1){
    GEN D1=bqf_checkdisc(q1);
    GEN D2=bqf_checkdisc(q2);
    if(!equalii(D1, D2)) pari_err_TYPE("discriminants not equal",q1);
    if(signe(D1)==1) return gerepileupto(top, bqf_comp_red(q1, q2, gsqrt(D1, prec), 1));
    avma=top;
    return bqf_comp_red(q1, q2, NULL, -1);
  }
  bqf_check(q1);
  bqf_check(q2);
  return bqf_comp(q1, q2);
}

//Returns the identity element
GEN bqf_idelt(GEN D){
  pari_sp top = avma;
  if(smodis(D, 2)==0){//D even
    GEN nD=negi(D);
    GEN q=cgetg(4, t_VEC);
    gel(q, 1)=gen_1;
    gel(q, 2)=gen_0;
    gel(q, 3)=shifti(nD, -2);//-D/4
    return gerepileupto(top, q);
  }//Now D is odd
  GEN omD=subsi(1, D);
  GEN q=cgetg(4, t_VEC);
  gel(q, 1)=gen_1;
  gel(q, 2)=gen_1;
  gel(q, 3)=shifti(omD, -2);//(1-D)/4
  return gerepileupto(top, q);
}

//Identifies q in terms of the given basis. Must pass in the lexicographic ncgp. If doing this on a large set, should sort this first, and use a set search.
GEN bqf_identify(GEN ncgp_lexic, GEN q, GEN rootD, int Dsign){
  pari_sp top=avma;
  long ind=itos(bqf_isequiv_set(q, gel(ncgp_lexic, 3), rootD, Dsign, 0));
  avma=top;
  if(ind==-1) pari_err_TYPE("Error, did not find q. Must have supplied the wrong data", q);
  return bqf_lexicind_tobasis(gel(ncgp_lexic, 2), ind);
}

//bqf_identify with type checking.
GEN bqf_identify_tc(GEN ncgp_lexic, GEN q, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)==-1) return gerepileupto(top, bqf_identify(ncgp_lexic, q, NULL, -1));
  GEN rootD=gsqrt(D, prec);
  return gerepileupto(top, bqf_identify(ncgp_lexic, q, rootD, 1));
}

//Given the VECSMALL of orders of elements in the narrow class group and an index ind, this returns the tuple [e1,...,er] where cgp[3][ind]=g1^e1*...*gr^er
GEN bqf_lexicind_tobasis(GEN orders, long ind){
  pari_sp top=avma;
  long lord=lg(orders);
  GEN es=cgetg(lord, t_VECSMALL);
  long u=ind-1, n, v;//we start with ind-1
  for(long i=lord-1;i>0;i--){//Go backwards.
    n=orders[i];
    u=sdivss_rem(u, n, &v);//u_{i+1}=u_i*n_i+v_i; u_{r+1}=i-1.
    es[i]=v;
  }
  return gerepilecopy(top, es);
}

//This computes the narrow class group
GEN bqf_ncgp(GEN D, long prec){
  pari_sp top = avma;
  GEN fd=coredisc(D);
  if(equalii(fd, D)){
    GEN field=Buchall(gsub(gsqr(pol_x(0)), D), 0, prec);//The field Q(sqrt(D))
    GEN gpclgp=bnfnarrow(field);//The narrow class group of the maximal order
    if(equali1(gel(gpclgp, 1))){//Class number 1, may as well do this separately
      GEN rvec=cgetg(4, t_VEC);
      gel(rvec, 1)=gen_1;
      gel(rvec, 2)=mkvecsmall(1);
      gel(rvec, 3)=cgetg(2, t_VEC);
      gmael(rvec, 3, 1)=bqf_idelt(D);
      return gerepileupto(top, rvec);
    }
    long j=lg(gel(gpclgp, 3));
    long ngens=j-1, lx;//Number of generators
    GEN idtobqf=cgetg_copy(gel(gpclgp,3), &lx);
    for(long i=1;i<=ngens;i++) gel(idtobqf, i)=ideal_tobqf(field, gmael(gpclgp, 3, i));//Making the generators into quadratic forms
    GEN rootD=gsqrt(D, prec);
    GEN rvec=cgetg(4, t_VEC);
    gel(rvec, 1)=icopy(gel(gpclgp, 1));//clgp size
    gel(rvec, 2)=cgetg(lg(gel(gpclgp, 3)), t_VECSMALL);//Shape of group
    gel(rvec, 3)=cgetg_copy(gel(gpclgp, 3), &lx);//Generators
    for(long i=1;i<=ngens;i++){//Inputting backwards, as bnfnarrow does the reverse order to what we want.
      j--;
      gel(rvec, 2)[i]=itos(gmael(gpclgp, 2, j));
      gmael(rvec, 3, i)=bqf_red(gel(idtobqf, j), rootD, signe(D), 0);
    }
    return gerepileupto(top, rvec);
  }//OK, so now we have a non-fundamental discrimiant
  GEN FULLgp=quadclassunit0(D, 0, NULL, prec);
  GEN rootD=gsqrt(D, prec);
  if(signe(D)==1 && gequal1(quadnorm(quadunit0(D, -1)))) return gerepileupto(top, bqf_ncgp_nonfundnarrow(FULLgp, D, rootD));//Must double in size
  GEN temp, newgens;//No modification needed, clgp=narrow clgp in this case
  long lx;
  if(equali1(gel(FULLgp, 1))){//Class number 1, may as well do this separately
    GEN rvec=cgetg(4, t_VEC);
    gel(rvec, 1)=gen_1;
    gel(rvec, 2)=mkvecsmall(1);//Universal object, okay
    gel(rvec, 3)=cgetg(2, t_VEC);
    gmael(rvec, 3, 1)=bqf_idelt(D);
    return gerepileupto(top, rvec);
  }
  newgens=cgetg_copy(gel(FULLgp, 3), &lx);
  for(long i=1;i<lg(gel(FULLgp, 3));i++){
    temp=gtovec(gmael(FULLgp, 3, i));
    setlg(temp, 4);//Truncate
    gel(newgens, i)=bqf_red(temp, rootD, signe(D), 0);
  }
  GEN rvec=cgetg(4, t_VEC);
  gel(rvec, 1)=icopy(gel(FULLgp, 1));//Size
  gel(rvec, 2)=cgetg(lg(gel(FULLgp, 2)), t_VECSMALL);
  gel(rvec, 3)=cgetg_copy(gel(FULLgp, 3), &lx);
  long j=lx;
  for(long i=1;i<lx;i++){
    j--;
    gel(rvec, 2)[j]=itos(gmael(FULLgp, 2, i));
    gmael(rvec, 3, j)=ZV_copy(gel(newgens, i));
  }
  return gerepileupto(top,rvec);
  //GEN modulus;
  //if(signe(D)==1) modulus=mkvec2(sqrti(diviiexact(D,fd)),mkvec2s(1,1));
  //else modulus=sqrti(diviiexact(D,fd));
  //GEN rayfield=bnrinit0(field, modulus, 1);//Initialize the ray class group
  //gpclgp=member_clgp(rayfield);
}

//D>0 is non-fundamental and the fundamental unit has norm 1. cgp is the current class group, the output of quadclassunit0. This modifies cgp to be the FULL narrow class group.
static GEN bqf_ncgp_nonfundnarrow(GEN cgp, GEN D, GEN rootD){
  pari_sp top=avma;
  GEN minus1=cgetg(4, t_VEC);
  gel(minus1, 1)=gen_m1;
  if(smodis(D, 2)==0) gel(minus1, 2)=gen_0;
  else gel(minus1, 2)=gen_1;
  gel(minus1, 3)=shifti(subii(D, gel(minus1, 2)), -2);//[-1,0/1, (D-0/1)/4]
  if(equali1(gel(cgp, 1))){//Class number 1=narrow class 2, may as well do this separately
    GEN rvec=cgetg(4, t_VEC);
    gel(rvec, 1)=gen_2;
    gel(rvec, 2)=mkvecsmall(2);
    gel(rvec, 3)=cgetg(2, t_VEC);
    gmael(rvec, 3, 1)=ZV_copy(minus1);
    return gerepileupto(top, rvec);
  }
  //Now we power through
  GEN neworder=ZV_copy(gel(cgp, 2));
  long lx;
  GEN newgens=cgetg_copy(gel(cgp, 3), &lx), temp;
  for(long i=1;i<lx;i++){
    temp=gtovec(gmael(cgp, 3, i));
    setlg(temp, 4);//Truncate
    gel(newgens, i)=temp;
  }//Just making the initial conversion to BQFs on our form.
  long firstm1=-1;
  for(long i=1;i<lx;i++){
    temp=bqf_pow_red(gel(newgens, i), gel(neworder, i), rootD, 1);
    if(equali1(ibqf_isequiv(minus1, temp, rootD))){firstm1=i;break;}
  }
  if(firstm1==-1){//Just add in new Z/2Z copy
    long lastodd=-1;
    for(long i=1;i<lx;i++){if(smodis(gel(neworder, i), 2)==1){lastodd=i;break;}}
    if(lastodd==-1){//All even, add a new Z/2Z
      GEN rvec=cgetg(4,t_VEC);
      long j=lx;
      gel(rvec, 1)=shifti(gel(cgp, 1), 1);
      gel(rvec, 2)=cgetg(j+1, t_VECSMALL);
      gel(rvec, 3)=cgetg(j+1, t_VEC);
      for(long i=1;i<lx;i++){
        gel(rvec, 2)[j]=itos(gel(neworder, i));
        gmael(rvec, 3, j)=ZV_copy(gel(newgens, i));
        j--;
      }
      gel(rvec, 2)[1]=2;
      gmael(rvec, 3, 1)=ibqf_red(minus1, rootD);
      return gerepileupto(top, rvec);
    }//OK, now there is an odd number.
    GEN rvec=cgetg(4, t_VEC);
    long j=lx-1;
    gel(rvec, 1)=shifti(gel(cgp, 1), 1);
    gel(rvec, 2)=cgetg(lx, t_VECSMALL);
    gel(rvec, 3)=cgetg_copy(newgens, &lx);
    for(long i=1;i<lastodd;i++){
      gel(rvec, 2)[j]=itos(gel(neworder, i));
      gmael(rvec, 3, j)=ZV_copy(gel(newgens, i));
      j--;
    }
    gel(rvec, 2)[j]=2*itos(gel(neworder, lastodd));
    gmael(rvec, 3, j)=bqf_comp_red(minus1, gel(newgens, lastodd), rootD, 1);
    j--;
    for(long i=lastodd+1;i<lx;i++){
      gel(rvec, 2)[j]=itos(gel(neworder, i));
      gmael(rvec, 3, j)=ZV_copy(gel(newgens, i));
      j--;
    }
    return gerepileupto(top, rvec);
  }//OK SO now we have an element powering to -1 and not 1. We must fix the previous and subsequent entries.
  GEN g=gel(newgens, firstm1), gpow=g, aj=gel(neworder, firstm1);//gpow represents g^(.5order(g)/aj)
  for(long j=firstm1+1;j<lx;j++){//Fix the later ones by modifying by gpow^(a_j/a_g) if they power to -1 instead of 1
    temp=bqf_pow_red(gel(newgens, j), gel(neworder, j), rootD, 1);
    if(equali1(ibqf_isequiv(minus1, temp, rootD))){//If =0 then we have the identity and there is no need to fix it.
        gpow=bqf_pow_red(gpow, diviiexact(aj, gel(neworder,j)), rootD, 1);
        aj=gel(neworder, j);//Updating
        gel(newgens, j)=bqf_comp_red(gel(newgens, j), gpow, rootD, 1);//Fixing it
    }
  }//Now we are done fixing the generators past firstm1; we shift focus to those before it.  
  long optplace=1;//The place where we should be inserting the 2x in the class group
  long lastv=vali(gel(neworder, firstm1));//2-adic valuation  
  for(long i=firstm1-1;i>0;i--){
    if(vali(gel(neworder, i))>lastv){optplace=i+1;break;}
  }//If this for method never triggers, we are left with optplace=1=the start, as desired.
  //Now we must modify the elements between firstm1 and lastv to only double the order of lastv and retain the other orders. In fact, we only need to modify these two elements
  if(firstm1!=optplace){//If equal, there is nothing left to do! If !=, we must shift the elements appropriately
    if(lastv==0){//The elements between firstm1 and optplace are all ODD powered, so we can adjust firstm1 and lastv by minus1
      gel(newgens, firstm1)=bqf_comp_red(minus1, gel(newgens, firstm1), rootD, 1);
      gel(newgens, optplace)=bqf_comp_red(minus1, gel(newgens, optplace), rootD, 1);
    }
    else{//firstm1 is g1 and has order 2^(v+1)h1, optplace is g2 and has order g^vh2 with h1, h2 odd and h1|h2. Replace g1 by g1^(2^(v+1))*g2^h2 which has order 2^vh1, and replace g2 by g2^(2^v)*g1^h1 which has order 2^(v+1)h2. These elements generate the same subgroup as well since 2h2 is coprime to 2^(2v+1)-h1h2.
      GEN twovali=shifti(gen_1, lastv);
      GEN h1=diviiexact(gel(neworder, firstm1), twovali);
      GEN h2=diviiexact(gel(neworder, optplace), twovali);
      gel(newgens, firstm1)=bqf_comp_red(bqf_square_red(bqf_pow_red(g, twovali, rootD, 1), rootD, 1), bqf_pow_red(gel(newgens, optplace), h2, rootD, 1), rootD, 1);
      gel(newgens, optplace)=bqf_comp_red(bqf_pow_red(g, h1, rootD, 1), bqf_pow_red(gel(newgens, optplace), twovali, rootD, 1), rootD, 1);
    }
  }
  gel(neworder, optplace)=shifti(gel(neworder, optplace), 1);//Doubling the correct place
  //Ok, now we just need to reverse the vectors and compile it into rvec and return it.
  GEN rvec=cgetg(4, t_VEC);
  long j=lx-1;
  gel(rvec, 1)=shifti(gel(cgp, 1), 1);
  gel(rvec, 2)=cgetg(lx, t_VECSMALL);
  gel(rvec, 3)=cgetg_copy(newgens, &lx);
  for(long i=1;i<lx;i++){
    gel(rvec, 2)[j]=itos(gel(neworder, i));
    gmael(rvec, 3, j)=ZV_copy(gel(newgens, i));
    j--;
  }
  return gerepileupto(top, rvec);
}

//Same as bqf_ncgp, but the third element is a lexicographic ordering of the elements of the ncgp.
GEN bqf_ncgp_lexic(GEN D, long prec){
  pari_sp top = avma;
  GEN ncgp=bqf_ncgp(D, prec);//Get the ncgp
  GEN rootD=gsqrt(D, prec);
  long nclno=itos(gel(ncgp, 1));//Class number as a long
  GEN rvec=cgetg(4, t_VEC);
  gel(rvec, 1)=icopy(gel(ncgp, 1));
  gel(rvec, 2)=zv_copy(gel(ncgp, 2));//Copy the first two elements
  gel(rvec, 3)=cgetg(nclno+1, t_VEC);//The forms vector
  gmael(rvec, 3, 1)=bqf_red(bqf_idelt(D), rootD, signe(D), 0);//Start with Id
  long f1, f2=1, pow, j, k=2;
  for(long i=lg(gel(ncgp, 2))-1;i>=1;i--){
    f1=1;
    f2=k-1;
    for(pow=1;pow<=gel(ncgp, 2)[i]-1;pow++){
      k=f2+1;
      for(j=f1;j<=f2;j++){
        gmael(rvec, 3, k)=bqf_comp_red(gmael(rvec, 3, j), gmael(ncgp, 3, i), rootD, signe(D));
        k++;
      }
      f1=f2+1;
      f2=k-1;//update f1, f2
    }
  }
  return gerepileupto(top,rvec);
}

//q^n without reduction
GEN bqf_pow(GEN q, GEN n){
  pari_sp top=avma;
  if(gequal0(n)){GEN D=bqf_disc(q);return gerepileupto(top, bqf_idelt(D));}
  GEN qpow=ZV_copy(q);
  if(signe(n)==-1){
    togglesign_safe(&gel(qpow, 2));//Taking the inverse
    if(equalis(n, -1)) return qpow;
    q=qpow;//Need to repoint q to here
    n=negi(n);
  }
  else if(equali1(n)) return qpow;
  GEN nbin=binary_zv(n);
  for(long i=2;i<lg(nbin);i++){
    qpow=bqf_square(qpow);
    if(nbin[i]==1) qpow=bqf_comp(qpow, q);
  }
  return gerepileupto(top, qpow);
}

//q^n with reduction
GEN bqf_pow_red(GEN q, GEN n, GEN rootD, int Dsign){
  pari_sp top=avma;
  if(gequal0(n)){GEN D=bqf_disc(q);return gerepileupto(top, bqf_idelt(D));}
  GEN qpow=ZV_copy(q);
  if(signe(n)==-1){
    togglesign_safe(&gel(qpow, 2));//Taking the inverse
    if(equalis(n, -1)) return gerepileupto(top, bqf_red(qpow, rootD, Dsign, 0));
    q=qpow;//Need to repoint q to here
    n=negi(n);
  }
  else if(equali1(n)) return gerepileupto(top, bqf_red(qpow, rootD, Dsign, 0));
  GEN nbin=binary_zv(n);
  for(long i=2;i<lg(nbin);i++){
    qpow=bqf_square_red(qpow, rootD, Dsign);
    if(nbin[i]==1) qpow=bqf_comp_red(qpow, q, rootD, Dsign);
  }
  return gerepileupto(top, qpow);
}

//bqf_pow_red with typecheck
GEN bqf_pow_tc(GEN q, GEN n, int tored, long prec){
  pari_sp top=avma;
  if(typ(n)!=t_INT) pari_err_TYPE("Please enter an integer n", n);
  if(tored==1){
    GEN D=bqf_checkdisc(q);
    if(signe(D)==1) return gerepileupto(top, bqf_pow_red(q, n, gsqrt(D, prec), 1));
    avma=top;
    return bqf_pow_red(q, n, NULL, -1);
  }//Now no reduction
  bqf_check(q);
  return bqf_pow(q, n);
}

//Squares the bqf q; code adapted from the PARI function "static void qfb_sqr(GEN z, GEN x)" in Qfb.c
GEN bqf_square(GEN q){
  pari_sp top=avma;
  GEN c, d1, x2, v1, v2, c3, m, p1, r;
  d1=bezout(gel(q,2), gel(q,1), &x2, NULL);//d1=gcd(A,B)
  c=gel(q, 3);
  m=mulii(c, x2);
  if(equali1(d1)) v1=v2=gel(q, 1);
  else{
    v1=diviiexact(gel(q,1), d1);
    v2=mulii(v1, gcdii(d1,c)); // = v1 iff q primitive 
    c=mulii(c, d1);
  }
  togglesign(m);
  r = modii(m, v2);
  p1 = mulii(r, v1);
  c3 = addii(c, mulii(r, addii(gel(q, 2), p1)));
  p1=shifti(p1, 1);
  GEN z=cgetg(4, t_VEC);
  gel(z,1) = mulii(v1, v2);
  gel(z,2) = addii(gel(q, 2), p1);
  gel(z,3) = diviiexact(c3, v2);
  return gerepileupto(top, z);
}

//bqf_square with reduction. If D<0, can pass rootD as NULL
GEN bqf_square_red(GEN q, GEN rootD, int Dsign){
  pari_sp top=avma;
  GEN qsqr=bqf_square(q);
  return gerepileupto(top, bqf_red(qsqr, rootD, Dsign, 0));
}

//bqf_square_red with typecheck
GEN bqf_square_tc(GEN q, int tored, long prec){
  pari_sp top=avma;
  if(tored==1){
    GEN D=bqf_checkdisc(q);
    if(signe(D)==1) return gerepileupto(top, bqf_square_red(q, gsqrt(D,prec), 1));
    avma=top;
    return bqf_square_red(q, NULL, -1);
  }
  bqf_check(q);
  return bqf_square(q);
}

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


