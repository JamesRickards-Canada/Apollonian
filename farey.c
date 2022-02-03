//Farey methods

//INCLUSIONS

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "apol.h"
#endif

//Returns the multiset of depths of fractions with reduced denominator n.
GEN fareydenom_depths(long n){
  pari_sp top=avma;
  long l=eulerphiu(n);
  GEN v=cgetg(l+1, t_VECSMALL);
  long ind=1;
  for(long num=1;num<n;num++){
	if(ugcd(num, n)>1) continue;
	v[ind]=fareydepth(mkfracss(num, n));
	ind++;
  }
  return gerepilecopy(top, v);
}

//Returns the most common depth
long fareydenom_dmode(long n){
  pari_sp top=avma;
  GEN v=fareydenom_depths(n);
  GEN s=vecsmallcount(v);
  long ind=vecindexmax(gel(s, 2));//Mode
  return gc_long(top, gel(s, 1)[ind]);
}

//Returns the depth of r, i.e. the number of times we must go "up" until we hit [0, 1]. If [a, b] is a pair with mediant r, then [a, b] has depth d whilst r has depth d+1.
long fareydepth(GEN r){
  pari_sp top=avma;
  long d=0;
  if(typ(r)!=t_VEC){
    r=fareyup(r);
	d++;
  }
  while(!gequal0(gel(r, 1)) || !gequal1(gel(r, 2))){
	r=fareyup(r);
	d++;
  }
  return gc_long(top, d);
}

//Returns the Farey pair that generates r. r can either be a rational or a Farey pair (not checked, so be careful).
GEN fareyup(GEN r){
  pari_sp top=avma;
  if(typ(r)==t_VEC){
	GEN n1=numer_i(gel(r, 1)), d1=denom_i(gel(r, 1));
	GEN n2=numer_i(gel(r, 2)), d2=denom_i(gel(r, 2));
	if(cmpii(d2, d1)>0){//Go right
	  GEN newd=subii(d2, d1);
	  if(equali1(newd)) return gerepilecopy(top, mkvec2(gel(r, 1), gen_1)); //Must be 1
	  return gerepilecopy(top, mkvec2(gel(r, 1), mkfrac(subii(n2, n1), newd)));
	}
	GEN newd=subii(d1, d2);
	if(equali1(newd)) return gerepilecopy(top, mkvec2(gen_0, gel(r, 2)));//Must be 0
	return gerepilecopy(top, mkvec2(mkfrac(subii(n1, n2), newd), gel(r, 2)));
  }
  //Now rational
  GEN num=numer_i(r);
  GEN den=denom_i(r);
  if(equali1(den)) return gerepileupto(top, mkvec2(gen_0, gen_1));//0 or 1, stop there.
  GEN bez=gcdext0(num, den);
  GEN a=gel(bez, 1), b=gel(bez, 2), r1, r2;//a*num+b*den=1
  if(signe(a)==1){//a>0, b<0, -b/a<r
    if(equali1(a)) r1=gen_0;
    else r1=mkfrac(negi(b), a);
	GEN newd=subii(den, a);
	if(equali1(newd)) r2=gen_1;
	else r2=mkfrac(addii(num, b), newd);
  }
  else{//a<0, b>0, r<b/-a
    if(gequalm1(a)) r2=gen_1;
    else r2=mkfrac(b, negi(a));
	GEN newd=addii(den, a);
	if(equali1(newd)) r1=gen_0;
    else r1=mkfrac(subii(num, b), newd);
  }
  return gerepilecopy(top, mkvec2(r1, r2));
}


