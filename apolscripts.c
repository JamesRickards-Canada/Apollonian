#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "apol.h"
#endif

//STATIC METHOD DECLARATIONS


//Returns the polynomials generating abelian subfields of a specified degree. G can be either a galoisinit or a list of polynomials. We output the polynomials and the groups.
GEN abelianfields(GEN G, long deg){
  pari_sp top=avma;
  long lG;
  if(lg(G)>2 && typ(gel(G, 2))==t_VEC){//We have a Galois init
    GEN GG=galoissubfields(G, 0, -1);
	G=cgetg_copy(GG, &lG);
	for(long i=1;i<lG;i++) gel(G, i)=gmael(GG, i, 1);//Copy the polynomials
  }
  else lG=lg(G);
  GEN pols=vectrunc_init(lG);
  GEN gps=vectrunc_init(lG);
  for(long i=1;i<lG;i++){
	GEN pol=gel(G, i);
	if(deg>0 && poldegree(pol, -1)!=deg) continue;//Wrong degree
	GEN G2=galoisinit(pol, NULL);
	if(gequal0(G2)) continue;//Not Galois
	GEN gab=galoisisabelian(G2, 0);
	if(gequal0(gab) && typ(gab)!=t_MAT) continue;//Not abelian. We include the type check as Q outputs [;]=0.
	vectrunc_append(pols, pol);
	vectrunc_append(gps, gab);
  }
  return gerepilecopy(top, mkvec2(pols, gps));
}

//Returns the polynomials generating galois subfields of a specified degree. G can be either a galoisinit or a list of polynomials.
GEN galoisfields(GEN G, long deg){
  pari_sp top=avma;
  long lG;
  if(lg(G)>2 && typ(gel(G, 2))==t_VEC){//We have a Galois init
    GEN GG=galoissubfields(G, 0, -1);
	G=cgetg_copy(GG, &lG);
	for(long i=1;i<lG;i++) gel(G, i)=gmael(GG, i, 1);//Copy the polynomials
  }
  else lG=lg(G);
  GEN pols=vectrunc_init(lG);
  for(long i=1;i<lG;i++){
	GEN pol=gel(G, i);
	if(deg>0 && poldegree(pol, -1)!=deg) continue;//Wrong degree
	GEN G2=galoisinit(pol, NULL);
	if(!gequal0(G2)) vectrunc_append(pols, pol);
  }
  return gerepilecopy(top, pols);
}

//Given a galoisinit, this checks if G is dihedral. If not, returns 0. If so, 1 if it's (Z/2Z)^n, and [sgp, disc] otherwise, where sgp is the maximal abelian subgroup and disc is the discriminant of the fixed field of sgp.
GEN galoisisdihedral(GEN G){
  pari_sp top=avma;
  GEN isab=galoisisabelian(G, 0);
  if(!gequal0(isab)){//Abelian
	GEN twogrpmat=gmul(matid(lg(isab)-1), gen_2);
	if(gequal(isab, twogrpmat)) return gc_const(top, gen_1);//Abelian and a 2-group
    return gc_const(top, gen_0);//Abelian, not a 2-group
  }
  GEN elts=gal_get_group(G);//There must be an element of order >2, else we have an abelian 2-group
  long le=lg(elts), ind;
  for(long i=1;i<le;i++){
	if(perm_orderu(gel(elts, i))<=2) continue;
	ind=i;//Found element of order >2, must be in the abelian part.
	break;
  }
  GEN abpart=vectrunc_init(le);
  GEN inv=NULL;
  for(long i=1;i<le;i++){//Finding the abelian part.
	if(perm_commute(gel(elts, i), gel(elts, ind))) vectrunc_append(abpart, gel(elts, i));
	else if(!inv && perm_orderu(gel(elts, i))==2) inv=gel(elts, i);
  }
  long lab=lg(abpart);
  if(2*lab-1!=le) return gc_const(top, gen_0);//Half of the elements aren't abelian
  GEN iperm=identity_perm(le-1);
  for(long i=1;i<lab;i++){
	if(!gequal(iperm, perm_mul(gel(abpart, i), perm_conj(inv, gel(abpart, i))))) return gc_const(top, gen_0);//Doesn't conjugate correctly
  }
  //Ok now it is dihedral.
  GEN quad=galoisfixedfield(G, abpart, 1, -1);//Must be quadratic
  GEN a=polcoef_i(quad, 2, 0);
  GEN b=polcoef_i(quad, 1, 0);
  GEN c=polcoef_i(quad, 0, 0);
  GEN D=subii(sqri(b), shifti(mulii(a, c), 2));
  D=coredisc(D);
  return gerepilecopy(top, mkvec2(abpart, D));
}

//Returns the polynomial generating the Hilbert class field for Q(sqrt(D)).
GEN hilbertpoly(GEN D, int redbest, long prec){
  pari_sp top=avma;
  GEN pol1=quadhilbert(D, prec);
  GEN pol2=mkpoln(3, gen_1, gen_0, negi(D));
  GEN pol=polcompositum0(pol1, pol2, 2);
  if(redbest) pol=polredbest(pol, 0);
  return gerepileupto(top, pol);
}

//Assume pol generates a Galois extension of Q, this returns the discriminants of quadratic subfields.
GEN quadsubfields(GEN pol, long prec){
  pari_sp top=avma;
  GEN G=galoisinit(pol, NULL);
  GEN sf=galoissubfields(G, 0, 0);
  long lsf=lg(sf);
  GEN v=vectrunc_init(lsf);
  for(long i=1;i<lsf;i++){
	GEN pol=gmael(sf, i, 1);
	if(degree(pol)!=2) continue;
	GEN b=polcoef_i(pol, 1, 0);
	GEN c=polcoef_i(pol, 0, 0);
	vectrunc_append(v, coredisc(subii(sqri(b), shifti(c, 2))));
  }
  return gerepileupto(top, ZV_sort(v));
}

//Returns the ring class field associated to the order of discriminant D. If check=1, we check it has the required ramification. Otherwise, we don't, so the answer might not be correct if there was not enough precision! We can supply a flag, which has the same meaning as for algdep.
GEN ringpoly(GEN D, long precinc, long flag, int check, long prec){
  pari_sp top=avma;
  long newprec=prec-1+precinc;
  long deg=itos(gel(quadclassunit0(D, 0, NULL, prec), 1));//Degree of the extension
  GEN pol;
  do{
	newprec++;
	GEN elt=gdivgs(gaddsg(smodis(D, 2), gsqrt(D, newprec)), 2);//(D%2+sqrt(D))/2
	GEN j=jell(elt, newprec);//j(...)
	if(gsigne(j)<0) j=gneg(j);//Want it positive for the cube root.
	GEN gener=gsqrtn(j, stoi(3), NULL, newprec);//Generates the ring class field
	pol=polredbest(algdep0(gener, deg, flag), 0);
	if(check){//Now we check it. TO BE ADDED LATER
	  pari_warn(warner,"Checking functionality not yet added! Result may be wrong!\n");
	  check=0;
	}
	else{
	  pari_warn(warner,"Result may not be correct if there is not enough precision\n");
	}
  }while(check);
  return gerepilecopy(top, pol);
}



