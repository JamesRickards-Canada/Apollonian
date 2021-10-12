//Apollonian circle packing methods

//INCLUSIONS

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "apol.h"
#endif

//STATIC DECLARATIONS
static long ZV_maxind(GEN v);

//Checks if v gives 4 circles that generate an integral Apollonian packing, returning 1 if so and 0 else
int apol_check(GEN v){
  pari_sp top=avma;
  if(typ(v)!=t_VEC || lg(v)!=5) pari_err_TYPE("Must be a length 4 integer vector", v);
  for(int i=1;i<=4;i++) if(typ(gel(v, i))!=t_INT) pari_err_TYPE("Entries must be integral", gel(v, i));
  GEN L=gen_0;
  for(int i=1;i<=4;i++) L=addii(L, sqri(gel(v, i)));
  L=shifti(L, 1);
  GEN R=gen_0;
  for(int i=1;i<=4;i++) R=addii(R, gel(v, i));
  R=sqri(R);
  return gc_int(top, equalii(L, R)? 1:0);
}

//Returns all primitive Apollonian root quadruples using the construction from x^2+m^2=d_1d_2 (page 19 of GLMWY Number Theory). This has first entry x=-n.
GEN apol_make(GEN n, GEN m, int red){
  pari_sp top=avma;
  if(typ(n)!=t_INT || typ(m)!=t_INT) pari_err_TYPE("Please input two integers", mkvec2(n, m));
  if(red==1 && signe(n)!=1) return cgetg(1, t_VEC);//Require n>0 for reduced.
  GEN d1d2=addii(sqri(n), sqri(m)), mn=negi(n), twom=shifti(m, 1);
  GEN divs=divisors(d1d2);//x^2+m^2=d_1d_2, and a=x, b=d1-x, c=d2-x, d=-2m+d1+d2-x.
  long ndivs=lg(divs);
  GEN apols=vectrunc_init(ndivs);
  for(long i=1;i<=(ndivs-1)/2;i++){
	GEN d1=gel(divs, i);
	if(red==1 && cmpii(twom, d1)==1) continue;//For root quadruple, we require -n<0<=2m<=d1<=d2
	GEN d2=gel(divs, ndivs-i);
	if(!equali1(gcdii(gcdii(d1, n), d2))) continue;
	GEN d1px=addii(d1, n);
	GEN v=mkvec4(mn, d1px, addii(d2, n), addii(d1px, subii(d2, twom)));
	if(red==2) vectrunc_append(apols, apol_reduce(v, 0));
	else vectrunc_append(apols, v);
  }
  return gerepilecopy(top, apols);
}

//Returns the set of four curvatures when we replace circle i.
GEN apol_move(GEN v, int ind){
  pari_sp top=avma;
  GEN rep=vecsmall_ei(4, ind);
  GEN S=gen_0;
  for(int i=1;i<=4;i++) if(!rep[i]) S=addii(S, gel(v, i));
  S=shifti(S, 1);
  long lv;//lv=5
  GEN newv=cgetg_copy(v, &lv);
  for(int i=1;i<=4;i++) gel(newv, i)=rep[i]? subii(S, gel(v, ind)):icopy(gel(v, i));
  return gerepileupto(top, newv);
}

//Returns a sorted list of curvatures of circles, where we go to depth depth, i.e. we do up to depth circle replacements. I may have to do some careful garbage collection for high depths.
GEN apol_orbit(GEN v, int depth){
  pari_sp top=avma;
  int ind=1;//We reuse ind to track which depth we are going towards.
  GEN W=zerovec(depth);//Tracks the sequence of APC's; W[ind] is at ind-1 depth
  gel(W, 1)=v;//The first one.
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1;//Tracks if we are going forward or not.
  long Nreps=2*(itos(powuu(3, depth))+1)+1, ireps=5;
  GEN reps=cgetg(Nreps, t_VEC);
  for(int i=1;i<=4;i++) gel(reps, i)=gel(v, i);//First 4 reps
  do{//1<=ind<=depth is assumed.
    I[ind]=forward? 1:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
	if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
	//At this point, we can go on with valid and new inputs
	GEN newv=apol_move(gel(W, ind), I[ind]);//Make the move
	gel(reps, ireps)=gel(newv, I[ind]);//Add the new element
	ireps++;//Increment ireps
	if(ind==depth) forward=0;
	else{//We can keep going forward
      ind++;
	  gel(W, ind)=newv;
	  forward=1;
	}
  }while(ind>0);
  return gerepileupto(top, ZV_sort(reps));
}

//Returns a sorted list of curvatures of circles surrounding v[ind]. We go to depth depth, i.e. we do up to depth circle replacements. I may have to do some careful garbage collection for high depths.
GEN apol_orbit_1(GEN v, int depth, int ind){
  pari_sp top=avma;
  GEN v1;
  if(ind==1) v1=v;
  else{//Just shifting it so we are replacing v[1]
    long l;
	v1=cgetg_copy(v, &l);//l=5 now
	gel(v1, 1)=gel(v, ind);
	gel(v1, ind)=gel(v, 1);
	for(int i=2;i<=4;i++) if(i!=ind) gel(v1, i)=gel(v, i);//We don't copy, so v1 is not a safe vector.
  }
  ind=1;//We reuse ind to track which depth we are going towards.
  GEN W=zerovec(depth);//Tracks the sequence of APC's; W[ind] is at ind-1 depth
  gel(W, 1)=v1;//The first one.
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1;//Tracks if we are going forward or not.
  long Nreps=3*itos(int2n(depth))+1, ireps=4;
  GEN reps=cgetg(Nreps, t_VEC);
  for(int i=2;i<=4;i++) gel(reps, i-1)=gel(v, i);//First 3 reps
  do{//1<=ind<=depth is assumed.
    I[ind]=forward? 2:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
	if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
	//At this point, we can go on with valid and new inputs
	GEN newv=apol_move(gel(W, ind), I[ind]);//Make the move
	gel(reps, ireps)=gel(newv, I[ind]);//Add the new element
	ireps++;//Increment ireps
	if(ind==depth) forward=0;
	else{//We can keep going forward
      ind++;
	  gel(W, ind)=newv;
	  forward=1;
	}
  }while(ind>0);
  return gerepileupto(top, ZV_sort(reps));
}

//Returns the quadratic form whose primitive values are the curvatures touching circle v[ind]. The formula is [a+b,a+b+c-d, a+c] if v=[a,b,c,d] and ind=1.
GEN apol_qf(GEN v, int ind){
  pari_sp top=avma;
  GEN is=cgetg(5, t_VECSMALL);
  is[1]=ind;
  for(int i=2;i<=4;i++) is[i]=is[i-1]%4+1;
  GEN apb=addii(gel(v, is[1]), gel(v, is[2]));//a+b
  GEN apbpc=addii(apb, gel(v, is[3]));//a+b+c
  GEN q=cgetg(4, t_VEC);
  gel(q, 1)=icopy(apb);//a+b
  gel(q, 2)=subii(apbpc, gel(v, is[4]));//a+b+c-d
  gel(q, 3)=addii(gel(v, is[1]), gel(v, is[3]));//a+c
  return gerepileupto(top, q);
}

//Returns the reduction of v. If seq=1, also returns a VECSMALL of the sequence of indices swapped to reduce.
GEN apol_reduce(GEN v, int seq){
  pari_sp top=avma;
  long ind;
  GEN dold;
  if(!seq){
	do{
	  ind=ZV_maxind(v);
	  dold=gel(v, ind);
	  v=apol_move(v, ind);
	}while(cmpii(gel(v, ind), dold)==-1);
    return gerepileupto(top, apol_move(v, ind));//Must go back one!
  }
  //Now keep track of the sequence
  llist *S=NULL;
  long len=-1;//We don't count the last move.
  do{
	len++;
	ind=ZV_maxind(v);
    dold=gel(v, ind);
    v=apol_move(v, ind);
	llist_putstart(&S, ind);
  }while(cmpii(gel(v, ind), dold)==-1);
  llist_pop(&S);//Remove last move.
  return gerepilecopy(top, mkvec2(apol_move(v, ind), llist_tovecsmall(S, len, -1)));
}

//Returns the maximum index of the ZV v. Returns the first such index if a tie.
static long ZV_maxind(GEN v){
  pari_sp top=avma;
  long cmax=1;
  GEN max=gel(v, 1);
  for(long i=2;i<lg(v);i++) if(cmpii(gel(v, i), max)==1){cmax=i;max=gel(v, i);}
  avma=top;
  return cmax;
}

