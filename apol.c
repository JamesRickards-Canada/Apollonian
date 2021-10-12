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
GEN apol_1orbit(GEN v, int depth, int ind){
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

