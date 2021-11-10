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
static long ZV_minind(GEN v);

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

//Returns [a, b, r], where the depth pairing corresponding to L is given by the circle (x-a)^2+(y-b)^2=r^2. If L is an integer, this corresponds to (Id, L). If L is a vecsmall/vector, this corresponds to (S_L[1]*...*S_L[n], L[1]). If a=r=oo, this corresponds to the line y=b.
GEN apol_dpair_circle(GEN L){
  pari_sp top=avma;
  if(typ(L)==t_INT){
	int sL=itos(L);
	switch(sL){
	  case 1:
		return gerepilecopy(top, mkvec3(mkoo(), gen_0, mkoo()));
	  case 2:
		return gerepilecopy(top, mkvec3(mkoo(), gen_1, mkoo()));
	  case 3:
		return mkvec3(gen_0, ghalf, ghalf);
	  default://i.e. case 4
		return mkvec3(gen_m1, ghalf, ghalf);
	}
  }
  if(typ(L)==t_VEC) L=gtovecsmall(L);
  GEN M=apol_getmatrices();
  GEN W=matid(4);
  for(long i=1;i<lg(L);i++) W=ZM_mul(W, gel(M, L[i]));
  W=ZM_mul(W, gel(M, 5));//Times k at the end
  GEN mtuvw=row(W, L[1]);//[-t,u,v,w]
  GEN twow=shifti(gel(mtuvw, 4), 1);//2w
  GEN a=Qdivii(gel(mtuvw, 3), gel(mtuvw, 4));
  GEN b=Qdivii(negi(gel(mtuvw, 1)), twow);
  GEN r=Qdivii(gen_1, twow);
  return gerepilecopy(top, mkvec3(a, b, r));
}

//Returns [S1, S2, S3, S4, K], where Si generate the Apollonian group, and K*[n,A,B,C]~=theta([A, B, C]).
GEN apol_getmatrices(){
  pari_sp top=avma;
  GEN S1=mkmat4(mkcol4s(-1, 0, 0, 0), mkcol4s(2, 1, 0, 0), mkcol4s(2, 0, 1, 0), mkcol4s(2, 0, 0, 1));
  GEN S2=mkmat4(mkcol4s(1, 2, 0, 0), mkcol4s(0, -1, 0, 0), mkcol4s(0, 2, 1 ,0), mkcol4s(0, 2, 0, 1));
  GEN S3=mkmat4(mkcol4s(1, 0, 2, 0), mkcol4s(0, 1, 2, 0), mkcol4s(0, 0, -1, 0), mkcol4s(0, 0, 2, 1));
  GEN S4=mkmat4(mkcol4s(1, 0, 0, 2), mkcol4s(0, 1, 0, 2), mkcol4s(0, 0, 1, 2), mkcol4s(0, 0, 0, -1));
  GEN K=mkmat4(mkcol4s(1, -1, -1, -1), mkcol4s(0, 1, 0, 1), mkcol4s(0, 0, 0, -1), mkcol4s(0, 0, 1, 1));
  return gerepilecopy(top, mkvec5(S1, S2, S3, S4, K));
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
	if(red==2) vectrunc_append(apols, apol_red(v, 0));
	else vectrunc_append(apols, v);
  }
  return gerepilecopy(top, apols);
}

//Given a bqf q, this gives the corresponding root quadruple. if pos=-1, we give the quadruple with -n, and if pos=1, we give the quadruple with +n. If red=1 we reduce, else we don't.
GEN apol_make_fromqf(GEN q, int pos, int red){
  pari_sp top=avma;
  GEN D=bqf_disc(q);
  if(signe(D)!=-1) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  long rem;
  GEN nsqr=divis_rem(D, -4, &rem);
  if(rem!=0) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  GEN sqrtrem;
  GEN n=sqrtremi(nsqr, &sqrtrem);
  if(!gequal0(sqrtrem)) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  //Now D=-4n^2, a=+/-n (- if pos=0), and q=[A, B, C]->[a, A-a, C-a, A+C-a-B]
  GEN a= pos? n:negi(n);//+/-n
  GEN Ama=subii(gel(q, 1), a);
  GEN v=mkvec4(a, Ama, subii(gel(q, 3), a), addii(Ama, subii(gel(q, 3), gel(q, 2))));//The APC
  if(red) v=apol_red(v, 0);
  return gerepileupto(top, ZV_sort(v));
}


//Returns the set of admissible residues modulo 24. There are 38 possible primitive sets: [[0, 1, 3, 4, 6, 9, 10, 12, 16, 18, 19, 22], [0, 1, 4, 6, 7, 9, 10, 12, 15, 16, 18, 22], [0, 1, 4, 9, 12, 13, 16, 21], [0, 1, 4, 9, 12, 16], [0, 2, 3, 5, 6, 8, 11, 12, 14, 18, 20, 21], [0, 2, 3, 6, 8, 9, 11, 12, 14, 17, 18, 20], [0, 2, 5, 6, 8, 12, 14, 15, 18, 20, 21, 23], [0, 2, 6, 8, 9, 12, 14, 15, 17, 18, 20, 23], [0, 3, 4, 6, 10, 12, 13, 16, 18, 19, 21, 22], [0, 3, 4, 7, 12, 15, 16, 19], [0, 3, 4, 12, 16, 19], [0, 3, 8, 11, 12, 15, 20, 23], [0, 3, 8, 11, 12, 20], [0, 4, 6, 7, 10, 12, 13, 15, 16, 18, 21, 22], [0, 4, 7, 12, 15, 16], [0, 4, 12, 13, 16, 21], [0, 5, 8, 9, 12, 17, 20, 21], [0, 5, 8, 12, 20, 21], [0, 8, 9, 12, 17, 20], [0, 8, 12, 15, 20, 23], [1, 3, 7, 9, 13, 15, 19, 21], [1, 6, 9, 10, 13, 18, 21, 22], [1, 6, 9, 10, 18, 22], [1, 9, 13, 21], [2, 3, 6, 11, 14, 15, 18, 23], [2, 3, 6, 11, 14, 18], [2, 5, 6, 9, 14, 17, 18, 21], [2, 5, 6, 14, 18, 21], [2, 6, 9, 14, 17, 18], [2, 6, 14, 15, 18, 23], [3, 5, 9, 11, 15, 17, 21, 23], [3, 6, 7, 10, 15, 18, 19, 22], [3, 6, 10, 18, 19, 22], [3, 7, 15, 19], [3, 11, 15, 23], [5, 9, 17, 21], [6, 7, 10, 15, 18, 22], [6, 10, 13, 18, 21, 22]]
//The lengths are 4 (4x), 6 (16x), 8 (10x), and 12 (8x). You ONLY need to go to depth 3 to find which class we are (proven by brute force check).
GEN apol_mod24(GEN v){
  pari_sp top=avma;
  long lv;
  GEN v24=cgetg_copy(v, &lv), tw4=stoi(24);//lv=5
  for(long i=1;i<lv;i++) gel(v24, i)=Fp_red(gel(v, i), tw4);
  GEN orb=apol_orbit(v24, 3, gen_0);//Only need depth 3
  for(long i=1;i<lg(orb);i++) gel(orb, i)=Fp_red(gel(orb, i), tw4);//Reduce modulo 24.
  return gerepileupto(top, ZV_sort_uniq(orb));//Sort the result.
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

//Computes apol_ncgpforms, and returns the depths of the corresponding circles of curvature n.
GEN apol_ncgp_depths(GEN n, long prec){
  pari_sp top=avma;
  GEN forms=apol_ncgp_forms(n, 1, 0, prec);//The forms
  long maxdepth=0, lf=lg(forms);
  GEN depths=cgetg(lf, t_VECSMALL);
  for(long i=1;i<lf;i++){//Computing depths
	long d=apol_quaddepth(gel(forms, i));
	depths[i]=d;
	if(d>maxdepth) maxdepth=d;
  }
  GEN dcount=zero_zv(maxdepth+1);//Depths 0 through d.
  for(long i=1;i<lf;i++) dcount[depths[i]+1]++;//Counting.
  return gerepileupto(top, dcount);
}

//Computes ncgp(-4n^2), and output the ACP's created from these quadratic forms with apol_make_qf. We only count ambiguous forms ONCE.
GEN apol_ncgp_forms(GEN n, int pos, int red, long prec){
  pari_sp top=avma;
  GEN D=mulis(sqri(n), -4);
  GEN U=bqf_ncgp_lexic(D, prec);
  GEN forms=gel(U, 3);
  long lf=lg(forms);
  GEN quads=vectrunc_init(lf);
  for(long i=1;i<lf;i++){//If we have [A, B, C] with B<0 we do not count it.
    GEN q=gel(forms, i);
	if(signe(gel(q, 2))==-1) continue;
	vectrunc_append(quads, apol_make_fromqf(q, pos, red));
  }
  return gerepileupto(top, quads);
}

//Computes apol_ncgpforms, and returns the sorted vector of smallest curvature for each example. We do not remove repeats, and output the negative of the curvatures (as they are all negative). If red=0, we actually just take the data as is, and output the smallest curvature in each of the quadruples (no negation).
GEN apol_ncgp_smallcurve(GEN n, int red, long prec){
  pari_sp top=avma;
  GEN forms=apol_ncgp_forms(n, 1, red, prec);
  long lf;
  GEN curves=cgetg_copy(forms, &lf);
  for(long i=1;i<lf;i++){
	gel(curves, i)=gmael(forms, i, 1);
	if(red) togglesign_safe(&gel(curves, i));
  }
  return gerepileupto(top, ZV_sort(curves));
}

//apol_ncgp_smallcurve, but we only reduce maxsteps steps.
GEN apol_ncgp_smallcurve_bsteps(GEN n, long maxsteps, long prec){
  pari_sp top=avma;
  GEN forms=apol_ncgp_forms(n, 1, 0, prec);
  long lf;
  GEN curves=cgetg_copy(forms, &lf);
  for(long i=1;i<lf;i++){
	GEN q=apol_red_bsteps(gel(forms, i), maxsteps);
	gel(curves, i)=gel(q, ZV_minind(q));
  }
  return gerepileupto(top, ZV_sort(curves));
}

//Returns a sorted list of curvatures of circles, where we go to depth depth, i.e. we do up to depth circle replacements. We also only retrieve curvatures <=bound, if this is passed in as non-zero.
GEN apol_orbit(GEN v, int depth, GEN bound){
  pari_sp top=avma;
  int ind=1;//We reuse ind to track which depth we are going towards.
  GEN W=zerovec(depth);//Tracks the sequence of APC's; W[ind] is at ind-1 depth
  gel(W, 1)=v;//The first one.
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1, usebound=1-gequal0(bound);//Tracks if we are going forward or not.
  long Nreps;
  if(gequal0(bound)) Nreps=2*(itos(powuu(3, depth))+1)+1;
  else Nreps=itos(bound);//If bound!=0, and depth is large, the previous Nreps definition may be too large
  GEN reps=vectrunc_init(Nreps);
  if(usebound){
	for(int i=1;i<=4;i++) if(cmpii(gel(v, i), bound)<=0) vectrunc_append(reps, gel(v, i));//First 4 reps
  }
  else{
    for(int i=1;i<=4;i++) vectrunc_append(reps, gel(v, i));//First 4 reps
  }
  do{//1<=ind<=depth is assumed.
    if(lg(reps)==Nreps){//We don't have enough space! Double the possible length of reps.
	  long newNreps=2*Nreps-1;
	  GEN newreps=vectrunc_init(newNreps);
	  for(long i=1;i<Nreps;i++) vectrunc_append(newreps, gel(reps, i));//Append the old list
	  Nreps=newNreps;
	  reps=newreps;
	}
    I[ind]=forward? 1:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
	if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
	//At this point, we can go on with valid and new inputs
	GEN newv=apol_move(gel(W, ind), I[ind]);//Make the move
	GEN elt=gel(newv, I[ind]);//The new element
	if(usebound && cmpii(elt, bound)==1) forward=0;//Must go back, elt too big
	else{
	  vectrunc_append(reps, elt);//Add the new element
	  if(ind==depth) forward=0;
	  else{//We can keep going forward
        ind++;
	    gel(W, ind)=newv;
	    forward=1;
	  }
    }
  }while(ind>0);
  return gerepileupto(top, ZV_sort(reps));
}

//Returns a sorted list of curvatures of circles surrounding v[ind]. We go to depth depth, i.e. we do up to depth circle replacements. We also only retrieve curvatures <=bound, if this is passed in as non-zero.
GEN apol_orbit_1(GEN v, int ind, int depth, GEN bound){
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
  int forward=1, usebound=1-gequal0(bound);//Tracks if we are going forward or not.
  long Nreps;
  if(gequal0(bound)) Nreps=3*itos(int2n(depth))+1;
  else Nreps=itos(bound);//If bound!=0, and depth is large, the previous Nreps definition may be too large
  GEN reps=vectrunc_init(Nreps);
  if(usebound){
	for(int i=2;i<=4;i++) if(cmpii(gel(v1, i), bound)<=0) vectrunc_append(reps, gel(v1, i));//First 3 reps
  }
  else{
    for(int i=2;i<=4;i++) vectrunc_append(reps, gel(v1, i));//First 3 reps
  }
  do{//1<=ind<=depth is assumed.
	if(lg(reps)==Nreps){//We don't have enough space! Double the possible length of reps.
	  long newNreps=2*Nreps-1;
	  GEN newreps=vectrunc_init(newNreps);
	  for(long i=1;i<Nreps;i++) vectrunc_append(newreps, gel(reps, i));//Append the old list
	  Nreps=newNreps;
	  reps=newreps;
	}
    I[ind]=forward? 2:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
	if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
	//At this point, we can go on with valid and new inputs
	GEN newv=apol_move(gel(W, ind), I[ind]);//Make the move
	GEN elt=gel(newv, I[ind]);//The new element
	if(usebound && cmpii(elt, bound)==1) forward=0;//Must go back, elt too big
	else{
	  vectrunc_append(reps, elt);//Add the new element
	  if(ind==depth) forward=0;
	  else{//We can keep going forward
        ind++;
	    gel(W, ind)=newv;
	    forward=1;
	  }
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

//Returns the depth of v, i.e. the minimal number of swaps requried to reach a quadruple with negative curvature.
long apol_quaddepth(GEN v){
  pari_sp top=avma;
  long ind, step=0;
  ind=ZV_minind(v);
  if(signe(gel(v, ind))!=1) return gc_long(top, step);//Start <0
  for(;;){
    step++;
    ind=ZV_maxind(v);
    v=apol_move(v, ind);
	ind=ZV_minind(v);
	if(signe(gel(v, ind))!=1) return gc_long(top, step);//Start <0
  }
}

//Returns the reduction of v. If seq=1, also returns a VECSMALL of the sequence of indices swapped to reduce.
GEN apol_red(GEN v, int seq){
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

//Reduces v, where we ONLY at most maxsteps steps.
GEN apol_red_bsteps(GEN v, long maxsteps){
  pari_sp top=avma;
  if(maxsteps==0) return ZV_copy(v);
  long ind, step=0;
  int mstepsreached=1;
  GEN dold;
  do{
    step++;
    ind=ZV_maxind(v);
    dold=gel(v, ind);
    v=apol_move(v, ind);
	if(cmpii(gel(v, ind), dold)!=-1){mstepsreached=0;break;}//We are reduced if we go back one.
  }while(step<maxsteps);
  if(mstepsreached) return gerepilecopy(top, v);
  else return gerepileupto(top, apol_move(v, ind));//Must go back one!
}

//Search for circles of curvature N up to depth depth. Returns the corresponding quadruples. Adops the code of apol_orbit. Returns the qf's if rqf=1, and both if rqf=2 (each entry is [ACP, qf]).
GEN apol_search(GEN v, GEN N, int depth, int rqf){
  pari_sp top=avma;
  glist *S=NULL;//Stores the ACP's.
  int ind=1;//We reuse ind to track which depth we are going towards.
  GEN W=zerovec(depth);//Tracks the sequence of APC's; W[ind] is at ind-1 depth
  gel(W, 1)=v;//The first one.
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1;//Tracks if we are going forward or not.
  long Nfound=0;
  for(int i=1;i<=4;i++){
	if(equalii(N, gel(v, i))){//Found one!
      if(rqf==0) glist_putstart(&S, gcopy(v));
	  else{
		GEN q=apol_qf(v, i);
		if(rqf==1) glist_putstart(&S, q);
		else glist_putstart(&S, mkvec2(gcopy(v), q));
	  }
	  Nfound++;
	  break;
	}
  }
  do{//1<=ind<=depth is assumed.
    I[ind]=forward? 1:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
	if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
	//At this point, we can go on with valid and new inputs
	GEN newv=apol_move(gel(W, ind), I[ind]);//Make the move
	int comp=cmpii(N, gel(newv, I[ind]));//Comparing the new element to N.
	if(comp==0){//Found it!
	  if(rqf==0) glist_putstart(&S, gcopy(newv));
	  else{
		GEN q=apol_qf(newv, I[ind]);
		if(rqf==1) glist_putstart(&S, q);
		else glist_putstart(&S, mkvec2(newv, q));
	  }
	  Nfound++;
	}
	if(ind==depth || comp<=0) forward=0;//Max depth OR the number is too big; once we reach or pass N, we cannot get N anymore.
	else{//We can keep going forward
      ind++;
	  gel(W, ind)=newv;
	  forward=1;
	}
  }while(ind>0);
  return gerepileupto(top, glist_togvec(S, Nfound, -1));
}

//Given a sorted ZV, counts how many entries are non-positive.
long ZV_countnonpos(GEN v){
  long i1=1, i2=lg(v)-1, i=0;
  if(signe(gel(v, i1))==1) return 0;//None <=0
  if(signe(gel(v, i2))!=1) return i2;//All <=0
  while(i1+1<i2){
	i=(i1+i2)/2;//Floor
	if(signe(gel(v, i))==1) i2=i;//v[i]>0
	else i1=i;//v[i]<=0
  }
  return i1;
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

//Returns the minimal index of the ZV v. Returns the first such index if a tie.
static long ZV_minind(GEN v){
  pari_sp top=avma;
  long cmin=1;
  GEN min=gel(v, 1);
  for(long i=2;i<lg(v);i++) if(cmpii(gel(v, i), min)==-1){cmin=i;min=gel(v, i);}
  avma=top;
  return cmin;
}

