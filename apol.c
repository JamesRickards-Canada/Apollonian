//(integral) Apollonian circle packing methods

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
static GEN apol_make_n(GEN q, GEN n, int red);

static GEN apol_search_bound(GEN v, GEN bound, int countsymm, GEN info, GEN (*getdata)(GEN, int, GEN, GEN*, int), GEN (*nextquad)(GEN, int, GEN), GEN (*retquad)(GEN));
static GEN apol_search_depth(GEN v, int depth, GEN bound, GEN info, GEN (*getdata)(GEN, int, GEN, GEN*, int), GEN (*nextquad)(GEN, int, GEN), GEN (*retquad)(GEN));
static GEN apol_circles_getdata(GEN vdat, int ind, GEN reps, GEN *nul, int state);
static GEN apol_thirdtangent(GEN circ1, GEN circ2, GEN c3, GEN c4, int right);
static GEN apol_circles_nextquad(GEN vdat, int ind, GEN reps);
static GEN apol_circles_retquad(GEN vdat);
static GEN apol_generic_nextquad(GEN vdat, int ind, GEN reps);
static GEN apol_generic_retquad(GEN vdat);
static GEN apol_curvatures_getdata(GEN vdat, int ind, GEN reps, GEN *nul, int state);
static GEN apol_find_getdata(GEN vdat, int ind, GEN reps, GEN *N, int state);



//BASIC METHODS


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

//Returns the external depth of v, i.e. the minimal number of swaps required to reach a quadruple with negative curvature.
long apol_extdepth(GEN v){
  pari_sp top=avma;
  long ind, step=0;
  for(;;){
    ind=vecindexmin(v);
    if(signe(gel(v, ind))!=1) return gc_long(top, step);//An index is <=0.
    step++;
    ind=vecindexmax(v);
    v=apol_move_1(v, ind);
  }
}

//Returns [S1, S2, S3, S4, K], where Si generate the Apollonian group, and K*[n,A,B,C]~=theta([A, B, C]) (see Staircase paper for theta description)
GEN apol_getmatrices(){
  pari_sp top=avma;
  GEN S1=mkmat4(mkcol4s(-1, 0, 0, 0), mkcol4s(2, 1, 0, 0), mkcol4s(2, 0, 1, 0), mkcol4s(2, 0, 0, 1));
  GEN S2=mkmat4(mkcol4s(1, 2, 0, 0), mkcol4s(0, -1, 0, 0), mkcol4s(0, 2, 1 ,0), mkcol4s(0, 2, 0, 1));
  GEN S3=mkmat4(mkcol4s(1, 0, 2, 0), mkcol4s(0, 1, 2, 0), mkcol4s(0, 0, -1, 0), mkcol4s(0, 0, 2, 1));
  GEN S4=mkmat4(mkcol4s(1, 0, 0, 2), mkcol4s(0, 1, 0, 2), mkcol4s(0, 0, 1, 2), mkcol4s(0, 0, 0, -1));
  GEN K=mkmat4(mkcol4s(1, -1, -1, -1), mkcol4s(0, 1, 0, 1), mkcol4s(0, 0, 0, -1), mkcol4s(0, 0, 1, 1));
  return gerepilecopy(top, mkvec5(S1, S2, S3, S4, K));
}

//Returns the possible obstructions modulo 24 of a primitive ACP, sorted lexicographically. They are stored in obstructions.dat
GEN apol_getobstructions(){
  GEN ret=cgetg(7, t_VEC);
  gel(ret, 1)=mkvecsmalln(6, 0L, 1L, 4L, 9L, 12L, 16L);
  gel(ret, 2)=mkvecsmalln(6, 0L, 4L, 12L, 13L, 16L, 21L);
  gel(ret, 3)=mkvecsmalln(6, 0L, 5L, 8L, 12L, 20L, 21L);
  gel(ret, 4)=mkvecsmalln(6, 0L, 8L, 9L, 12L, 17L, 20L);
  gel(ret, 5)=mkvecsmalln(8, 2L, 3L, 6L, 11L, 14L, 15L, 18L, 23L);
  gel(ret, 6)=mkvecsmalln(8, 3L, 6L, 7L, 10L, 15L, 18L, 19L, 22L);
  return ret;
}

/*Returns the set of admissible residues modulo 24. There are 6 possible primitive sets: 
[0, 1, 4, 9, 12, 16]; primes are 1 mod 24
[0, 4, 12, 13, 16, 21]; primes are 13 mod 24
[0, 5, 8, 12, 20, 21]; primes are 5 mod 24
[0, 8, 9, 12, 17, 20]; primes are 17 mod 24
[2, 3, 6, 11, 14, 15, 18, 23]; primes are 11, 23 mod 24
[3, 6, 7, 10, 15, 18, 19, 22]; primes are 7, 19 mod 24
You ONLY need to go to depth 3 to find which class we are in (proven by brute force check).*/
GEN apol_mod24(GEN v){
  pari_sp top=avma;
  long lv;
  GEN v24=cgetg_copy(v, &lv), tw4=stoi(24);//lv=5
  for(long i=1;i<lv;i++) gel(v24, i)=Fp_red(gel(v, i), tw4);//Reduce v modulo 24
  GEN orb=apol_curvatures_depth(v24, 3, gen_0);//Only need depth 3
  for(long i=1;i<lg(orb);i++) gel(orb, i)=Fp_red(gel(orb, i), tw4);//Reduce modulo 24.
  return gerepileupto(top, ZV_sort_uniq(orb));//Sort the result.
}

//Calls apol_move_1 or apol_move_batch, depending on if command is an integer or a vector. This is intended for use in gp. Use apl_move_1 or apol_move_batch with PARI.
GEN apol_move(GEN v, GEN command){
  long t=typ(command);
  switch(t){
    case t_INT: return apol_move_1(v, itos(command));//itos does not clutter the stack.
    case t_VECSMALL: return apol_move_batch(v, command);
  }
  pari_sp top=avma;
  command=gtovecsmall(command);
  return gerepileupto(top, apol_move_batch(v, command));
}

//Returns the set of four curvatures when we replace circle ind.
GEN apol_move_1(GEN v, int ind){
  if(ind<=0 || ind>=5) return v;//Do nothing.
  pari_sp top=avma;
  GEN rep=vecsmall_ei(4, ind);
  GEN S=gen_0;
  for(int i=1;i<=4;i++) if(!rep[i]) S=addii(S, gel(v, i));
  S=shifti(S, 1);//2(b+c+d)
  long lv;
  GEN newv=cgetg_copy(v, &lv);//lv=5
  for(int i=1;i<=4;i++) gel(newv, i)=rep[i]? subii(S, gel(v, ind)):icopy(gel(v, i));
  return gerepileupto(top, newv);
}

//Does apol_move_1 for v and bat[1], bat[2], ... The input bat needs to be a Vecsmall.
GEN apol_move_batch(GEN v, GEN bat){
  pari_sp top=avma;
  GEN newv=ZV_copy(v);
  for(long bind=1;bind<lg(bat);bind++){
    int ind=bat[bind];
    GEN S=gen_0;
    for(int i=1;i<=4;i++) if(i!=ind) S=addii(S, gel(newv, i));
    S=shifti(S, 1);//2(b+c+d)
    gel(newv, ind)=subii(S, gel(newv, ind));
  }
  return gerepilecopy(top, newv);
}

//Returns the quadratic form whose primitive values are v[ind]+curvatures touching the ind^th circle. The formula is [a+b,a+b+c-d, a+c] if v=[a,b,c,d] and ind=1.
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
GEN apol_red(GEN v, int seq){
  pari_sp top=avma;
  long ind;
  GEN dold;
  if(!seq){
    do{
      ind=vecindexmax(v);
      dold=gel(v, ind);
      v=apol_move_1(v, ind);
    }while(cmpii(gel(v, ind), dold)==-1);
    return gerepileupto(top, apol_move_1(v, ind));//Must go back one!
  }
  //Now keep track of the sequence
  long Sind=0, Slen=10;
  GEN S=cgetg(Slen+1, t_VECSMALL);//Initialize S
  do{
    ind=vecindexmax(v);
    dold=gel(v, ind);
    v=apol_move_1(v, ind);
    S=vecsmalllist_append(S, &Sind, &Slen, ind);
  }while(cmpii(gel(v, ind), dold)==-1);
  S=vecsmall_shorten(S, Sind-1);//Remove last move.
  return gerepilecopy(top, mkvec2(apol_move_1(v, ind), S));
}

//Reduces v, where we go ONLY at most maxsteps steps.
GEN apol_red_partial(GEN v, long maxsteps){
  pari_sp top=avma;
  if(maxsteps==0) return ZV_copy(v);
  long ind, step=0;
  int mstepsreached=1;
  GEN dold;
  do{
    step++;
    ind=vecindexmax(v);
    dold=gel(v, ind);
    v=apol_move_1(v, ind);
    if(cmpii(gel(v, ind), dold)!=-1){mstepsreached=0;break;}//We are reduced if we go back one.
  }while(step<maxsteps);
  if(mstepsreached) return gerepilecopy(top, v);
  else return gerepileupto(top, apol_move_1(v, ind));//Must go back one!
}



//CREATION OF ACPS


//Given a bqf q of discriminant -4n^2, this gives the corresponding Descartes quadruple. if pos=-1, we give the quadruple with -n, and if pos=1, we give the quadruple with +n. If red=1 we reduce, else we don't.
GEN apol_make(GEN q, int pos, int red){
  pari_sp top=avma;
  GEN D=bqf_disc(q);
  if(signe(D)!=-1) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  long rem;
  GEN nsqr=divis_rem(D, -4, &rem);
  if(rem!=0) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  GEN sqrtrem;
  GEN n=sqrtremi(nsqr, &sqrtrem);
  if(!gequal0(sqrtrem)) pari_err_TYPE("q must have discriminant -4n^2 for some integer n.", D);
  n= pos? n:negi(n);//+/-n
  return gerepileupto(top, apol_make_n(q, n, red));
}

//apol_make, but we are given the value of n, where disc(q)=-4n^2.
static GEN apol_make_n(GEN q, GEN n, int red){
  pari_sp top=avma;//q=[A, B, C]->[n, A-n, C-n, A+C-n-B]
  GEN Amn=subii(gel(q, 1), n);
  GEN v=mkvec4(n, Amn, subii(gel(q, 3), n), addii(Amn, subii(gel(q, 3), gel(q, 2))));//The ACP
  if(red) v=apol_red(v, 0);
  return gerepileupto(top, ZV_sort(v));
}

//Computes ncgp(-4n^2), and output the ACP's created from these quadratic forms with apol_make. We only count ambiguous forms ONCE, and we take pos=sign(n).
GEN apol_makeall(GEN n, int red, long prec){
  pari_sp top=avma;
  GEN D=mulis(sqri(n), -4);
  GEN U=bqf_ncgp_lexic(D, prec);
  GEN forms=gel(U, 3);
  long lf=lg(forms);
  GEN quads=vectrunc_init(lf);
  for(long i=1;i<lf;i++){//If we have [A, B, C] with B<0 we do not count it.
    GEN q=gel(forms, i);
    if(signe(gel(q, 2))==-1) continue;
    vectrunc_append(quads, apol_make_n(q, n, red));
  }
  return gerepileupto(top, quads);
}



//SEARCHING FOR CURVATURES


//Write one method for up to a bound, one method by depth, and one by layers?
//Allow input of a pointer to a function, as well as keeping track of some vector?


/*
The outline of the searching methods is:
	We use depth first search on the Apollonian graph.
	The current partial path is stored in the variable W. Each element of W is normally a Descartes quadruple, but there may be extra information.
	We have functions to:
		Retrieve the data we want from a new element (the curvature, entire quadruple, etc.), or returning NULL if we don't want to count it (e.g. only counting quadruples with certain properties)
		Move 1 layer deeper in the tree, and return the new element to be stored in W (again, normally a Descartes quadruple
		Retrieve the Descartes quadruple from the current element of W
	These functions are (respectively):
		static GEN getdata(GEN vdat, int ind, GEN reps, GEN *info, int state)
		static GEN nextquad(GEN vdat, int ind, GEN reps)
		static GEN retquad(GEN vdat)
		
	getdata:
		vdat is the newest node, where ind is the index of the circle we last swapped (stored in the format of v)
		reps passes in the current set of data found. This is used if W stores indices of elements in reps, and is also used at the end to clean up the data (we may want to sort the final answer, or do something else, etc.)
		info tracks extra information you may want, and may be modified by getdata. It is not touched in apol_search_(type). If not needed, pass it as NULL.
		state tracks what we are doing with this method.
			state=1 means we are adding a new piece of information potentially. Your method should compute the data, and return it if you want it added, and return NULL if we don't want to count it.
			state=2 is for the initial 4 circles, and you can do something special there if you wish. We always run the method on these first 4 circles in order atthe start.
			state=0 happens at the end, and is used as a wrap-up. You should return a modified version of reps that is gerepileupto compatible, e.g. sort it, just copy it, etc.
	
	nextquad:
		essentially does apol_move_1(vdat, ind), but formats it correctly in case vdat is not just the Descartes quadruple. Pass in ind=0 for the initial setup.
  
	retquad:
		Returns the Descartes quadruple from vdat. Normally this is just {return vdat;}, but if we keep track of more it is not.
*/


//Starting at the integral Descartes quadruple v, we search through the circle packing, finding all circles with curvature <=bound. For each such circle we compute some data (from getdata), and return the final set. If countsymm=0, we do not double count when there is symmetry, otherwise we do. For the strip packing, we force countsymm=0, since otherwise it would be infinite. It can be shown that all symmetries are realized for the reduced form (i.e. [a, b, c, d] -> [a, b, c, d] or [a, b, c, c], if this happens in a packing, it happens for the reduced form, so we can just check it there).
static GEN apol_search_bound(GEN v, GEN bound, int countsymm, GEN info, GEN (*getdata)(GEN, int, GEN, GEN*, int), GEN (*nextquad)(GEN, int, GEN), GEN (*retquad)(GEN)){
  pari_sp top=avma, mid;
  v=apol_red(v, 0);//Start by reducing v, since we want curvatures up to a given bound.
  ZV_sort_inplace(v);//May as well make the minimal curvature first.
  int ind=1;//We reuse ind to track which depth we are going towards.
  int depth=10;//Initially we go up to depth 10, but this may change if we need to go deeper.
  GEN W=zerovec(depth);//Tracks the sequence of Descartes quadruples; W[ind] is at ind-1 depth
  GEN vdat=nextquad(v, 0, NULL);//The first one. Stores the quadruple + data
  gel(W, 1)=vdat;
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of indices of the replacements
  int forward=1;//Tracks if we are going forward or not.
  long Nreps=400;//Initial size of the return.
  GEN reps=vectrunc_init(Nreps);//We may need to extend the length of reps later on.
  if(countsymm){
    if(gequal0(gel(v, 1))){
      countsymm=0;
      pari_warn(warner, "Strip packing, must not count symmetry");
    }
    else{
      for(int i=1;i<=4;i++){
        if(cmpii(gel(v, i), bound)<=0){//First 4 reps
          GEN dat=getdata(vdat, i, reps, &info, 2);
          if(dat) vectrunc_append(reps, dat);//We have data to store.
        }
      }
    }
  }
  if(!countsymm){//Must ignore symmetry. v is already sorted.
    for(int i=1;i<=4;i++){
      if(cmpii(gel(v, i), bound)<=0 && (i==1 || !equalii(gel(v, i-1), gel(v, i)))){//First 4 reps
        GEN dat=getdata(vdat, i, reps, &info, 2);
        if(dat) vectrunc_append(reps, dat);//We have data to store.
      }
    }
  }
  do{//1<=ind<=depth is assumed.
    if(gc_needed(top, 2)){//Garbage day!
      GEN oldreps=reps;
      mid=avma;
      W=gcopy(W);
      I=gcopy(I);
      info=gcopy(info);
      reps=vectrunc_init(Nreps);
      for(long i=1;i<lg(oldreps);i++) vectrunc_append(reps, gcopy(gel(oldreps, i)));//Copying reps; can't use vectrunc_append_batch.
      gerepileallsp(top, mid, 4, &W, &I, &info, &reps);
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
    if(ind==1 && !countsymm){//Make sure we haven't seen this element before. If this happens at some point, then it ALSO happens at the reduced form!
      int goback=0;
      for(int j=1;j<I[ind];j++){
        if(equalii(gel(v, j), gel(v, I[ind]))){goback=1;break;}
      }
      if(goback){forward=0;continue;}//We've already done this element in a previous move.
    }
    GEN newvdat=nextquad(gel(W, ind), I[ind], reps);//The data for the move.
	GEN newcurv=gel(retquad(newvdat), I[ind]);
    if(cmpii(newcurv, bound)>0){forward=0;continue;}//Must go back, elt too big
    if(ind==1 && !countsymm){//The replacement being the same can ONLY happen for the reduced form.
      if(equalii(gel(v, I[ind]), newcurv)){forward=0;continue;}//Must go back, same thing (e.g. the strip packing)
    }
    GEN dat=getdata(newvdat, I[ind], reps, &info, 1);//Retrieve the data
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
static GEN apol_search_depth(GEN v, int depth, GEN bound, GEN info, GEN (*getdata)(GEN, int, GEN, GEN*, int), GEN (*nextquad)(GEN, int, GEN), GEN (*retquad)(GEN)){
  pari_sp top=avma, mid;
  int ind=1;//We reuse ind to track which depth we are going towards.
  GEN W=zerovec(depth);//Tracks the sequence of Descartes quadruples; W[ind] is at ind-1 depth
  GEN vdat=nextquad(v, 0, NULL);//The first one. Stores the quadruple + data
  gel(W, 1)=vdat;
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1, usebound=1-gequal0(bound);;//Tracks if we are going forward or not, and whether or not to use the bound.
  long Nreps=400;
  GEN reps=vectrunc_init(Nreps);//We may need to extend the length of reps later on.
  for(int i=1;i<=4;i++){
    if(!usebound || cmpii(gel(v, i), bound)<=0){//First 4 reps
      GEN dat=getdata(vdat, i, reps, &info, 2);
      if(dat) vectrunc_append(reps, dat);//We have data to store.
    }
  }  
  do{//1<=ind<=depth is assumed.
    if(gc_needed(top, 2)){//Garbage day!
      GEN oldreps=reps;
      mid=avma;
      W=gcopy(W);
      I=gcopy(I);
      info=gcopy(info);
      reps=vectrunc_init(Nreps);
      for(long i=1;i<lg(oldreps);i++) vectrunc_append(reps, gcopy(gel(oldreps, i)));//Copying reps; can't use vectrunc_append_batch.
      gerepileallsp(top, mid, 4, &W, &I, &info, &reps);
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
    GEN newvdat=nextquad(gel(W, ind), I[ind], reps);//Make the move
    GEN newcurv=gel(retquad(newvdat), I[ind]);//The new element
    if(usebound && cmpii(newcurv, bound)>0){forward=0;continue;}//Must go back, new curvature too big
    GEN dat=getdata(newvdat, I[ind], reps, &info, 1);//Retrieve the data
    if(dat) vectrunc_append(reps, dat);//Add the new piece of data if we are allowed.
    ind++;
    if(ind>depth){forward=0;continue;}//Reached maximum depth, must go back
    gel(W, ind)=newvdat;
    forward=1;
  }while(ind>0);
  return gerepileupto(top, getdata(NULL, 0, reps, NULL, 0));
}


//vdat=[v, [indices in reps of the circles]]. We want to return the circle, i.e. [curvature, radius, x, y].
static GEN apol_circles_getdata(GEN vdat, int ind, GEN reps, GEN *nul, int state){
  if(state==0) return gcopy(reps);
  GEN v=gel(vdat, 1);//The new quadruple
  if(state==2){//Initial 4 circles
	switch(ind){
	  case 1:;
	    GEN c1=gel(v, 1);//First curvature
		return mkvec4(c1, Qdivii(gen_1, c1), gen_0, gen_0);//Outer circle
	  case 2:;
	    GEN c2=gel(v, 2);
		GEN r2=Qdivii(gen_1, c2);
		return mkvec4(c2, r2, gen_0, gneg(gadd(gmael(reps, 1, 2), r2)));//first inner circle, placed vertically at the top. Note gmael(reps, 1, 2)=r1<0.
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
  if(gequal(newcirc, prevcirc)) newcirc=apol_thirdtangent(oldcirc1, oldcirc2, newc, gel(v, is[2]), 0);//If the two curvatures were the same, this could trigger. If we lack oo precision, this could not work, and must be changed slightly.
  else{//This block must also be updated if there is not oo precision.
    GEN oldcirc3=gel(reps, prevind[is[2]]);//The unused old circle. Our newcirc must be tangent to it.
    GEN rsums=gsqr(gadd(gel(oldcirc3, 2), gel(newcirc, 2)));//(r1+r2)^2
    GEN dcentres=gadd(gsqr(gsub(gel(oldcirc3, 3), gel(newcirc, 3))), gsqr(gsub(gel(oldcirc3, 4), gel(newcirc, 4))));//dist(centres)^2
    if(!gequal(rsums, dcentres)) newcirc=apol_thirdtangent(oldcirc1, oldcirc2, newc, gel(v, is[2]), 0);//Must be the other side.
  }
  prevind[ind]=lg(reps);//Updating the location of the circle.
  return newcirc;
}

//Store a circle as [curvature, radius, x, y]. Given two tangent circles and a third curvature, this finds this third circle that is tangent to the first two. For internal tangency, we need a negative radius & curvature. There are always 2 places to put the circle: left or right of the line from circ1 to circ2. If right=1, we put it right, else we put it left. c4 is one of the curvatures to complete an Apollonian quadruple (supplying it allows us to always work with exact numbers in the case of integral ACPs).
static GEN apol_thirdtangent(GEN circ1, GEN circ2, GEN c3, GEN c4, int right){
  pari_sp top=avma;
  long prec=3;//Does not matter, things here are exact.
  //The centres form a triangle with sides r1+r2, r1+r3, r2+r3, or -r1-r2, -r1-r3, r2+r3 (if internal tangency). Let theta be the angle at the centre of c1.
  GEN c1=gel(circ1, 1), c2=gel(circ2, 1);//Curvatures
  GEN c1pc2=gadd(c1, c2), c1pc3=gadd(c1, c3);
  GEN denom=gmul(c1pc2, c1pc3);
  GEN costheta=gsubsg(1, gdiv(gmulsg(2, gsqr(c1)), denom));//1-2c1^2/((c1+c2)(c1+c3))
  GEN sintheta=gdiv(gabs(gmul(c1, gadd(c1pc2, gsub(c3, c4))), prec), denom);//|c1(c1+c2+c3-c4)|/((c1+c2)(c1+c3))
  GEN r1=gel(circ1, 2), r2=gel(circ2, 2), r3=gdivsg(1, c3);//The radii
  GEN x1=gel(circ1, 3), y1=gel(circ1, 4), x2=gel(circ2, 3), y2=gel(circ2, 4);
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
  return gerepilecopy(top, mkvec4(c3, r3, x, y));
}

//vdat=[v, [indices in reps of previous circles]], unless ind=0, when vdat=v
static GEN apol_circles_nextquad(GEN vdat, int ind, GEN reps){
  if(ind==0) return mkvec2(vdat, mkvecsmall4(1, 2, 3, 4));//The initial vdat.
  GEN v=gel(vdat, 1);//The actual quadruple
  GEN newv=apol_move_1(v, ind);//The new v
  return mkvec2(newv, zv_copy(gel(vdat, 2)));//We do NOT update the index in ind, which will update to lg(reps). This will be done in the getdata method, since we still need the old one for now.
}

//vdat=[v, [indices in reps of previous circles]]. This returns the first element of vdat
static GEN apol_circles_retquad(GEN vdat){return gel(vdat, 1);}

//Given a bounded integral Descartes quadruple, this returns equations for the circles with curvatures <=maxcurv at depth<=depth. An equation takes the form [curvature, radius, x, y]. The outer circle has centre (0, 0).
GEN apol_circles(GEN v, GEN maxcurv){
  return apol_search_bound(v, maxcurv, 1, NULL, &apol_circles_getdata, &apol_circles_nextquad, &apol_circles_retquad);
}

//vdat=v=Descartes quadruple. This returns the next one
static GEN apol_generic_nextquad(GEN vdat, int ind, GEN reps){return apol_move_1(vdat, ind);}

//vdat=v=Descartes quadruple. This returns it
static GEN apol_generic_retquad(GEN vdat){return vdat;}

//Helper method for apol_curvatures, to feed into apol_search_bound. vdat=v=Descartes quadruple
static GEN apol_curvatures_getdata(GEN vdat, int ind, GEN reps, GEN *nul, int state){
  if(state==0) return ZV_sort(reps);//Sort it and return.
  return gel(vdat, ind);//state=1/2, we want the new curvature
}

//Returns the curvatures in the packing v up to bound.
GEN apol_curvatures(GEN v, GEN bound, int countsymm){
  return apol_search_bound(v, bound, countsymm, NULL, &apol_curvatures_getdata, &apol_generic_nextquad, &apol_generic_retquad);
}

//Returns the curvatures up to depth depth from v. If bound>0, we only count those at most bound.
GEN apol_curvatures_depth(GEN v, int depth, GEN bound){
  if(depth<=0) depth=1;//To avoid errors from bad input.
  return apol_search_depth(v, depth, bound, NULL, &apol_curvatures_getdata, &apol_generic_nextquad, &apol_generic_retquad);
}

//Helper method for apol_find, to feed into apol_search_bound.
static GEN apol_find_getdata(GEN vdat, int ind, GEN reps, GEN *N, int state){
  if(state==0) return gcopy(reps);//Nothing to do.
  if(equalii(gel(vdat, ind), *N)) return vdat;//We have found N!
  return NULL;//This was not N, do not return anything.
}

//Searches for all circles with curvature N, and returns the corresponding quadruples. If countsymm=1, we may have repeats coming from the symmetry.
GEN apol_find(GEN v, GEN N, int countsymm){
  return apol_search_bound(v, N, countsymm, N, &apol_find_getdata, &apol_generic_nextquad, &apol_generic_retquad);
}





//Returns a sorted list of curvatures of circles that are maxlayers layers in from v[1]. Thus maxlayers=1 means that we only consider circles tangent to v[1] (along with v[1] itself). We only retrieve curvatures up to the given bound.
GEN apol_orbit_layers(GEN v, int maxlayers, GEN bound){
  pari_sp top=avma, mid;
  int ind=1;//We reuse ind to track which depth we are going towards.
  int depth=10;//Initially we go up to depth 10, but this may change if we need to go deeper.
  GEN W=zerovec(depth);//Tracks the sequence of APC's; W[ind] is at ind-1 depth
  GEN Wlayers=zerovec(depth);//Tracks the layers.
  gel(W, 1)=v;//The first one.
  gel(Wlayers, 1)=mkvecsmalln(6, 0L, 1L, 1L, 1L, 0L, 1L);//First circles are in layers 0, 1, 1, 1; min 0 1 time (format is layer 1, layer 2, layer 3, layer 4, min, freq of min.
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1;//Tracks if we are going forward or not.
  long Nreps=itos(bound);
  GEN reps=vectrunc_init(Nreps);//We may need to extend the length of reps later on.
  for(int i=1;i<=4;i++) if(cmpii(gel(v, i), bound)<=0) vectrunc_append(reps, gel(v, i));//First 4 reps
  do{//1<=ind<=depth is assumed.
    if(gc_needed(top, 2)){//Garbage day!
      GEN oldreps=reps;
      mid=avma;
      W=gcopy(W);
      Wlayers=gcopy(Wlayers);
      I=gcopy(I);
      reps=vectrunc_init(Nreps);
      for(long i=1;i<lg(oldreps);i++) vectrunc_append(reps, gcopy(gel(oldreps, i)));//Copying reps
      gerepileallsp(top, mid, 4, &W, &Wlayers, &I, &reps);
    }
    if(lg(reps)==Nreps){//We don't have enough space! Double the possible length of reps.
      Nreps=2*Nreps-1;
      GEN newreps=vectrunc_init(Nreps);
      vectrunc_append_batch(newreps, reps);//Append the old list
      reps=newreps;
    }
    I[ind]=forward? 1:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
    if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
    //At this point, we can go on with valid and new inputs. Start with checking the layer.
    GEN curlayer=gel(Wlayers, ind), nextlayer=vecsmall_copy(curlayer);
    if(curlayer[I[ind]]==curlayer[5]){//Currently min. If max, it MUST be repeated, and does not change anything.
      if(curlayer[6]==1){nextlayer[I[ind]]=nextlayer[I[ind]]+2;nextlayer[5]++;nextlayer[6]=3;}//go from x, x+1, x+1, x+1 to x+2, x+1, x+1, x+1.
      else{nextlayer[I[ind]]++;nextlayer[6]--;}//Repeated minimum just moves up 1.
    }
    if(nextlayer[I[ind]]>maxlayers){forward=0;continue;}//Too many layers deep.
    GEN newv=apol_move_1(gel(W, ind), I[ind]);//Make the move
    GEN elt=gel(newv, I[ind]);//The new element
    if(cmpii(elt, bound)==1 || ZV_equal(gel(W, ind), newv)) forward=0;//Must go back, elt too big OR same thing (e.g. happens with strip packing)
    else{
      vectrunc_append(reps, elt);//Add the new element
      ind++;
      if(ind==depth){
        int newdepth=2*depth;
        GEN newW=zerovec(newdepth), newWlayers=zerovec(newdepth), newI=vecsmall_ei(newdepth, 1);
        for(long i=1;i<=depth;i++){gel(newW, i)=gel(W, i);gel(newWlayers, i)=gel(Wlayers, i);newI[i]=I[i];}//Copy them over.
        W=newW;
        Wlayers=newWlayers;
        I=newI;
        depth=newdepth;
      }
      gel(W, ind)=newv;
      gel(Wlayers, ind)=nextlayer;
      forward=1;
    }
  }while(ind>0);
  return gerepileupto(top, ZV_sort(reps));
}

//Does apol_orbit_layers, but only retrieves the primes that are found (no repeats), NOT counting 2 or 3.
GEN apol_orbit_primes(GEN v, int maxlayers, GEN bound){
  pari_sp top=avma, mid;
  int ind=1;//We reuse ind to track which depth we are going towards.
  int depth=10;//Initially we go up to depth 10, but this may change if we need to go deeper.
  GEN W=zerovec(depth);//Tracks the sequence of APC's; W[ind] is at ind-1 depth
  GEN Wlayers=zerovec(depth);//Tracks the layers.
  gel(W, 1)=v;//The first one.
  gel(Wlayers, 1)=mkvecsmalln(6, 0L, 1L, 1L, 1L, 0L, 1L);//First circles are in layers 0, 1, 1, 1; min 0 1 time (format is layer 1, layer 2, layer 3, layer 4, min, freq of min.
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1;//Tracks if we are going forward or not.
  GEN m24=apol_mod24(v);
  GEN rposs;//Residue classes possible.
  switch(itos(gel(m24, 2))){//Uniquely determines the residue class
    case 1:
      rposs=mkvecsmall(1);break;
    case 4:
      rposs=mkvecsmall(13);break;
    case 5:
      rposs=mkvecsmall(5);break;
    case 8:
      rposs=mkvecsmall(17);break;
    case 3:
      rposs=mkvecsmall2(11, 23);break;
    default://Case 6
      rposs=mkvecsmall2(7, 19);
  }
  
  long Nreps=itos(bound), lgrposs=lg(rposs), rem;
  GEN reps=vectrunc_init(Nreps);//We may need to extend the length of reps later on.
  for(int i=1;i<=4;i++){
    if(cmpii(gel(v, i), bound)<=0){//First 4 reps
      rem=smodis(gel(v, i), 24);
      if((rem==rposs[1] || (lgrposs==3 && rem==rposs[2])) && isprime(gel(v, i))) vectrunc_append(reps, gel(v, i));
    }
  }
  do{//1<=ind<=depth is assumed.
    if(gc_needed(top, 2)){//Garbage day!
      GEN oldreps=reps;
      mid=avma;
      W=gcopy(W);
      Wlayers=gcopy(Wlayers);
      I=gcopy(I);
      reps=vectrunc_init(Nreps);
      for(long i=1;i<lg(oldreps);i++) vectrunc_append(reps, gcopy(gel(oldreps, i)));//Copying reps
      gerepileallsp(top, mid, 4, &W, &Wlayers, &I, &reps);
    }
    if(lg(reps)==Nreps){//We don't have enough space! Double the possible length of reps.
      Nreps=2*Nreps-1;
      GEN newreps=vectrunc_init(Nreps);
      vectrunc_append_batch(newreps, reps);//Append the old list
      reps=newreps;
    }
    I[ind]=forward? 1:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
    if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
    //At this point, we can go on with valid and new inputs. Start with checking the layer.
    GEN curlayer=gel(Wlayers, ind), nextlayer=vecsmall_copy(curlayer);
    if(curlayer[I[ind]]==curlayer[5]){//Currently min. If max, it MUST be repeated, and does not change anything.
      if(curlayer[6]==1){nextlayer[I[ind]]=nextlayer[I[ind]]+2;nextlayer[5]++;nextlayer[6]=3;}//go from x, x+1, x+1, x+1 to x+2, x+1, x+1, x+1.
      else{nextlayer[I[ind]]++;nextlayer[6]--;}//Repeated minimum just moves up 1.
    }
    if(nextlayer[I[ind]]>maxlayers){forward=0;continue;}//Too many layers deep.
    GEN newv=apol_move_1(gel(W, ind), I[ind]);//Make the move
    GEN elt=gel(newv, I[ind]);//The new element
    if(cmpii(elt, bound)==1 || ZV_equal(gel(W, ind), newv)) forward=0;//Must go back, elt too big OR same thing (e.g. happens with strip packing)
    else{
      rem=smodis(elt, 24);
      if((rem==rposs[1] || (lgrposs==3 && rem==rposs[2])) && isprime(elt)) vectrunc_append(reps, elt);//Add the new element
      ind++;
      if(ind==depth){
        int newdepth=2*depth;
        GEN newW=zerovec(newdepth), newWlayers=zerovec(newdepth), newI=vecsmall_ei(newdepth, 1);
        for(long i=1;i<=depth;i++){gel(newW, i)=gel(W, i);gel(newWlayers, i)=gel(Wlayers, i);newI[i]=I[i];}//Copy them over.
        W=newW;
        Wlayers=newWlayers;
        I=newI;
        depth=newdepth;
      }
      gel(W, ind)=newv;
      gel(Wlayers, ind)=nextlayer;
      forward=1;
    }
  }while(ind>0);
  return gerepileupto(top, ZV_sort_uniq(reps));
}




//STRIP PACKING METHODS


//Returns [curvature, r, a, b], where the depth pairing corresponding to L is given by the circle (x-a)^2+(y-b)^2=r^2. If L is an integer, this corresponds to (Id, L). If L is a vecsmall/vector, this corresponds to (S_L[1]*...*S_L[n], L[1]). If a=r=oo, this corresponds to the line y=b.
GEN apol_dpair_circle(GEN L){
  pari_sp top=avma;
  int t=typ(L);
  switch(t){
    case t_INT:;//Required to stop error with int sL
      int sL=itos(L);
      switch(sL){
        case 1:
          return gerepilecopy(top, mkvec4(gen_0, mkoo(), mkoo(), gen_0));
        case 2:
          return gerepilecopy(top, mkvec4(gen_0, mkoo(), mkoo(), gen_1));
        case 3:
          return mkvec4(gen_2, ghalf, gen_0, ghalf);
        default://i.e. case 4
          return mkvec4(gen_2, ghalf, gen_m1, ghalf);
      }
    case t_VEC:
      L=gtovecsmall(L);//Make it a vecsmall
    case t_VECSMALL:
      break;
    default:
      pari_err_TYPE("L needs to be an integer from 1 to 4, or a vecsmall/vector of such integers", L);
  }
  GEN M=apol_getmatrices();
  GEN W=matid(4);
  for(long i=1;i<lg(L);i++) W=ZM_mul(W, gel(M, L[i]));
  W=ZM_mul(W, gel(M, 5));//Times k at the end
  GEN mtuvw=row(W, L[1]);//[-t, u, v, w]
  GEN twow=shifti(gel(mtuvw, 4), 1);//2w=curvature
  GEN a=Qdivii(gel(mtuvw, 3), gel(mtuvw, 4));
  GEN b=Qdivii(negi(gel(mtuvw, 1)), twow);
  GEN r=Qdivii(gen_1, twow);
  return gerepilecopy(top, mkvec4(twow, r, a, b));
}

//Returns the quadratic form corresponding to the circle in the stip packing designated by L. If L is an integer, this corresponds to Id_L. If L is a vecsmall/vector, this corresponds to S_L[1]*...*S_L[n]. We can't have L=1 or 2, this doesn't give a circle.
GEN apol_strip_qf(GEN L, int red){
  pari_sp top=avma;
  GEN c=apol_dpair_circle(L);//The corresponding circle, but it's scaled by 1/2, so need to scale it by 2.
  if(typ(gel(c, 2))==t_INFINITY) pari_err_TYPE("Please give a circle instead, i.e. don't input L=1 or 2", L);
  GEN C=shifti(gel(c, 1), -1);//C=curvature
  GEN B=gmulgs(gmul(gel(c, 1), gel(c, 3)), 2);//B=2*curvature*Real(centre)
  GEN A=gsub(gmul(gadd(gsqr(gel(c, 3)), gsqr(gel(c, 4))), shifti(C, 2)), gmulgs(gel(c, 2), 2));//A=cocurvature=|centre|^2*curvature-radius
  GEN q=mkvec3(A, B, C);
  if(red) return gerepileupto(top, dbqf_red(q));
  return gerepilecopy(top, q);
}



//VISUALIZATION


//Given a list of circles, this prints them to the screen in a format suitable for Desmos.
void printcircles_desmos(GEN c){
  for(long i=1;i<lg(c);i++){
    if(gequal0(gmael(c, i, 1))) pari_printf("y=%Ps\n", gmael(c, i, 4));//Horizontal line
    else pari_printf("(x-%Ps)^2+(y-%Ps)^2=1/(%Ps)^2\n", gmael(c, i, 3), gmael(c, i, 4), gmael(c, i, 1));
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
  if(gcmpgs(gmael(c, 1, 2), 0)<0){//First circle neg curvature, so everything is a circle.
    GEN scalingfactor=gdiv(largestcirc, gneg(gmael(c, 1, 2)));//Scaling factor.
    gel(cscale, 1)=mkvec3(largestcirc, gen_0, gen_0);
    for(long i=2;i<lc;i++){
      gel(cscale, i)=mkvec3(gmul(gmael(c, i, 2), scalingfactor), gmul(gmael(c, i, 3), scalingfactor), gmul(gmael(c, i, 4), scalingfactor));//r, x, y
    }//Circles have been scaled!
  }
  else{//Strip packing, OR the largest curvature does not come first. The width should be at least as much as the height.
    GEN minx=mkoo(), maxx=mkmoo();
    for(long i=1;i<lc;i++){
      if(gequal0(gmael(c, i, 1))) continue;//Horizontal line.
      GEN circxmin=gsub(gmael(c, i, 3), gmael(c, i, 2));
      GEN circxmax=gadd(gmael(c, i, 3), gmael(c, i, 2));
      minx=gmin_shallow(minx, circxmin);
      maxx=gmax_shallow(maxx, circxmax);
    }
    if(typ(minx)==t_INFINITY || typ(maxx)==t_INFINITY) pari_err_TYPE("You don't have any circles", c);
    GEN scalingfactor=gmulgs(gdiv(largestcirc, gsub(maxx, minx)), 2);
    GEN horizshift=gdivgs(gadd(maxx, minx), 2);
    for(long i=1;i<lc;i++){
      if(gequal0(gmael(c, i, 1))){//Horizontal line
        gel(cscale, i)=mkvec3(mkoo(), largestcirc, gmul(gmael(c, i, 4), scalingfactor));//oo, "radius", y-intercept ("radius"=r means going from [-r, y] to [r, y])
        continue;
      }
      gel(cscale, i)=mkvec3(gmul(gmael(c, i, 2), scalingfactor), gmul(gsub(gmael(c, i, 3), horizshift), scalingfactor), gmul(gmael(c, i, 4), scalingfactor));//r, x, y
    }
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
      if(typ(gmael(cscale, i, 1))==t_INFINITY) pari_fprintf(f, "  \\draw[ultra thin, fill=col%d] (%P.10fin, %P.10fin) -- (%P.10fin, %P.10fin);\n", smodis(gmael(c, i, 1), modcolours), gneg(gmael(cscale, i, 2)), gmael(cscale, i, 3), gmael(cscale, i, 2), gmael(cscale, i, 3));
      else{
        if(signe(gmael(c, i, 1))==-1) pari_fprintf(f, "  \\draw[ultra thin] (%P.10fin, %P.10fin) circle (%P.10fin);\n", gmael(cscale, i, 2), gmael(cscale, i, 3), gmael(cscale, i, 1));
        else pari_fprintf(f, "  \\draw[ultra thin, fill=col%d] (%P.10fin, %P.10fin) circle (%P.10fin);\n", smodis(gmael(c, i, 1), modcolours), gmael(cscale, i, 2), gmael(cscale, i, 3), gmael(cscale, i, 1));
      }
    }
  }
  GEN ten=stoi(10);
  if(addnumbers){
      char *options;
      if(modcolours>0) options=", white";
      else options="";
    for(long i=1;i<lc;i++){
      GEN curv=gmael(c, i, 1), scaleby;
      if(gsigne(curv)!=1) continue;//Don't add numbers for curvatures <=0.
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



//SPECIALIZED METHODS
//These are scripts that are useful to me, but not likely useful in general.



//Computes apol_makeall, and returns the external depths of the corresponding circles of curvature n. This works because we use the reduced quadratic forms: A-n<=C-n<=A+C-B-n. Swapping A+C-B-n gives A+C+B-n, so this cannot decrease it. Thus either we are already reduced (and A-n<=0 necessarily), or it is not the largest term (which thus must be n). Thus n is immediately swapped if depth>0, and this gives the circle depth.
GEN apol_ncgp_depths(GEN n, long prec){
  pari_sp top=avma;
  GEN forms=apol_makeall(n, 0, prec);//The forms
  long maxdepth=0, lf=lg(forms);
  GEN depths=cgetg(lf, t_VECSMALL);
  for(long i=1;i<lf;i++){//Computing depths
    long d=apol_extdepth(gel(forms, i));
    depths[i]=d;
    if(d>maxdepth) maxdepth=d;
  }
  GEN dcount=zero_zv(maxdepth+1);//Depths 0 through d.
  for(long i=1;i<lf;i++) dcount[depths[i]+1]++;//Counting.
  return gerepileupto(top, dcount);
}

//Computes apol_ncgpforms, and returns the sorted vector of smallest curvature for each example. We do not remove repeats, and output the negative of the curvatures (as they are all negative). If red=0, we actually just take the data as is, and output the smallest curvature in each of the quadruples (no negation).
GEN apol_ncgp_smallcurve(GEN n, int red, long prec){
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

//apol_ncgp_smallcurve, but we only reduce maxsteps steps.
GEN apol_ncgp_smallcurve_bsteps(GEN n, long maxsteps, long prec){
  pari_sp top=avma;
  GEN forms=apol_makeall(n, 0, prec);
  long lf;
  GEN curves=cgetg_copy(forms, &lf);
  for(long i=1;i<lf;i++){
    GEN q=apol_red_partial(gel(forms, i), maxsteps);
    gel(curves, i)=gel(q, vecindexmin(q));
  }
  return gerepileupto(top, ZV_sort(curves));
}


