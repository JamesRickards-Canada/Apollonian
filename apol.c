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
static GEN thirdtangent(GEN circ1, GEN circ2, GEN c3, GEN c4, int right, long prec);



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

//Returns the depth of v, i.e. the minimal number of swaps required to reach a quadruple with negative curvature.
long apol_depth(GEN v){
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
  GEN orb=apol_orbit(v24, 3, gen_0);//Only need depth 3
  for(long i=1;i<lg(orb);i++) gel(orb, i)=Fp_red(gel(orb, i), tw4);//Reduce modulo 24.
  return gerepileupto(top, ZV_sort_uniq(orb));//Sort the result.
}

//Calls apol_move_1 or apol_move_batch, depending on if command is an integer or a vector. This is intended for use in gp.
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

//Returns the set of four curvatures when we replace circle i.
GEN apol_move_1(GEN v, int ind){
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


//Given a bqf q, this gives the corresponding Descartes quadruple. if pos=-1, we give the quadruple with -n, and if pos=1, we give the quadruple with +n. If red=1 we reduce, else we don't.
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
  //Now D=-4n^2, a=+/-n (- if pos=0), and q=[A, B, C]->[a, A-a, C-a, A+C-a-B]
  GEN a= pos? n:negi(n);//+/-n
  GEN Ama=subii(gel(q, 1), a);
  GEN v=mkvec4(a, Ama, subii(gel(q, 3), a), addii(Ama, subii(gel(q, 3), gel(q, 2))));//The APC
  if(red) v=apol_red(v, 0);
  return gerepileupto(top, ZV_sort(v));
}

//Computes apol_ncgpforms, and returns the depths of the corresponding circles of curvature n.
GEN apol_ncgp_depths(GEN n, long prec){
  pari_sp top=avma;
  GEN forms=apol_ncgp_forms(n, 1, 0, prec);//The forms
  long maxdepth=0, lf=lg(forms);
  GEN depths=cgetg(lf, t_VECSMALL);
  for(long i=1;i<lf;i++){//Computing depths
    long d=apol_depth(gel(forms, i));
    depths[i]=d;
    if(d>maxdepth) maxdepth=d;
  }
  GEN dcount=zero_zv(maxdepth+1);//Depths 0 through d.
  for(long i=1;i<lf;i++) dcount[depths[i]+1]++;//Counting.
  return gerepileupto(top, dcount);
}

//Computes ncgp(-4n^2), and output the ACP's created from these quadratic forms with apol_make. We only count ambiguous forms ONCE.
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
    vectrunc_append(quads, apol_make(q, pos, red));
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
    GEN q=apol_red_partial(gel(forms, i), maxsteps);
    gel(curves, i)=gel(q, vecindexmin(q));
  }
  return gerepileupto(top, ZV_sort(curves));
}



//SEARCHING FOR CURVATURES


//Given a BOUNDED integral ACP, this returns equations for the circles with curvatures <=maxcurv at depth<=depth. An equation takes the form [curvature, radius, x, y].
GEN apol_circles(GEN v, GEN maxcurv, int depth, long prec){
  pari_sp top=avma;
  GEN vred=apol_red(v, 0);
  ZV_sort_inplace(vred);//So the minimal curvature occurs first.
  if(gequal0(gel(vred, 1))) pari_err_TYPE("Cannot be the strip packing", v);
  long maxcircs=10;//Stores the maximal number of circles+1, double it every time we try to go over.
  GEN clist=vectrunc_init(maxcircs);//Stores the cicles. Each entry is [[curvature, radius, x, y], previous indices], where previous index is the 3 previous circles it is tangent to (after the first four).
  GEN c1=gel(v, 1), c2=gel(v, 2), c3=gel(v, 3), c4=gel(v, 4);//Starting curvatures.
  GEN r1=Qdivii(gen_1, c1), r2=Qdivii(gen_1, c2);
  GEN circ1=mkvec4(c1, r1, gen_0, gen_0);//Outer circle
  vectrunc_append(clist, circ1);
  GEN circ2=mkvec4(c2, r2, gen_0, gneg(gadd(r1, r2)));//first inner circle, placed vertically at the top.
  vectrunc_append(clist, circ2);
  GEN circ3=thirdtangent(circ1, circ2, c3, c4, 0, prec);//Third circle goes left.
  vectrunc_append(clist, circ3);
  GEN circ4=thirdtangent(circ2, circ3, c4, c1, 0, prec);//Fourth circle is left of circ2 ->circ3.
  vectrunc_append(clist, circ4);
  long clistind=5;
  
  //Now we go down! Adopts the code of apol_search
  int ind=1;//We ind to track which depth we are going towards.
  GEN W=zerovec(depth);//Tracks the sequence of ACP's; W[ind] is at ind-1 depth
  GEN Winds=zerovec(depth);//Tracks the corresponding indices in clist
  gel(W, 1)=v;//The first one.
  gel(Winds, 1)=mkvecsmall4(1, 2, 3, 4);
  GEN I=vecsmall_ei(depth, 1);//Tracks the sequence of replacements
  int forward=1;//Tracks if we are going forward or not.
  do{//1<=ind<=depth is assumed.
    I[ind]=forward? 1:I[ind]+1;
    if(ind>1 && I[ind-1]==I[ind]) I[ind]++;//Don't repeat
    if(I[ind]>4){ind--;continue;}//Go back. Forward already must =0, so no need to update.
    //At this point, we can go on with valid and new inputs
    GEN newv=apol_move_1(gel(W, ind), I[ind]);//Make the move
    GEN newc=gel(newv, I[ind]);
    GEN newvecsmall=gen_0;//Need this to be accessible to the next if/else block
    int comp=cmpii(maxcurv, newc);//Comparing the new element to maxcurv.
    if(comp>=0){//Small enough!
      int is[3]={0, 0, 0};
      int isind=0;
      for(int i=1;i<=4;i++){
        if(I[ind]==i) continue;
        is[isind]=i;
        isind++;
      }//is are the three non-I[ind] indices in {1, 2, 3, 4}.
      GEN oldcirc1=gel(clist, gel(Winds, ind)[is[0]]);//One of the old circles
      GEN oldcirc2=gel(clist, gel(Winds, ind)[is[1]]);//One of the old circles
      GEN newcirc=thirdtangent(oldcirc1, oldcirc2, newc, gel(newv, is[2]), 1, prec);//The new circle, if it is to the right of newcirc1 ->newcirc2.
      GEN prevcirc=gel(clist, gel(Winds, ind)[I[ind]]);//The circle we are "replacing"
      if(gequal(newcirc, prevcirc)) newcirc=thirdtangent(oldcirc1, oldcirc2, newc, gel(newv, is[2]), 0, prec);//If the two curvatures were the same, this could trigger. If we lack oo precision, this could not work, and must be changed slightly.
      else{//This block must also be updated if there is not oo precision.
        GEN oldcirc3=gel(clist, gel(Winds, ind)[is[2]]);//The unused old circle. Our newcirc must be tangent to it.
        GEN rsums=gsqr(gadd(gel(oldcirc3, 2), gel(newcirc, 2)));//(r1+r2)^2
        GEN dcentres=gadd(gsqr(gsub(gel(oldcirc3, 3), gel(newcirc, 3))), gsqr(gsub(gel(oldcirc3, 4), gel(newcirc, 4))));//dist(centres)^2
        if(!gequal(rsums, dcentres)) newcirc=thirdtangent(oldcirc1, oldcirc2, newc, gel(newv, is[2]), 0, prec);//Must be the other side.
      }
      //Now we update things in clist.
      if(clistind==maxcircs){//Double the size
        maxcircs=2*maxcircs;
        GEN oldclist=clist;
        clist=vectrunc_init(maxcircs);
        vectrunc_append_batch(clist, oldclist);//Put the old circles back.
      }
      newvecsmall=cgetg(5, t_VECSMALL);
      for(int i=1;i<=4;i++){
        if(I[ind]==i) newvecsmall[i]=clistind;
        else newvecsmall[i]=gel(Winds, ind)[i];
      }
      vectrunc_append(clist, newcirc);
      clistind++;
    }
    if(ind==depth || comp<=0) forward=0;//Max depth OR the number is too big; once we reach or pass N, we cannot get N anymore.
    else{//We can keep going forward
      ind++;
      gel(W, ind)=newv;
      gel(Winds, ind)=newvecsmall;
      forward=1;
    }
  }while(ind>0);
  return gerepilecopy(top, clist);  
}

//Store a circle as [curvature, radius, x, y]. Given two tangent circles and a third curvature, this finds this third circle that is tangent to the first two. For internal tangency, we need a negative radius & curvature. There are always 2 places to put the circle: left or right of the line from circ1 to circ2. If right=1, we put it right, else we put it left. c4 is one of the curvatures to complete an Apollonian quadruple (supplying it allows us to always work with exact numbers in the case of integral ACPs).
static GEN thirdtangent(GEN circ1, GEN circ2, GEN c3, GEN c4, int right, long prec){
  pari_sp top=avma;
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
    GEN newv=apol_move_1(gel(W, ind), I[ind]);//Make the move
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

//Search for circles of curvature N up to depth depth. Returns the corresponding quadruples. Adopts the code of apol_orbit. Returns the qf's if rqf=1, and both if rqf=2 (each entry is [ACP, qf]).
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
    GEN newv=apol_move_1(gel(W, ind), I[ind]);//Make the move
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

