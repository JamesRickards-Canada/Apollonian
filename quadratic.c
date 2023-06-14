/*Methods to deal with quadratic forms.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"




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

//Generate the list of primes dividing D for which D/p^2 is a discriminant, can pass in facs=factorization of D
GEN discprimeindex(GEN D, GEN facs){
  pari_sp top=avma;
  if(!isdisc(D)) pari_err_TYPE("Not a discriminant", D);
  if(gequal0(facs)) facs=Z_factor(D);
  long numprimes=itos(gel(matsize(facs), 1));
  GEN curp=gen_0, curexp=gen_0;
  GEN plist=vectrunc_init(numprimes+1);
  for(long i=1;i<=numprimes;i++){//Run through the primes, test if we can divide out or not.
    curp=gcoeff(facs, i, 1);
    curexp=gcoeff(facs, i, 2);
    if(cmpis(curexp, 2)>=0){//Exponent needs to be at least 2
      if(cmpis(curp, 2)>0) vectrunc_append(plist, curp);//If p>2 we are good
      else if(cmpis(curexp, 4)>=0) vectrunc_append(plist, curp);//Power of 2 is at least 4=good
      else if(equalis(curexp, 2) && equali1(modis(shifti(D, -2), 4))) vectrunc_append(plist, curp);//4*1 mod 4
    }
  }
  return gerepilecopy(top, plist);
}

//Returns the set of divisors of D that are discriminants and whose quotient to D is a square.
GEN discsuperorders(GEN D){
  pari_sp top=avma;
  if(!isdisc(D)) pari_err_TYPE("D needs to be a discriminant", D);
  GEN divs=divisors(D);
  long ld=lg(divs), j=ld, dsig=signe(D);
  GEN list=vectrunc_init(ld);
  for(long i=1;i<ld;i++){
    j--;//i+j=ld
    if(Z_issquare(gel(divs, j))){
      GEN D;
      if(dsig==1) D=gel(divs, i);
      else D=negi(gel(divs, i));//Making the appropriate sign change.
      if(isdisc(D)) vectrunc_append(list, D);
    }
  }
  return gerepilecopy(top, list);
}

//Returns 1 if discriminant and 0 if not.
int isdisc(GEN D){
  if(typ(D)!=t_INT) return 0;//Checking integrality
  if(smodis(D, 4)<2 && !Z_issquare(D)) return 1;//0,1 mod 4 and not a square.
  return 0;
}






























//STATIC DECLARATIONS

//BASIC OPERATIONS ON BINARY QUADRATIC FORMS
static int bqf_compare(void *data, GEN q1, GEN q2);
static int bqf_compare_tmat(void *data, GEN d1, GEN d2);

//COMPOSITION/CLASS GROUP
static GEN bqf_ncgp_nonfundnarrow(GEN cgp, GEN D, GEN rootD);

//TUPLE REPS OF PRIMES
static int mod_collapse(GEN L);
static int bqf_tuplevalid_cmp(void *data, GEN x, GEN y);





//TEMPORARY FOR NOW
static long gen_search_temp(GEN T, GEN x, void *data, int (*cmp)(void*,GEN,GEN));

//Simply copied from bibli2.c, used in case we switch between 2.13.4 and 2.14
static long gen_search_temp(GEN T, GEN x, void *data, int (*cmp)(void*,GEN,GEN)){
  long u = lg(T)-1, i, l, s;

  if (!u) return -1;
  l = 1;
  do
  {
    i = (l+u) >> 1; s = cmp(data, x, gel(T,i));
    if (!s) return i;
    if (s < 0) u = i-1; else l = i+1;
  } while (u >= l);
  return -((s < 0)? i: i+1);
}


//BASIC OPERATIONS ON BINARY QUADRATIC FORMS



//Returns the generator of the automorphism group of q in PSL(2,Z) (which is Z if D>0, Z/2Z if D=-4, Z/3Z if D=-3, and 1 if D<-4
GEN bqf_automorph(GEN q){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)==-1) return gerepileupto(top, dbqf_automorph(q, D));
  return gerepileupto(top, ibqf_automorph_D(q, D));
}

//Compares two bqfs based on A, then B, then C. It is safe to modify this method, as sorting/searching will call the pointer to this function
static int bqf_compare(void *data, GEN q1, GEN q2){
  for(long j=1;j<=3;++j){
    switch(cmpii(gel(q1, j), gel(q2, j))){
      case -1:
        return -1;
      case 1:
        return 1;
    }
  }
  return 0;
}

//bqf_compare, but assumes the inputs are d1=[q1, mat1] and d2=[q2, mat2], and compares q1 and q2. This returns 0 if q1=q2, even if mat1!=mat2 (the mats don't matter)
static int bqf_compare_tmat(void *data, GEN d1, GEN d2){
  return bqf_compare(data, gel(d1, 1), gel(d2, 1));
}

//Returns the discriminant of q, which must be length 3 vector and integral.
GEN bqf_disc(GEN q){
  pari_sp top = avma;
  return gerepileupto(top, subii(sqri(gel(q, 2)), shifti(mulii(gel(q, 1), gel(q, 3)), 2)));
}

//Returns 1 if the BQF q (indefinite/positive definite) of discriminant D is reduced, and 0 if not. D needs to only have the correct sign.
int bqf_isreduced(GEN q, GEN D){
  pari_sp top=avma;
  if(!D) D=bqf_disc(q);
  int Dsign=signe(D);
  if(Dsign==1){//D>0, so AC<0 and B>|A+C|
    if(signe(gel(q, 1))==signe(gel(q, 3))) return gc_int(top, 0);//AC>0, not reduced.
    GEN s=addii(gel(q, 1),gel(q, 3));
    if(cmpii(gel(q, 2), absi(s))>0) return gc_int(top, 1);//Reduced!
	return gc_int(top, 0);//B<=|A+C|, not reduced.
  }
  else{//D<0, so |B|<=A<=C and if A=|B| or A=C then B>=0
    int ACcomp=cmpii(gel(q, 1), gel(q, 3));
    if(ACcomp>0) return gc_int(top, 0);//A>C, not reduced
    GEN babs=absi(gel(q, 2));//|B|
	int Bcomp=cmpii(babs, gel(q, 1));
    if(Bcomp>0) return gc_int(top, 0);//|B|>A, not reduced.
    if(Bcomp<0 && ACcomp<0) return gc_int(top, 1);//|B|<A<C
    if(signe(gel(q, 2))>=0) return gc_int(top, 1);//A=C or |B|=A, so need B>=0
	return gc_int(top, 0);//A=C or |B|=A, and B<0.
  }
  return gc_int(top, -1);//Some error with the inputs
}

//Generates a random proper bqf with max coefficient maxc. If type=1 it will be indefinite, type=-1 positive definite, type=0 either. primitive=1 means primitive, =0 means don't care. This is not designed for efficiency
GEN bqf_random(GEN maxc, int type, int primitive){
  pari_sp top=avma;
  GEN rmax=addis(shifti(maxc, 1), 1);//2*maxc+1
  for(;;){
	GEN q=mkvec3(subii(randomi(rmax), maxc), subii(randomi(rmax), maxc), subii(randomi(rmax), maxc));//Random numbers between -maxc and maxc
    GEN D=bqf_disc(q);
	if(primitive==1){if(!equali1(ZV_content(q))) continue;}// Not primitive
    if(type==1){if(signe(D)==1 && isdisc(D)==1) return gerepilecopy(top, q);}//Correct sign!
    else if(type==-1){
      if(signe(D)==-1){
        if(signe(gel(q, 1))==-1) ZV_togglesign(q);//Making positive definite
        return gerepilecopy(top, q);
      }
    }
    else{
      if(signe(D)==-1){
        if(signe(gel(q, 1))==-1) ZV_togglesign(q);//Making positive definite
        return gerepilecopy(top,q);
      }
      if(isdisc(D)==1) return gerepilecopy(top,q);
    }
  }
}

//Generates a random bqf of discriminant D with |B|<=2maxc and primitive. This is not designed for efficiency
GEN bqf_random_D(GEN maxc, GEN D){
  pari_sp top=avma;
  if(!isdisc(D)) return gen_0;
  GEN B=shifti(addis(randomi(maxc), 1), 1);//Random even number in [2, 2*maxc]
  if(smodis(D, 2)) B=subis(B, 1);//Random odd number in [1, 2*maxc-1]
  GEN AC=shifti(subii(sqri(B), D), -2);//A*C
  if(random_Fl(2)) B=negi(B);
  GEN g=gcdii(D, B);
  GEN Aposs=divisors(AC), A, C;
  long r, lx=lg(Aposs);
  for(;;){
	r=1+random_Fl(lx-1);
    A=gel(Aposs, r);
    C=gel(Aposs, lx-r);
    if(equali1(gcdii(gcdii(A, C), g))) break;
  }
  if(signe(D)==1){
	if(signe(AC)==-1){
      if(random_Fl(2)) A=negi(A);
      else C=negi(C);
	}
	else{
      if(random_Fl(2)){A=negi(A);C=negi(C);}
    }
  }
  return gerepilecopy(top, mkvec3(A, B, C));
}

//Reduce the bqf q of discriminant D with sign(D)=Dsign, rootD=sqrt(D) (can pass as NULL if D<0 as it is not needed), tmat=1 if we return the transition matrix and 0 else.
GEN bqf_red(GEN q, GEN rootD, int Dsign, int tmat){
  if(Dsign==1){
    if(tmat==0) return ibqf_red(q,rootD);
    return ibqf_red_tmat(q,rootD);
  }
  if(tmat==0) return dbqf_red(q);
  return dbqf_red_tmat(q);
}

//bqf_red with typecheck
GEN bqf_red_tc(GEN q, int tmat, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);//No squares allowed
  return gerepileupto(top,bqf_red(q,gsqrt(absi(D),prec),signe(D),tmat));
}

//Returns the roots of q in order.
GEN bqf_roots(GEN q, GEN D, GEN w){
  pari_sp top = avma;
  if(Z_issquare(D)){
    if(signe(gel(q, 1))==0){//q[1]=0
      if(signe(gel(q,2))==0){//Both roots oo
        GEN rvec=cgetg(3,t_VEC);//Return vector
        gel(rvec, 1) = mkoo();
        gel(rvec, 2) = mkoo();
        return rvec;//Nothing unnecessary is added to the stack
      }
      GEN a=Qdivii(negi(gel(q,3)),gel(q,2));//One of the roots
      GEN rvec=cgetg(3,t_VEC);//Return vector
      if(signe(gel(q,2))==1){//Sign is 1
        gel(rvec,1)=gcopy(a);//copy a into rvec[1]
        gel(rvec, 2) = mkoo();
        return gerepileupto(top,rvec);
      }
      else{//Sign is -1
        gel(rvec, 2)=gcopy(a);//copy a into rvec[2]
        gel(rvec, 1) = mkoo();
        return gerepileupto(top,rvec);
      }
    }//Now the roots are distinct, finite, and rational.
    GEN den=shifti(gel(q,1),1);
    GEN num1=addii(negi(gel(q,2)),w);//first root is (-q[2]+sqrt(D))/(2q[1])
    GEN num2=subii(negi(num1),shifti(gel(q,2),1));//Second root is (-q[2]-sqrt(D))/(2q[1])
    GEN rvec=cgetg(3,t_VEC);//Return vector
    gel(rvec,1)=Qdivii(num1,den);
    gel(rvec,2)=Qdivii(num2,den);
    return gerepileupto(top,rvec);
  }//Now the roots are as above, but need to use gen functions instead
  GEN den=shifti(gel(q,1),1);
  GEN num1=gadd(negi(gel(q,2)),w);//first root is (-q[2]+sqrt(D))/(2q[1])
  GEN num2=gsub(gneg(num1),shifti(gel(q,2),1));//Second root is (-q[2]-sqrt(D))/(2q[1])
  GEN rvec=cgetg(3,t_VEC);//Return vector
  gel(rvec,1)=gdiv(num1,den);
  gel(rvec,2)=gdiv(num2,den);
  return gerepileupto(top,rvec);
}

//bqf_roots, where we only submit q and perform checks
GEN bqf_roots_tc(GEN q){
  pari_sp top=avma;
  bqf_check(q);//We allow square discriminants
  GEN D=bqf_disc(q);
  GEN w=gen_0;
  if(Z_issquare(D)) w=sqrti(D);
  else w=quadroot(D);
  return gerepileupto(top,bqf_roots(q,D,w));
}

//takes in BQF q and matrix mtx, and outputs mtx*q.
GEN bqf_trans(GEN q, GEN mtx){
  pari_sp top = avma;//q=[A,B,C], mtx=[a,b;c,d], rvec=[Aa^2+Bac+Cc^2,2Aab+Bbc+Bad+2Ccd,Ab^2+Bbd+Cd^2]
  GEN Aa=mulii(gel(q,1),gcoeff(mtx,1,1));
  GEN Bc=mulii(gel(q,2),gcoeff(mtx,2,1));
  GEN Cc=mulii(gel(q,3),gcoeff(mtx,2,1));
  GEN AapBc=addii(Aa,Bc);
  GEN coef1=addii(mulii(AapBc,gcoeff(mtx,1,1)),mulii(Cc,gcoeff(mtx,2,1)));
  GEN coef2=addii(addii(mulii(addii(AapBc,Aa),gcoeff(mtx,1,2)),shifti(mulii(Cc,gcoeff(mtx,2,2)),1)),mulii(gel(q,2),mulii(gcoeff(mtx,1,1),gcoeff(mtx,2,2))));
  GEN coef3=addii(mulii(addii(mulii(gel(q,1),gcoeff(mtx,1,2)),mulii(gel(q,2),gcoeff(mtx,2,2))),gcoeff(mtx,1,2)),mulii(gel(q,3),sqri(gcoeff(mtx,2,2))));
  GEN rvec=cgetg(4,t_VEC);
  gel(rvec,1)=icopy(coef1);
  gel(rvec,2)=icopy(coef2);
  gel(rvec,3)=icopy(coef3);
  return gerepileupto(top,rvec);
}

//bqf_trans with typecheck.
GEN bqf_trans_tc(GEN q, GEN mtx){
  bqf_check(q);
  intmatrix_check(mtx);
  return bqf_trans(q,mtx);
}

//bqf_trans by L^n
GEN bqf_transL(GEN q, GEN n){
  GEN qnew=cgetg(4,t_VEC);
  pari_sp top=avma;
  GEN An=mulii(gel(q,1),n);
  GEN AnpB=addii(An,gel(q,2));
  GEN An2pBn=mulii(AnpB,n);
  pari_sp bot=avma;
  GEN qnewA=icopy(gel(q,1));
  GEN qnewB=addii(AnpB,An);
  GEN qnewC=addii(An2pBn,gel(q,3));
  gerepileallsp(top,bot,3,&qnewA,&qnewB,&qnewC);
  gel(qnew,1)=qnewA;
  gel(qnew,2)=qnewB;
  gel(qnew,3)=qnewC;
  return qnew;
}

//bqf_trans by R^n
GEN bqf_transR(GEN q, GEN n){
  GEN qnew=cgetg(4,t_VEC);
  pari_sp top=avma;
  GEN Cn=mulii(gel(q,3),n);
  GEN CnpB=addii(Cn,gel(q,2));
  GEN Cn2pBn=mulii(CnpB,n);
  pari_sp bot=avma;
  GEN qnewC=icopy(gel(q,3));
  GEN qnewB=addii(CnpB,Cn);
  GEN qnewA=addii(Cn2pBn,gel(q,1));
  gerepileallsp(top,bot,3,&qnewA,&qnewB,&qnewC);
  gel(qnew,1)=qnewA;
  gel(qnew,2)=qnewB;
  gel(qnew,3)=qnewC;
  return qnew;
}

//bqf_trans by S, i.e. [A,B,C]->[C,-B,A]
GEN bqf_transS(GEN q){
  long lx;
  GEN qnew=cgetg_copy(q,&lx);
  gel(qnew,1)=icopy(gel(q,3));
  gel(qnew,2)=icopy(gel(q,2));
  togglesign_safe(&gel(qnew,2));
  gel(qnew,3)=icopy(gel(q,1));
  return qnew;
}

//Outputs a form similar to q whose first coefficient is coprime to n. Useful for passing between ideals and quadratic forms, where we want to control the prime factors of the ideal. q must be definite/indefinite, and primitive (at least gcd(q,n)=1 is necessary). Note that this method is not very efficient or very good, but it works.
GEN bqf_trans_coprime(GEN q, GEN n){
  pari_sp top=avma;
  q=ZV_copy(q);
  GEN r, x=stoi(3);
  while(!equali1(gcdii(gel(q,1),n))){
    r=randomi(x);
    if(signe(r)==0) q=bqf_transL(q, gen_1);
    else if(equali1(r)) q=bqf_transR(q, gen_1);
    else q=bqf_transS(q);
  }
  return gerepileupto(top,q);
}

//bqf_trans_coprime with typechecking
GEN bqf_trans_coprime_tc(GEN q, GEN n){
  pari_sp top=avma;
  bqf_check(q);//Checking the disc
  if(typ(n)!=t_INT || gequal0(n)) pari_err_TYPE("Please enter a non-zero integer n", n);
  if(!equali1(gcdii(ZV_content(q),n))) pari_err_TYPE("gcd(q) must be coprime to n", n);
  avma=top;
  return bqf_trans_coprime(q, n);
}



//BASIC METHODS FOR NEGATIVE DISCRIMINANTS



//bqf_automorph for negative discriminants
GEN dbqf_automorph(GEN q, GEN D){
  pari_sp top=avma;
  int c=cmpis(D,-4);
  if(c<0) return mkmat22(gen_1,gen_0,gen_0,gen_1);//D<-4, no automorph.
  GEN M=gel(dbqf_red_tmat(q),2);
  if(c==0){//D=-4, and the autom is M[0,1;-1,0]M^-1
    return gerepileupto(top,ZM_mul(ZM_mul(M,mkmat22(gen_0,gen_1,gen_m1,gen_0)),ZM_inv(M,NULL)));
  }
  //Now D=-3, autom is M[0,1;-1,-1]M^-1
  return gerepileupto(top,ZM_mul(ZM_mul(M,mkmat22(gen_0,gen_1,gen_m1,gen_m1)),ZM_inv(M,NULL)));
}

//bqf_red for negative discriminants
GEN dbqf_red(GEN q){
  if(signe(gel(q,1))==-1){//Negating q if if negative definite
    pari_sp top=avma;
    GEN qnew=cgetg(3,t_VEC);
    for(int i=1;i<=3;++i) gel(qnew,i)=negi(gel(q,i));
    GEN negans=dbqf_red(qnew);
    return gerepileupto(top,gneg(negans));
  }//Now q is positive definite.
  pari_sp top=avma;
  if(bqf_isreduced(q, gen_m1)) return ZV_copy(q);//No garbage needed
  //Now q is not reduced, so the pointer to it will change and not point at the original q.
  GEN n=gen_0;//n represents if we are doing L^n or R^n
  if(signe(gel(q,2))==1) q=bqf_trans(q,mkmat22s(0,1,-1,0));//Pointing wrong way
  if(cmpii(gel(q,1),gel(q,3))<=0){//A<=C, so start by going L
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,1),1)));//n=floor(-B/2A)
    q=bqf_transL(q,n);
  }//Now we can loop doing R then L
  for(;;){
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,3),1)));//n=floor(-B/2C)
    if(equalii(n,gen_0)) break;//Can't go on!
    q=bqf_transR(q,n);
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,1),1)));//n=floor(-B/2A)
    if(equalii(n,gen_0)) break;//Can't go on!
    q=bqf_transL(q,n);
    if(gc_needed(top,1)) gerepileupto(top,q);//Garbage collection if memory needed.
  }
  //At this point we have B<0 and B+2A>0, B+2C>0. One of q, Sq, Lq, Rq, LSq, RSq is now reduced.
  int arelation=cmpii(gel(q,1),gel(q,3));
  GEN mB=negi(gel(q,2));
  GEN Lmat=gen_0;
  switch(arelation){
    case -1://A<C
      if(cmpii(mB,gel(q,1))==-1) Lmat=mkmat22s(1,0,0,1);//Already reduced
      else if(cmpii(mB,gel(q,3))<=0) Lmat=mkmat22s(1,1,0,1);//A<=-B<=C and A!=C-> L
      else Lmat=mkmat22s(-1,1,-1,0);//A<C<-B ->LS
      break;
    case 0://A=C
      if(cmpii(mB,gel(q,1))<=0) Lmat=mkmat22s(0,1,-1,0);//-B<=A=C ->S
      else Lmat=mkmat22s(1,0,1,1); //A=C<-B ->R
      break;
    case 1://A>C
      if(cmpii(mB,gel(q,3))<=0) Lmat=mkmat22s(0,1,-1,0);//-B<=C<A ->S
      else if(cmpii(mB,gel(q,1))==-1) Lmat=mkmat22s(0,1,-1,1);//C<-B<A ->RS
      else Lmat=mkmat22s(1,0,1,1);//C<A<=-B ->R
  }
  return gerepileupto(top,bqf_trans(q,Lmat));
}

//bqf_red for negative discriminants, also returns transition matrix
GEN dbqf_red_tmat(GEN q){
  if(signe(gel(q,1))==-1){//Negating q if if negative definite
    pari_sp top=avma;
    GEN qnew=cgetg(3,t_VEC);
    for(int i=1;i<=3;++i) gel(qnew,i)=negi(gel(q,i));
    GEN negans=dbqf_red_tmat(qnew);
    GEN rvec=cgetg(3,t_VEC);
    gel(rvec,1)=gneg(gel(negans,1));
    gel(rvec,2)=ZM_copy(gel(negans,2));
    return gerepileupto(top,rvec);
  }//Now q is positive definite.
  GEN rvec=cgetg(3,t_VEC);
  pari_sp top=avma;
  GEN Lmat=mkmat22s(1,0,0,1);//Will represent L^n
  if(bqf_isreduced(q, gen_m1)){
    gel(rvec,1)=ZV_copy(q);
    gel(rvec,2)=Lmat;//Just reusing the variable
    return rvec;//No garbage needed
  }//Now q is not reduced, so the pointer to it will change and not point at the original q.
  GEN Rmat=mkmat22s(1,0,0,1);//Will represent R^n
  pari_sp mid=avma, bot=avma;
  GEN n=gen_0;//n represents if we are doing L^n or R^n
  GEN mat=gen_0;//Will represent transition matrix
  if(signe(gel(q,2))==1){//Pointing wrong way
    mat=mkmat22s(0,1,-1,0);
    q=bqf_trans(q,mat);
  }
  else mat=mkmat22s(1,0,0,1);
  if(cmpii(gel(q,1),gel(q,3))<=0){//A<=C, so start by going L
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,1),1)));//n=floor(-B/2A)
    gcoeff(Lmat,1,2)=n;
    bot=avma;
    q=bqf_transL(q,n);
    mat=ZM_mul(mat,Lmat);
  }//Now we can loop doing R then L
  for(;;){
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,3),1)));//n=floor(-B/2C)
    if(equalii(n,gen_0)) break;//Can't go on!
    gcoeff(Rmat,2,1)=n;
    bot=avma;//If we exit loop after L, this is required for the final gerepileallsp
    q=bqf_transR(q,n);
    mat=ZM_mul(mat,Rmat);
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,1),1)));//n=floor(-B/2A)
    if(equalii(n,gen_0)) break;//Can't go on!
    gcoeff(Lmat,1,2)=n;
    bot=avma;
    q=bqf_transL(q,n);
    mat=ZM_mul(mat,Lmat);
    if(gc_needed(mid,1)) gerepileallsp(mid,bot,2,&q,&mat);//Garbage collection if memory needed.
  }
  //At this point we have B<0 and B+2A>0, B+2C>0. One of q, Sq, Lq, Rq, LSq, RSq is now reduced.
  int arelation=cmpii(gel(q,1),gel(q,3));
  GEN mB=negi(gel(q,2));
  switch(arelation){
    case -1://A<C
      if(cmpii(mB,gel(q,1))==-1) Lmat=mkmat22s(1,0,0,1);//Just reusing Lmat. Already reduced
      else if(cmpii(mB,gel(q,3))<=0) Lmat=mkmat22s(1,1,0,1);//A<=-B<=C and A!=C-> L
      else Lmat=mkmat22s(-1,1,-1,0);//A<C<-B ->LS
      break;
    case 0://A=C
      if(cmpii(mB,gel(q,1))<=0) Lmat=mkmat22s(0,1,-1,0);//-B<=A=C ->S
      else Lmat=mkmat22s(1,0,1,1); //A=C<-B ->R
      break;
    case 1://A>C
      if(cmpii(mB,gel(q,3))<=0) Lmat=mkmat22s(0,1,-1,0);//-B<=C<A ->S
      else if(cmpii(mB,gel(q,1))==-1) Lmat=mkmat22s(0,1,-1,1);//C<-B<A ->RS
      else Lmat=mkmat22s(1,0,1,1);//C<A<=-B ->R
  }
  bot=avma;
  q=bqf_trans(q,Lmat);
  mat=ZM_mul(mat,Lmat);
  gerepileallsp(top,bot,2,&q,&mat);//Garbage
  gel(rvec,1)=q;
  gel(rvec,2)=mat;
  return rvec;
}



//BASIC OPERATIONS SPECIFIC TO INDEFINITE FORMS/POSITIVE DISCRIMINANTS



//ibqf_automorph, but we pass in q and D and don't check. This is useful in a situation where we don't want to care about pell(D) in the ambient method.
GEN ibqf_automorph_D(GEN q, GEN D){
  pari_sp top=avma;
  return gerepileupto(top,ibqf_automorph_pell(q,pell(D)));
}

//Returns the invariant automorph of the PIBQF q. If q is not primitive, the output may be wrong.
GEN ibqf_automorph_pell(GEN q, GEN qpell){
  pari_sp top = avma;
  GEN a=subii(gel(qpell,1),mulii(gel(q,2),gel(qpell,2)));//Top left is a/2
  GEN d=addii(negi(a),shifti(gel(qpell,1),1));//Bot right is d/2
  GEN b=mulii(gel(q,3),gel(qpell,2));//top right is -b
  GEN rvec=cgetg(3,t_MAT);
  gel(rvec,1)=cgetg(3,t_COL);
  gel(rvec,2)=cgetg(3,t_COL);
  gcoeff(rvec,1,1)=shifti(a,-1);//Divide by 2
  gcoeff(rvec,1,2)=negi(b);
  gcoeff(rvec,2,1)=mulii(gel(q,1),gel(qpell,2));
  gcoeff(rvec,2,2)=shifti(d,-1);//Divide by 2
  return gerepileupto(top,rvec);
}

//Returns 1 if q is reciprocal
int ibqf_isrecip(GEN q, GEN rootD){
  pari_sp top=avma;
  long lx;
  GEN q2=cgetg_copy(q,&lx);
  for(int i=1;i<4;i++) gel(q2,i)=negi(gel(q,i));
  GEN z=ibqf_isequiv(q,q2,rootD);
  if(equali1(z)){avma=top;return 1;}
  avma=top;
  return 0;
}

//ibqf_isrecip with typechecking
int ibqf_isrecip_tc(GEN q, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)==-1) pari_err_TYPE("Please enter an indefinite binary quadratic form q", q);
  int i=ibqf_isrecip(q, gsqrt(D,prec));
  avma=top;
  return i;
}

//Given river form q and rootD=sqrt(D), finds the left neighbour of q and returns it
GEN ibqf_leftnbr(GEN q, GEN rootD){
  pari_sp top = avma;
  if(signe(gel(q,1))==1){//A>0
    GEN AmB=subii(gel(q,1),gel(q,2));
    GEN AmBpC=addii(AmB,gel(q,3));
    if(signe(AmBpC)==1){//Apply (0 1;-1 0), then L maximally, then (0 1;-1 0)
      avma=top;
      GEN delta=floorr(divri(subir(gel(q,2),rootD),shifti(gel(q,3),1)));//floor((B-sqrt(D))/(2C))
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_m1;
      gcoeff(tmat,1,2)=gen_0;
      gcoeff(tmat,2,1)=delta;
      gcoeff(tmat,2,2)=gen_m1;
      return gerepileupto(top,bqf_trans(q,tmat));
    }
    else{//Apply (0 1;-1 0), then R maximally
      avma=top;
      GEN delta=floorr(divri(addri(rootD,gel(q,2)),shifti(gel(q,1),1)));//floor((sqrt(D)+B)/(2A)
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=delta;
      gcoeff(tmat,1,2)=gen_1;
      gcoeff(tmat,2,1)=gen_m1;
      gcoeff(tmat,2,2)=gen_0;
      return gerepileupto(top,bqf_trans(q,tmat));
    }
  }
  else{
    GEN ApB=addii(gel(q,1),gel(q,2));
    GEN ApBpC=addii(ApB,gel(q,3));
    if(signe(ApBpC)==1){//Apply L maximally, then (0 1;-1 0)
      avma=top;
      GEN delta=floorr(divri(negi(addri(rootD,gel(q,2))),shifti(gel(q,1),1)));//floor((-sqrt(D)-B)/2A)
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=delta;
      gcoeff(tmat,1,2)=gen_m1;
      gcoeff(tmat,2,1)=gen_1;
      gcoeff(tmat,2,2)=gen_0;
      return gerepileupto(top,bqf_trans(q,tmat));
    }
    else{//Apply then R maximally
      avma=top;
      GEN delta=floorr(divri(subri(rootD,gel(q,2)),shifti(gel(q,3),1)));//floor((sqrt(D)-B)/2C)
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_1;
      gcoeff(tmat,1,2)=gen_0;
      gcoeff(tmat,2,1)=icopy(delta);
      gcoeff(tmat,2,2)=gen_1;
      return gerepileupto(top,bqf_trans(q,tmat));
    }
  } 
}

//Given river form q and rootD=sqrt(D), finds the left neighbour of q and the transition matrix
GEN ibqf_leftnbr_tmat(GEN q, GEN rootD){
  pari_sp top = avma;
  if(signe(gel(q,1))==1){//A>0
    GEN AmB=subii(gel(q,1),gel(q,2));
    GEN AmBpC=addii(AmB,gel(q,3));
    if(signe(AmBpC)==1){//Apply (0 1;-1 0), then L maximally, then (0 1;-1 0)
      avma=top;
      GEN delta=floorr(divri(subir(gel(q,2),rootD),shifti(gel(q,3),1)));//floor((B-sqrt(D))/(2C))
      GEN rvec=cgetg(3,t_VEC);
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_m1;
      gcoeff(tmat,1,2)=gen_0;
      gcoeff(tmat,2,1)=icopy(delta);
      gcoeff(tmat,2,2)=gen_m1;
      gel(rvec,2)=tmat;
      gel(rvec,1)=bqf_trans(q,tmat);
      return gerepileupto(top,rvec);
    }
    else{//Apply (0 1;-1 0), then R maximally
      avma=top;
      GEN delta=floorr(divri(addri(rootD,gel(q,2)),shifti(gel(q,1),1)));//floor((sqrt(D)+B)/(2A)
      GEN rvec=cgetg(3,t_VEC);
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=icopy(delta);
      gcoeff(tmat,1,2)=gen_1;
      gcoeff(tmat,2,1)=gen_m1;
      gcoeff(tmat,2,2)=gen_0;
      gel(rvec,2)=tmat;
      gel(rvec,1)=bqf_trans(q,tmat);
      return gerepileupto(top,rvec);
    }
  }
  else{
    GEN ApB=addii(gel(q,1),gel(q,2));
    GEN ApBpC=addii(ApB,gel(q,3));
    if(signe(ApBpC)==1){//Apply L maximally, then (0 1;-1 0)
      avma=top;
      GEN delta=floorr(divri(negi(addri(rootD,gel(q,2))),shifti(gel(q,1),1)));//floor((-sqrt(D)-B)/2A)
      GEN rvec=cgetg(3,t_VEC);
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=icopy(delta);
      gcoeff(tmat,1,2)=gen_m1;
      gcoeff(tmat,2,1)=gen_1;
      gcoeff(tmat,2,2)=gen_0;
      gel(rvec,2)=tmat;
      gel(rvec,1)=bqf_trans(q,tmat);
      return gerepileupto(top,rvec);
    }
    else{//Apply then R maximally
      avma=top;
      GEN delta=floorr(divri(subri(rootD,gel(q,2)),shifti(gel(q,3),1)));//floor((sqrt(D)-B)/2C)
      GEN rvec=cgetg(3,t_VEC);
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_1;
      gcoeff(tmat,1,2)=gen_0;
      gcoeff(tmat,2,1)=icopy(delta);
      gcoeff(tmat,2,2)=gen_1;
      gel(rvec,2)=tmat;
      gel(rvec,1)=bqf_trans(q,tmat);
      return gerepileupto(top,rvec);
    }
  } 
}

//ibqf_leftnbr, updating the transition matrix
GEN ibqf_leftnbr_update(GEN qvec, GEN rootD){
  pari_sp top=avma;
  GEN rn=ibqf_leftnbr_tmat(gel(qvec,1),rootD);
  long lx;
  GEN rvec=cgetg_copy(rn,&lx);
  gel(rvec,1)=ZV_copy(gel(rn,1));
  gel(rvec,2)=ZM_mul(gel(qvec,2),gel(rn,2));
  return gerepileupto(top,rvec);
}

//ibqf_leftnbr but with type checking
GEN ibqf_leftnbr_tc(GEN q, int tmat, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q);
  if(signe(gel(q,1))==signe(gel(q,3))) pari_warn(warner,"form is not on the river!");
  if(tmat==0) return gerepileupto(top,ibqf_leftnbr(q,gsqrt(D,prec)));
  return gerepileupto(top,ibqf_leftnbr_tmat(q,gsqrt(D,prec)));
}

//bqf_red for positive discriminants
GEN ibqf_red(GEN q, GEN rootD){
  pari_sp top=avma;
  GEN toriv=ibqf_toriver(q,rootD);//Now the form is on the river
  if(bqf_isreduced(toriv, gen_1)) return toriv;//No garbage!
  return gerepilecopy(top,ibqf_rightnbr(toriv,rootD));//If not reduced, take the right neighbour
}

//bqf_red for positive discriminants, also returns transition matrix.
GEN ibqf_red_tmat(GEN q, GEN rootD){
  pari_sp top=avma;
  GEN toriv=ibqf_toriver_tmat(q,rootD);//Now the form is on the river
  if(bqf_isreduced(gel(toriv, 1), gen_1)) return toriv;
  return gerepileupto(top,ibqf_rightnbr_update(toriv,rootD));
}

//Reduces q to a form with A>0
GEN ibqf_red_pos(GEN q, GEN rootD){
  pari_sp top=avma;
  GEN red1=ibqf_red(q,rootD);
  if(signe(gel(red1,1))==-1) red1=ibqf_leftnbr(red1,rootD);
  return gerepileupto(top,red1);
}

//ibqf_red_pos plus, transition matrix.
GEN ibqf_red_pos_tmat(GEN q, GEN rootD){
  pari_sp top=avma;
  GEN red1=ibqf_red_tmat(q,rootD);
  if(signe(gel(gel(red1,1),1))==1) return red1;//Done!
  return gerepileupto(top,ibqf_leftnbr_update(red1,rootD));//Go left to make A>0
}

//Returns the reduced orbit of q
GEN ibqf_redorbit(GEN q, GEN rootD){
  pari_sp top = avma;
  GEN qred=ibqf_red(q,rootD);
  glist *orbit=NULL;//orbit list
  GEN qseek=qred;
  long nforms=0;
  do{
    glist_putstart(&orbit,qseek);
    qseek=ibqf_rightnbr(qseek,rootD);
    nforms++;
  }
  while(!ZV_equal(qred,qseek));
  return gerepileupto(top,glist_togvec(orbit, nforms, -1));
}

//Returns the reduced orbit of q along with the corresponding transition matrices
GEN ibqf_redorbit_tmat(GEN q, GEN rootD){
  pari_sp top = avma;
  GEN qred=ibqf_red_tmat(q,rootD);
  clist *origin=NULL;//orbit list
  clist_putbefore(&origin,qred);
  clist *left=origin, *right=origin;
  GEN qseekL=qred, qseekR=qred;
  long nforms=1;
  for(;;){
    qseekR=ibqf_rightnbr_update(qseekR,rootD);//No need to check as this will NEVER equal qseekL; the signs of the first coefficient will always be distinct here.
    clist_putafter(&right,qseekR);
    nforms++;
    qseekL=ibqf_leftnbr_update(qseekL,rootD);
    if(ZV_equal(gel(qseekL,1),gel(qseekR,1))) break;//Done!
    clist_putbefore(&left,qseekL);
    nforms++;
  }
  return gerepileupto(top,clist_togvec(origin,nforms,1));
}

//Returns the forms with A>0 in the reduced orbit of q
GEN ibqf_redorbit_posonly(GEN q, GEN rootD){
  pari_sp top = avma;
  GEN qred=ibqf_red_pos(q,rootD);
  glist *orbit=NULL;//orbit list
  GEN qseek=qred;
  long nforms=0;
  do{
    glist_putstart(&orbit,qseek);
    qseek=ibqf_rightnbr(ibqf_rightnbr(qseek,rootD),rootD);
    nforms++;
  }
  while(!ZV_equal(qred,qseek));
  return gerepileupto(top,glist_togvec(orbit, nforms, -1));
}

//Returns the forms with A>0 in the reduced orbit of q along with the transition matrices
GEN ibqf_redorbit_posonly_tmat(GEN q, GEN rootD){
  pari_sp top = avma;
  GEN qred=ibqf_red_pos_tmat(q,rootD);
  clist *origin=NULL;//orbit list
  clist_putbefore(&origin,qred);
  clist *left=origin, *right=origin;
  GEN qseekL=qred, qseekR=qred;
  long nforms=1;
  for(;;){
    qseekR=ibqf_rightnbr_update(ibqf_rightnbr_update(qseekR,rootD),rootD);
    if(ZV_equal(gel(qseekL,1),gel(qseekR,1))) break;//Done!
    clist_putafter(&right,qseekR);
    nforms++;
    qseekL=ibqf_leftnbr_update(ibqf_leftnbr_update(qseekL,rootD),rootD);
    if(ZV_equal(gel(qseekL,1),gel(qseekR,1))) break;//Done!
    clist_putbefore(&left,qseekL);
    nforms++;
  }
  return gerepileupto(top,clist_togvec(origin,nforms,1));
}

//ibqf_redorbit with type checking
GEN ibqf_redorbit_tc(GEN q, int tmat, int posonly, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q);
  if(posonly==0){
    if(tmat==0) return gerepileupto(top,ibqf_redorbit(q,gsqrt(D,prec)));
    else return gerepileupto(top,ibqf_redorbit_tmat(q,gsqrt(D,prec)));
  }
  if(tmat==0) return gerepileupto(top,ibqf_redorbit_posonly(q,gsqrt(D,prec)));
  return gerepileupto(top,ibqf_redorbit_posonly_tmat(q,gsqrt(D,prec)));
}

//Given river form q and rootD=sqrt(D), finds the right neighbour of q and returns it
GEN ibqf_rightnbr(GEN q, GEN rootD){
  pari_sp top = avma;
  if(signe(gel(q,1))==1){//A>0
    GEN ApB=addii(gel(q,1),gel(q,2));
    GEN ApBpC=addii(ApB,gel(q,3));
    if(signe(ApBpC)==1){//Applying R maximally, then (0 1;-1 0)
      avma=top;
      GEN delta=floorr(divri(negr(addir(gel(q,2),rootD)),shifti(gel(q,3),1)));//floor((-B-sqrt(D))/(2C)
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_0;
      gcoeff(tmat,1,2)=gen_1;
      gcoeff(tmat,2,1)=gen_m1;
      gcoeff(tmat,2,2)=delta;
      return gerepileupto(top,bqf_trans(q,tmat));
    }
    else{//Applying L maximally
      avma=top;
      GEN delta=floorr(divri(subri(rootD,gel(q,2)),shifti(gel(q,1),1)));//floor((sqrt(D)-B)/(2A)
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_1;
      gcoeff(tmat,1,2)=delta;
      gcoeff(tmat,2,1)=gen_0;
      gcoeff(tmat,2,2)=gen_1;
      return gerepileupto(top,bqf_trans(q,tmat));
    }
  }
  else{
    GEN AmB=subii(gel(q,1),gel(q,2));
    GEN AmBpC=addii(AmB,gel(q,3));
    if(signe(AmBpC)==1){//Apply (0 1;-1 0), then R maximally, then (0 1;-1 0)
      avma=top;
      GEN delta=ceilr(divri(subri(rootD,gel(q,2)),shifti(gel(q,1),1)));//ceil((sqrt(D)-B)/2A)
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_1;
      gcoeff(tmat,1,2)=delta;
      gcoeff(tmat,2,1)=gen_0;
      gcoeff(tmat,2,2)=gen_1;
      return gerepileupto(top,bqf_trans(q,tmat));
    }
    else{//Apply (0 1;-1 0) then L maximally
      avma=top;
      GEN delta=floorr(divri(addri(rootD,gel(q,2)),shifti(gel(q,3),1)));//floor((B+sqrt(D)/2C)
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_0;
      gcoeff(tmat,1,2)=gen_m1;
      gcoeff(tmat,2,1)=gen_1;
      gcoeff(tmat,2,2)=delta;
      return gerepileupto(top,bqf_trans(q,tmat));
    }
  }
}

//Given river form q and rootD=sqrt(D), finds the right neighbour of q and the transition matrix
GEN ibqf_rightnbr_tmat(GEN q, GEN rootD){
  pari_sp top = avma;
  if(signe(gel(q,1))==1){//A>0
    GEN ApB=addii(gel(q,1),gel(q,2));
    GEN ApBpC=addii(ApB,gel(q,3));
    if(signe(ApBpC)==1){//Applying R maximally, then (0 1;-1 0)
      avma=top;
      GEN delta=floorr(divri(negr(addir(gel(q,2),rootD)),shifti(gel(q,3),1)));//floor((-B-sqrt(D))/(2C)
      GEN rvec=cgetg(3,t_VEC);
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_0;
      gcoeff(tmat,1,2)=gen_1;
      gcoeff(tmat,2,1)=gen_m1;
      gcoeff(tmat,2,2)=icopy(delta);
      gel(rvec,2)=tmat;
      gel(rvec,1)=bqf_trans(q,tmat);
      return gerepileupto(top,rvec);
    }
    else{//Applying L maximally
      avma=top;
      GEN delta=floorr(divri(subri(rootD,gel(q,2)),shifti(gel(q,1),1)));//floor((sqrt(D)-B)/(2A)
      GEN rvec=cgetg(3,t_VEC);
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_1;
      gcoeff(tmat,1,2)=icopy(delta);
      gcoeff(tmat,2,1)=gen_0;
      gcoeff(tmat,2,2)=gen_1;
      gel(rvec,2)=tmat;
      gel(rvec,1)=bqf_trans(q,tmat);
      return gerepileupto(top,rvec);
    }
  }
  else{
    GEN AmB=subii(gel(q,1),gel(q,2));
    GEN AmBpC=addii(AmB,gel(q,3));
    if(signe(AmBpC)==1){//Apply (0 1;-1 0), then R maximally, then (0 1;-1 0)
      avma=top;
      GEN delta=ceilr(divri(subri(rootD,gel(q,2)),shifti(gel(q,1),1)));//ceil((sqrt(D)-B)/2A)
      GEN rvec=cgetg(3,t_VEC);
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_1;
      gcoeff(tmat,1,2)=icopy(delta);
      gcoeff(tmat,2,1)=gen_0;
      gcoeff(tmat,2,2)=gen_1;
      gel(rvec,2)=tmat;
      gel(rvec,1)=bqf_trans(q,tmat);
      return gerepileupto(top,rvec);
    }
    else{//Apply (0 1;-1 0) then L maximally
      avma=top;
      GEN delta=floorr(divri(addri(rootD,gel(q,2)),shifti(gel(q,3),1)));//floor((B+sqrt(D)/2C)
      GEN rvec=cgetg(3,t_VEC);
      GEN tmat=cgetg(3,t_MAT);
      gel(tmat,1)=cgetg(3, t_COL);
      gel(tmat,2)=cgetg(3, t_COL);
      gcoeff(tmat,1,1)=gen_0;
      gcoeff(tmat,1,2)=gen_m1;
      gcoeff(tmat,2,1)=gen_1;
      gcoeff(tmat,2,2)=icopy(delta);
      gel(rvec,2)=tmat;
      gel(rvec,1)=bqf_trans(q,tmat);
      return gerepileupto(top,rvec);
    }
  }
}

//ibqf_rightnbr, updating the transition matrix
GEN ibqf_rightnbr_update(GEN qvec, GEN rootD){
  pari_sp top=avma;
  GEN rn=ibqf_rightnbr_tmat(gel(qvec,1),rootD);
  long lx;
  GEN rvec=cgetg_copy(rn,&lx);
  gel(rvec,1)=ZV_copy(gel(rn,1));
  gel(rvec,2)=ZM_mul(gel(qvec,2),gel(rn,2));
  return gerepileupto(top,rvec);
}

//ibqf_rightnbr but with type checking
GEN ibqf_rightnbr_tc(GEN q, int tmat, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q);
  if(signe(gel(q,1))==signe(gel(q,3))) pari_warn(warner,"form is not on the river!");
  if(tmat==0) return gerepileupto(top,ibqf_rightnbr(q,gsqrt(D,prec)));
  return gerepileupto(top,ibqf_rightnbr_tmat(q,gsqrt(D,prec)));
}

//Returns the river sequence corresponding to q. If q is already reduced, the sequence will start at q. Otherwise it will start at bqf_red(q). 1=R and 0=L
GEN ibqf_river(GEN q, GEN rootD){
  pari_sp top = avma;
  GEN qriv=ibqf_red(q,rootD);
  int dir=2;
  if(signe(gel(qriv,1))==1) dir=1;//which direction we start going.
  llist *seq=NULL;//The sequence of moves.
  pari_sp mid=avma;
  GEN qpos=qriv;
  GEN trans=gen_0;
  long rlen=0;
  long change;
  do{
    trans=ibqf_rightnbr_tmat(qpos,rootD);
    qpos=gel(trans,1);
    change=itos(gcoeff(gel(trans,2),2,2));
    llist_putstart(&seq,change);
    rlen=rlen+change;
    if(gc_needed(mid,1)) gerepileupto(mid,qpos);//Garbage collection.
  }
  while(!ZV_equal(qriv,qpos));//Generate the sequence
  avma=top;//Everything necessary is now stored in longs.
  GEN rvec=cgetg(rlen+1,t_VEC);
  mid=avma;
  GEN elts=mkvec2(gen_0,gen_1);
  long pos=rlen;
  while(seq!=NULL){
    for(change=1;change<=seq->data;++change){
      gel(rvec,pos)=gel(elts,dir);//Okay as using universal objects
      pos--;
    }
    dir=3-dir;//Swapping directions
    seq=seq->next;//Moving on
  }
  avma=mid;//Okay as rvec just stores gen_0 and gen_1's
  return rvec;
}

//Returns [VECSMALL of indices going L, VECSMALL of indices going R, VECSMALL of ibqf_river]
GEN ibqf_river_positions(GEN q, GEN rootD){
  pari_sp top = avma;
  GEN qriv=ibqf_red(q,rootD);
  int dir=2;
  if(signe(gel(qriv,1))==1) dir=1;//which direction we start going.
  llist *L=NULL;//Seq of L's
  llist *R=NULL;//Seq of R's
  long ind=1, llen=0, rlen=0, change;//Current index, Length of L's, length of R's, the length of the next part
  int curdir=dir;
  pari_sp mid=avma;
  GEN qpos=qriv;
  GEN trans=gen_0;
  do{
    trans=ibqf_rightnbr_tmat(qpos,rootD);
    qpos=gel(trans,1);
    change=itos(gcoeff(gel(trans,2),2,2));
    if(curdir==1){
      for(long i=1;i<=change;i++){llist_putstart(&R,ind);ind++;}
      rlen=rlen+change;curdir=2;
    }
    else{
      for(long i=1;i<=change;i++){llist_putstart(&L,ind);ind++;}
      llen=llen+change;curdir=1;
    }
    if(gc_needed(mid,1)) gerepileupto(mid,qpos);//Garbage collection.
  }
  while(!ZV_equal(qriv,qpos));//Generate the sequence
  avma=top;//Everything necessary is now stored in longs.
  GEN rvec=cgetg(4,t_VEC);
  gel(rvec,1)=llist_tovecsmall(L,llen,-1);
  gel(rvec,2)=llist_tovecsmall(R,rlen,-1);
  gel(rvec,3)=cgetg(rlen+llen+1,t_VECSMALL);
  for(long i=1;i<=llen;i++) gel(rvec,3)[gel(rvec,1)[i]]=0;
  for(long i=1;i<=rlen;i++) gel(rvec,3)[gel(rvec,2)[i]]=1;
  return rvec;
}

//Returns [VECSMALL of indices going L, VECSMALL of indices going R, VECSMALL of ibqf_river, vector of corresponding forms]
GEN ibqf_river_positions_forms(GEN q, GEN rootD){
  GEN rvec=cgetg(5,t_VEC);
  pari_sp top = avma;
  GEN qstart=ibqf_red(q,rootD);
  if(signe(gel(qstart,1))==-1) qstart=bqf_trans(qstart,mkmat22s(0,1,-1,0));//Making A>0
  glist *orbit=NULL;//orbit list
  llist *ldir=NULL;//List of indices going left
  llist *rdir=NULL;//List of indices going right
  GEN qseek=qstart;
  long nlforms=0;
  long nrforms=0;
  long nforms=0;
  GEN ApB=gen_0;
  pari_sp mid;
  do{
    glist_putstart(&orbit,qseek);
    mid=avma;
    nforms++;
    ApB=addii(gel(qseek,1),gel(qseek,2));
    if(cmpii(ApB,negi(gel(qseek,3)))==1){//A+B+C>0, go right 1
      avma=mid;
      qseek=bqf_transR(qseek,gen_1);
      llist_putstart(&rdir,nforms);
      nrforms++;
    }
    else{//A+B+C<0, go left 1
      avma=mid;
      qseek=bqf_transL(qseek,gen_1);
      llist_putstart(&ldir,nforms);
      nlforms++;
    }
  }
  while(!ZV_equal(qstart,qseek));
  gel(rvec,4)=gerepileupto(top,glist_togvec(orbit, nforms, -1));//Now memory is totally clean.
  gel(rvec,1)=llist_tovecsmall(ldir, nlforms, -1);
  gel(rvec,2)=llist_tovecsmall(rdir, nrforms, -1);
  gel(rvec,3)=cgetg(nforms+1,t_VECSMALL);
  for(long i=1;i<=nlforms;i++) gel(rvec,3)[gel(rvec,1)[i]]=0;
  for(long i=1;i<=nrforms;i++) gel(rvec,3)[gel(rvec,2)[i]]=1;
  return rvec;
}

//ibqf_river with typechecking
GEN ibqf_river_tc(GEN q, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q);
  return gerepileupto(top,ibqf_river(q,gsqrt(D,prec)));
}

//Returns the forms on the river of q with first coefficient positive
GEN ibqf_riverforms(GEN q, GEN rootD){
  pari_sp top = avma;
  GEN qstart=ibqf_red(q,rootD);
  if(signe(gel(qstart,1))==-1) qstart=gerepileupto(top,bqf_trans(qstart,mkmat22s(0,1,-1,0)));//Making A>0
  glist *orbit=NULL;//orbit list
  GEN qseek=qstart;
  long nforms=0;
  GEN ApB=gen_0;
  pari_sp mid;
  do{
    glist_putstart(&orbit,qseek);
    mid=avma;
    ApB=addii(gel(qseek,1),gel(qseek,2));
    if(cmpii(ApB,negi(gel(qseek,3)))==1){//A+B+C>0, go right 1
      avma=mid;
      qseek=bqf_transR(qseek,gen_1);
    }
    else{//A+B+C<0, go left 1
      avma=mid;
      qseek=bqf_transL(qseek,gen_1);
    }   
    nforms++;
  }
  while(!ZV_equal(qstart,qseek));
  return gerepileupto(top,glist_togvec(orbit, nforms, -1));
}

//ibqf_riverforms with typechecking.
GEN ibqf_riverforms_tc(GEN q, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q);
  return gerepileupto(top,ibqf_riverforms(q,gsqrt(D,prec)));
}

//Outputs [z, gamma_q(z)] where z is on the root geodesic of q and these points are symmetric on the circle. The equation for z is x+iy with (x,y)=(-B/2A-DU/2AT,sqrt(D)/(|A|T))
GEN ibqf_symmetricarc(GEN q, GEN D, GEN rootD, GEN qpell, long prec){
  pari_sp top=avma;
  GEN AT=absi(mulii(gel(q,1),gel(qpell,1)));
  GEN x=gneg(gadd(Qdivii(gel(q,2),shifti(gel(q,1),1)),Qdivii(mulii(D,gel(qpell,2)),shifti(AT,1))));//the x val
  GEN mid=Qdivii(negi(gel(q,2)),gel(q,1));//-B/A
  GEN rvec=cgetg(3,t_VEC);
  GEN z=cgetc(prec);
  gel(z,1)=gtofp(x,prec);
  gel(z,2)=gdiv(rootD,AT);
  gel(rvec,1)=z;
  gel(rvec,2)=cgetc(prec);
  gel(gel(rvec,2),1)=gsub(mid,gel(z,1));//The real part is reflected
  gel(gel(rvec,2),2)=gcopy(gel(z,2));//Same imaginary part
  return gerepileupto(top,rvec);
}

//ibqf_symmetricarc with typechecking
GEN ibqf_symmetricarc_tc(GEN q, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)!=1) pari_err_TYPE("Please enter a PIBQF q", q);
  GEN qpell=pell(D);
  return gerepileupto(top,ibqf_symmetricarc(q, D, gsqrt(D,prec), qpell, prec));
}

//Reduces a PIBQF to the river, and not necessarily a reduced form.
GEN ibqf_toriver(GEN q, GEN rootD){
  if(signe(gel(q,1))!=signe(gel(q,3))) return ZV_copy(q);//Already on river!
  pari_sp top=avma;
  GEN n=gen_0;//n represents if we are doing L^n or R^n
  if(signe(gel(q,1))==signe(gel(q,2))) q=bqf_trans(q,mkmat22s(0,1,-1,0));//Pointing wrong way
  if((signe(gel(q,1))==1 && cmpii(gel(q,1),gel(q,3))<=0) || (signe(gel(q,1))==-1 && cmpii(gel(q,1),gel(q,3))>=0)){//|A|<=|C|, so go with L
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,1),1)));//n=floor(-B/2A)
    q=bqf_transL(q,n);
  }//Now we can loop doing R then L
  while(signe(gel(q,1))==signe(gel(q,3))){//A,C still positive
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,3),1)));//n=floor(-B/2C)
    q=bqf_transR(q,n);
    if(signe(gel(q,1))!=signe(gel(q,3))) break;//Done! Don't go L anymore.
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,1),1)));//n=floor(-B/2A)
    q=bqf_transL(q,n);
    if(gc_needed(top,1)) gerepileupto(top,q);//Garbage collection if memory needed.
  }
  return gerepileupto(top,q);
}

//Reduces a PIBQF to the river, and not necessarily a reduced form. Also returns the transition matrix.
GEN ibqf_toriver_tmat(GEN q, GEN rootD){
  GEN rvec=cgetg(3,t_VEC);
  pari_sp top=avma;
  if(signe(gel(q,1))!=signe(gel(q,3))){//Already on river!
    gel(rvec,1)=ZV_copy(q);
    gel(rvec,2)=mkmat22s(1,0,0,1);
    return rvec;//No garbage
  }
  GEN Lmat=mkmat22s(1,0,0,1);//Will represent L^n
  GEN Rmat=mkmat22s(1,0,0,1);//Will represent R^n
  pari_sp mid=avma, bot=avma;
  GEN n=gen_0;//n represents if we are doing L^n or R^n
  GEN mat=gen_0;//Will represent transition matrix
  if(signe(gel(q,1))==signe(gel(q,2))){//Pointing wrong way
    mat=mkmat22s(0,1,-1,0);
    q=bqf_trans(q,mat);
  }
  else mat=mkmat22s(1,0,0,1);
  if((signe(gel(q,1))==1 && cmpii(gel(q,1),gel(q,3))<=0) || (signe(gel(q,1))==-1 && cmpii(gel(q,1),gel(q,3))>=0)){//|A|<=|C|, so go with L
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,1),1)));//n=floor(-B/2A)
    gcoeff(Lmat,1,2)=n;
    bot=avma;
    q=bqf_transL(q,n);
    mat=ZM_mul(mat,Lmat);
  }//Now we can loop doing R then L
  while(signe(gel(q,1))==signe(gel(q,3))){//A,C still positive
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,3),1)));//n=floor(-B/2C)
    gcoeff(Rmat,2,1)=n;
    bot=avma;//If we exit loop after L, this is required for the final gerepileallsp
    q=bqf_transR(q,n);
    mat=ZM_mul(mat,Rmat);
    if(signe(gel(q,1))!=signe(gel(q,3))) break;//Done! Don't go L anymore.
    n=gfloor(Qdivii(negi(gel(q,2)),shifti(gel(q,1),1)));//n=floor(-B/2A)
    gcoeff(Lmat,1,2)=n;
    bot=avma;
    q=bqf_transL(q,n);
    mat=ZM_mul(mat,Lmat);
    if(gc_needed(mid,1)) gerepileallsp(mid,bot,2,&q,&mat);//Garbage collection if memory needed.
  }
  gerepileallsp(top,bot,2,&q,&mat);//Garbage
  gel(rvec,1)=q;
  gel(rvec,2)=mat;
  return rvec;
}

//Given an integral 2x2 matrix mtx, outputs the primitive quadratic form corresponding to mtx(x)=x. Typically used when matdet(mtx)=1 and mtx is hyperbolic with positive trace.
GEN mat_toibqf(GEN mtx){
  pari_sp top = avma;
  GEN B=subii(gcoeff(mtx,2,2),gcoeff(mtx,1,1));//d-a if mtx=[a,b;c,d]
  GEN gc=gcdii(gcdii(gcoeff(mtx,2,1),B),gcoeff(mtx,1,2));//gc=gcd(c,d-a,b)
  GEN C=negi(gcoeff(mtx,1,2));//-b
  GEN rvec=cgetg(4,t_VEC);
  gel(rvec,1)=Qdivii(gcoeff(mtx,2,1),gc);
  gel(rvec,2)=Qdivii(B,gc);
  gel(rvec,3)=Qdivii(C,gc);
  return gerepileupto(top,rvec);
}

//mat_toibqf but checks for 2x2 integral matrix. Does not check for det=1 or hyperbolic
GEN mat_toibqf_tc(GEN mtx){
  intmatrix_check(mtx);
  return mat_toibqf(mtx);//No garbage
}



//EQUIVALENCE OF BQFs



//Given forms q1 q2 of the same discriminant D of sign Dsign and rootD=sqrt(D), finds if they are equivalent, and returns a transition matrix if tmat=1. If Dsign=-1, can pass rootD as null
GEN bqf_isequiv(GEN q1, GEN q2, GEN rootD, int Dsign, int tmat){
  if(Dsign==1){
    if(tmat==0) return ibqf_isequiv(q1,q2,rootD);
    return ibqf_isequiv_tmat(q1,q2,rootD);
  }
  if(tmat==0) return dbqf_isequiv(q1,q2);//negdisc no tmat
  return dbqf_isequiv_tmat(q1,q2);//negdisc and tmat
}

//The collector method for determining if q is equivalent to an element of S. Returns -1 if no, and the first index for which it is otherwise. Pass tmat=1 to also return the transition matrix as the second output. If D<0, rootD is not required and can be passed in as gen_0. Does NOT check that disc(w)=disc(q) for w in S, so if matrices in S may not obey this then you should really eliminate them first for optimal efficiency
GEN bqf_isequiv_set(GEN q, GEN S, GEN rootD, int Dsign, int tmat){
  if(Dsign==-1){
    if(tmat==1) return dbqf_isequiv_set_tmat(q,S);
    return stoi(dbqf_isequiv_set(q,S));
  }//Now D>0
  if(tmat==1) return ibqf_isequiv_set_byS_tmat(q,S,rootD);
  return stoi(ibqf_isequiv_set_byS(q,S,rootD));
}

//bqf_isequiv(_set) with typechecking.
GEN bqf_isequiv_tc(GEN q1, GEN q2, int tmat, long prec){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q1);
  if(typ(q2)!=t_VEC)  pari_err_TYPE("please supply a BQF or a set of BQFs",q2);
  if(typ(gel(q2,1))==t_VEC){//q2 is a set of BQFs
    GEN D2;
    long ncount=1;//Let's eliminate the non-D-discriminants
    int keep[lg(q2)-1];//To keep or not; initially all zeroes
    for(long i=1;i<lg(q2);++i){
      D2=bqf_checkdisc(gel(q2,i));
      if(equalii(D,D2)){ncount++;keep[i]=1;}
      else{keep[i]=0;}
    }
    if(ncount==1){avma=top;return gen_m1;}
    GEN S=cgetg(ncount,t_VEC);
    long pos=1;
    long translate[ncount];//Translate back to the originial indices
    for(long i=1;i<lg(q2);++i){//Putting the correct discriminant entries back; note S is NOT gerepile safe
      if(keep[i]==1){translate[pos]=i;gel(S,pos)=gel(q2,i);pos++;}
    }
    GEN result=bqf_isequiv_set(q1,S,gsqrt(D,prec),signe(D),tmat);
    if(gequal(result,gen_m1)){avma=top; return gen_m1;}//Nope
    if(tmat==0) return gerepileupto(top,stoi(translate[itos(result)]));
    GEN rvec=cgetg(3,t_VEC);
    gel(rvec,1)=stoi(translate[itos(gel(result,1))]);
    gel(rvec,2)=ZM_copy(gel(result,2));
    return gerepileupto(top,rvec);
  }
  GEN D2=bqf_checkdisc(q2);
  if(!equalii(D,D2)) return gen_0;//Not equivalent
  return gerepileupto(top,bqf_isequiv(q1,q2,gsqrt(D,prec),signe(D),tmat));
}

//bqf_isequiv for negative discriminants
GEN dbqf_isequiv(GEN q1, GEN q2){
  pari_sp top = avma;
  GEN q1red=dbqf_red(q1);
  GEN q2red=dbqf_red(q2);
  if(ZV_equal(q1red,q2red)){
    avma=top;
    return gen_1;//Equal
  }
  avma=top;
  return gen_0;//Not equal
}

//bqf_isequiv for negative discriminants, also returns transition matrix
GEN dbqf_isequiv_tmat(GEN q1, GEN q2){
  pari_sp top = avma;
  GEN q1red=dbqf_red_tmat(q1);
  GEN q2red=dbqf_red_tmat(q2);
  if(ZV_equal(gel(q1red,1),gel(q2red,1))) return gerepileupto(top,ZM_mul(gel(q1red,2),ZM_inv(gel(q2red,2),NULL)));//The det is 1, no need to track d.
  avma=top;//Not equal
  return gen_0;
}

//Returns the first index for which q is equivalent to, and -1 if no such index
long dbqf_isequiv_set(GEN q, GEN S){
  pari_sp top=avma;
  GEN qred=dbqf_red(q);
  GEN sred;
  for(long i=1;i<lg(S);++i){
    sred=dbqf_red(gel(S,i));
    if(ZV_equal(qred,sred)){avma=top;return i;}//Done
    cgiv(sred);//Don't need this memory slot
  }
  avma=top;
  return -1;
}

//Returns a length 2 vector with the first index being the first index for which qred (reduced) is equivalent to, and the second entry a transition matrix. Returns -1 if no such index
GEN dbqf_isequiv_set_tmat(GEN q, GEN S){
  pari_sp top=avma;
  GEN qred=dbqf_red_tmat(q);
  GEN sred;
  for(long i=1;i<lg(S);++i){
    sred=dbqf_red_tmat(gel(S,i));
    if(ZV_equal(gel(qred,1),gel(sred,1))){
      GEN smatinv=ZM_inv(gel(sred,2),NULL);
      GEN rvec=cgetg(3,t_VEC);
      gel(rvec,1)=stoi(i);
      gel(rvec,2)=ZM_mul(gel(qred,2),smatinv);
      return gerepileupto(top,rvec);
    }//Done
    cgiv(sred);//Don't need this memory slot
  }
  avma=top;
  return gen_m1;
}

//bqf_isequiv for positive discriminants
GEN ibqf_isequiv(GEN q1, GEN q2, GEN rootD){
  pari_sp top = avma, mid, bot;
  GEN q1red=ibqf_red(q1,rootD);
  GEN q2red=ibqf_red(q2,rootD);
  mid=avma;//For gc_needed
  GEN q1seek=q1red;
  GEN q2seek=q2red;
  if(ZV_equal(q1seek,q2seek)){//Equal
    avma=top;
    return gen_1;
  }
  for(;;){//seek varaiables get hit with ibqf_rightnbr until loop or hit q1red,q2red. Should have approximately the same average time usage as only doing this with q1, has lower variance, more memory needed, but (crucially) will produce matrices with smaller coefficients, which is why we use this approach.
    bot=avma;
    q1seek=ibqf_rightnbr(q1seek,rootD);
    if(ZV_equal(q1seek,q1red)){//Loop completed, not equal :(
      avma=top;
      return gen_0;
    }
    if(ZV_equal(q1seek,q2red)){//Loop completed, equal
      avma=top;
      return gen_1;
    }
    q2seek=ibqf_rightnbr(q2seek,rootD);//Now doing q2
    if(ZV_equal(q2seek,q2red)){//Loop completed, not equal :(
      avma=top;
      return gen_0;
    }
    if(ZV_equal(q2seek,q1red)){//Loop completed, equal
      avma=top;
      return gen_1;
    }
    if(gc_needed(mid,1)) gerepileallsp(mid,bot,2,&q1seek,&q2seek);//Garbage collection if memory needed.
    }
}

//bqf_isequiv for positive discriminants, also returns transition matrix
GEN ibqf_isequiv_tmat(GEN q1, GEN q2, GEN rootD){
  pari_sp top = avma;
  GEN q1red=ibqf_red_pos_tmat(q1,rootD);
  GEN q2red=ibqf_red_pos_tmat(q2,rootD);
  pari_sp mid=avma, bot;//For gc_needed
  GEN q1R=q1red, q1L=q1red;//Going left and right from q1
  if(ZV_equal(gel(q1red,1),gel(q2red,1))) return gerepileupto(top,ZM_mul(gel(q1red,2),ZM_inv(gel(q2red,2),NULL)));
  for(;;){//go right with q1R and left with q1L until loop or hit q2red. Should have approximately the same average time usage as only going right and take more memory, but has lower variance, and (crucially) will produce matrices with smaller coefficients, which is why we use this approach.
    bot=avma;
    q1R=ibqf_rightnbr_update(ibqf_rightnbr_update(q1R,rootD),rootD);//Right twice
    if(ZV_equal(gel(q1R,1),gel(q2red,1))) return gerepileupto(top,ZM_mul(gel(q1R,2),ZM_inv(gel(q2red,2),NULL)));//Done!
    if(ZV_equal(gel(q1R,1),gel(q1L,1))){avma=top;return gen_0;}//Loop completed, not equal :(
    q1L=ibqf_leftnbr_update(ibqf_leftnbr_update(q1L,rootD),rootD);//Left twice
    if(ZV_equal(gel(q1L,1),gel(q2red,1))) return gerepileupto(top,ZM_mul(gel(q1L,2),ZM_inv(gel(q2red,2),NULL)));//Done!
    if(ZV_equal(gel(q1R,1),gel(q1L,1))){avma=top;return gen_0;}//Loop completed, not equal :(
    if(gc_needed(mid,1)) gerepileallsp(mid,bot,2,&q1R,&q1L);//Garbage collection if memory needed.
  }
}

//ibqf_isequiv_set where we find the reduced orbit of q, and search the reductions of S in it.
long ibqf_isequiv_set_byq(GEN q, GEN S, GEN rootD){
  pari_sp top=avma;
  GEN qredorb=ibqf_redorbit_posonly(q,rootD);
  gen_sort_inplace(qredorb,NULL,&bqf_compare,NULL);//Sort
  long rind=ibqf_isequiv_set_byq_presorted(qredorb, S, rootD);
  avma=top;
  return rind;
}

//ibqf_isequiv_set with qcode coming presorted.
long ibqf_isequiv_set_byq_presorted(GEN qredsorted, GEN S, GEN rootD){
  pari_sp top=avma;
  GEN Sred;
  long ind;
  for(long i=1;i<lg(S);++i){
    Sred=ibqf_red_pos(gel(S, i), rootD);
    ind=gen_search_temp(qredsorted, Sred, NULL, &bqf_compare);
    if(ind>0){avma=top;return i;}
    cgiv(Sred);
  }
  return -1;//No garbage!
}

//ibqf_isequiv_set while also returning the tmat.
GEN ibqf_isequiv_set_byq_tmat(GEN q, GEN S, GEN rootD){
  pari_sp top=avma;
  GEN qredorb=ibqf_redorbit_posonly_tmat(q,rootD);
  gen_sort_inplace(qredorb,NULL,&bqf_compare_tmat,NULL);
  return gerepileupto(top,ibqf_isequiv_set_byq_tmat_presorted(qredorb, S, rootD));
}

//ibqf_isequiv_set_byq_tmat with presorting.
GEN ibqf_isequiv_set_byq_tmat_presorted(GEN qredsorted, GEN S, GEN rootD){
  pari_sp top=avma;
  GEN Sred;
  long ind;
  for(long i=1;i<lg(S);++i){
    Sred=ibqf_red_pos_tmat(gel(S, i), rootD);
    ind=gen_search_temp(qredsorted, Sred, NULL, &bqf_compare_tmat);
    if(ind>0){
      GEN Sredinv=ZM_inv(gel(Sred, 2), NULL);
      GEN rvec=cgetg(3, t_VEC);
      gel(rvec, 1)=stoi(i);
      gel(rvec, 2)=ZM_mul(gel(gel(qredsorted, ind), 2), Sredinv);
      return gerepileupto(top, rvec);
    }
    cgiv(Sred);
  }
  return gen_m1;
}

//ibqf_isequiv_set where we reduce S and run through neighbours of q checking for equivalence.
long ibqf_isequiv_set_byS(GEN q, GEN S, GEN rootD){
  pari_sp top=avma;
  long lx;
  GEN Scodes=cgetg_copy(S,&lx);
  for(long i=1;i<lx;++i) gel(Scodes,i)=ibqf_red_pos(gel(S,i),rootD);
  GEN perm;
  gen_sort_inplace(Scodes,NULL,&bqf_compare,&perm);//Sort
  long rind=ibqf_isequiv_set_byS_presorted(q, Scodes, perm, rootD);
  avma=top;
  return rind;
}

//ibqf_isequiv_set_byS with Scodes sorted by the permutation perm^(-1)
long ibqf_isequiv_set_byS_presorted(GEN q, GEN Sreds, GEN perm, GEN rootD){
  pari_sp top=avma;
  GEN qred=ibqf_red_pos(q, rootD);
  GEN qseek=qred;
  long ind;
  do{
    ind=gen_search_temp(Sreds, qseek, NULL, &bqf_compare);
    if(ind>0){avma=top; return perm[ind];}//Done!
    qseek=ibqf_rightnbr(ibqf_rightnbr(qseek, rootD), rootD);
  }
  while(!ZV_equal(qred, qseek));
  avma=top;
  return -1;//No solution
}

//ibqf_isequiv_set_byS with tmat
GEN ibqf_isequiv_set_byS_tmat(GEN q, GEN S, GEN rootD){
  pari_sp top=avma;
  long lx;
  GEN Sreds=cgetg_copy(S,&lx);
  for(long i=1;i<lx;++i) gel(Sreds,i)=ibqf_red_pos_tmat(gel(S,i),rootD);
  GEN perm;
  gen_sort_inplace(Sreds,NULL,&bqf_compare_tmat,&perm);//Sort
  return gerepileupto(top,ibqf_isequiv_set_byS_tmat_presorted(q, Sreds, perm, rootD));
}

//ibqf_isequiv_set_byS_tmat with presorting
GEN ibqf_isequiv_set_byS_tmat_presorted(GEN q, GEN Sreds, GEN perm, GEN rootD){
  pari_sp top=avma;
  GEN qred=ibqf_red_pos_tmat(q,rootD);
  GEN q1L=qred, q1R=qred;
  long ind;
  int isdone=0;
  ind=gen_search_temp(Sreds,q1R,NULL,&bqf_compare_tmat);
  if(ind>0) isdone=1;
  else{
    for(;;){//At start of each loop, both q1L and q1R have been checked already
      q1R=ibqf_rightnbr_update(ibqf_rightnbr_update(q1R,rootD),rootD);//Update q1R
      if(ZV_equal(gel(q1L,1),gel(q1R,1))) break;//Done, no solution
      ind=gen_search_temp(Sreds,q1R,NULL,&bqf_compare_tmat);
      if(ind>0){isdone=1;break;}//Finished with R
      q1L=ibqf_leftnbr_update(ibqf_leftnbr_update(q1L,rootD),rootD);//Update q1L
      if(ZV_equal(gel(q1L,1),gel(q1R,1))) break;//Done, no solution
      ind=gen_search_temp(Sreds,q1L,NULL,&bqf_compare_tmat);
      if(ind>0){isdone=-1;break;}//Finished with L
    }
  }
  if(isdone==0){avma=top;return gen_m1;}//No solution
  GEN M=ZM_inv(gel(gel(Sreds,ind),2),NULL);
  GEN rvec=cgetg(3,t_VEC);
  gel(rvec,1)=stoi(perm[ind]);
  if(isdone==1) gel(rvec,2)=ZM_mul(gel(q1R,2),M);//Finished with R
  else gel(rvec,2)=ZM_mul(gel(q1L,2),M);//Finished with L
  return gerepileupto(top,rvec);
} 



//CLASS GROUPS AND COMPOSITION OF FORMS



//Computes all equiv classes forms of disc D. If prim=0, we include non-primitive forms. If GL=1, we take up to GL-equivalence. This relies on bqf_ncgp_lexic.
GEN bqf_allforms(GEN D, int prim, int GL, long prec){
  pari_sp top=avma;
  GEN forms;
  if(prim) forms=gel(bqf_ncgp_lexic(D, prec), 3);//Just primitive
  else{
    GEN discs=discsuperorders(D);//Possible discriminants
    long ld;
    forms=cgetg_copy(discs, &ld);
    for(long i=1;i<ld;i++){
      gel(forms, i)=gel(bqf_ncgp_lexic(gel(discs, i), prec), 3);//The forms.
      GEN scale=sqrti(diviiexact(D, gel(discs, i)));//What we must scale the forms by.
      if(!equali1(scale)){//If scale=1, no scaling!
        for(long j=1;j<lg(gel(forms, i));j++){
          for(int k=1;k<=3;k++) gmael3(forms, i, j, k)=mulii(gmael3(forms, i, j, k), scale);//Scale up!
        }
      }
    }
    forms=shallowconcat1(forms);
  }
  GEN newforms=forms;
  if(GL){
    if(signe(D)==1) pari_warn(warner, "GL not yet implemented for D>0");
    else{
      long lf=lg(forms);
      newforms=vectrunc_init(lf);
      for(long i=1;i<lf;i++){
        if(signe(gmael(forms, i, 2))!=-1) vectrunc_append(newforms, gel(forms, i));//Only keep B>=0.
      }
    }
  }
  return gerepilecopy(top, newforms);
}

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




//GENERAL CHECKING METHODS



//Checks that an input is an integral bqf (does NOT check discriminant!=square), and produces a pari_error if not. Used for not getting segmentation faults from a GP interface user. 
void bqf_check(GEN q){
  if(typ(q)!=t_VEC) pari_err_TYPE("Please input a length 3 integral VECTOR", q);
  if(lg(q)!=4) pari_err_TYPE("Please input a LENGTH 3 integral vector", q);
  for(int i=1;i<=3;i++) if(typ(gel(q, i))!=t_INT) pari_err_TYPE("Please input a length 3 INTEGRAL vector", q);
}

//Checks that an input is an integral bqf with disc!=square and produces a pari_error if not. Returns the discriminant
GEN bqf_checkdisc(GEN q){
  bqf_check(q);
  GEN D=bqf_disc(q);
  if(isdisc(D)) return D;
  pari_err(e_TYPE,"bqf: please input a bqf with non-square discriminant", q);
  return gen_0;
}

//Checks that an input is an integral 2x2 matrix, and produces a pari_error if not. Used for not getting segmentation faults from a GP interface user. 
void intmatrix_check(GEN mtx){
  if(typ(mtx)!=t_MAT) pari_err(e_TYPE, "Please input a 2x2 intgral MATRIX", mtx);
  if(lg(mtx)!=3) pari_err(e_TYPE, "Please input a 2x2 intgral matrix", mtx);
  if(lg(gel(mtx, 1))!=3) pari_err(e_TYPE, "Please input a 2x2 intgral matrix", mtx);
  for(int i=1;i<=2;i++) for(int j=1;j<=2;j++) if(typ(gcoeff(mtx, i, j))!=t_INT) pari_err(e_TYPE, "Please input a 2x2 INTEGRAL matrix", mtx);
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


