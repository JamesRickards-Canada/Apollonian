//IMPORTED FROM Q-QUADRATIC
//This is a collection of methods dealing with primitive integral binary quadratic forms.

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "apol.h"
#endif

//STATIC DECLARATIONS

//COMPOSITION/CLASS GROUP
static GEN bqf_ncgp_nonfundnarrow(GEN cgp, GEN D, GEN rootD);

//REPRESENTATIONS OF NUMBERS BY BQFs
static GEN bqf_reps_all(GEN n);
static GEN bqf_reps_trivial(void);
static void dbqf_reps_proper(GEN qred, GEN D, GEN n, glist **sols, long *nsols, GEN f, int *terminate);
static void ibqf_reps_proper(GEN qorb, GEN D, GEN rootD, GEN n, glist **sols, long *nsols, GEN f, int *terminate);
static GEN bqf_reps_creatervec(glist *sols, glist *scale, llist *nsolslist, long *totnsols, long *count, int half);
static GEN bqf_reps_creatervec_proper(glist *sols, long nsols, int half);
static GEN bqf_reps_makeprimitive(GEN q, GEN *n);
static void bqf_reps_updatesolutions(glist **sols, long *nsols, GEN *a, GEN *b);

//MORE REPRESENTATION OF NUMBERS
static GEN bqf_bigreps_creatervecfin(GEN newsols, GEN a, GEN b, GEN disc);
static GEN bqf_bigreps_creatervecpos(GEN newsols, GEN a, GEN b, GEN disc);
static GEN bqf_bigreps_createrveclin(GEN newsols, GEN a, GEN b, GEN disc);
static GEN zbqf_bigreps(GEN q, GEN n);
static GEN bqf_linearsolve_zall(GEN yzsols, GEN n2, GEN Minv);
static GEN bqf_linearsolve_zfin(GEN yzsols, GEN n2, GEN Minv);
static GEN bqf_linearsolve_zlin(GEN yzsols, GEN n2, GEN Minv);
static GEN bqf_linearsolve_zpos(GEN yzsols, GEN n2, GEN Minv, GEN M);
static GEN bqf_linearsolve_zquad(GEN yzsols, GEN n2, GEN Minv);



//DISCRIMINANT METHODS



//Generates list of discriminants from D1 to D2, can specify if they are fundamental and coprime to a given input.
GEN disclist(GEN D1, GEN D2, int fund, GEN cop){
  pari_sp ltop = avma;
  if (typ(D1) != t_INT) pari_err_TYPE("disclist",D1);
  if (typ(D2) != t_INT) pari_err_TYPE("disclist",D2);
  if (typ(cop) != t_INT) pari_err_TYPE("disclist",cop);
  glist *S=NULL;//Pointer to the list start
  long count=0;//Counts how many
  if(fund==0){
    if(gequal0(cop)){
      GEN D = gen_0; //int
      for(D=icopy(D1); cmpii(D, D2) <= 0; D = addis(D, 1)){//Need to icopy as we give back the memory space for D
        if(isdisc(D)){
          glist_putstart(&S,D);
          count++;
        }
        else cgiv(D);//This is OK as addis(D,1) will still work, it will just overwrite the memory location which is OK.
      }
    }
    else{
      GEN D = gen_0;//int
      GEN gc=gen_0;
      for(D=icopy(D1); cmpii(D, D2) <= 0; D = addis(D, 1)){
        gc=gcdii(cop, D);
        if(equali1(gc) && isdisc(D)){
          cgiv(gc);
          glist_putstart(&S,D);
          count++;
        }
        else{
          cgiv(gc);
          cgiv(D);
        }
      }
    }
  }
  else{
    if(gequal0(cop)){
      GEN D=gen_0;
      GEN fD=gen_0;
      for(D=icopy(D1); cmpii(D, D2) <= 0; D = addis(D, 1)){
		if(isdisc(D)==0){cgiv(D);continue;}
        fD=coredisc(D);
        if(equalii(fD, D)){
          cgiv(fD);
          glist_putstart(&S,D);
          count++;
        }
        else{
          cgiv(fD);
          cgiv(D);
        }
      }
    }
    else{
      GEN D= gen_0;
      GEN fD=gen_0;
      GEN gc=gen_0;
      for(D=icopy(D1); cmpii(D, D2) <= 0; D = addis(D, 1)){
        gc=gcdii(cop,D);
        if(equali1(gc)){
            cgiv(gc);
			if(isdisc(D)==0){cgiv(D);continue;}
            fD=coredisc(D);
            if(equalii(fD, D)){
              cgiv(fD);
              glist_putstart(&S,D);
              count++;
            }
            else{
              cgiv(fD);
              cgiv(D);
            }
        }
        else{
          cgiv(gc);
          cgiv(D);
        }
      }
    }
  }
  GEN Svec=glist_togvec(S, count, -1);
  Svec = gerepileupto(ltop, Svec);
  return Svec;
}

//Generate the list of primes dividing D for which D/p^2 is a discriminant, can pass in facs=factorization of D
GEN discprimeindex(GEN D, GEN facs){
  pari_sp top=avma;
  if(!isdisc(D)) pari_err_TYPE("discprimeindex, not a discriminant", D);
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

//Returns minimal positive solution [T,U] to T^2-DU^2=4
GEN pell(GEN D){
  pari_sp top=avma;
  if(!isdisc(D)) pari_err_TYPE("Please input a positive discriminant", D);
  if(signe(D)==-1) pari_err_TYPE("Please input a positive discriminant", D);
  GEN u=quadunit0(D, -1);
  if(gequalm1(quadnorm(u))) u=gsqr(u);//We need the smallest unit with norm 1, not -1
  if(smodis(D, 2)==1){//D is odd
    GEN a=shifti(real_i(u), 1);//a,b are guarenteed integers since the input is a t_QUAD, and a/2+b*omega=u, omega=(D%2+sqrt(D))/2
    GEN rvec=cgetg(3, t_VEC);
    gel(rvec, 2)=gimag(u);//Must copy
    gel(rvec, 1)=addii(a, gel(rvec, 2));
    return gerepileupto(top, rvec);
  }
  else{//D is even
    GEN a=real_i(u);
    GEN rvec=cgetg(3, t_VEC);
    gel(rvec, 1)=shifti(a, 1);
    gel(rvec, 2)=gimag(u);//Must copy
    return gerepileupto(top, rvec);
  }
}

//Returns log(epsilon(D)), where epsilon(D) is the fundamental unit of positive norm.
GEN posreg(GEN D, long prec){
  pari_sp top = avma;
  if(!isdisc(D)) pari_err_TYPE("posreg, not a positive discriminant", D);
  if(signe(D)==-1) pari_err_TYPE("posreg, not a positive discriminant", D);
  if(gequalm1(quadnorm(quadunit(D)))) return gerepileupto(top, gmulsg(2, quadregulator(D, prec)));
  return gerepileupto(top, quadregulator(D, prec));
}

//Returns sqrt(D) of type t_QUAD
GEN quadroot(GEN D){
  pari_sp top = avma;
  if(typ(D)!=t_INT) pari_err_TYPE("Please enter a non-square INTEGER", D);
  if(Z_issquare(D)) pari_err_TYPE("Please enter a NON-SQUARE integer", D);
  return gerepileupto(top, quadgen0(shifti(D, 2), -1));//Just do quadgen on 4D
}



//BASIC OPERATIONS ON BINARY QUADRATIC FORMS



//Returns the generator of the automorphism group of q in PSL(2,Z) (which is Z if D>0, Z/2Z if D=-4, Z/3Z if D=-3, and 1 if D<-4
GEN bqf_automorph_tc(GEN q){
  pari_sp top=avma;
  GEN D=bqf_checkdisc(q);
  if(signe(D)==-1) return gerepileupto(top, dbqf_automorph(q, D));
  return gerepileupto(top, ibqf_automorph_D(q, D));
}

//Compares two bqfs based on A, then B, then C. It is safe to modify this method, as sorting/searching will call the pointer to this function
int bqf_compare(void *data, GEN q1, GEN q2){
  int i;
  for(long j=1;j<=3;++j){
    i=cmpii(gel(q1,j),gel(q2,j));
    switch(i){
	  case -1:
	    return -1;
	  case 1:
	    return 1;
    }
  }
  return 0;
}

//bqf_compare, but assumes the inputs are d1=[q1,mat1] and d2=[q2,mat2], and compares q1 and q2. This returns 0 if q1=q2, even if mat1!=mat2 (the mats don't matter)
int bqf_compare_tmat(void *data, GEN d1, GEN d2){
  return bqf_compare(data,gel(d1,1),gel(d2,1));
}

//Returns the discriminant of q, which must be length 3 vector and integral.
GEN bqf_disc(GEN q){
  pari_sp top = avma;
  return gerepileupto(top,subii(sqri(gel(q,2)),shifti(mulii(gel(q,1),gel(q,3)),2)));
}

//Returns the discriminant of q after performing a type check.
GEN bqf_disc_tc(GEN q){
  bqf_check(q);
  return bqf_disc(q);//No avma necessary
}

//Returns 1 if the BQF q (indefinite/positive definite) of discriminant D with sign Dsign is reduced, and 0 if not.
int bqf_isreduced(GEN q, int Dsign){
  pari_sp top=avma;
  int answer=0;
  if(Dsign==1){//D>0, so AC<0 and B>|A+C|
    if(signe(gel(q,1))==signe(gel(q,3))) return 0;
    GEN s=addii(gel(q,1),gel(q,3));
    GEN abss=absi(s);
    if(cmpii(gel(q,2),abss)==1) answer=1;//Only need to update it if this happens
    avma=top;
    return answer;
  }
  else{//D<0, so |B|<=A<=C and if A=|B| or A=C then B>=0
    if(cmpii(gel(q,1),gel(q,3))<=0){
      GEN babs=absi(gel(q,2));
      if(cmpii(babs,gel(q,1))<=0){
        if(equalii(gel(q,1),babs) || equalii(gel(q,1),gel(q,3))){
          if(cmpis(gel(q,2),0)>=0) answer=1;
        }
        else answer=1;
      }
    }
    avma=top;
    return answer;
  }
}

//bqf_isreduced with typecheck
int bqf_isreduced_tc(GEN q){
  pari_sp top=avma;
  bqf_check(q);
  GEN D=bqf_disc(q);
  int ans=bqf_isreduced(q,signe(D));
  avma=top;
  return ans;
}

//Generates a random proper bqf with max coefficient maxc. If type=1 it will be indefinite, type=-1 positive definite, type=0 either. primitive=1 means primitive, =0 means don't care. This is not designed for efficiency
GEN bqf_random(GEN maxc, int type, int primitive){
  pari_sp top=avma;
  setrand(getwalltime());
  if(typ(maxc)!=t_INT) return gen_0;
  if(signe(maxc)!=1) return gen_0;//Just making sure maxc is a positive integer.
  GEN q=cgetg(4,t_VEC);
  GEN A, B, C, D;
  for(;;){
    A=randomi(maxc);
	if(signe(randomi(gen_2))==0) togglesign_safe(&A);
	B=randomi(maxc);
	if(signe(randomi(gen_2))==0) togglesign_safe(&B);
	C=randomi(maxc);
	if(signe(randomi(gen_2))==0) togglesign_safe(&C);
	gel(q,1)=A;
	gel(q,2)=B;
	gel(q,3)=C;
	D=bqf_disc(q);
	if(type==1){
	  if(primitive==1){if(!equali1(ZV_content(q))) continue;}// Not primitive
	  if(signe(D)==1 && isdisc(D)==1) return gerepilecopy(top,q);
	}
	else if(type==-1){
	  if(primitive==1){if(!equali1(ZV_content(q))) continue;}// Not primitive
	  if(signe(D)==-1){
		if(signe(A)==-1) ZV_togglesign(q);//Making positive definite
		return gerepilecopy(top,q);
	  }
	}
	else{
	  if(primitive==1){if(!equali1(ZV_content(q))) continue;}// Not primitive
	  if(signe(D)==-1){
		if(signe(A)==-1) ZV_togglesign(q);//Making positive definite
		return gerepilecopy(top,q);
	  }
	  if(isdisc(D)==1) return gerepilecopy(top,q);
	}
  }
}

//Generates a random bqf of discriminant D with |B|<=2maxc and primitive. This is not designed for efficiency
GEN bqf_random_D(GEN maxc, GEN D){
  pari_sp top=avma;
  //setrand(getwalltime());
  if(typ(maxc)!=t_INT) return gen_0;
  if(signe(maxc)!=1) return gen_0;//Just making sure maxc is a positive integer.
  if(!isdisc(D)) return gen_0;
  GEN B=randomi(maxc);
  B=shifti(B,1);
  if(smodis(D,2)==1) B=addii(B,gen_1);
  GEN AC=shifti(subii(sqri(B),D),-2);
  if(signe(randomi(gen_2))==0) B=negi(B);
  GEN g=ggcd(D,B);
  GEN Aposs=divisors(AC), A, C;
  long r, lx=lg(Aposs);
  for(;;){
    r=itos(randomi(stoi(lx-1)))+1;
	A=gel(Aposs,r);
	C=gel(Aposs,lx-r);
	if(equali1(ggcd(ggcd(A,C),g))) break;
  }
  if(signe(D)==1 && signe(AC)==-1){
	if(signe(randomi(gen_2))==0) A=negi(A);
	else togglesign_safe(&C);
  }
  else if(signe(D)==1){
	if(signe(randomi(gen_2))==0){A=negi(A);C=negi(C);}
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
  if(bqf_isreduced(q,-1)) return ZV_copy(q);//No garbage needed
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
  if(bqf_isreduced(q,-1)){
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
  if(bqf_isreduced(toriv,1)) return toriv;//No garbage!
  return gerepilecopy(top,ibqf_rightnbr(toriv,rootD));//If not reduced, take the right neighbour
}

//bqf_red for positive discriminants, also returns transition matrix.
GEN ibqf_red_tmat(GEN q, GEN rootD){
  pari_sp top=avma;
  GEN toriv=ibqf_toriver_tmat(q,rootD);//Now the form is on the river
  if(bqf_isreduced(gel(toriv,1),1)) return toriv;
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
	Sred=ibqf_red_pos(gel(S,i), rootD);
	ind=gen_search(qredsorted,Sred,0,NULL,&bqf_compare);
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
	Sred=ibqf_red_pos_tmat(gel(S,i), rootD);
	ind=gen_search(qredsorted,Sred,0,NULL,&bqf_compare_tmat);
	if(ind>0){
	  GEN Sredinv=ZM_inv(gel(Sred,2),NULL);
	  GEN rvec=cgetg(3,t_VEC);
	  gel(rvec,1)=stoi(i);
	  gel(rvec,2)=ZM_mul(gel(gel(qredsorted,ind),2),Sredinv);
	  return gerepileupto(top,rvec);
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
  GEN qred=ibqf_red_pos(q,rootD);
  GEN qseek=qred;
  long ind;
  do{
	ind=gen_search(Sreds,qseek,0,NULL,&bqf_compare);
	if(ind>0){avma=top; return perm[ind];}//Done!
	qseek=ibqf_rightnbr(ibqf_rightnbr(qseek,rootD),rootD);
  }
  while(!ZV_equal(qred,qseek));
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
  ind=gen_search(Sreds,q1R,0,NULL,&bqf_compare_tmat);
  if(ind>0) isdone=1;
  else{
    for(;;){//At start of each loop, both q1L and q1R have been checked already
	  q1R=ibqf_rightnbr_update(ibqf_rightnbr_update(q1R,rootD),rootD);//Update q1R
	  if(ZV_equal(gel(q1L,1),gel(q1R,1))) break;//Done, no solution
	  ind=gen_search(Sreds,q1R,0,NULL,&bqf_compare_tmat);
	  if(ind>0){isdone=1;break;}//Finished with R
	  q1L=ibqf_leftnbr_update(ibqf_leftnbr_update(q1L,rootD),rootD);//Update q1L
	  if(ZV_equal(gel(q1L,1),gel(q1R,1))) break;//Done, no solution
	  ind=gen_search(Sreds,q1L,0,NULL,&bqf_compare_tmat);
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



//Composes q1, q2. Code adapted from "static void qfb_comp(GEN z, GEN x, GEN y)" in Qfb.c
GEN bqf_comp(GEN q1, GEN q2){
  pari_sp top=avma;
  GEN n, c, d, y1, v1, v2, c3, m, p1, r;
  if(ZV_equal(q1,q2)) return bqf_square(q1);
  n=shifti(subii(gel(q2,2),gel(q1,2)), -1);
  v1=gel(q1,1);
  v2=gel(q2,1);
  c =gel(q2,3);
  d =bezout(v2,v1,&y1,NULL);
  if (equali1(d)) m = mulii(y1,n);
  else{
    GEN s = subii(gel(q2,2), n);
    GEN x2, y2, d1 = bezout(s,d,&x2,&y2); /* x2 s + y2 (x1 v1 + y1 v2) = d1 */
    if(!equali1(d1)){
      v1 = diviiexact(v1,d1);
      v2 = diviiexact(v2,d1); // gcd = 1 iff q1 or q2 primitive
      v1 = mulii(v1, gcdii(c,gcdii(gel(q1,3),gcdii(d1,n))));
      c = mulii(c, d1);
    }
    m = addii(mulii(mulii(y1,y2),n), mulii(gel(q2,3),x2));
  }
  togglesign(m);
  r = modii(m, v1);
  p1 = mulii(r, v2);
  c3 = addii(c, mulii(r,addii(gel(q2,2),p1)));
  p1=shifti(p1,1);
  GEN z=cgetg(4,t_VEC);
  gel(z,1) = mulii(v1,v2);
  gel(z,2) = addii(gel(q2,2), p1);
  gel(z,3) = diviiexact(c3,v1); 
  return gerepileupto(top,z);
}

//bqf_comp with reduction. If D<0, can pass rootD as NULL
GEN bqf_comp_red(GEN q1, GEN q2, GEN rootD, int Dsign){
  pari_sp top=avma;
  GEN qcomp=bqf_comp(q1,q2);
  return gerepileupto(top,bqf_red(qcomp,rootD,Dsign,0));
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
    GEN q=cgetg(4,t_VEC);
    gel(q,1)=gen_1;
    gel(q,2)=gen_0;
    gel(q,3)=shifti(nD,-2);//-D/4
    return gerepileupto(top,q);
  }//Now D is odd
  GEN omD=subii(gen_1,D);
  GEN q=cgetg(4,t_VEC);
  gel(q,1)=gen_1;
  gel(q,2)=gen_1;
  gel(q,3)=shifti(omD,-2);//(1-D)/4
  return gerepileupto(top,q);
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
      GEN rvec=cgetg(4,t_VEC);
      gel(rvec,1)=gen_1;
      gel(rvec,2)=mkvecsmall(1);
      gel(rvec,3)=cgetg(2,t_VEC);
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
    gel(rvec, 3)=cgetg_copy(gel(gpclgp,3),&lx);//Generators
    for(long i=1;i<=ngens;++i){//Inputting backwards, as bnfnarrow does the reverse order to what we want.
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
    gel(rvec, 3)=cgetg(2,t_VEC);
    gmael(rvec, 3, 1)=bqf_idelt(D);
    return gerepileupto(top, rvec);
  }
  newgens=cgetg_copy(gel(FULLgp, 3), &lx);
  for(long i=1;i<lg(gel(FULLgp, 3));i++){
	temp=gtovec(gel(gel(FULLgp, 3),i));
	setlg(temp, 4);//Truncate
	gel(newgens, i)=bqf_red(temp, rootD, signe(D), 0);
  }
  GEN rvec=cgetg(4, t_VEC);
  gel(rvec, 1)=icopy(gel(FULLgp,1));//Size
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
    gel(rvec, 3)=cgetg(2,t_VEC);
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
	gel(rvec, 1)=shifti(gel(cgp, 1),1);
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
  gel(rvec, 3)=cgetg_copy(newgens,&lx);
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
  if(gequal0(n)){GEN D=bqf_disc(q);return gerepileupto(top,bqf_idelt(D));}
  GEN qpow=ZV_copy(q);
  if(signe(n)==-1){
	togglesign_safe(&gel(qpow,2));//Taking the inverse
	if(equalis(n,-1)) return qpow;
	q=qpow;//Need to repoint q to here
	n=negi(n);
  }
  else if(equali1(n)) return qpow;
  GEN nbin=binary_zv(n);
  for(long i=2;i<lg(nbin);i++){
	qpow=bqf_square(qpow);
	if(nbin[i]==1) qpow=bqf_comp(qpow,q);
  }
  return gerepileupto(top,qpow);
}

//q^n with reduction
GEN bqf_pow_red(GEN q, GEN n, GEN rootD, int Dsign){
  pari_sp top=avma;
  if(gequal0(n)){GEN D=bqf_disc(q);return gerepileupto(top,bqf_idelt(D));}
  GEN qpow=ZV_copy(q);
  if(signe(n)==-1){
	togglesign_safe(&gel(qpow,2));//Taking the inverse
	if(equalis(n,-1)) return gerepileupto(top,bqf_red(qpow, rootD, Dsign, 0));
	q=qpow;//Need to repoint q to here
	n=negi(n);
  }
  else if(equali1(n)) return gerepileupto(top,bqf_red(qpow, rootD, Dsign, 0));
  GEN nbin=binary_zv(n);
  for(long i=2;i<lg(nbin);i++){
	qpow=bqf_square_red(qpow, rootD, Dsign);
	if(nbin[i]==1) qpow=bqf_comp_red(qpow, q, rootD, Dsign);
  }
  return gerepileupto(top,qpow);
}

//bqf_pow_red with typecheck
GEN bqf_pow_tc(GEN q, GEN n, int tored, long prec){
  pari_sp top=avma;
  if(typ(n)!=t_INT) pari_err_TYPE("Please enter an integer n", n);
  if(tored==1){
    GEN D=bqf_checkdisc(q);
    if(signe(D)==1) return gerepileupto(top,bqf_pow_red(q, n, gsqrt(D,prec),1));
    avma=top;
    return bqf_pow_red(q, n, NULL, -1);
  }//Now no reduction
  bqf_check(q);
  return bqf_pow(q,n);
}

//Squares the bqf q; code adapted from the PARI function "static void qfb_sqr(GEN z, GEN x)" in Qfb.c
GEN bqf_square(GEN q){
  pari_sp top=avma;
  GEN c, d1, x2, v1, v2, c3, m, p1, r;
  d1=bezout(gel(q,2),gel(q,1),&x2, NULL);//d1=gcd(A,B)
  c=gel(q,3);
  m=mulii(c,x2);
  if(equali1(d1)) v1=v2=gel(q,1);
  else{
    v1=diviiexact(gel(q,1),d1);
    v2=mulii(v1,gcdii(d1,c)); // = v1 iff q primitive 
    c=mulii(c,d1);
  }
  togglesign(m);
  r = modii(m,v2);
  p1 = mulii(r, v1);
  c3 = addii(c, mulii(r,addii(gel(q,2),p1)));
  p1=shifti(p1,1);
  GEN z=cgetg(4,t_VEC);
  gel(z,1) = mulii(v1,v2);
  gel(z,2) = addii(gel(q,2), p1);
  gel(z,3) = diviiexact(c3,v2);
  return gerepileupto(top,z);
}

//bqf_square with reduction. If D<0, can pass rootD as NULL
GEN bqf_square_red(GEN q, GEN rootD, int Dsign){
  pari_sp top=avma;
  GEN qsqr=bqf_square(q);
  return gerepileupto(top,bqf_red(qsqr,rootD,Dsign,0));
}

//bqf_square_red with typecheck
GEN bqf_square_tc(GEN q, int tored, long prec){
  pari_sp top=avma;
  if(tored==1){
    GEN D=bqf_checkdisc(q);
    if(signe(D)==1) return gerepileupto(top,bqf_square_red(q,gsqrt(D,prec),1));
    avma=top;
    return bqf_square_red(q,NULL,-1);
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
    A=nfnorm(numf,alph2);
    C=nfnorm(numf,alph1);
  }
  else{
    A=nfnorm(numf,alph1);
    C=nfnorm(numf,alph2);
  }
  GEN B=gmul(polcoef_i(a1,0,varno),gen_2);//B/2 may be a half integer, so can't use shifti or mulii sadly
  togglesign_safe(&B);//Negate B in place
  GEN d=gcdii(gcdii(A,B),C);
  GEN Q=cgetg(4,t_VEC);
  if(equali1(d)){
    gel(Q,1)=icopy(A);
    gel(Q,2)=icopy(B);
    gel(Q,3)=icopy(C);
  }
  else{
    gel(Q,1)=diviiexact(A,d);
    gel(Q,2)=diviiexact(B,d);
    gel(Q,3)=diviiexact(C,d);
  }
  return gerepileupto(top,Q);
}



//REPRESENTATION OF NUMBERS BY BQFs



//Solves q(x,y)=n for x,y integers (inputs are integral quadratic form and integral n). proper=1 means gcd(x,y)=1 (only when disc!=square), and proper=0 means no such restriction. This method just calls the appropriate subroutine based on D. half=1 means only output one family between (x,y) and (-x,-y), and =0 means output both.
GEN bqf_reps(GEN q, GEN n, int proper, int half, long prec){
  pari_sp top=avma;
  if(isexactzero(q)) return bqf_reps_all(n);//q=0
  n=icopy(n);//Copy n
  q=bqf_reps_makeprimitive(q,&n);//Updating q and n; this also makes a copy of q so as to not affect it
  if(q==NULL){avma=top; return gen_0;}//No solution; gcd(q) did not divide n.
  GEN D=bqf_disc(q);
  if(signe(D)==-1){//Negative disc
    if(signe(gel(q,1))==-1){ZV_togglesign(q);togglesign_safe(&n);}//Negating q and n to make positive definite.
    GEN qred=dbqf_red_tmat(q);//The reduction of q and the transition matrix
	return gerepileupto(top,dbqf_reps(qred,D,n,proper,half));
  }
  if(signe(D)==0){//Disc=0
    if(signe(gel(q,1))==-1 || signe(gel(q,3))==-1){ZV_togglesign(q);togglesign_safe(&n);}//Making positive
	GEN A=sqrti(gel(q,1));
	GEN B=sqrti(gel(q,3));//must be squares
	if(signe(gel(q,2))==-1) togglesign_safe(&B);//Now (Ax+By)^2=n is to be solved
	return gerepileupto(top,zbqf_reps(A,B,n,half));
  }//Now D>0
  GEN rootD;
  long issquare=Z_issquareall(D,&rootD);//sets rootD to sqrt(D) if it is a square
  if(issquare) return gerepileupto(top,sbqf_reps(q,D,rootD,n,half));
  //Now D>0 is not a square
  rootD=gsqrt(D,prec);
  GEN qorb=ibqf_redorbit_posonly_tmat(q,rootD);
  gen_sort_inplace(qorb,NULL,&bqf_compare_tmat,NULL);//Sort qorb
  GEN qautom=ibqf_automorph_D(q,D);
  return gerepileupto(top,ibqf_reps(qorb,qautom,D,rootD,n,proper,half));
}

//bqf_reps with typechecking.
GEN bqf_reps_tc(GEN q, GEN n, int proper, int half, long prec){
  bqf_check(q);
  if(typ(n)!=t_INT) pari_err_TYPE("bqf_reps: please input an integer n",n);
  return bqf_reps(q,n,proper,half,prec);
}

//qred=dbqf_red_tmat(q), q has discriminant D, and is primitive AND positive definite. If you do not have qred/have not checked primitiveness, please use bqf_reps instead.
GEN dbqf_reps(GEN qred, GEN D, GEN n, int proper, int half){
  pari_sp top=avma;
  if(isexactzero(n)) return bqf_reps_trivial();//Checking n=0
  if(signe(n)==-1) return gen_0;//No solution when n<0 SINCE we are assumeing positive definite
  glist *sols=NULL;//The pairs of solutions
  long nsols=0;//Number of solutions
  GEN f=cgetg(4,t_VEC);//f will store [d,b,c] during our search, where d|n and n/d is a square (d=n only if proper=1)
  int terminate=0;//Only used when proper=0
  if(proper){//Only seek proper
    dbqf_reps_proper(qred,D,n,&sols,&nsols,f,&terminate);
    if(sols==NULL){avma=top;return gen_0;}//No solution
    GEN rvec=bqf_reps_creatervec_proper(sols, nsols, half);//Make return vector minus the first component.
    gel(rvec,1)=cgetg(2,t_VEC);
    gel(gel(rvec,1),1)=gen_0;//Finite
    return gerepileupto(top,rvec);
  }//Now we seek all solutions, not just proper ones
  GEN divs=divisors(n);//Divisors of n
  llist *nsolslist=NULL;//List of how many solutions we get each time
  glist *scales=NULL;//List of how much we must scale each solution by
  long oppi=lg(divs);
  long count=0;//Number of divisors which we try
  long totnsols=0;//Total number of solutions overall
  GEN thescale;
  for(long i=1;i<lg(divs);++i){
	oppi--;
	if(!Z_issquareall(gel(divs, oppi), &thescale)) continue;//Continue on
	count++;
	glist_putstart(&scales, thescale);//The scale
	nsols=0;//Resetting solution count
	dbqf_reps_proper(qred, D, gel(divs,i), &sols, &nsols, f, &terminate);//Doing the proper solve
	if(count==1 && terminate==1){avma=top;return gen_0;}//No solutions to the residue for the fundamental disc, so there will be none EVER!
	llist_putstart(&nsolslist, nsols);
	totnsols=totnsols+nsols;
  }//Now we have our solutions
  if(sols==NULL){avma=top;return gen_0;}//No solution
  GEN rvec=bqf_reps_creatervec(sols, scales, nsolslist, &totnsols, &count, half);//Initialize the solution vector
  gel(rvec,1)=cgetg(2, t_VEC);
  gel(gel(rvec,1),1)=gen_0;//Finite
  return gerepileupto(top,rvec);
}

//Assume q is primitive, updates sols and nsols with the solutions. Will set terminate to 0 if the residue test failed. nsols should be passed in as pointing to 0
static void dbqf_reps_proper(GEN qred, GEN D, GEN n, glist **sols, long *nsols, GEN f, int *terminate){
  GEN fourn=shifti(n, 2);//=4n
  GEN resclass=Zn_quad_roots(fourn, gen_0, gneg(D));//square roots of D modulo 4n
  if(resclass==NULL){*terminate=1;return;}//No solution.
  //We have residue classes, so SOME form of discriminant D represents n. We now have to check if such forms can be q.
  GEN modu=gel(resclass, 1);//The modulus
  //q properly represents n <=> q is similar to [n,b,c] for 0<=b<2|n|. For this to happen, b^2-4nc=D -> D==b^2 mod 4n, hence the residue condition. gel(resclass, 1)=modulus necessarily divides 2n.
  GEN B, fred, transmat;//fred will be the reduction of f, transmat circ f=[n,b,c]
  gel(f,1)=n;
  long kmax=itos(diviiexact(shifti(absi(n),1),modu))-1;//2|n|/modu-1
  for(long j=1;j<lg(gel(resclass, 2));++j){//Trying each residue
    B=gel(gel(resclass, 2),j);//Initially B=residue
    gel(f, 2)=B;
    gel(f, 3)=diviiexact(subii(sqri(B), D),fourn);//(B^2-D)/(4n)
    fred=dbqf_red_tmat(f);
    if(ZV_equal(gel(qred,1),gel(fred,1))){//Solution!
      transmat=ZM_mul(gel(qred,2),ZM_inv(gel(fred,2),NULL));//transmat circ f=[n,b,c]
	  bqf_reps_updatesolutions(sols, nsols, &gcoeff(transmat,1,1), &gcoeff(transmat,2,1));//Update!
    }
	else cgiv(fred);//May as well discard
    for(long k=1;k<=kmax;++k){//Go through the other possibilities
	  gel(f,3)=addii(gel(f,3),diviiexact(mulii(modu,addii(modu,shifti(B,1))),fourn));//f[3]=f[3]+modu(modu+2B)/(4n)
	  B=addii(B,modu);//So B=k*modu+resclass[2][j] in the end
	  gel(f,2)=B;
	  fred=dbqf_red_tmat(f);
      if(ZV_equal(gel(qred,1),gel(fred,1))){//Solution!
        transmat=ZM_mul(gel(qred,2),ZM_inv(gel(fred,2),NULL));//transmat circ f=[n,b,c]
	    bqf_reps_updatesolutions(sols, nsols, &gcoeff(transmat,1,1), &gcoeff(transmat,2,1));//Update!
      }
	  else cgiv(fred);//May as well discard
	}
  }
  if(equalis(D,-4)){//Must account for extra automorph
	GEN mat=ZM_mul(gel(qred,2),ZM_mul(mkmat22(gen_0,gen_1,gen_m1,gen_0),ZM_inv(gel(qred,2),NULL)));//Conjugating [0,1;-1,0] by gel(qred,2) to get the non-trivial autormorph of q.
	B=gcoeff(mat,1,2);
	gcoeff(mat,1,2)=gcoeff(mat,2,1);
	gcoeff(mat,2,1)=B;//Transposing mat
	glist *seeker=*sols;//This runs through sols, and we must add sol*mat each time.
	for(long i=1;i<=*nsols;i++){
	  transmat=ZV_ZM_mul(seeker->data,mat);
	  glist_putstart(sols,transmat);
      seeker=seeker->next;
	}
	*nsols=*nsols*2;//This doubled the solutions
  }
  else if(equalis(D,-3)){//Must account for extra automorphs (two of them)
	GEN mat=ZM_mul(gel(qred,2),ZM_mul(mkmat22(gen_1,gen_1,gen_m1,gen_0),ZM_inv(gel(qred,2),NULL)));//Conjugating [1,1;-1,0] by gel(qred,2) to get a non-trivial autormorph of q.
	B=gcoeff(mat,1,2);
	gcoeff(mat,1,2)=gcoeff(mat,2,1);
	gcoeff(mat,2,1)=B;//Transposing mat
	glist *seeker=*sols;//This runs through sols, and we must add sol*mat each time.
	for(long i=1;i<=*nsols;i++){
	  transmat=ZV_ZM_mul(seeker->data,mat);
	  glist_putstart(sols,transmat);
	  transmat=ZV_ZM_mul(transmat,mat);//The square must also be added
	  glist_putstart(sols,transmat);
      seeker=seeker->next;
	}
	*nsols=*nsols*3;//This tripled the solutions
  }
}

//qorb the sorted output of iqbf_redorbit_posonly_tmat, q has discriminant D, and is primitive. If you do not have qorb/have not checked primitiveness, please use bqf_reps instead. You also need to pass qautom as the automorph of q, though this method will just copy it for the return (it does not use it in any way) so if you don't really care you can just pass in [;], the 0 matrix
GEN ibqf_reps(GEN qorb, GEN qautom, GEN D, GEN rootD, GEN n, int proper, int half){
  pari_sp top=avma;
  if(isexactzero(n)) return bqf_reps_trivial();//Checking n=0
  glist *sols=NULL;//The pairs of solutions
  long nsols=0;//Number of solutions
  GEN f=cgetg(4,t_VEC);//f will store [d,b,c] during our search, where d|n and n/d is a square (d=n only if proper=1)
  int terminate=0;//Only used when proper=0
  if(proper){//Only seek proper
    ibqf_reps_proper(qorb,D,rootD,n,&sols,&nsols,f,&terminate);
    if(sols==NULL){avma=top;return gen_0;}//No solution
    GEN rvec=bqf_reps_creatervec_proper(sols, nsols, half);//Make return vector minus the first component.
    gel(rvec,1)=cgetg(3,t_VEC);
    gel(gel(rvec,1),1)=gen_1;//Positive
	gel(gel(rvec,1),2)=ZM_copy(qautom);
    return gerepileupto(top,rvec);
  }//Now we seek all solutions, not just proper ones
  GEN divs=divisors(n);//Divisors of n
  llist *nsolslist=NULL;//List of how many solutions we get each time
  glist *scales=NULL;//List of how much we must scale each solution by
  long oppi=lg(divs);
  long count=0;//Number of divisors which we try
  long totnsols=0;//Total number of solutions overall
  for(long i=1;i<lg(divs);++i){
	oppi--;
	if(!Z_issquare(gel(divs,oppi))) continue;//Continue on
	count++;
	glist_putstart(&scales,sqrti(gel(divs,oppi)));//The scale
	nsols=0;//Resetting solution count
	if(signe(n)==-1) ibqf_reps_proper(qorb,D,rootD,negi(gel(divs,i)),&sols,&nsols,f,&terminate);
	else ibqf_reps_proper(qorb,D,rootD,gel(divs,i),&sols,&nsols,f,&terminate);//Doing the proper solve
	if(count==1 && terminate==1){avma=top;return gen_0;}//No solutions to the residue for the fundamental disc, so there will be none EVER!
	llist_putstart(&nsolslist,nsols);
	totnsols=totnsols+nsols;
  }//Now we have our solutions
  if(sols==NULL){avma=top;return gen_0;}//No solution
  GEN rvec=bqf_reps_creatervec(sols,scales,nsolslist,&totnsols,&count,half);//Initialize the solution vector
  gel(rvec,1)=cgetg(3,t_VEC);
  gel(gel(rvec,1),1)=gen_1;//Positive
  gel(gel(rvec,1),2)=ZM_copy(qautom);
  return gerepileupto(top,rvec);
}

//Assume q is primitive, updates sols and nsols with the solutions. Will set terminate to 0 if the residue test failed.
static void ibqf_reps_proper(GEN qorb, GEN D, GEN rootD, GEN n, glist **sols, long *nsols, GEN f, int *terminate){
  GEN fourn=shifti(n, 2);//=4n
  GEN resclass=Zn_quad_roots(fourn, gen_0, gneg(D));//square roots of D modulo 4n
  if(resclass==NULL){*terminate=1;return;}//No solution.
  //We have residue classes, so SOME form of discriminant D represents n. We now have to check if such forms can be q.
  GEN modu=gel(resclass, 1);//The modulus
  //q properly represents n <=> q is similar to [n,b,c] for 0<=b<2|n|. For this to happen, b^2-4nc=D -> D==b^2 mod 4n, hence the residue condition. gel(resclass, 1)=modulus necessarily divides 2n.
  GEN B, fred, transmat;//fred will be the reduction of f, transmat circ f=[n,b,c]
  gel(f,1)=n;
  long kmax=itos(diviiexact(shifti(absi(n),1),modu))-1;//2|n|/modu-1
  long ind;
  for(long j=1;j<lg(gel(resclass, 2));++j){//Trying each residue
    B=gel(gel(resclass, 2), j);//Initially B=residue
    gel(f, 2)=B;
    gel(f, 3)=diviiexact(subii(sqri(B), D),fourn);//(B^2-D)/(4n)
    fred=ibqf_red_pos_tmat(f, rootD);
	ind=gen_search(qorb, fred, 0, NULL, &bqf_compare_tmat);
    if(ind>0){//Solution!
      transmat=ZM_mul(gel(gel(qorb,ind),2),ZM_inv(gel(fred,2),NULL));//transmat circ f=[n,b,c]
	  bqf_reps_updatesolutions(sols, nsols, &gcoeff(transmat,1,1), &gcoeff(transmat,2,1));//Update!
    }
	else cgiv(fred);//May as well discard
    for(long k=1;k<=kmax;++k){//Go through the other possibilities
	  gel(f,3)=addii(gel(f,3),diviiexact(mulii(modu,addii(modu,shifti(B,1))),fourn));//f[3]=f[3]+modu(modu+2B)/(4n)
	  B=addii(B,modu);//So B=k*modu+resclass[2][j] in the end
	  gel(f,2)=B;
	  fred=ibqf_red_pos_tmat(f,rootD);
	  ind=gen_search(qorb,fred,0,NULL,&bqf_compare_tmat);
      if(ind>0){//Solution!
        transmat=ZM_mul(gel(gel(qorb,ind),2),ZM_inv(gel(fred,2),NULL));//transmat circ f=[n,b,c]
	    bqf_reps_updatesolutions(sols, nsols, &gcoeff(transmat,1,1), &gcoeff(transmat,2,1));//Update!
      }
	  else cgiv(fred);//May as well discard
	}
  }
}

//Assume q is primitive, rootD is the positive square root of D (a positive square).
GEN sbqf_reps(GEN q, GEN D, GEN rootD, GEN n, int half){
  pari_sp top=avma;
  //We want to factorize q as (ax+by)(cx+dy)
  GEN a,b,c,d;
  if(signe(gel(q,1))!=0){
    d=shifti(subii(gel(q,2),rootD),-1);//(B-sqrt(D))/2
    b=addii(d,rootD);//(B+sqrt(D))/2
    a=gcdii(d,gel(q,1));//Call a=g
    c=diviiexact(gel(q,1),a);//c=A/g
    b=diviiexact(b,c);//b=(b+sqrt(D))/(2A/g)
    d=diviiexact(d,a);//d=d/g
  }
  else{//A=0
	if(signe(gel(q,2))==1){//We want to make sure det(a,b;c,d)=-sqrt(D) to be consistent with the A!=0 case
	  a=gen_0;
	  b=gen_1;
	  c=icopy(gel(q,2));
	  d=icopy(gel(q,3));
	}
	else{
		a=icopy(gel(q,2));
		b=icopy(gel(q,3));
		c=gen_0;
		d=gen_1;
	}
  }
  //Now we are factored as (gx+(B+sqrt(D))/2(A/g))(A/gx+(B-sqrt(D))/(2g)) with corresponding a,b,c,d for coefficients.
  if(isexactzero(n)){//Two linears through 0,0
	GEN rvec=cgetg(4,t_VEC);
	gel(rvec,1)=mkvec(gen_2);//Linear
	togglesign_safe(&a);
	gel(rvec,2)=cgetg(3,t_VEC);
	gel(gel(rvec,2),1)=mkvec2copy(b, a);
	gel(gel(rvec,2),2)=mkvec2(gen_0, gen_0);
	togglesign_safe(&c);
	long lx;
	gel(rvec,3)=cgetg_copy(gel(rvec,2),&lx);
	gel(gel(rvec,3),1)=mkvec2copy(d, c);
	gel(gel(rvec,3),2)=mkvec2(gen_0, gen_0);
	return gerepileupto(top,rvec);
  }
  GEN den;
  if(signe(n)==-1){togglesign_safe(&a);togglesign_safe(&b);den=rootD;}//Toggling a and b; there is no need to update n to -n as divisors(n) spits out positive numbers
  else den=negi(rootD);//den represents the determinant of M=[a,b;c,d]. We must solve M[x;y]=[u;n/u] for all divisors of n. Finitely many solutions now
  togglesign_safe(&b);
  togglesign_safe(&c);//M^(-1)=[d,b;c,a]/den now.
  GEN x, y, r;
  GEN divs=divisors(n);
  glist *sols=NULL;
  long nsols=0;
  long mi=lg(divs);
  for(long i=1;i<lg(divs);++i){//Run through the divisors
	mi--;//i+mi=lg(divs), so divs[i]*divs[mi]=|n|
	x=dvmdii(addii(mulii(d,gel(divs,i)),mulii(b,gel(divs,mi))),den,&r);
	if(!isexactzero(r)){cgiv(r);cgiv(x);continue;}//Does not divide exactly!
	cgiv(r);
	y=dvmdii(addii(mulii(c,gel(divs,i)),mulii(a,gel(divs,mi))),den,&r);
	if(!isexactzero(r)){cgiv(r);cgiv(y);cgiv(x);continue;}//Does not divide exactly!
	cgiv(r);
	bqf_reps_updatesolutions(&sols,&nsols,&x,&y);
  }
  GEN rvec=bqf_reps_creatervec_proper(sols,nsols,half);
  gel(rvec,1)=cgetg(2,t_VEC);
  gel(gel(rvec,1),1)=gen_0;//Finite
  return gerepileupto(top,rvec);
}

//Please input the form as (Ax+By)^2=n where A,B are coprime. If you do not have this form, please use bqf_reps.
GEN zbqf_reps(GEN A, GEN B, GEN n, int half){
  pari_sp top=avma;
  long lx;
  if(isexactzero(n)){//n=0
	GEN rvec=cgetg(3,t_VEC);
	gel(rvec,1)=mkvec(gen_2);//Linear
	gel(rvec,2)=cgetg_copy(rvec,&lx);
	gel(gel(rvec,2),1)=cgetg_copy(rvec,&lx);
	gel(gel(gel(rvec,2),1),1)=icopy(B);
	gel(gel(gel(rvec,2),1),2)=negi(A);
	gel(gel(rvec,2),2)=mkvec2(gen_0,gen_0);
	return rvec;	  
  }
  GEN rootn;
  long sq=Z_issquareall(n,&rootn);
  if(sq==0){avma=top;return gen_0;}
  GEN S=lin_intsolve(A,B,rootn);
  if(half==1){
	GEN rvec=cgetg(3,t_VEC);
	gel(rvec,1)=mkvec(gen_2);//Linear
	gel(rvec,2)=cgetg_copy(rvec,&lx);
	gel(gel(rvec,2),1)=ZV_copy(gel(S,1));
	gel(gel(rvec,2),2)=ZV_copy(gel(S,2));
	return gerepileupto(top,rvec);
  }
  GEN rvec=cgetg(4,t_VEC);
  gel(rvec,1)=mkvec(gen_2);//Linear
  gel(rvec,2)=cgetg(3,t_VEC);
  gel(gel(rvec,2),1)=ZV_copy(gel(S,1));
  gel(gel(rvec,2),2)=ZV_copy(gel(S,2));
  gel(rvec,3)=cgetg_copy(gel(rvec,2),&lx);
  gel(gel(rvec,3),1)=ZV_copy(gel(S,1));
  gel(gel(rvec,3),2)=ZV_copy(gel(S,2));
  ZV_togglesign(gel(gel(rvec,3),2));
  return gerepileupto(top,rvec);
}

//bqf_reps for q=0
static GEN bqf_reps_all(GEN n){
  if(signe(n)!=0) return gen_0;
  GEN rvec=cgetg(2,t_VEC);
  gel(rvec,1)=cgetg(2,t_VEC);
  gel(gel(rvec,1),1)=gen_m1;//All
  return rvec;//"A" stands for all
}

//Returns the trivial solution set.
static GEN bqf_reps_trivial(void){
  GEN rvec=cgetg(3,t_VEC);
  gel(rvec,1)=cgetg(2,t_VEC);
  gel(rvec,2)=cgetg(3,t_VEC);
  gel(gel(rvec,1),1)=gen_0;//Finite
  gel(gel(rvec,2),1)=gen_0;
  gel(gel(rvec,2),2)=gen_0;
  return rvec;
}

//Creates the return vector given the list of solutions, the number of solutions, and the integer half. This does NOT initialize the first component, must be done in the individual method.
static GEN bqf_reps_creatervec(glist *sols, glist *scales, llist *nsolslist, long *totnsols, long *count, int half){
  if(half==0){
	long lind1=2*(*totnsols)+1;
	long lind2=*totnsols+1, nsols;
	GEN rvec=cgetg(lind1+1,t_VEC);
	GEN scale;
	for(long i=1;i<=*count;++i){
	  scale=scales->data;
	  nsols=nsolslist->data;
	  for(long j=1;j<=nsols;++j){
        gel(rvec,lind2)=ZV_Z_mul(sols->data,scale);
		gel(rvec,lind1)=ZV_copy(gel(rvec,lind2));
		ZV_togglesign(gel(rvec,lind1));//Negating it
	    sols=sols->next;
	    lind1--;
		lind2--;
	  }
	  scales=scales->next;
	  nsolslist=nsolslist->next;
    }
	return rvec;	  
  }
  GEN rvec=cgetg(*totnsols+2,t_VEC);
  long lind=*totnsols+1, nsols;
  GEN scale;
  for(long i=1;i<=*count;++i){
	scale=scales->data;
	nsols=nsolslist->data;
	for(long j=1;j<=nsols;++j){
      gel(rvec,lind)=ZV_Z_mul(sols->data,scale);
	  sols=sols->next;
	  lind--;
	}
	scales=scales->next;
	nsolslist=nsolslist->next;
  }
  return rvec;
}

//Creates the return vector given the list of solutions, the number of solutions, and the integer half. This does NOT initialize the first component, must be done in the individual method.
static GEN bqf_reps_creatervec_proper(glist *sols, long nsols, int half){
  if(half==0){
    long lind1=2*nsols+1;
	long lind2=nsols+1;
	GEN rvec=cgetg(lind1+1,t_VEC);
    while(lind2>1){
	  gel(rvec,lind2)=ZV_copy(sols->data);
	  gel(rvec,lind1)=ZV_copy(sols->data);
	  ZV_togglesign(gel(rvec,lind1));//Negating it
	  sols=sols->next;
	  lind1--;
	  lind2--;
    }
	return rvec;
  }
  GEN rvec=cgetg(nsols+2,t_VEC);
  long lind=nsols+1;
  while(lind>1){
	gel(rvec,lind)=ZV_copy(sols->data);
	sols=sols->next;
	lind--;
  }
  return rvec;
}

//Returns q/gcd(q), and sets n to n/gcd(q). If gcd(q) does not divide n, returns NULL. If gcd(q)=1 does not clutter stack (gcd=gen_1), but otherwise this does clutter the stack
static GEN bqf_reps_makeprimitive(GEN q, GEN *n){
  GEN qgcd=ZV_content(q);
  if(equali1(qgcd)) return gcopy(q);
  //Now q not primitive; now qgcd=gcd(q[1],q[2],q[3])
  GEN remainder;
  *n=dvmdii(*n,qgcd,&remainder);//Replace n by n/qgcd
  if(!isexactzero(remainder)) return NULL;
  return ZV_Z_divexact(q,qgcd);
}

//sols stores the list of solutions, of length nsols, this method adds [a,b] to the list and updates the number of solutions. Does not clutter stack, but not gerepile safe as the components of solution are created before solution.
static void bqf_reps_updatesolutions(glist **sols, long *nsols, GEN *a, GEN *b){
  GEN solution=cgetg(3,t_VEC);
  gel(solution,1)=*a;
  gel(solution,2)=*b;
  glist_putstart(sols,solution);
  *nsols=*nsols+1;
}



//MORE REPRESENTATION OF NUMBERS



//Solves Ax^2+Bxy+Cy^2+Dx+Ey=n, where q=[A,B,C,D,E]. Can save computation time by using the more precice dbqf/ibqf_reps methods, but this is not important right now...
GEN bqf_bigreps(GEN q, GEN n, long prec){
  pari_sp top=avma;
  if(isexactzero(q)) return bqf_reps_all(n);//q=0
  n=icopy(n);//Copying n
  q=bqf_reps_makeprimitive(q,&n);//Updating q and n. This works even though q now has length 5
  if(q==NULL){avma=top; return gen_0;}//No solution; gcd(q) did not divide n.
  GEN disc=bqf_disc(q);//Still works as B^2-4AC even though q has length 5
  if(gequal0(disc)){//Disc 0
	if(signe(gel(q,1))==-1 || signe(gel(q,3))==-1){ZV_togglesign(q);togglesign_safe(&n);}//Making positive
	return gerepileupto(top,zbqf_bigreps(q, n));
  }
  GEN newq=cgetg(4,t_VEC);
  for(int i=1;i<=3;i++) gel(newq,i)=icopy(gel(q,i));//newq=[A,B,C]
  GEN a=subii(shifti(mulii(gel(q,3),gel(q,4)),1),mulii(gel(q,2),gel(q,5)));//a=2CD-BE
  GEN b=subii(shifti(mulii(gel(q,1),gel(q,5)),1),mulii(gel(q,2),gel(q,4)));//b=2AE-BD
  pari_sp mid=avma;
  GEN Aa2pBab=addii(mulii(gel(q,1),sqri(a)),mulii(gel(q,2),mulii(a,b)));//A*a^2+B*a*b
  GEN Cb2pDdisca=addii(mulii(gel(q,3),sqri(b)),mulii(gel(q,4),mulii(disc,a)));//C*b^2+D*disc*a
  GEN Ediscb=mulii(gel(q,5),mulii(disc,b));//E*disc*b
  GEN newn=gerepileupto(mid,subii(mulii(sqri(disc),n),addii(Aa2pBab,addii(Cb2pDdisca,Ediscb))));//newn=disc^2*n-(A*a^2+B*a*b+C*b^2+D*disc*a+E*disc*b)
  GEN newsols=bqf_reps(newq, newn, 0, 0, prec);//Solving the bqf_translated form. We have x=(X+a)/disc, y=(Y+b)/disc with (X,Y) solutions to the new equation. Thus we only take solutions with X==-a mod disc and Y==-b mod disc
  if(gequal0(newsols)){avma=top;return gen_0;}//No solutions
  GEN first=gel(gel(newsols,1),1);
  if(gequal0(first)) return gerepileupto(top,bqf_bigreps_creatervecfin(newsols, a, b, disc));//For Finite solutions, i.e. disc<0, disc=square!=0 and n!=0, disc!=square and n=0.
  if(equali1(first)) return gerepileupto(top,bqf_bigreps_creatervecpos(newsols, a, b, disc));//For disc>0 non-square and n!=0
  return gerepileupto(top,bqf_bigreps_createrveclin(newsols, a, b, disc));//For disc!=square and n=0
}

//bqf_bigreps with typecheck
GEN bqf_bigreps_tc(GEN q, GEN n, long prec){
  if(typ(q)!=t_VEC) pari_err_TYPE("bqf_bigreps: please enter a length five integral VECTOR",q);
  if(lg(q)!=6) pari_err_TYPE("bqf_bigreps: please enter a length FIVE integral vector",q);
  for(int i=1;i<6;i++) if(typ(gel(q,i))!=t_INT) pari_err_TYPE("bqf_bigreps: please enter a length five INTEGRAL vector",q);
  if(typ(n)!=t_INT) pari_err_TYPE("bqf_bigreps: please input an integer n",n);
  return bqf_bigreps(q, n, prec);
}

//bqf_bigreps for discriminant 0; q must be primitive and have A,C>=0 and Q!=0
static GEN zbqf_bigreps(GEN q, GEN n){
  pari_sp top=avma, mid;
  if(gequal0(gel(q,4)) && gequal0(gel(q,5))){//D=E=0, this is just bqf_reps
    GEN A=sqrti(gel(q,1));
	GEN B=sqrti(gel(q,3));//must be squares
	if(signe(gel(q,2))==-1) togglesign_safe(&B);//Now (Ax+By)^2=n is to be solved
	return gerepileupto(top,zbqf_reps(A,B,n,0));
  }
  if(gequal0(gel(q,1)) && gequal0(gel(q,3))){//A=B=C=0, just a linear equation
	GEN rvec=cgetg(3,t_VEC);
	gel(rvec,1)=cgetg(2,t_VEC);
	gel(gel(rvec,1),1)=gen_2;//Linear
	gel(rvec,2)=lin_intsolve(gel(q,4), gel(q,5), n);//As q is primitive, this is guarenteed to have solutions
	return rvec;//No garbage!!
  }
  GEN g=gcdii(gcdii(gel(q,1),gel(q,2)),gel(q,3));
  GEN A=sqrti(diviiexact(gel(q,1),g)), B;
  if(!gequal0(A)) B=diviiexact(gel(q,2),shifti(mulii(g,A),1));//Now we need to solve g(Ax+By)^2+Dx+Ey=n. So let u=Ax+By, then Dx+Ey=n-gu^2 (and note gcd(A,B)=1 since A=B=0 has already been dealt with).
  else B=gen_1;//If A=0, must have B=1
  GEN det=subii(mulii(A,gel(q,5)),mulii(B,gel(q,4)));//The determinant of [A,B;D,E]
  if(!gequal0(det)){//Invertible
    glist *S=NULL;
    long nsols=1;//Only because we have the first element
	glist_putstart(&S, mkvec(gen_m2));//Type
	GEN basex=mkvec3(mulii(B, g), gel(q, 5), negi(mulii(B, n)));
	GEN basey=mkvec3(negi(mulii(A, g)), negi(gel(q,4)), mulii(A, n));
	GEN gnew=gcdii(ZV_content(basex), gcdii(det, ZV_content(basey)));//Removing trivial factors
	det=diviiexact(det,gnew);
	basex=ZV_Z_divexact(basex, gnew);
	basey=ZV_Z_divexact(basey, gnew);
	//Now det*x=basex[1]u^2+basex[2]u+basex[3], similarly for det*y with basey instead. So we just need to solve this quadratic to be ==0 mod det.
	GEN u=gen_0, dx, dy, x, y, r, v;
	for(;cmpii(u,absi(det))<0;u=addis(u,1)){
	  mid=avma;
	  dx=addii(mulii(addii(mulii(gel(basex,1),u),gel(basex,2)),u),gel(basex,3));//Should really do this for u=0 and add the difference each time, maybe I'll do this later
	  x=dvmdii(dx,det,&r);
	  if(!gequal0(r)){avma=mid;continue;}//No good, continue on
	  dy=addii(mulii(addii(mulii(gel(basey,1),u),gel(basey,2)),u),gel(basey,3));//Same comment as above
	  y=dvmdii(dy,det,&r);
	  if(!gequal0(r)){avma=mid;continue;}//No good, continue on
	  v=mkvec2(mkvec3(mulii(gel(basex,1),det),addii(shifti(mulii(gel(basex,1),u),1),gel(basex,2)),x), mkvec3(mulii(gel(basey,1),det),addii(shifti(mulii(gel(basey,1),u),1),gel(basey,2)),y));//Substituting u=det*X+u and dividing by det.
	  glist_putstart(&S, v);
	  nsols++;
	}
	if(nsols==1){glist_free(S);avma=top;return gen_0;}//No solutions
    return gerepileupto(top,glist_togvec(S, nsols, -1));//Solutions!
  }//Not invertible. We also know that [A,B] and [D,E] are both not [0,0] as we already checked for this.
  GEN w;
  if(gequal0(B)) w=diviiexact(gel(q,4),A);//Must divide exactly as A,B are coprime and [A,B] and [D,E] are scalar multiples
  else w=diviiexact(gel(q,5),B);//w represents the scale, so [D,E]=w[A,B] and w!=0 or oo. Thus we get gu^2+wu-n=0 where Ax+By=u.
  GEN urootpart=addii(sqri(w),shifti(mulii(g,n),2)), r;//w^2+4gn
  urootpart=sqrtremi(urootpart,&r);
  if(!gequal0(r)){avma=top;return gen_0;}//Irrational, no bueno!
  GEN u1=Qdivii(negi(addii(w,urootpart)),shifti(g,1));//(-w-urootpart)/(2g), one of the roots for u
  if(gequal0(urootpart)){//Only one solution
    if(typ(u1)!=t_INT){avma=top;return gen_0;}//Not integral, no solution.
	GEN rvec=cgetg(3,t_VEC);
	gel(rvec,1)=cgetg(2,t_VEC);
	gel(gel(rvec,1),1)=gen_2;//Linear
	gel(rvec,2)=lin_intsolve(A, B, u1);
	return gerepileupto(top,rvec);
  }//Two solutions potentially
  GEN u2=gadd(u1,Qdivii(urootpart,g));//The other root.
  if(typ(u1)==t_INT){
    if(typ(u2)==t_INT){
      GEN rvec=cgetg(4,t_VEC);gel(rvec,1)=cgetg(2,t_VEC);gel(gel(rvec,1),1)=gen_2;//Linear
	  gel(rvec,2)=lin_intsolve(A,B,u1);gel(rvec,3)=lin_intsolve(A,B,u2);
	  return gerepileupto(top,rvec);
	}
    else{
      GEN rvec=cgetg(3,t_VEC);gel(rvec,1)=cgetg(2,t_VEC);gel(gel(rvec,1),1)=gen_2;//Linear
	  gel(rvec,2)=lin_intsolve(A,B,u1);
	  return gerepileupto(top,rvec);
    }	
  }
  else{
	if(typ(u2)==t_INT){
	  GEN rvec=cgetg(3,t_VEC);gel(rvec,1)=cgetg(2,t_VEC);gel(gel(rvec,1),1)=gen_2;//Linear
	  gel(rvec,2)=lin_intsolve(A,B,u2);
	  return gerepileupto(top,rvec);
	}
	else{avma=top;return gen_0;}//No solution
  }
}

//For finite solutions, i.e. disc<0, disc=square!=0 and n!=0, disc!=square and n=0.
static GEN bqf_bigreps_creatervecfin(GEN newsols, GEN a, GEN b, GEN disc){
  pari_sp top=avma, mid;
  glist *S=NULL;
  long nsols=1;//Only because we have the first element
  glist_putstart(&S, mkvec(gen_0));//Finite solutions
  GEN xpa, ypb, newx, newy, r;
  for(long i=2;i<lg(newsols);i++){
	mid=avma;
	xpa=addii(gel(gel(newsols,i),1),a);
	newx=dvmdii(xpa,disc,&r);
	if(!gequal0(r)){avma=mid;continue;}//Did not satisfy congruence
	ypb=addii(gel(gel(newsols,i),2),b);
	newy=dvmdii(ypb,disc,&r);
	if(!gequal0(r)){avma=mid;continue;}//Did not satisfy congruence
	glist_putstart(&S, mkvec2(newx, newy));//Not clean, but we copy it later with glist_togvec so OK
	nsols++;
  }
  if(nsols==1){glist_free(S);avma=top;return gen_0;}//No solutions
  return gerepileupto(top,glist_togvec(S, nsols, -1));//Solutions!
}

//For disc!=square and n=0
static GEN bqf_bigreps_createrveclin(GEN newsols, GEN a, GEN b, GEN disc){
  pari_sp top=avma, mid;
  //For the solution [s1,s2] and [x0,y0] we get X=(Us1+x0+a)/disc, Y=(Us2+y0+b)/disc. Rearrange to Us_1==-x_0-a and Us_2==-x_0-b both modulo disc. So solve this by writing 1=g_1s_1+g_2s_2. Since we could only have x0=y0=0, this is even easier.
  glist *S=NULL;
  long nsols=1;//Only because we have the first element
  glist_putstart(&S, mkvec(gen_2));//Linear solutions
  GEN slope,g, u, v, res, amod=Fp_red(a, disc), bmod=Fp_red(b, disc), l1, l2;
  for(long i=2;i<lg(newsols);i++){//Each solution is [slope,[0,0]]
	mid=avma;
	slope=gel(gel(newsols,i),1);
	g=bezout(gel(slope,1),gel(slope,2),&u,&v);//us1+vs2=g=1 since s1 and s2 are coprime
	u=Fp_red(u, disc);v=Fp_red(v, disc);
	res=Fp_neg(Fp_add(Fp_mul(u,amod,disc),Fp_mul(v,bmod,disc), disc), disc);//res=-ua-vb; now N===res modulo disc
	l1=addii(mulii(res,gel(slope,1)),a);
	l1=dvmdii(l1,disc,&g);
	if(signe(g)!=0){avma=mid;continue;}//No solution
	l2=addii(mulii(res,gel(slope,2)),b);
	l2=dvmdii(l2,disc,&g);
	if(signe(g)!=0){avma=mid;continue;}//No solution
	//Now we have Us1+l1, Us2+l2
	glist_putstart(&S, mkvec2(slope,mkvec2(l1,l2)));//Very unclean, but everything is copied later.
	nsols++;
  }
  if(nsols==1){glist_free(S);avma=top;return gen_0;}//No solutions
  return gerepileupto(top,glist_togvec(S, nsols, -1));//Solutions!
}

//For disc>0 non-square and n!=0
static GEN bqf_bigreps_creatervecpos(GEN newsols, GEN a, GEN b, GEN disc){
  pari_sp top=avma;
  glist *S=NULL;
  long nsols=1;//Only because we have the first element
  GEN firstelt=mkvec3(gen_1, gel(gel(newsols,1),2),mkvec2(Qdivii(a,disc),Qdivii(b,disc)));
  glist_putstart(&S, firstelt);//Entering the first elt
  GEN Mbase=shallowtrans(ZM_copy(gel(gel(newsols,1),2))), M=mkmat22s(1,0,0,1);//Invariant automorph transposed, and M represents (Mbase^pow)^transposed  
  GEN Mmod=FpM_red(ginv(gel(gel(newsols,1),2)),disc);//inv autom^-1 mod disc
  GEN baseres=FpC_red(mkcol2(negi(a),negi(b)),disc);//[-a;-b] modulo disc
  GEN newsolsmod=cgetg(lg(newsols)-1,t_VEC);
  for(long i=1;i<lg(newsols)-1;i++) gel(newsolsmod,i)=mkcol2(Fp_red(gel(gel(newsols,i+1),1), disc), Fp_red(gel(gel(newsols,i+1),2), disc));//Making [x;y] mod disc for all solutions [x,y]
  GEN sol=baseres;//Represents M^(-pow)[a;b] mod disc. Any time this appears in newsols, we append (M^pow[x;y]+[a;b])/disc to the solution set
  do{
	for(long i=1;i<lg(newsolsmod);i++){
	  if(gequal(sol,gel(newsolsmod,i))){//Can make slightly? faster with Qdivii on a vector
		glist_putstart(&S, gdiv(ZV_ZM_mul(gel(newsols,i+1),M),disc));//We use M^T so that ([x,y]*M^T)^T=M*[x;y] as desired.
		nsols++;
	  }
	}
    M=ZM_mul(M, Mbase);//Updating M
    sol=FpM_FpC_mul(Mmod, sol, disc);//Updating sol
  }
  while(!gequal(sol,baseres));
  if(nsols==1){glist_free(S);avma=top;return gen_0;}//No solutions
  M=shallowtrans(M);//Transpose back
  gel(firstelt,2)=M;//Update M
  return gerepileupto(top,glist_togvec(S, nsols, -1));//Solutions!
}

//Solves AX^2+BY^2+CZ^2+DXY+EXZ+FYZ=n1 and UX+VY+WZ=n2.
GEN bqf_linearsolve(GEN q, GEN n1, GEN lin, GEN n2, long prec){
  pari_sp top=avma;
  GEN g=ZV_content(q), r;
  if(!equali1(g)){
	n1=dvmdii(n1, g, &r);
	if(!gequal0(r)){avma=top;return gen_0;}//Instant fail
	q=ZV_Z_divexact(q,g);
  }
  g=ZV_content(lin);
  if(!equali1(g)){
	n2=dvmdii(n2, g, &r);
	if(!gequal0(r)){avma=top;return gen_0;}//Instant fail
	lin=ZV_Z_divexact(lin,g);
  }
  GEN Mshift=mat3_complete(gel(lin,1),gel(lin,2),gel(lin,3));//M[x;y;z]=[x';y';z'], and lin is now x'=n2
  GEN Mshiftinv=ZM_inv(Mshift, NULL);
  GEN qM=mkmat3(mkcol3(gel(q,1),Qdivii(gel(q,4),gen_2),Qdivii(gel(q,5),gen_2)), mkcol3(Qdivii(gel(q,4),gen_2),gel(q,2),Qdivii(gel(q,6),gen_2)), mkcol3(Qdivii(gel(q,5),gen_2), Qdivii(gel(q,6),gen_2), gel(q,3)));//[x,y,z]qM[x;y;z]=n1
  qM=gmul(gtrans(Mshiftinv), gmul(qM, Mshiftinv));//Now [x',y',z']qM[x';y';z']=n1
  GEN yzeqn=cgetg(6,t_VEC);
  gel(yzeqn,1)=gcoeff(qM,2,2);//y'^2
  gel(yzeqn,2)=gmulgs(gcoeff(qM,2,3),2);//y'z'
  gel(yzeqn,3)=gcoeff(qM,3,3);//z'^2
  gel(yzeqn,4)=mulii(gmulgs(gcoeff(qM,1,2),2),n2);//y'
  gel(yzeqn,5)=mulii(gmulgs(gcoeff(qM,1,3),2),n2);//z'
  GEN newn1=subii(n1,mulii(gcoeff(qM,1,1),sqri(n2)));
  GEN yzsols=bqf_bigreps(yzeqn, newn1, prec);
  if(typ(yzsols)!=t_VEC){avma=top;return gen_0;}//No solutions
  if(gequal0(gel(gel(yzsols,1),1))) return gerepileupto(top,bqf_linearsolve_zfin(yzsols, n2, Mshiftinv));//Finite
  if(equali1(gel(gel(yzsols,1),1))) return gerepileupto(top,bqf_linearsolve_zpos(yzsols, n2, Mshiftinv, Mshift));//Positive
  if(equalii(gel(gel(yzsols,1),1),gen_2)) return gerepileupto(top,bqf_linearsolve_zlin(yzsols, n2, Mshiftinv));//Linear
  if(equalim1(gel(gel(yzsols,1),1))) return gerepileupto(top,bqf_linearsolve_zall(yzsols, n2, Mshiftinv));//ALL, so hyperplane
  return gerepileupto(top,bqf_linearsolve_zquad(yzsols, n2, Mshiftinv));//Quadratic
}

//bqf_linearsolve with typecheck
GEN bqf_linearsolve_tc(GEN q, GEN n1, GEN lin, GEN n2, long prec){
  if(typ(q)!=t_VEC) pari_err_TYPE("Please enter a length 6 integral VECTOR", q);
  if(lg(q)!=7) pari_err_TYPE("Please enter a length SIX integral vector", q);
  for(long i=1;i<7;i++) {if(typ(gel(q,i))!=t_INT) pari_err_TYPE("Please enter a length 6 INTEGRAL vector", q);}
  if(typ(n1)!=t_INT) pari_err_TYPE("Please enter an integer", n1);
  if(gequal0(lin)) pari_err_TYPE("Please enter a NON-ZERO length 3 integral vector", lin);
  if(typ(lin)!=t_VEC) pari_err_TYPE("Please enter a non-zero length 3 integral VECTOR", lin);
  if(lg(lin)!=4) pari_err_TYPE("Please enter a non-zero length THREE integral vector", lin);
  for(long i=1;i<4;i++) {if(typ(gel(lin,i))!=t_INT) pari_err_TYPE("Please enter a non-zero length 3 INTEGRAL vector", lin);}
  if(typ(n2)!=t_INT) pari_err_TYPE("Please enter an integer", n2);
  return bqf_linearsolve(q, n1, lin, n2, prec);
}

//Given lin, n2, and the y',z' solutions, this boosts them up to x,y,z solutions when y' and z' can be ANYTHING
static GEN bqf_linearsolve_zall(GEN yzsols, GEN n2, GEN Minv){//We get a hyperplane
  pari_sp top=avma;
  GEN v1=mkcol3(n2, gen_0, gen_0);
  GEN v2=mkcol3(gen_0, gen_1, gen_0);
  GEN v3=mkcol3(gen_0, gen_0, gen_1);//[x',y',z']=[n1,U,V] for U,V anything
  GEN sols=cgetg(3,t_VEC);
  gel(sols,1)=mkvec(gen_m1);//Hyperplane
  gel(sols,2)=cgetg(4,t_VEC);
  gel(gel(sols,2),1)=ZM_ZC_mul(Minv,v2);settyp(gel(gel(sols,2),1), t_VEC);
  gel(gel(sols,2),2)=ZM_ZC_mul(Minv,v3);settyp(gel(gel(sols,2),2), t_VEC);
  gel(gel(sols,2),3)=ZM_ZC_mul(Minv,v1);settyp(gel(gel(sols,2),3), t_VEC);
  return gerepileupto(top, sols);
}

//Given lin, n2, and the y',z' solutions, this boosts them up to x,y,z solutions when this is finite
static GEN bqf_linearsolve_zfin(GEN yzsols, GEN n2, GEN Minv){
  pari_sp top=avma;
  long lx;
  GEN v=mkcol3(n2, gen_0, gen_0);
  GEN sols=cgetg_copy(yzsols, &lx);
  gel(sols,1)=mkvec(gen_0);//Finite
  for(long i=2;i<lg(yzsols);i++){
	gel(v,2)=gel(gel(yzsols,i),1);
	gel(v,3)=gel(gel(yzsols,i),2);
	gel(sols,i)=ZM_ZC_mul(Minv, v);
	settyp(gel(sols,i), t_VEC);
  }
  return gerepileupto(top, sols);
}

//Given lin, n2, and the y',z' solutions, this boosts them up to x,y,z solutions when this is linear
static GEN bqf_linearsolve_zlin(GEN yzsols, GEN n2, GEN Minv){
  pari_sp top=avma;
  long lx;
  GEN v1=mkcol3(n2, gen_0, gen_0);
  GEN v2=mkcol3(gen_0, gen_0, gen_0);//Will be [x';y';z']=v1+Uv2, so v1=[n2;y0;z0] and v2=[0;s;t]
  GEN sols=cgetg_copy(yzsols, &lx);
  gel(sols,1)=mkvec(gen_2);//Linear
  for(long i=2;i<lg(yzsols);i++){
	gel(v1,2)=gel(gel(gel(yzsols,i),2),1);//y0
	gel(v1,3)=gel(gel(gel(yzsols,i),2),2);//z0
	gel(v2,2)=gel(gel(gel(yzsols,i),1),1);//s
	gel(v2,3)=gel(gel(gel(yzsols,i),1),2);//t
	gel(sols,i)=cgetg(3,t_VEC);
	gel(gel(sols,i),1)=ZM_ZC_mul(Minv, v2);
	gel(gel(sols,i),2)=ZM_ZC_mul(Minv, v1);
	settyp(gel(gel(sols,i),1), t_VEC);
	settyp(gel(gel(sols,i),2), t_VEC);
  }
  return gerepileupto(top, sols);
}

//Given lin, n2, and the y',z' solutions, this boosts them up to x,y,z solutions when this is positive
static GEN bqf_linearsolve_zpos(GEN yzsols, GEN n2, GEN Minv, GEN M){
  pari_sp top=avma;
  GEN tmat=gel(gel(yzsols,1),2);
  GEN Mext=zeromatcopy(3,3);
  gcoeff(Mext,1,1)=gen_1;
  gcoeff(Mext,2,2)=gcoeff(tmat,1,1);
  gcoeff(Mext,2,3)=gcoeff(tmat,1,2);
  gcoeff(Mext,3,2)=gcoeff(tmat,2,1);
  gcoeff(Mext,3,3)=gcoeff(tmat,2,2);//[1,0,0;0,tmat;0,tmat]
  GEN shiftext=mkcol3(n2,gel(gel(gel(yzsols,1),3),1),gel(gel(gel(yzsols,1),3),2));//[y';z']=Mext^k[0;y0;z0]+shiftext for [y0,z0] in yzsols[2:]
  long lx;
  GEN v=mkcol3(gen_0, gen_0, gen_0);//v represents the base solution, i.e. [0;y0;z0]
  GEN M1=ZM_mul(Minv,Mext);
  GEN sols=cgetg_copy(yzsols, &lx);
  gel(sols,1)=cgetg(4,t_VEC);
  gel(gel(sols,1),1)=gen_1;//Positive
  gel(gel(sols,1),2)=ZM_mul(M1,M);//Transition matrix
  gel(gel(sols,1),3)=gmul(Minv, shiftext);//New shift
  settyp(gel(gel(sols,1),3), t_VEC);
  for(long i=2;i<lg(yzsols);i++){
	gel(v,2)=gel(gel(yzsols,i),1);
	gel(v,3)=gel(gel(yzsols,i),2);
	gel(sols,i)=gmul(Minv, v);
	settyp(gel(sols,i), t_VEC);
  }
  return gerepileupto(top, sols);
}

//Given lin, n2, and the y',z' solutions, this boosts them up to x,y,z solutions when this is quadratic
static GEN bqf_linearsolve_zquad(GEN yzsols, GEN n2, GEN Minv){
  pari_sp top=avma;
  long lx;
  GEN v1=mkcol3(gen_0, gen_0, gen_0);
  GEN v2=mkcol3(gen_0, gen_0, gen_0);
  GEN v3=mkcol3(n2, gen_0, gen_0);//[n2;y';z'=]U^2v1+Uv2+v3 is the idea
  GEN t1, t2, t3;
  GEN sols=cgetg_copy(yzsols, &lx);
  gel(sols,1)=mkvec(gen_m2);//quadratic
  for(long i=2;i<lg(yzsols);i++){
	gel(v1,2)=gel(gel(gel(yzsols,i),1),1);
	gel(v1,3)=gel(gel(gel(yzsols,i),2),1);
	gel(v2,2)=gel(gel(gel(yzsols,i),1),2);
	gel(v2,3)=gel(gel(gel(yzsols,i),2),2);
	gel(v3,2)=gel(gel(gel(yzsols,i),1),3);
	gel(v3,3)=gel(gel(gel(yzsols,i),2),3);
	t1=ZM_ZC_mul(Minv, v1);
	t2=ZM_ZC_mul(Minv, v2);
	t3=ZM_ZC_mul(Minv, v3);
	gel(sols,i)=cgetg(4, t_VEC);
	for(int j=1;j<4;j++) gel(gel(sols,i),j)=mkvec3(gel(t1,j),gel(t2,j),gel(t3,j));
  }
  return gerepilecopy(top, sols);//Must copy
}



//GENERAL CHECKING METHODS



//Checks that an input is an integral bqf (does NOT check discriminant!=square), and produces a pari_error if not. Used for not getting segmentation faults from a GP interface user. 
void bqf_check(GEN q){
  if(typ(q)!=t_VEC) pari_err(e_TYPE,"bqf: please input a length 3 integral VECTOR",q);
  if(lg(q)!=4) pari_err(e_TYPE,"bqf: please input a LENGTH 3 integral vector",q);
  if(typ(gel(q,1))!=t_INT) pari_err(e_TYPE,"bqf: please input a length 3 INTEGRAL vector",q);
  if(typ(gel(q,2))!=t_INT) pari_err(e_TYPE,"bqf: please input a length 3 INTEGRAL vector",q);
  if(typ(gel(q,3))!=t_INT) pari_err(e_TYPE,"bqf: please input a length 3 INTEGRAL vector",q);
}

//Checks that an input is an integral bqf with disc!=square and produces a pari_error if not. Returns the discriminant
GEN bqf_checkdisc(GEN q){
  bqf_check(q);
  GEN D=bqf_disc(q);
  if(isdisc(D)) return D;
  pari_err(e_TYPE,"bqf: please input a bqf with non-square discriminant",q);
  return gen_0;
}

//Checks that an input is an integral 2x2 matrix, and produces a pari_error if not. Used for not getting segmentation faults from a GP interface user. 
void intmatrix_check(GEN mtx){
  if(typ(mtx)!=t_MAT) pari_err(e_TYPE,"mat: please input a hyberbolic two x two intgral MATRIX",mtx);
  if(lg(mtx)!=3) pari_err(e_TYPE,"mat: please input a hyberbolic two x TWO intgral matrix",mtx);
  if(lg(gel(mtx,1))!=3) pari_err(e_TYPE,"mat: please input a hyberbolic TWO x two intgral matrix",mtx);
  if(typ(gcoeff(mtx,1,1))!=t_INT) pari_err(e_TYPE,"mat: please input a hyberbolic two x two INTEGRAL matrix",mtx);
  if(typ(gcoeff(mtx,1,2))!=t_INT) pari_err(e_TYPE,"mat: please input a hyberbolic two x two INTEGRAL matrix",mtx);
  if(typ(gcoeff(mtx,2,1))!=t_INT) pari_err(e_TYPE,"mat: please input a hyberbolic two x two INTEGRAL matrix",mtx);
  if(typ(gcoeff(mtx,2,2))!=t_INT) pari_err(e_TYPE,"mat: please input a hyberbolic two x two INTEGRAL matrix",mtx);
}

