//IMPORTED FROM Q-QUADRATIC
//This is a collection of miscellaneous methods that may be useful in a variety of settings, and not just for the programs they were originally created for.

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "apol.h"
#endif


//STATIC METHOD DECLARATIONS


//INFINITY 


//Adds a,b, and allows for oo
GEN addoo(GEN a, GEN b){//No need to do garbage collection
  if(typ(a)==t_INFINITY){
    if(inf_get_sign(a)==1) return mkoo();
	else return mkmoo();
  }
  if(typ(b)==t_INFINITY){
    if(inf_get_sign(b)==1) return mkoo();
	else return mkmoo();
  }
  return gadd(a,b);
}

//Divides a and b, and allows for oo and division by 0. Returns oo for 0/0.
GEN divoo(GEN a, GEN b){//No garbage collection necessary
  if(gequal0(b)){//b=0
    if(gcmpgs(a,0)>=0) return mkoo();
    return mkmoo();
  }
  if(typ(a)==t_INFINITY){//a=+/-oo
	if(gsigne(a)==gsigne(b)) return mkoo();
	return mkmoo();
  }
  if(typ(b)==t_INFINITY) return gen_0;
  return gdiv(a,b);
}



//INTEGER VECTORS


//Copies an integer vector
GEN ZV_copy(GEN v){
  long len=lg(v);
  GEN rvec=cgetg(len,t_VEC);
  for(long i=1;i<len;++i) gel(rvec,i)=icopy(gel(v,i));
  return rvec;
}

//Checks if the integral vectors v1, v2 are equal, returns 1 if so and 0 else. Segmentation faults will occur if the entries are not integer.
int ZV_equal(GEN v1, GEN v2){
  long len=lg(v1);
  if(lg(v2)!=len) return 0;//Different length
  for(long i=1;i<len;++i){if(!equalii(gel(v1,i),gel(v2,i))) return 0;}
  return 1;
}

//v is a Z-vector, divides v by y, assuming y divides each term exactly
GEN ZV_Z_divexact(GEN v, GEN y){
  long lx;
  GEN rvec=cgetg_copy(v, &lx);
  for(long i=1;i<lx;i++) gel(rvec, i)=diviiexact(gel(v, i), y);
  return rvec;
  
}

//v is a Z-vector, multiplies v by the integer x
GEN ZV_Z_mul(GEN v, GEN x){
  long lx;
  GEN rvec=cgetg_copy(v, &lx);
  for(long i=1;i<lx;i++) gel(rvec, i)=mulii(gel(v, i), x);
  return rvec; 
}



//LINEAR ALGEBRA


//Returns the vector of [eigenvalue, [eigenvectors]] of an FpM, where p is PRIME
GEN FpM_eigenvecs(GEN M, GEN p){
  pari_sp top=avma;
  GEN pol=FpM_charpoly(M, p);
  GEN rts=FpX_roots(pol, p);
  long nrts=lg(rts);
  if(nrts==1){avma=top;return cgetg(1, t_VEC);}//No roots
  GEN shiftM=cgetg(nrts, t_VEC), id=matid(lg(M)-1);
  for(long i=1;i<nrts;i++) gel(shiftM, i)=FpM_sub(M, FpM_Fp_mul(id, gel(rts, i), p), p);//Stores M-eval*Id
  GEN ret=cgetg(nrts, t_VEC);
  for(long i=1;i<nrts;i++){
	gel(ret, i)=cgetg(3, t_VEC);
	gel(gel(ret, i), 1)=icopy(gel(rts, i));//Eigenvalue
	gel(gel(ret, i), 2)=FpM_ker(gel(shiftM, i), p);//Eigenvectors
  }
  return gerepileupto(top, ret);
}

//Solves Ax+By=n
GEN lin_intsolve(GEN A, GEN B, GEN n){
  pari_sp top = avma;
  if(gequal0(A) && gequal0(B)) return gen_0;
  GEN g,X,Y;
  g=gbezout(A,B,&X,&Y);//XA+YB=g now
  GEN scale=Qdivii(n,g);
  if(typ(scale)!=t_INT){//g does not divide n
	avma=top;
	return gen_0;
  }
  GEN rvec=cgetg(3,t_VEC);//Creating the return vector
  gel(rvec,1)=cgetg(3,t_VEC);
  gel(rvec,2)=cgetg(3,t_VEC);
  gel(gel(rvec,2),1)=mulii(X,scale);
  gel(gel(rvec,2),2)=mulii(Y,scale);//Second part initialized
  if(gequal0(A)){
    gel(gel(rvec,1),1)=gen_1;
	gel(gel(rvec,1),2)=gen_0;
  }
  else if(gequal0(B)){
	gel(gel(rvec,1),1)=gen_0;
	gel(gel(rvec,1),2)=gen_1;
  }
  else{
	gel(gel(rvec,1),1)=diviiexact(negi(B),g);
	gel(gel(rvec,1),2)=diviiexact(A,g);
  }
  return(gerepilecopy(top, rvec));
}

//lin_intsolve with typecheck
GEN lin_intsolve_tc(GEN A, GEN B, GEN n){
  if (typ(A)!=t_INT) pari_err_TYPE("lin_intsolve. A should be an integer",A);
  if (typ(B)!=t_INT) pari_err_TYPE("lin_intsolve. B should be an integer",B);
  if (typ(n)!=t_INT) pari_err_TYPE("lin_intsolve. n should be an integer",n);
  return lin_intsolve(A,B,n);
}

//Returns a 3x3 ZM with top row A, B, C. Assumes gcd(A,B,C)=1
GEN mat3_complete(GEN A, GEN B, GEN C){
  pari_sp top=avma;
  GEN u, v, w, x;
  GEN g=bezout(A, B, &u, &v);
  bezout(g,C, &x, &w);
  u=mulii(u,x);
  v=mulii(v,x);
  togglesign_safe(&v);//Now uA-vB+wC=g2=1 since gcd(A,B,C)=1. Write the bottom two rows as [d,f,g;g,h,i]
  GEN m;//We want to now solve ei-fh=u; di-fg=v; dh-eg=w.
  if(gequal0(u) && gequal0(v)){
	m=zeromatcopy(3,3);
	gcoeff(m,2,1)=gen_1;gcoeff(m,2,2)=gen_0;gcoeff(m,2,3)=gen_0;
    gcoeff(m,3,1)=gen_0;gcoeff(m,3,3)=gen_0;
	if(equali1(w)) gcoeff(m,3,2)=gen_1;
	else gcoeff(m,3,2)=gen_m1;//Making det 1
  }
  else{
	GEN p, q;
    g=bezout(u, v, &q, &p);
    GEN H=diviiexact(u, g);togglesign_safe(&H);//We use captials to represent the matrix elements now
    GEN G=diviiexact(v, g);togglesign_safe(&G);
    GEN E=mulii(p, w);
    GEN D=mulii(q, w);togglesign_safe(&D);//DH-EG=w, F=g and I=0
    m=zeromatcopy(3,3);
	gcoeff(m,2,1)=icopy(D);gcoeff(m,2,2)=icopy(E);gcoeff(m,2,3)=icopy(g);
    gcoeff(m,3,1)=icopy(G);gcoeff(m,3,2)=icopy(H);gcoeff(m,3,3)=gen_0;
  }
  gcoeff(m,1,1)=icopy(A);gcoeff(m,1,2)=icopy(B);gcoeff(m,1,3)=icopy(C);
  return gerepileupto(top, m);
}

//Returns a 3x3 ZM with top row A, B, C. Assumes gcd(A,B,C)=1
GEN mat3_complete_tc(GEN A, GEN B, GEN C){
  pari_sp top=avma;
  if(typ(A)!=t_INT) pari_err_TYPE("Please enter three integers with gcd 1", A);
  if(typ(B)!=t_INT) pari_err_TYPE("Please enter three integers with gcd 1", B);
  if(typ(C)!=t_INT) pari_err_TYPE("Please enter three integers with gcd 1", C);
  GEN g=gcdii(A, B);
  GEN g1=gcdii(g,C);
  if(!equali1(g1)) pari_err_TYPE("GCD is not equal to 1", mkvec3(A, B, C));
  avma=top;
  return mat3_complete(A, B, C);
}



//LISTS


//Appends x to v, returning v, and updating vind to vind++. If vind++>vlen, then we double the length of v as well. If this happens, the resulting vector is not suitable for gerepileupto; this must be done at the end (necessary anyway since it's likely we have to call vec_shorten at some point).
GEN veclist_append(GEN v, long *vind, long *vlen, GEN x){
  if(*vind==*vlen){//Need to lengthen!
    *vlen=2**vlen;
    v=vec_lengthen(v, *vlen);
  }
  *vind=*vind+1;
  gel(v, *vind)=x;
  return v;
}

//Appends x to v, returning v, and updating vind to vind++. If vind++>vlen, then we double the length of v as well. Don't forget to call vec_shorten at the end, since some positions are uninitialized.
GEN vecsmalllist_append(GEN v, long *vind, long *vlen, long x){
  if(*vind==*vlen){//Need to lengthen!
    *vlen=2**vlen;
    v=vecsmall_lengthen(v, *vlen);
  }
  *vind=*vind+1;
  v[*vind]=x;
  return v;
}



//PRIMES


//Returns the set of primes in the given range that lie in one of the given residue classes.
GEN primes_mod(GEN range, GEN residues, long modulus){
  pari_sp top=avma;
  if(typ(range)!=t_VEC || lg(range)<3) pari_err_TYPE("range should be a vector of length 2", range);
  if(typ(residues)!=t_VEC){
	if(typ(residues)==t_INT) residues=mkvec(residues);
	else pari_err_TYPE("residues should be a vector of integers, or a single integer", residues);
  }//Checking inputs, making a single residue a vector.
  long nres=lg(residues);
  GEN mods=cgetg(nres, t_VECSMALL);
  for(long i=1;i<nres;i++) mods[i]=smodis(gel(residues, i), modulus);
  mods=vecsmall_uniq(mods);//mods is the sorted list, duplicates removed.
  GEN plist=primes_interval(gel(range, 1), gel(range, 2));//Initial list of all primes.
  long np=lg(plist);
  GEN pmodlist=vectrunc_init(np);//List of primes in the correct residue classes.
  for(long i=1;i<np;i++){
	GEN p=gel(plist, i);
	long rem=smodis(p, modulus);
	if(zv_search(mods, rem)) vectrunc_append(pmodlist, p);//non-zero, so the index was found!
  }
  return gerepilecopy(top, pmodlist);
}



//RANDOM


//Returns a random element from the vector v. NOT STACK CLEAN
GEN rand_elt(GEN v){
  if(typ(v)!=t_VEC) pari_err_TYPE("Not a vector", v);
  long i=rand_l(lg(v)-1);
  return gel(v,i);
}

//Returns random long from 1 to len
long rand_l(long len){
  pari_sp top=avma;
  GEN ind=randomi(stoi(len));
  long ret=itos(ind)+1;
  avma=top;
  return ret;
}



//LISTS


//List of GENs

//Frees the memory pari_malloc'ed by glist
void glist_free(glist *l){
  glist *temp=l;
  while(l!=NULL){
	temp=l;
	l=l->next;
	pari_free(temp);
  }
}

//Removes the last element of the glist and returns it without copying
GEN glist_pop(glist **head_ref){
  glist *temp=*head_ref;
  GEN x=temp->data;
  *head_ref=temp->next;
  pari_free(temp);
  return x;
}

//Put an element at the start of the glist
void glist_putstart(glist **head_ref, GEN new_data){
  glist *new_elt = (glist*)pari_malloc(sizeof(glist)); 
  new_elt->data = new_data; 
  new_elt->next = *head_ref; 
  *head_ref = new_elt; 

}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a vector, makes a clean copy. This also frees the list, but we also need to clean up the list data at the list creation location. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN glist_togvec(glist *l, long length, int dir){
	glist *lcopy=l;
	GEN rvec=cgetg(length+1, t_VEC);
	if(dir==1){
	  long lind=1;
	  while(l!=NULL && lind<=length){
		  gel(rvec,lind)=gcopy(l->data);
		  l=l->next;
		  lind++;
	  }
	  if(lind<=length){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	else{
      long lind=length;
	  while(l!=NULL && lind>0){
		gel(rvec,lind)=gcopy(l->data);
		l=l->next;
		lind--;
	  }
	  if(lind>0){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	glist_free(lcopy);
	return rvec;
}

//Appends l to the end of v, returning a clean copy. dir=-1 means forward, dir=-1 backward. This also frees the list, but we also need to clean up the list data at the list creation location. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN glist_togvec_append(glist *l, GEN v, long length, int dir){
	glist *lcopy=l;
	long vlen=lg(v), rveclen=length+vlen;
	GEN rvec=cgetg(rveclen, t_VEC);
	for(long i=1;i<vlen;i++) gel(rvec, i)=gcopy(gel(v, i));//Copying v
	if(dir==1){
	  long lind=vlen;
	  while(l!=NULL && lind<rveclen){
		  gel(rvec,lind)=gcopy(l->data);
		  l=l->next;
		  lind++;
	  }
	  if(lind<rveclen){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	else{
      long lind=rveclen-1;
	  while(l!=NULL && lind>=vlen){
		gel(rvec,lind)=gcopy(l->data);
		l=l->next;
		lind--;
	  }
	  if(lind>=vlen){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	glist_free(lcopy);
	return rvec;
}


//List of longs

//Frees the memory pari_malloc'ed by llist
void llist_free(llist *l){
  llist *temp=l;
  while(l!=NULL){
	temp=l;
	l=l->next;
	pari_free(temp);
  }
}

//Removes the last element of the llist and returns it
long llist_pop(llist **head_ref){
  llist *temp=*head_ref;
  long x=temp->data;
  *head_ref=temp->next;
  pari_free(temp);
  return x;
}

//Put an element at the start of the llist
void llist_putstart(llist **head_ref, long new_data){
  llist *new_elt = (llist*)pari_malloc(sizeof(llist)); 
  new_elt->data = new_data; 
  new_elt->next = *head_ref; 
  *head_ref = new_elt; 
}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a vector. This also frees the list. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN llist_togvec(llist *l, long length, int dir){//No garbage collection necessary with longs!
  llist *lcopy=l;
  GEN rvec=cgetg(length+1, t_VEC);
  if(dir==1){
    long lind=1;
    while(l!=NULL && lind<=length){
	  gel(rvec,lind)=stoi(l->data);
	  l=l->next;
	  lind++;
    }
    if(lind<=length){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  else{
    long lind=length;
    while(l!=NULL && lind>0){
      gel(rvec,lind)=stoi(l->data);
	  l=l->next;
	  lind--;
    }
    if(lind>0){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  llist_free(lcopy);
  return(rvec);
}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a VECSMALL. This also frees the list. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN llist_tovecsmall(llist *l, long length, int dir){//No garbage collection necessary with longs!
  llist *lcopy=l;
  GEN rvec=cgetg(length+1, t_VECSMALL);
  if(dir==1){
    long lind=1;
    while(l!=NULL && lind<=length){
	  rvec[lind]=l->data;
	  l=l->next;
	  lind++;
    }
    if(lind<=length){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  else{
    long lind=length;
    while(l!=NULL && lind>0){
      rvec[lind]=l->data;
	  l=l->next;
	  lind--;
    }
    if(lind>0){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  llist_free(lcopy);
  return(rvec);
}
