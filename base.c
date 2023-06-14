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



//VECTORS


//Clean version of vecreverse
GEN vecrev(GEN v){
  long lv;
  GEN w=cgetg_copy(v, &lv);
  for(long i=1;i<lv;i++) gel(w, i)=gcopy(gel(v, lv-i));
  return w;
}

//Copies an integer vector
GEN ZV_copy(GEN v){
  long len=lg(v);
  GEN rvec=cgetg(len, t_VEC);
  for(long i=1;i<len;i++) gel(rvec, i)=icopy(gel(v, i));
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



//LISTS OF VARIABLE LENGTH


/* Sample use of veclist_append:
pari_sp top=avma;
long vind=0, vlen=10;
GEN v=zerovec(vlen);//So that we can garbage collect.
...
GEN x=...;
v=veclist_append(v, &vind, &vlen, x);
...
v=vec_shorten(v, vind);
return gerepilecopy(top, v);
*/


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



//MODS

//Returns the squares mod n. If cop=1, only returns coprime squares
GEN modsquares(GEN n, long cop){
  pari_sp top=avma;
  GEN no2=truedivis(n, 2);//floor(n/2)
  GEN L=vectrunc_init(itos(no2)+1);
  for(GEN y=gen_0;cmpii(y, no2)<=0;y=addis(y, 1)){//Only need to go 
	if(cop && !equali1(gcdii(y, n))) continue;//cop=1 and gcd(y, n)>1
	vectrunc_append(L, Fp_sqr(y, n));//Append y^2 mod n
  }
  return gerepileupto(top, ZV_sort_uniq(L));
}

//Given a bunch of residues modulo n, this breaks them down into classes mod q^e, for all prime powers q^e||n. The format is a vector, where each entry corresponds to a prime, and has the format [residues, [q, e, q^e]]. Assume n>1
GEN mod_breakdown(GEN res, GEN n){
  pari_sp top=avma;
  n=absi(n);
  GEN fact=factor(n);
  long lgprimes=lg(gel(fact, 1));
  GEN ret=cgetg(lgprimes, t_VEC);
  for(long i=1;i<lgprimes;i++){
	GEN p=gcoeff(fact, i, 1);//prime
	GEN e=gcoeff(fact, i, 2);//exponent
	GEN pe=powii(p, e);//p^e
	long lgres;
	GEN ppart=cgetg_copy(res, &lgres);
	for(long j=1;j<lgres;j++) gel(ppart, j)=Fp_red(gel(res, j), pe);
	gel(ret, i)=mkvec2(ZV_sort_uniq(ppart), mkvec3(p, e, pe));
  }
  return gerepilecopy(top, ret);
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

//Returns the vector of prime factors of N.
GEN primefactors(GEN N){
  pari_sp top=avma;
  if(gequal0(N)) return cgetg(1, t_VEC);
  GEN ps=gel(absZ_factor(N), 1);
  return gerepilecopy(top, shallowtrans(ps));
}


//RANDOM


//Returns a random element from the vector v. NOT STACK CLEAN
GEN rand_elt(GEN v){return gel(v, 1+random_Fl(lg(v)-1));}



//LISTS


//Circular list of GENs

//Frees the memory pari_malloc'ed by clist
void clist_free(clist *l, long length){
  clist *temp=l;
  long i=1;
  while(i<length){
	temp=l;
	l=l->next;
	pari_free(temp);
	i++;
  }
  pari_free(l);
}

//Put an element before *head_ref, and update *head_ref to point there
void clist_putbefore(clist **head_ref, GEN new_data){
  clist *new_elt = (clist*)pari_malloc(sizeof(clist)); 
  new_elt->data = new_data;
  if(*head_ref!=NULL){
    new_elt->next = *head_ref; 
    new_elt->prev = (*head_ref)->prev;
    (*head_ref)->prev = new_elt;
	(new_elt->prev)->next=new_elt;
  }
  else{
    new_elt->next = new_elt; 
    new_elt->prev = new_elt;
  }
  *head_ref = new_elt;
}

//Put an element after *head_ref, and update *head_ref to point there
void clist_putafter(clist **head_ref, GEN new_data){
  clist *new_elt = (clist*)pari_malloc(sizeof(clist)); 
  new_elt->data = new_data;
  if(*head_ref!=NULL){
    new_elt->prev = *head_ref; 
    new_elt->next = (*head_ref)->next;
    (*head_ref)->next = new_elt;
	(new_elt->next)->prev=new_elt;
  }
  else{
    new_elt->next = new_elt; 
    new_elt->prev = new_elt;
  }
  *head_ref = new_elt;
}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a vector, and makes a clean copy. This also frees the list, but we also need to clean up the list data at the list creation location. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN clist_togvec(clist *l, long length, int dir){
  if(l==NULL){//Empty list, return the empty vector.
    GEN rvec=cgetg(1,t_VEC);
    return(rvec);	  
  }
  GEN rvec=cgetg(length+1, t_VEC);
  long lind=1;
  if(dir==1){
    while(lind<=length){
	  gel(rvec,lind)=gcopy(l->data);
	  l=l->next;
	  lind++;
    }
  }
  else{
    while(lind<=length){
	  gel(rvec,lind)=gcopy(l->data);
	  l=l->prev;
	  lind++;
    }
  }
  clist_free(l,length);
  return rvec;
}


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
