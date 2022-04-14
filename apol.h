//STRUCTURES

typedef struct listtype1{//A circular list of GENs, stores data, next, and previous terms
  GEN data; 
  struct listtype1 *next;
  struct listtype1 *prev;
}clist;

typedef struct listtype2{//A generic linked list of GENs, stores data and next term
  GEN data; 
  struct listtype2 *next;
}glist;

typedef struct listtype3{//A generic linked list of longs, stores data and next term
  long data; 
  struct listtype3 *next;
}llist;

//METHODS


//apol.c


int apol_check(GEN v);
GEN apol_circles(GEN v, GEN maxcurv, int depth, long prec);
void printcircles_desmos(GEN c);
GEN printcircles_tex(GEN c, char *imagename, int addnumbers, int compile, int open, long prec);
GEN apol_dpair_circle(GEN L);
GEN apol_getmatrices();
GEN apol_getobstructions();
GEN apol_make(GEN q, int pos, int red);
GEN apol_mod24(GEN v);
GEN apol_move(GEN v, int ind);
GEN apol_ncgp_depths(GEN n, long prec);
GEN apol_ncgp_forms(GEN n, int pos, int red, long prec);
GEN apol_ncgp_smallcurve(GEN n, int red, long prec);
GEN apol_ncgp_smallcurve_bsteps(GEN n, long maxsteps, long prec);
GEN apol_orbit(GEN v, int depth, GEN bound);
GEN apol_orbit_layers(GEN v, int maxlayers, GEN bound);
GEN apol_orbit_primes(GEN v, int maxlayers, GEN bound);
GEN apol_qf(GEN v, int ind);
long apol_quaddepth(GEN v);
GEN apol_red(GEN v, int seq);
GEN apol_red_bsteps(GEN v, long maxsteps);
GEN apol_search(GEN v, GEN N, int depth, int rqf);
GEN apol_strip_qf(GEN L, int red);
long mod24_search(GEN L, GEN v);
long ZV_countnonpos(GEN v);


//base.c


//INFINITY 
GEN addoo(GEN a, GEN b);
GEN divoo(GEN a, GEN b);

//INTEGER VECTORS
GEN ZV_copy(GEN v);
int ZV_equal(GEN v1, GEN v2);
GEN ZV_Z_divexact(GEN v, GEN y);
GEN ZV_Z_mul(GEN v, GEN x);

//LINEAR ALGEBRA
GEN FpM_eigenvecs(GEN M, GEN p);
GEN lin_intsolve(GEN A, GEN B, GEN n);
GEN lin_intsolve_tc(GEN A, GEN B, GEN n);
GEN mat3_complete(GEN A, GEN B, GEN C);
GEN mat3_complete_tc(GEN A, GEN B, GEN C);

//LISTS
GEN veclist_append(GEN v, long *vind, long *vlen, GEN x);
GEN vecsmalllist_append(GEN v, long *vind, long *vlen, long x);

//PRIMES
GEN primes_mod(GEN range, GEN residues, long modulus);

//RANDOM
GEN rand_elt(GEN v);

//LISTS
void clist_free(clist *l, long length);
void clist_putbefore(clist **head_ref, GEN new_data);
void clist_putafter(clist **head_ref, GEN new_data);
GEN clist_togvec(clist *l, long length, int dir);
void glist_free(glist *l);
GEN glist_pop(glist **head_ref);
void glist_putstart(glist **head_ref, GEN new_data);
GEN glist_togvec(glist *l, long length, int dir);
GEN glist_togvec_append(glist *l, GEN v, long length, int dir);
void llist_free(llist *l);
long llist_pop(llist **head_ref);
void llist_putstart(llist **head_ref, long new_data);
GEN llist_togvec(llist *l, long length, int dir);
GEN llist_tovecsmall(llist *l, long length, int dir);


//bqf.c


//DISCRIMINANT METHODS
GEN disclist(GEN D1, GEN D2, int fund, GEN cop);
GEN discprimeindex(GEN D, GEN facs);
GEN discsuperorders(GEN D);
int isdisc(GEN D);
GEN pell(GEN D);
GEN posreg(GEN D, long prec);
GEN quadroot(GEN D);

//BASIC OPERATIONS ON BINARY QUADRATIC FORMS
GEN bqf_automorph_tc(GEN q);
int bqf_compare(void *data, GEN q1, GEN q2);
int bqf_compare_tmat(void *data, GEN d1, GEN d2);
GEN bqf_disc(GEN q);
GEN bqf_disc_tc(GEN q);
int bqf_isreduced(GEN q, int Dsign);
int bqf_isreduced_tc(GEN q);
GEN bqf_random(GEN maxc, int type, int primitive);
GEN bqf_random_D(GEN maxc, GEN D);
GEN bqf_red(GEN q, GEN rootD, int Dsign, int tmat);
GEN bqf_red_tc(GEN q, int tmat, long prec);
GEN bqf_roots(GEN q, GEN D, GEN w);
GEN bqf_roots_tc(GEN q);
GEN bqf_trans(GEN q, GEN mtx);
GEN bqf_trans_tc(GEN q, GEN mtx);
GEN bqf_transL(GEN q, GEN n);
GEN bqf_transR(GEN q, GEN n);
GEN bqf_transS(GEN q);
GEN bqf_trans_coprime(GEN q, GEN n);
GEN bqf_trans_coprime_tc(GEN q, GEN n);

//BASIC METHODS FOR NEGATIVE DISCRIMINANTS
GEN dbqf_automorph(GEN q, GEN D);
GEN dbqf_red(GEN q);
GEN dbqf_red_tmat(GEN q);

//BASIC OPERATIONS SPECIFIC TO INDEFINITE FORMS/POSITIVE DISCRIMINANTS
GEN ibqf_automorph_D(GEN q, GEN D);
GEN ibqf_automorph_pell(GEN q, GEN qpell);
int ibqf_isrecip(GEN q, GEN rootD);
int ibqf_isrecip_tc(GEN q, long prec);
GEN ibqf_leftnbr(GEN q, GEN rootD);
GEN ibqf_leftnbr_tmat(GEN q, GEN rootD);
GEN ibqf_leftnbr_tc(GEN q, int tmat, long prec);
GEN ibqf_leftnbr_update(GEN qvec, GEN rootD);
GEN ibqf_red(GEN q, GEN rootD);
GEN ibqf_red_tmat(GEN q, GEN rootD);
GEN ibqf_red_pos(GEN q, GEN rootD);
GEN ibqf_red_pos_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit(GEN q, GEN rootD);
GEN ibqf_redorbit_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit_posonly(GEN q, GEN rootD);
GEN ibqf_redorbit_posonly_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit_tc(GEN q, int tmat, int posonly, long prec);
GEN ibqf_rightnbr(GEN q, GEN rootD);
GEN ibqf_rightnbr_tmat(GEN q, GEN rootD);
GEN ibqf_rightnbr_update(GEN qvec, GEN rootD);
GEN ibqf_rightnbr_tc(GEN q, int tmat, long prec);
GEN ibqf_river(GEN q, GEN rootD);
GEN ibqf_river_positions(GEN q, GEN rootD);
GEN ibqf_river_positions_forms(GEN q, GEN rootD);
GEN ibqf_river_tc(GEN q, long prec);
GEN ibqf_riverforms(GEN q, GEN rootD);
GEN ibqf_riverforms_tc(GEN q, long prec);
GEN ibqf_symmetricarc(GEN q, GEN D, GEN rootD, GEN qpell, long prec);
GEN ibqf_symmetricarc_tc(GEN q, long prec);
GEN ibqf_toriver(GEN q, GEN rootD);
GEN ibqf_toriver_tmat(GEN q, GEN rootD);
GEN mat_toibqf(GEN mtx);
GEN mat_toibqf_tc(GEN mtx);

//EQUIVALENCE OF BQFs
GEN bqf_isequiv(GEN q1, GEN q2, GEN rootD, int Dsign, int tmat);
GEN bqf_isequiv_set(GEN q, GEN S, GEN rootD, int Dsign, int tmat);
GEN bqf_isequiv_tc(GEN q1, GEN q2, int tmat, long prec);
GEN dbqf_isequiv(GEN q1, GEN q2);
GEN dbqf_isequiv_tmat(GEN q1, GEN q2);
long dbqf_isequiv_set(GEN q, GEN S);
GEN dbqf_isequiv_set_tmat(GEN q, GEN S);
GEN ibqf_isequiv(GEN q1, GEN q2, GEN rootD);
GEN ibqf_isequiv_tmat(GEN q1, GEN q2, GEN rootD);
long ibqf_isequiv_set_byq(GEN q, GEN S, GEN rootD);
long ibqf_isequiv_set_byq_presorted(GEN qredsorted, GEN S, GEN rootD);
GEN ibqf_isequiv_set_byq_tmat(GEN q, GEN S, GEN rootD);
GEN ibqf_isequiv_set_byq_tmat_presorted(GEN qredsorted, GEN S, GEN rootD);
long ibqf_isequiv_set_byS(GEN q, GEN S, GEN rootD);
long ibqf_isequiv_set_byS_presorted(GEN q, GEN Sreds, GEN perm, GEN rootD);
GEN ibqf_isequiv_set_byS_tmat(GEN q, GEN S, GEN rootD);
GEN ibqf_isequiv_set_byS_tmat_presorted(GEN q, GEN Sreds, GEN perm, GEN rootD);


//CLASS GROUPS AND COMPOSITION OF FORMS
GEN bqf_allforms(GEN D, int prim, int GL, long prec);
GEN bqf_comp(GEN q1, GEN q2);
GEN bqf_comp_red(GEN q1, GEN q2, GEN rootD, int Dsign);
GEN bqf_comp_tc(GEN q1, GEN q2, int tored, long prec);
GEN bqf_idelt(GEN D);
GEN bqf_identify(GEN ncgp_lexic, GEN q, GEN rootD, int Dsign);
GEN bqf_identify_tc(GEN ncgp_lexic, GEN q, long prec);
GEN bqf_lexicind_tobasis(GEN orders, long ind);
GEN bqf_ncgp(GEN D, long prec);
GEN bqf_ncgp_lexic(GEN D, long prec);
GEN bqf_pow(GEN q, GEN n);
GEN bqf_pow_red(GEN q, GEN n, GEN rootD, int Dsign);
GEN bqf_pow_tc(GEN q, GEN n, int tored, long prec);
GEN bqf_square(GEN q);
GEN bqf_square_red(GEN q, GEN rootD, int Dsign);
GEN bqf_square_tc(GEN q, int tored, long prec);
GEN ideal_tobqf(GEN numf, GEN ideal);

//REPRESENTATION OF NUMBERS BY BQFs
GEN bqf_reps(GEN q, GEN n, int proper, int half, long prec);
GEN bqf_reps_tc(GEN q, GEN n, int proper, int half, long prec);
GEN dbqf_reps(GEN qred, GEN D, GEN n, int proper, int half);
GEN ibqf_reps(GEN qorb, GEN qautom, GEN D, GEN rootD, GEN n, int proper, int half);
GEN sbqf_reps(GEN q, GEN D, GEN rootD, GEN n, int half);
GEN zbqf_reps(GEN A, GEN B, GEN n, int half);

//MORE REPRESENTATION OF NUMBERS
GEN bqf_bigreps(GEN q, GEN n, long prec);
GEN bqf_bigreps_tc(GEN q, GEN n, long prec);
GEN bqf_linearsolve(GEN q, GEN n1, GEN lin, GEN n2, long prec);
GEN bqf_linearsolve_tc(GEN q, GEN n1, GEN lin, GEN n2, long prec);

//GENERAL CHECKING METHODS
void bqf_check(GEN q);
GEN bqf_checkdisc(GEN q);
void intmatrix_check(GEN mtx);


//farey.c
GEN fareydenom_depths(long n);
long fareydenom_dmode(long n);
long fareydepth(GEN r);
GEN fareyup(GEN r);


//visual.c

//DATA
GEN integerbin(GEN v, GEN binlen, GEN binstart);
GEN integerbin_cumu(GEN v, GEN binlen, GEN binstart);
GEN veccount(GEN v);
GEN vecsmallcount(GEN v);

//HISTOGRAMS
void hist_autocompile(GEN minx, GEN maxx, char *imagename, char *plotoptions, int open);
GEN hist_make(GEN data, char *imagename, int compilenew, int open, char *plotoptions, long prec);
GEN hist_tobins(GEN data, GEN minx, GEN maxx, GEN nbins, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec);
GEN hist_tobins_defaultbins(GEN data, GEN minx, GEN maxx, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec);
GEN hist_rebin(GEN data, GEN histdata, GEN nbins, long prec);
GEN hist_rerange(GEN data, GEN histdata, GEN minx, GEN maxx, long prec);
GEN hist_rescale(GEN data, GEN histdata, int scale, long prec);

//REGRESSIONS
GEN OLS(GEN X, GEN y, int retrsqr);
GEN OLS_nointercept(GEN X, GEN y, int retrsqr);
GEN OLS_single(GEN x, GEN y, int retrsqr);
GEN rsquared(GEN X, GEN y, GEN fit);

//TEX
void tex_compile(char *imagename, int open);
void tex_recompile(GEN data);