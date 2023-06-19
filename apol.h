//apol.c

//BASIC METHODS
int apol_check(GEN v);
long apol_extdepth(GEN v);
GEN apol_getmatrices();
GEN apol_getobstructions();
GEN apol_mod24(GEN v);
GEN apol_move_1(GEN v, int ind);
GEN apol_move_1GEN(GEN v, int ind);
GEN apol_move_batch(GEN v, GEN bat);
GEN apol_move_batchGEN(GEN v, GEN bat);
GEN apol_qf(GEN v, int ind);
GEN apol_red(GEN v, int seq);
GEN apol_red_partial(GEN v, long maxsteps);

//CREATION OF ACPS
GEN apol_make(GEN q, int pos, int red);
GEN apol_makeall(GEN n, int red, long prec);

//SEARCHING FOR CURVATURES
GEN apol_circles(GEN v, GEN maxcurv);
GEN apol_circles_depth(GEN v, int depth, GEN maxcurv);
GEN apol_curvatures(GEN v, GEN bound, int countsymm);
GEN apol_curvatures_depth(GEN v, int depth, GEN bound);
GEN apol_curvatures_layer(GEN v, int maxlayers, GEN bound, int countsymm);
GEN apol_find(GEN v, GEN N, int countsymm);
GEN apol_primes(GEN v, GEN bound, int countsymm);
GEN apol_primes_layer(GEN v, int maxlayers, GEN bound, int countsymm);
GEN apol_thirdtangent(GEN circ1, GEN circ2, GEN c3, GEN c4, int right);


GEN apol_missing(GEN v, GEN B, int family, int load);
GEN apol_missing_load(GEN v, GEN B, int family);

//STRIP PACKING METHODS
GEN apol_depthelt_circle(GEN L);
GEN apol_farey_allqf(GEN q);
GEN apol_farey_qf(GEN p, GEN q);
GEN apol_stair(GEN L, int format, long prec);
GEN apol_stairs(GEN tmax);
GEN apol_strip_qf(GEN L, int red);

//VISUALIZATION
void printcircles_desmos(GEN c);
GEN printcircles_tex(GEN c, char *imagename, int addnumbers, int modcolours, int compile, int open, long prec);

//SUPPORTING METHODS
GEN apol_words(int d);

//LISTS OF VARIABLE LENGTH TO DELETE
GEN vecsmalllist_append(GEN v, long *vind, long *vlen, long x);

//SPECIALIZED METHODS
GEN apol_makeall_extdepths(GEN n, long prec);
GEN apol_makeall_small(GEN n, int red, long prec);
GEN apol_makeall_small_maxsteps(GEN n, long maxsteps, long prec);


/*data.c*/

/*SECTION 1: DATA*/
GEN integerbin(GEN v, GEN blen, GEN bstart);
GEN integerbin_cumu(GEN v, GEN blen, GEN bstart);
GEN vecreduce(GEN v);
GEN vecsmallreduce(GEN v);


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
GEN tex_makecolours(int ncol);
void tex_compile(char *imagename, int open);
void tex_recompile(GEN data);




/*geometry.c*/

/*SECTION 2: MOBIUS*/
GEN mobius(GEN M, GEN x, long prec);

/*SECTION 3: TOLERANCE*/
GEN deftol(long prec);
int toleq(GEN x, GEN y, GEN tol);

/*quadratic.c*/

/*SECTION 1: DISCRIMINANT METHODS*/
GEN disclist(GEN D1, GEN D2, int fund, GEN cop);
int isdisc(GEN D);

/*SECTION 2: BASIC QUADRATIC FORM METHODS*/
GEN qfbapplyL(GEN q, GEN n);
GEN qfbapplyR(GEN q, GEN n);
GEN qfbapplyS(GEN q);
GEN idealtoqfb(GEN nf, GEN x);
GEN qfbtoideal(GEN nf, GEN q);

/*SECTION 3: CLASS GROUP*/
GEN lexind(GEN v, long ind);
GEN qfbnarrow(GEN D, long prec);
GEN qfbnarrowlex(GEN D, long prec);
