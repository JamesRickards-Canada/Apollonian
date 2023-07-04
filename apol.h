/*apol.c*/

/*SECTION 1: BASIC METHODS*/
GEN apol_admissiblesets();
int apol_check(GEN v, long prec);
long apol_chi(GEN v);
GEN apol_chi4(GEN v);
GEN apol_complete(GEN a, GEN b, GEN c, long prec);
long apol_extdepth(GEN v, long prec);
GEN apol_matrices();
GEN apol_mod24(GEN v);
GEN apol_move(GEN v, GEN command, long prec);
GEN apol_qf(GEN v, int ind);
GEN apol_red(GEN v, int seq, long prec);
GEN apol_red_partial(GEN v, long maxsteps, long prec);
GEN apol_type(GEN v, int chi);

/*SECTION 2: CREATION OF ACPS*/
GEN apol_make(GEN q, int pos, int red);
GEN apol_makeall(GEN n, int red, long prec);

/*SECTION 3: COMPUTING THE CIRCLES*/
GEN apol_circles(GEN v, GEN B, long depth, long prec);
void printcircles_desmos(GEN c);
GEN printcircles_tex(GEN c, char *imagename, int addnumbers, int modcolours, int compile, int open, long prec);

/*SECTION 4: STRIP PACKING METHODS*/

/*SECTION 5: SUPPORTING METHODS*/
GEN apol_words(int d);
GEN quarticresidue(GEN x, GEN y);


/*apol_fast.c*/

/*SECTION 1: MISSING CURVATURES*/
/*1: GP ACCESS*/
GEN apol_missing(GEN v, GEN B, int family, int load);
GEN apol_missing_load(GEN v, GEN B, int family);
GEN apol_missingfamilies(GEN v);

/*SECTION 2: SEARCHING FOR CURVATURES*/
GEN apol_curvatures(GEN v, GEN B, int tofile);
GEN apol_find(GEN v, GEN c, int all);


/*data.c*/

/*SECTION 1: DATA*/
GEN integerbin(GEN v, GEN blen, GEN bstart);
GEN integerbin_cumu(GEN v, GEN blen, GEN bstart);
GEN vecreduce(GEN v);
GEN vecsmallreduce(GEN v);

/*SECTION 2: HISTOGRAMS*/
GEN hist_make(GEN v, char *imagename, int compilenew, int open, char *plotoptions, long prec);
GEN hist_rebin(GEN v, GEN histdata, GEN nbins, long prec);
GEN hist_rerange(GEN v, GEN histdata, GEN minx, GEN maxx, long prec);
GEN hist_rescale(GEN v, GEN histdata, int scale, long prec);

/*SECTION 3: LINEAR REGRESSION*/
GEN OLS(GEN X, GEN y, int retrsqr);
GEN OLS_nointercept(GEN x, GEN y, int retrsqr);
GEN OLS_single(GEN x, GEN y, int retrsqr);
GEN rsquared(GEN X, GEN y, GEN fit);

/*SECTION 4: LATEX*/
GEN tex_makecolours(int ncol);
void tex_compile(char *imagename, int open);
void tex_recompile(GEN data);


/*geometry.c*/

/*SECTION 2: MOBIUS*/
GEN mobius(GEN M, GEN x, long prec);

/*SECTION 3: TOLERANCE*/
GEN deftol(long prec);
int toleq(GEN x, GEN y, GEN tol);
int toleq0(GEN x, GEN tol);


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
GEN qfbsos(GEN q);
GEN qfbsos1(GEN q);

/*SECTION 3: CLASS GROUP*/
GEN lexind(GEN v, long ind);
GEN qfbnarrow(GEN D, long prec);
GEN qfbnarrowlex(GEN D, long prec);
