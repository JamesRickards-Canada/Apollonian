print("\n\nType '?apol' for help.\n\n");
addhelp(apol, "For each package P, call ?P to access a basic description and list of methods. Installed packages:\napollonian\ndata\n geometry\nquadratic");
parigp_version=version();
apol_library=strprintf("./libapol-%d-%d.so", parigp_version[1], parigp_version[2]);

/*apol.c*/
	addhelp(apollonian,"Basic methods:\n\tapol_admissiblesets, apol_check, apol_chi2, apol_chi4, apol_complete, apol_extdepth, apol_matrices, apol_mod24, apol_move, apol_qf, apol_red, apol_red_partial, apol_type.\n\nCreation of ACPs:\n\tapol_make, apol_makeall.\n\nSearching the ACP:\n\tapol_curvatures, apol_find, apol_missing, apol_missing_load, apol_missingfamilies.\n\nCircles and pictures:\n\tapol_circles, printcircles_desmos, printcircles_tex.\n\nStrip packing methods:\n\t\n\nSupporting methods:\n\tapol_words, quarticresidue.");
	
	/*SECTION 1: BASIC METHODS*/
		install(apol_admissiblesets,"",,apol_library);
		addhelp(apol_admissiblesets,"apol_admissiblesets(): returns the possible classes modulo 24 of an Apollonian circle packing.");
		install(apol_check,"iGp");
		addhelp(apol_check,"apol_check(v): retuns 1 if this is a Descartes quadruple, i.e. if 2(a^2+b^2+c^2+d^2)=(a+b+c+d)^2. If the terms are inexact, we only check up to tolerance (half of the precision).");
		install(apol_chi2,"lG");
		addhelp(apol_chi2,"apol_chi2(v): returns the chi_2 value of the packing, which determines which quadratic obstruction the packing has.");
		install(apol_chi4,"G");
		addhelp(apol_chi4,"apol_chi4(v): returns the chi_4 value of the packing, which determines which quartic obstructions the packing has. Only valid for types (6, 1) and (6, 17); otherwise there will be no meaning..");
		install(apol_complete,"GDGDGp");
		addhelp(apol_complete,"apol_complete(a, {b}, {c}): given three curvatures, returns the Descartes quadruple containing them, choosing the one with minimal curvature, and sorting the result. Can pass a as a length 3 vector, or as three separate inputs.");
		install(apol_extdepth,"lGp");
		addhelp(apol_extdepth,"apol_extdepth(v): returns the external depth of the quadruple v, i.e. the minimal number of swaps to reach a quadruple with nonpositive curvature. If the packing is full-plane, the output will be meaningless as eventually floating point errors will add up.");
		install(apol_matrices,"");
		addhelp(apol_matrices, "apol_matrices(): returns [S1, S2, S3, S4, K], where Si generate the Apollonian group, and K*[n,A,B,C]~=theta([A, B, C]) (see the 'Apollonian Staircase' paper for a description of K).");
		install(apol_mod24,"G");
		addhelp(apol_mod24,"apol_mod24(v): returns the set of curvatures modulo 24 possible in the necessarily primitive integral ACP. There are 6 possible sets, obtainable with apol_admissiblesets(). See also apol_type.");
		install(apol_move,"GGp");
		addhelp(apol_move, "apol_move(v, ind): returns the Descartes quadruple where we replace circle ind with the other possible circle. Can also supply v as a vector/vecsmall, where we apply the sequence of moves from left to right.");
		install(apol_qf,"GD1,L,");
		addhelp(apol_qf, "apol_qf(v, {ind=1}): returns a quadratic form q where the curvatures of circles tangent to v[ind] are the values q(x, y)-v[ind] with gcd(x, y)=1. If ind=1 and v=[a, b, c, d], the formula is q=[a+b, a+b+c-d, a+c].");
		install(apol_red,"GD0,L,p");
		addhelp(apol_red,"apol_red(v, {seq=0}): reduces v and returns the reduction. If seq=1, also returns a Vecsmall of the sequence of indices used to reach the reduced form. Note that this sequence reads from left to right, and is compatible with apol_move. If the packing is full-plane, the output will be meaningless as eventually floating point errors will add up.");
		install(apol_red_partial,"GLp");
		addhelp(apol_red_partial,"apol_red_partial(v, maxsteps): reduce v, doing at most maxsteps. In particular, the returned quadruple may not be reduced!");
		install(apol_type,"GD2,L,");
		addhelp(apol_type,"apol_type(v, {chi=2}): returns the type of the primitive Apollonian circle packing v, i.e. the pair (a, b) where a is how many residues modulo 24 are hit, and b is the smallest residue coprime to 6, which uniquely determines the admissible set. If chi>=1, the third entry is the chi_2 value (determines the quadratic obstructions). If chi=2, the fourth entry is the chi_4 value, which is 0 unless the packing has type (6, 1) or (6, 17). In those cases, it determines the quartic obstructions.");

	/*SECTION 2: CREATION OF ACPS*/
		install(apol_make,"GD1,L,D1,L,");
		addhelp(apol_make,"apol_make(q, {pos=1}, {red=1}); returns a (sorted) Descartes quadruple corresponding to q, which is a positive definite binary quadratic form of discriminant -4n^2. If pos=0, the root quadruple starts with -n<0. Otherwise, the quadruple has a circle of curvature n>0. If red=1 we reduce the form, which may eliminate the curvature n circle from the quadruple if pos=1.");
		install(apol_makeall,"GD1,L,p");
		addhelp(apol_makeall,"apol_makeall(n, {red=1}): returns all Descartes quadruples containing n. We reduce the quadruples iff red=1. The output has length h^{+/-}(-4n^2), and may contain Descartes quadruples in the same packing (if n appears multiple times in 'unique ways'.");

	/*SECTION 3: COMPUTING THE CIRCLES*/
		install(apol_circles,"GGD0,L,DGp");
		addhelp(apol_circles, "apol_circles(v, bound, {maxdepth=0}, {maxxval}): computes all circles with curvature <= bound in v, also limiting to depth at most maxdepth if set, and ensuring all points circles have |x|<=maxxval if set (useful for strip, half-plane, and full-plane packings).\n\nCan also pass bound=[Bmin, Bmax], where we make sure they have curvatures in this interval instead.\nIf v is not a half or full plane packing, ensure that it is reduced (otherwise incorrect output WILL occur).\n\nWe always compute the starting four circles, even if one curvature there is > bounds. Returns the vector of equations, where each element is of the form [centre, radius, curvature], representing the circle centred at (x, y) with given radius/curvature. Negative radius/curvature corresponds to the outermost circle. The outer circle is centred at 0, the next largest circle is tangent at the top, and the third circle is on the left of them. For a strip, half-plane, or full-plane packing, you must set the depth to something non-zero. The output of this method can be used in printcircles_tex to display it in LaTeX.");
		install(printcircles_desmos,"vG");
		addhelp(printcircles_desmos,"printcircles_desmos(c): prints to the screen the list of equations of the circles in c, suitable for copying and pasting into Desmos.");
		install(printcircles_tex,"GrD1,L,D0,L,DGD1,L,D1,L,p");
		addhelp(printcircles_tex,"printcircles_tex(c, imagename, {addnumbers=1}, {modcolours=0}, {outerrad=3}, {compile=1}, {open=1}): prints the circles in c to the tex file images/build/imagename_build.tex. If addnumbers=1, we add the curvatures to each circle (only valid for integral ones). If modcolours>=1, we colour the circles based on their remainders mod modcolours (only valid for integral). If outerrad is given, this sets the radius of the largest circle to this many inches (defualt 3). If compile=1 we compile the file and move the output to images/imagename.pdf, and if open=1 (only valid with WSL), we also open the resulting image. Returns [imagename, open].");

	/*SECTION 4: STRIP PACKING METHODS*/

	/*SECTION 5: SUPPORTING METHODS*/
		install(apol_words,"L");
		addhelp(apol_words,"apol_words(d): returns all 4*3^(d-1) reduced words of length d in the Apollonian group, as Vecsmalls of 1-4's with no consecutive repeats.");
		install(quarticresidue,"GG");
		addhelp(quarticresidue,"quarticresidue(x, y): returns the quartic residue symbol [x/y] for the coprime Gaussian integers x, y with y odd. Does not check that y is odd or that x and y are coprime.");


/*apol_fast.c*/

	/*SECTION 1: MISSING CURVATURES*/
		/*1: GP ACCESS*/
		install(apol_missing,"GGD1,L,D1,L,");
		addhelp(apol_missing,"apol_missing(v, B, {remobstruct=1}, {load=1}): saves all the curvatures up to B in the correct residue classes modulo 24 but are missing from the packing to a file in the folder 'missing'. B can be given either as a positive integer, or as [B1, B2] to only look in the range B1<=c<=B2. The curvatures are separated by residue, each on a different line. If remobstruct=1, we also remove obstructions of the form ax^2 and ax^4. We then print the admissible set modulo 24, the quadratic obstructions, and then the quartic obstructions on the first three lines. If load=1, we load this output into a vector in gp.");
		install(apol_missing_load,"GGD1,L,");
		addhelp(apol_missing_load,"apol_missing_load(v, B, {family=1}): loads the saved curvature data. If family=1, the first entry is the vector of modulo 24 residues, then the quadratic obstructions found in each class, then the quartic obstructions found in each class, then the missing curvatures that are not in one of the obstruction classes. Creates an error if the data does not exist.");
		install(apol_missingfamilies,"G");
		addhelp(apol_missingfamilies,"apol_missingfamilies(v): returns the quadratic and quartic obstructions for v, in order. Each element of the obstruction is the vector of u's such that ux^2 or ux^4 is completely absent from the corresponding residue class (where gcd(x, 1/2/3/6)=1 to ensure we land in that class).");
		
	/*SECTION 2: SEARCHING FOR CURVATURES*/	
		install(apol_find,"GGD1,L,");
		addhelp(apol_find,"apol_find(v, c, {all=1}): returns the Descartes quadruples in the ACP of v containing curvature c. Note that two distinct circles in the picture may correspond to the same Descartes quadruple (due to symmetries), and they are counted once only. The output does not necessarily preserve the ordering in v, i.e. apol_find([5, 45, 8, 8], 45, 0) returns [8, 45, 8, 5]. If all=1 returns the vector of all quadruples, and if all=0 we stop at the first one (returning 0 if none found). Both the packing and curvature c must be of long size, i.e. at most 2^63-1.");
		install(apol_curvatures,"GGD0,L,");
		addhelp(apol_curvatures,"apol_curvatures(v, B, {tofile=0}): finds all positive curvatures in the packing for v up to the bound B (counted up to symmetry). B can be given either as a positive integer, or as [B1, B2] to only look in the range B1<=c<=B2. If tofile = 1, results are output to the folder ./curv_freq, and are labelled by v, Bmin, Bmax, and the residue class (one file per residue class). Each line of the file contains the frequency of the corresponding curvature (runs from the smallest to largest admissible curvature in this residue class in [Bmin, Bmax]). If tofile = 0, instead returns [curvs, freqs], where curvs is a Vecsmall of the curvatures that do appear, and freqs is a Vecsmall of their frequency (all >= 1). If tofile = 2, do both. Note that all curvatures must fit into unsigned longs (2^64 - 1 on 64-bit hardware), and the frequencies must be at most unsigned int sized (limiting to 2^32 - 1 on 64-bit hardware).");
    install(apol_quadruples,"GG");
    addhelp(apol_quadruples,"apol_quadruples(v, B): finds all quadruples (up to symmetry) in the packing for v where all entries of the quadruples are bounded by B. These quadruples are all reachable from v without permuting the entries.");


/*data.c*/
	addhelp(data,"Data methods:\n\tintegerbin, integerbin_cumu, vecreduce, vecsmallreduce.\n\nHistograms:\n\thist_make, hist_rebin, hist_rereange, hist_rescale.\n\nLinear regression:\n\tOLS, OLS_nointercept, OLS_single, rsquared.\n\nLatex:\n\ttex_recompile.");

	/*SECTION 1: DATA*/
		install(integerbin,"GGD0,G,");
		addhelp(integerbin,"integerbin(v, blen, {bstart=0}): assumes v is a sorted vector of integers, and puts them into bins of length blen, starting with bstart (defaults to 0). Returns [bends, counts], with bends being the ends of each bin.");
		install(integerbin_cumu,"GGD0,G,");
		addhelp(integerbin_cumu,"integerbin_cumu(v, blen, {bstart=0}): assumes v is a sorted vector of integers, and puts them into bins of length blen, starting with bstart (defaults to 0). Returns [bends, counts], with bends being the ends of each bin. This is cumulative, so counts is nondecreasing.");
		install(vecreduce,"G");
		addhelp(vecreduce,"vecreduce(v): returns [uniq, count], where uniq is the sorted vector v with repeats removed, and count is the corresponding number of times they appear in v (as a Vecsmall).");
		install(vecsmallreduce,"G");
		addhelp(vecsmallreduce,"vecsmallreduce(v): returns [uniq, count], where uniq is the sorted Vecsmall v with repeats removed, and count is the corresponding number of times they appear in v (as a Vecsmall).");
		
	/*SECTION 2: HISTOGRAMS*/
		install(hist_make,"GrD0,L,D0,L,Drp");
		addhelp(hist_make,"hist_make(v, imagename, {compilenew=0}, {open=0}, {plotoptions=NULL}): assuming v is a sorted vector of real numbers, this bins the data and makes a histogram in LaTeX using tikz and externalize. The output is in the folder /images, with the build file (named imagename_build.tex) being in the folder /images/build. If compilenew=0, assumes the LaTeX document to compile the plot is pre-made, and otherwise this method automatically writes it. If additionally plotoptions!=NULL, this string is inserted in between \\begin{axis} and \\end{axis} in the LaTeX document (allowing one to customize how the histogram looks). If open=1, the pdf is automatically opened (only works with Linux subsystem for Windows). The returned value is used to modify the histogram, e.g. changing the bins, scaling it, and changing the range.");
		install(hist_rebin,"GGGp");
		addhelp(hist_rebin,"hist_rebin(v, histdata, nbins): rebins the data according to the new number of bins, and updates histdata (which is the output of any hist_ method).");
		install(hist_rerange,"GGGGp");
		addhelp(hist_rerange,"hist_rerange(v, histdata, minx, maxx): rebins v according to the new minimum and maximum value. This is useful when there are outliers that skew the look of the graph. Returns the updated histdata.");
		install(hist_rescale,"GGLp");
		addhelp(hist_rescale,"hist_rescale(v, histdata, scale): if scale=1 scales the data so the total area is 1, and if scale=0 uses the absolute count for the y-axis. Returns the updated histdata.");

	/*SECTION 3: LINEAR REGRESSION*/
		install(OLS,"GGD1,L,");
		addhelp(OLS,"OLS(X, y, {retrsqr=1}): performs ordinary least squares regression on the data, with the inputs being the n columns of the matrix X, and the outputs being the entries of the column vector y. We include a constant term, so the first row of X must be all 1's. If retrsqr=1, returns [params, R^2], and otherwise returns params, where params is the column vector of best fit parameters.");
		install(OLS_nointercept,"GGD1,L,");
		addhelp(OLS_nointercept,"OLS_nointercept(x, y, {retrsqr=1}): performs ordinary least squares regression assuming that y[i]=c*X[i] (where X and y are both (column) vectors), i.e. the y-intercept is 0. Returns c if retrsqr=0, or [c, R^2] otherwise.");
		install(OLS_single,"GGD1,L,");
		addhelp(OLS_single,"OLS_single(x, y, {retrsqr=1}): performs linear regression for a single variable (essentially a macro for OLS with y=mx+b). Requires y to be a column vector and x is either a vector or column vector.");
		install(rsquared,"GGG");
		addhelp(rsquared,"rsquared(X, y, fit): returns the R^2 value for the proposed linear regression, where the input is X, output is y, and fit is the proposed parameters.");

	/*SECTION 4: LATEX*/
		install(tex_recompile,"vG");
		addhelp(tex_recompile,"tex_recompile(data): given the output of a tex image creation call, Recompiles the image, returning nothing. This is used when you edit the LaTeX document by hand.");


/*geometry.c*/
	addhelp(geometry,"Lines:\n\t[slope, intercept]\n\ty=slope*x+intercept unless slope=oo, where it is x=intercept instead.\n\nCircles:\n\t[centre, radius, curvature]\n\tRadius > 0, but curvature can be < 0 to mean the 'outside' is the inside. This does not get preserved under the mobius function.\n\nGeometry methods:\n\tmobius.");

	/*SECTION 2: MOBIUS*/
		install(mobius,"GGp");
		addhelp(mobius,"mobius(M, x): returns the image of x under the Mobius map M, which is a 2x2 invertible matrix. Allows points (including oo), circles, and lines as input. If the return is a circle, the curvature will always be positive.");


/*quadratic.c*/
	addhelp(quadratic,"Discriminant methods:\n\tdisclist, isdisc.\n\nBasic quadratic form methods:\n\tqfbapply, qfbapplyL, qfbapplyR, qfbapplyS, idealtoqfb, qfbtoideal, qfbsos.\n\nClass groups:\n\tlexind, qfbnarrow, qfbnarrowlex.");

	/*SECTION 1: DISCRIMINANT METHODS*/
		install(disclist,"GGD0,L,D0,G,");
		addhelp(disclist, "disclist(D1, D2, {fund=0}, {cop=0}): returns the set of discriminants between D1 and D2, inclusive. If fund=1, only returns fundamental discriminants. If cop!=0, only returns discriminants coprime to cop.");
		install(isdisc,"iG");
		addhelp(isdisc,"isdisc(D): returns 1 if D is a discriminant, 0 else.");

	/*SECTION 2: BASIC QUADRATIC FORM METHODS*/
		install(qfb_apply_ZM,"GG",qfbapply);/*In PARI but not installed.*/
		addhelp(qfbapply,"qfbapply(q, g): returns the quadratic form formed by g acting on q, where g is a matrix with integral coefficients.");
		install(qfbapplyL,"GD1,G,");
		addhelp(qfbapplyL,"qfbapplyL(q, {n=1}): returns L^n acting on q, where L=[1, 1;0, 1].");
		install(qfbapplyR,"GD1,G,");
		addhelp(qfbapplyR,"qfbapplyR(q, {n=1}): returns R^n acting on q, where R=[1, 0;1, 1].");
		install(qfbapplyS,"G");
		addhelp(qfbapplyS,"qfbapplyS(q): returns S acting on q, where S=[0, 1;-1, 0].");
		install(idealtoqfb,"GG");
		addhelp(idealtoqfb,"idealtoqfb(nf, x): given a quadratic number field, returns the primitive integral binary quadratic form corresponding to the fractional ideal x, positive definite if the field is imaginary. We also assume that we are working in the maximal ideal.");
		install(qfbtoideal,"GG");
		addhelp(qfbtoideal,"qfbtoideal(nf, q): given a quadratic number field, returns the fractional ideal corresponding to the primitive integral binary quadratic form q (positive definite if the field is imaginary). We also assume that its discriminant is fundamental.");
		install(qfbsos,"G");
		addhelp(qfbsos,"qfbsos(q): given a quadratic form of discriminant -4n^2 < 0, writes it as a sum of two squares (ax+by)^2+(cx+dy)^2, returning [a, b;c, d].");
		
	/*SECTION 3: CLASS GROUP*/
		install(lexind,"GL");
		addhelp(lexind,"lexind(v, ind): returns the index ind output of forvec(a=vector(#v, i, [0, v[i]-1]), print(a)), i.e. finds the corresponding lexicographic ordering element.");
		install(qfbnarrow,"Gp");
		addhelp(qfbnarrow,"qfbnarrow(D): returns the narrow class group C=Cl^+(D) in terms of quadratic forms. C[1] is the class number, C[2] are the orders of generators as a Vecsmall (largest to smallest, with each term dividing the previous one), and C[3] are the corresponding generators. Note that class number 1 will return [1, Vecsmall([1]), [idelt]], not [1, [], []], and the second return element is always a Vecsmall, not a vector.");
		install(qfbnarrowlex,"Gp");
		addhelp(qfbnarrowlex,"qfbnarrowlex(D): does qfbnarrow, except the third entry is the lexicographic ordering of representatives of the class group (with respect to the generators and their orders). Can pass in qfbnarrow(D) for D.");

default(parisize, "4096M");\\Must come at the end