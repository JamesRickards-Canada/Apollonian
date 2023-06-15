print("\n\nType '?apol' for help.\n\n");
addhelp(apol, "For each package P, call ?P to access a basic description and list of methods. Installed packages:\n apollonian\n base\n geometry\n quadratic\n visual");
parigp_version=version();
apol_library=strprintf("./libapol-%d-%d.so", parigp_version[1], parigp_version[2]);

\\apol.c

\\ACP=Apollonian circle packaing
	\\BASIC METHODS
		install("apol_check","iG",,apol_library);
		addhelp(apol_check, "Input v, a length 4 integral vector.\n Retuns 1 if this is a Descartes quadruple, i.e. if 2(a^2+b^2+c^2+d^2)=(a+b+c+d)^2.");
		install("apol_extdepth","lG",,apol_library);
		addhelp(apol_extdepth,"Input v, a Descartes quadruple.\n Returns the external depth of the quadruple v, i.e. the minimal number of swaps to reach a quadruple with negative curvature.");
		install("apol_getmatrices","",,apol_library);
		addhelp(apol_getmatrices, "Returns [S1, S2, S3, S4, K], where Si generate the Apollonian group, and K*[n,A,B,C]~=theta([A, B, C]).");
		install("apol_getobstructions","",,apol_library);
		addhelp(apol_getobstructions,"Returns the possible classes modulo 24 of an Apollonian circle packing.");
		install("apol_mod24","G",,apol_library);
		addhelp(apol_mod24,"Input v, a Descartes quadruple.\n Returns the set of curvatures modulo 24 possible in the correponding ACP. There are 6 possible primitive sets.");
		install("apol_move","GG",,apol_library);
		addhelp(apol_move, "Inputs v, ind: Descartes quadruple v, 1<=ind<=4 or a vector/vecsmall of integers between 1 and 4.\n Returns the Descartes quadruple where we replace circle ind with the other possible circle (applying this left to right if ind is a vector/vecsmall).");
		install("apol_qf","GD1,L,",,apol_library);
		addhelp(apol_qf, "Inputs v, {ind=1}: Descartes quadruple v, 1<=ind<=4.\n Returns a quadratic form q where the integers primitively represented by q are a+the curvatures of the circles surrounding the circle with curvature a, a=v[ind].");
		install("apol_red","GD0,L,",,apol_library);
		addhelp(apol_red,"Inputs v, {seq=0}: Descartes quadruple v.\n Reduces v and returns the reduction. If seq=1, also returns a vecsmall of the sequence of indices used to reach the reduced form.");
		install("apol_red_partial","GL",,apol_library);
		addhelp(apol_red_partial,"Inputs: Descartes quadruple v and positive integer maxsteps.\n We reduce v, doing at most maxsteps. In particular, the returned quadruple may not be reduced!");
		install("apol_thirdtangent","GGGGD1,L,",,apol_library);
		addhelp(apol_thirdtangent,"Inputs: circ1, circ2, c3, c4, {right=1}.\n Given two tangent circles circ1 and circ2 (given by [centre, radius, curvature] and curvatures c3 and c4 completing a Descartes quadruple, this computes the equation of the circle of curvature c3. We place it to the right of the ray from circ1 to circ2 if and only if right=1.");

		addhelp(ap_basic,"Installed methods:\napol_check, apol_extdepth, apol_getmatrices, apol_getobstructions, apol_mod24, apol_move, apol_qf, apol_red, apol_red_partial, apol_thirdtangent.");

	\\CREATION OF ACPS
		install("apol_make","GD1,L,D1,L,", ,apol_library);
		addhelp(apol_make,"Input q, {pos=1}, {red=1}; q a quadratic form of discriminant -4n^2.\n Returns the Descartes quadruple it corresponds to. If pos=0, the root quadruple starts with -n<0. Else, the form has a +n circle. If red=1 we reduce the form, otherwise we don't.");
		install("apol_makeall","GD1,L,p",,apol_library);
		addhelp(apol_makeall,"Inputs n, {red=1}.\n Returns all Descartes quadruples containing n. We reduce the quadruples iff red=1. The output has length h^{+/-}(-4n^2), and may contain Descartes quadruples in the same packing (if n appears multiple times in ''unique ways''.");

		addhelp(ap_make,"Installed methods:\napol_make, apol_makeall.");

	\\SEARCHING FOR CURVATURES
		install("apol_circles","GG",,apol_library);
		addhelp(apol_circles, "Inputs v, maxcurv: Descartes quadruple v, positive integer maxcurv.\n Computes all circles with curvature <=maxcurv in v. Returns the list, where each element is of the form [centre, radius, curvature], representing the circle centred at (x, y) with given radius/curvature. Negative radius/curvature corresponds to the outermost circle. The outer circle is centred at 0, the next largest circle is tangent at the top, and the third circle is on the left of them.");
		install("apol_circles_depth","GLD0,G,",,apol_library);
		addhelp(apol_circles_depth,"Inputs v, depth, {maxcurv=0}: Descartes quadruple v, positive integer maxcurv.\n Computes all circles with curvature <=maxcurv and depth<=depth in v. Returns the list, where each element is of the form [centre, radius, curvature], representing the circle centred at (x, y) with given radius/curvature. Negative radius/curvature corresponds to the outermost circle. The outer circle is centred at 0, the next largest circle is tangent at the top, and the third circle is on the left of them.");
		install("apol_curvatures","GGD0,L,",,apol_library);
		addhelp(apol_curvatures,"Inputs v, bound, {countsymm=0}: Descartes quadruple v, bound>=0, countsymm=0, 1.\n Returns a sorted list of curvatures of circles in v at most bound. If countsymm=1, symmetries of the packing are counted with their multiplicity.");
		install("apol_curvatures_depth","GLD0,G,",,apol_library);
		addhelp(apol_curvatures_depth,"Inputs v, depth, {bound=0}: Descartes quadruple v, depth>=1, bound>=0.\n Returns a sorted list of curvatures of circles at most depth circle replacements away from v. If bound>0, we also only save the ones of size at most bound.");
		install("apol_curvatures_layer","GLGD0,L,",,apol_library);
		addhelp(apol_curvatures_layer,"Inputs v, maxlayers, bound, {countsymm=0}: Descartes quadruple v, maxlayers>=1 bound>=0, countsymm=0, 1. Returns the curvatures at most bound in the first maxlayers outer layers of the ACP corresponding to v. If countsymm=1, symmetries are counted with their multiplicity. If v=[0,0,1,1], this does not work correctly.");
		install("apol_find","GGD0,L,",,apol_library);
		addhelp(apol_find,"Inputs v, N, {countsymm=0}: Descartes quadruple v, positive integer N, countsymm=0, 1.\n Returns the Descartes quadruples in the ACP of v containing N. If countsymm=1, symmetries of the packing are counted with their multiplicity.");
		install("apol_primes","GGD0,L,",,apol_library);
		addhelp(apol_primes,"Inputs v, bound, {countsymm=0}: Descartes quadruple v, bound>=0, countsymm=0, 1.\n Returns the prime curvatures at most bound in the ACP corresponding to v. If countsymm=1, symmetries are counted with their multiplicities.");
		install("apol_primes_layer","GLGD0,L,",,apol_library);
		addhelp(apol_primes_layer,"Inptus v, maxlayer, bound, {countsymm=0}: Descartes quadruple v, maxlayer>=1, bound>=0, countsymm=0, 1.\n Returns the prime curvatures in layer at most maxlayer with curvature at most bound. If countsymm=1, symmetries are counted with their multiplicities.");

		addhelp(ap_search,"Installed methods:\napol_circles, apol_curvatures, apol_curvatures_depth, apol_curvatures_layer, apol_find, apol_primes, apol_primes_layer.");

	\\STRIP PACKING METHODS
		install("apol_depthelt_circle","G",,apol_library);
		addhelp(apol_depthelt_circle,"Input L, an integer between 1 and 4, or a vector/vecsmall of integers between 1 and 4.\n Returns the circle/line corresponding to the depth element L. If L is an integer, this corresponds to Id_L. If L is a vecsmall/vector, this corresponds to S_L[1]*...*S_L[n].");
		install("apol_farey_allqf","G",,apol_library);
		addhelp(apol_farey_allqf,"Input q, a positive integer.\n Returns the set of primitive quadratic forms (Kate's construction) corresponding to the upside down Farey circle at x=p/q, over all p. The forms are all non-equivalent.");
		install("apol_farey_qf","GG",,apol_library);
		addhelp(apol_farey_qf,"Inputs p and q, positive integers.\n Returns the quadratic forms (Kate's construction) corresponding to the upside down Farey circle at x=p/q.");
		install("apol_stair","GD1,L,p",,apol_library);
		addhelp(apol_stair, "Input L, {format=1}: L is an integer between 1 and 4 or a vector/vecsmall of integers between 1 and 4, and format=0,1.\n Returns the data for the stair corresponding to the depth element L. If format=1, we return [t, a_W] as in my paper. If we don't intersect the fundamental domain or L=2, then we return 0. If format=0, we return [cutoff, height], and return [0, 0] if we don't intersect the fundamental domain.");
		install("apol_stairs","G",,apol_library);
		addhelp(apol_stairs,"Input tmax, a positive integer.\n Returns the stairs in the packing to tmax. We use the format of apol_stair with format=1, i.e. [t, a_W], hence we skip the identity element. We also don't combine stairs of the same height.");
		install("apol_strip_qf","GD0,L,",,apol_library);
		addhelp(apol_strip_qf,"Inputs L, {red=0}: L an integer between 1 and 4 or a vector/vecsmall of integers between 1 and 4, red=0 or 1.\n Returns the quadratic form corresponding to this circle in the strip packing, i.e. generating the curvatures of PSL(2, Z) times this circle.");

		addhelp(ap_strip,"Installed methods:\napol_depthelt_circle, apol_farey_allqf, apol_farey_qf, apol_stair, apol_strip_qf.");

	\\VISUALIZATION
		install("printcircles_desmos","vG",,apol_library);
		addhelp(printcircles_desmos,"Input c, a list of circles.\n Prints to the screen the list of equations of the circles, suitable for copying and pasting into Desmos.");
		install("printcircles_tex","GrD1,L,D0,L,D1,L,D1,L,p",,apol_library);
		addhelp(printcircles_tex,"Input c, imagename, {addnumbers=1}, {modcolours=0}, {compile=1}, {open=1}: list of circles c, string imagename, addnumbers/compile/open =0, 1, modcolours>=0.\n Prints the circles in c to the tex file images/build/imagename_build.tex. If addnumbers=1, we add the curvatures to each circle. If modcolours>=1, we colour the circles based on their remainders mod modcolours. If compile=1 we compile the file and move the output to images/imagename.pdf, and if open=1 (only valid with WSL), we also open the resulting image. Returns [imagename, open].");

		addhelp(ap_visual,"Installed methods:\nprintcircles_desmos, printcircles_tex.");

	\\SUPPORTING METHODS
		install("apol_words","L",,apol_library);
		addhelp(apol_words,"Input d>0, the depth.\n Returns all reduced words of length d in the Apollonian group, as a Vecsmall of 1-4's (no consecutive repeats).");

		addhelp(ap_support,"Installed methods:\napol_words.");

	\\SPECIALIZED METHODS
		install("apol_makeall_extdepths","Gp",,apol_library);
		addhelp(apol_makeall_extdepths,"Input n.\n Computes apol_makeall(n), and returns [d0,d1,...,dk], where the number of forms at depth i is given by di.");
		install("apol_makeall_small","GD1,L,p",,apol_library);
		addhelp(apol_makeall_small,"Input n, {red=1}: n a positive integer.\n Computes apol_makeall(n, red), takes the corresponding smallest curvatures and returns the sorted list. If red=1, we reduce the quadruples, and negate the curvature (since it is <=0 always). If red=0, we do not do so, since it could be positive or negative.");
		install("apol_makeall_small_maxsteps","GLp",,apol_library);
		addhelp(apol_makeall_small_maxsteps,"Input n, maxsteps.\n Does apol_makeall_small, but does NOT reduce the forms; instead, we reduce them by at most maxsteps only. We then return the raw data of the smallest element (without negating the curvature).");
		
		addhelp(ap_special,"Installed methods:\napol_makeall_extdepths, apol_makeall_small, apol_makeall_small_maxsteps.");

	\\GENERAL HELP
		addhelp(apollonian,"This package is a collection of methods used to deal with integral Apollonian circle packings. Subtopics:\n Basic methods (?ap_basic)\nCreation of ACPs (?ap_make)\nSearching for curvatures (?ap_search)\nStrip packing methods (?ap_strip)\nVisualization (?ap_visual)\nSupporting methods (?ap_support)\nSpecialized methods (?ap_special)");

/*mobius.c*/
	addhelp(geometry,"Lines:\n\t[slope, intercept]\n\ty=slope*x+intercept unless slope=oo, where it is x=intercept instead.\n\nCircles:\n\t[centre, radius, curvature]\n\tRadius > 0, but curvature can be < 0 to mean the 'outside' is the inside. This does not get preserved under the mobius function.\n\nGeometry methods:\n\tmobius.");

	/*SECTION 2: MOBIUS*/
		install("mobius","GGp");
		addhelp(mobius,"mobius(M, x): returns the image of x under the Mobius map M, which is a 2x2 invertible matrix. Allows points (including oo), circles, and lines as input. If the return is a circle, the curvature will always be positive.");

/*quadratic.c*/
	addhelp(quadratic,"Discriminant methods:\n\tdisclist, isdisc.\n\nBasic quadratic form methods:\n\tqfbapply, qfbapplyL, qfbapplyR, qfbapplyS, idealtoqfb, qfbtoideal.\n\nClass groups:\n\tlexind, qfbnarrow, qfbnarrowlex.");

	/*SECTION 1: DISCRIMINANT METHODS*/
		install(disclist,"GGD0,L,D0,G,",,apol_library);
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
		
	/*SECTION 3: CLASS GROUP*/
		install(lexind,"GL");
		addhelp(lexind,"lexind(v, ind): returns the index ind output of forvec(a=vector(#v, i, [0, v[i]-1]), print(a)), i.e. finds the corresponding lexicographic ordering element.");
		install(qfbnarrow,"Gp");
		addhelp(qfbnarrow,"qfbnarrow(D): returns the narrow class group C=Cl^+(D) in terms of quadratic forms. C[1] is the class number, C[2] are the orders of generators as a Vecsmall (largest to smallest, with each term dividing the previous one), and C[3] are the corresponding generators. Note that class number 1 will return [1, Vecsmall([1]), [idelt]], not [1, [], []], and the second return element is always a Vecsmall, not a vector.");
		install(qfbnarrowlex,"Gp");
		addhelp(qfbnarrowlex,"qfbnarrowlex(D): does qfbnarrow, except the third entry is the lexicographic ordering of representatives of the class group (with respect to the generators and their orders). Can pass in qfbnarrow(D) for D.");


\\visual.c
	\\DATA
		install("integerbin","GGD0,G,",,apol_library);
		addhelp(integerbin,"Inputs v, binlen, {binstart=0}.\n Assumes v is a sorted list of integers, and puts them into bins of length binlen, starting with binstart (assumed to be 0). Returns [binends, counts], with binends being the last number in the bin.");
		install("integerbin_cumu","GGD0,G,",,apol_library);
		addhelp(integerbin_cumu,"Inputs v, binlen, {binstart=0}.\n Assumes v is a sorted list of integers, and puts them into bins of length binlen, starting with binstart (assumed to be 0). Returns [binends, counts], with binends being the last number in the bin. This is cumulative, so counts is increasing.");
		install("veccount","G",,apol_library);
		addhelp(veccount,"Input v, a vector.\n Returns [uniq, count], where uniq is the sorted vector v with repeats removed, and count is the corresponding number of times they appear in v.");
		install("vecsmallcount","G",,apol_library);
		addhelp(vecsmallcount,"Input v, a vecsmall.\n Returns [uniq, count], where uniq is the sorted vecsmall v with repeats removed, and count is the corresponding number of times they appear in v.");
		install("ZV_countnonpos","lG",,apol_library);
		addhelp(ZV_countnonpos,"Input v, a sorted vector of integers.\n Returns the number of entries that are nonpositive.");
		
	\\HISTOGRAMS
		install("hist_make","GrD0,L,D0,L,Drp",,apol_library);
		addhelp(hist_make,"Inputs data, imagename, {compilenew=0}, {open=0}, {plotoptions=NULL}: sorted list of real numbers data, name of the tikz picture, {compilenew=0, 1}, {open=0, 1}, {plotoptions=NULL, string}.\n Automatically bins the data, and creates a pdf of the histogram using tikz and externalize. The output is in the folder /images, with the build file (named imagename_build.tex) being in the folder /images/build. If compilenew=0 assumes the LaTeX document to compile the plot is pre-made, and otherwise this method automatically writes it. If additionally, plotoptions!=NULL, this string is inserted in between \\begin{axis} and \\end{axis} in the LaTeX document (allowing one to customize how the histogram looks). If open=1, the pdf is automatically opened (only works with Linux subsystem for Windows). The returned value is used to modify the histogram, e.g. changing the bins, scaling it, and changing the range.");
		install("hist_rebin","GGGp",,apol_library);
		addhelp(hist_rebin,"Inputs data, histdata, nbins: the sorted data, the length 7 vector output of a hist_ method, the number of bins.\n Rebins the data according to the new number of bins, and updates histdata.");
		install("hist_rerange","GGGGp",,apol_library);
		addhelp(hist_rerange,"Inputs data, histdata, minx, maxx: the sorted data, the length 7 vector output of a hist_ method, minimum x-value, maximum x-value.\n Rebins the data according to the new minimum and maximum value. This is useful when there are outliers that skew the look of the graph. Returns the updated histdata.");
		install("hist_rescale","GGLp",,apol_library);
		addhelp(hist_rescale,"Inputs data, histdata, scale: the sorted data, the length 7 vector output of a hist_ method, and scale=0, 1.\n If scale=1 scales the data so the total area is 1, and if scale=0 uses the absolute count for the y-axis. Returns the updated histdata.");

	\\REGRESSIONS & PLOTS
		install("OLS","GGD1,L,",,apol_library);
		addhelp(OLS,"Inputs X, y, {retrsqr=1}:  m*n matrix X with top row being all 1's, length n column vector y, retrsqr=0, 1.\n Performs ordinary least squares regression on the data, where the n inputs are the columns of X, and the outputs are the entries of y. We must include a constant term, hence why the first row of X must be all 1's. If retrsqr=1, returns [pararms, R^2], and otherwise returns params, where params is the length m column vector of best fit parameters.");
		install("OLS_nointercept","GGD1,L,",,apol_library);
		addhelp(OLS_nointercept,"Inputs X, y, {retrsqr=1}: vector X, column vector y (of same length), retrsqr=0, 1.\n Performs ordinary least squares regression on the data assuming that y[i]=c*X[i], i.e. the y-intercept is 0. Returns c if retrsqr=0, or [c, R^2] otherwise.");
		install("OLS_single","GGD1,L,",,apol_library);
		addhelp(OLS_single,"Inputs x, y, {retrsqr=1}: vector x, column vector y, retrsqr=0, 1. Performs linear regression for a single variable (essentially a macro for OLS with y=mx+b.");
		install("rsquared","GGG",,apol_library);
		addhelp(rsquared,"Inputs X, y, fit: X and y data supplied to OLS, and fit the proposed fit (a column vector of parameters). This returns the R^2 value for this proposal.");

	\\TEX
		install("tex_recompile","vG",,apol_library);
		addhelp(tex_recompile,"Input data, the output of a tex image creation call.\n Recompiles the image, returning nothing. This is used when you edit the LaTeX document by hand.");

	\\GENERAL HELP
		addhelp(visual,"This package deals with visualizing data. Subtopics:\n Data (data)\n Histograms (hist)\n Regressions/plots (reg)\n Tex (tex)");
		addhelp(data,"integerbin, veccount, vecsmallcount.");
		addhelp(hist,"Installed methods:\n hist_make, hist_rebin, hist_rerange, hist_rescale.");
		addhelp(reg,"OLS, OLS_nointercept, OLS_single, rsquared.");
		addhelp(tex,"tex_recompile");

default(parisize, "4096M");\\Must come at the end