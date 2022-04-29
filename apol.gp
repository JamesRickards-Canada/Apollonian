print("\n\nType '?apol' for help.\n\n");
addhelp(apol, "For each package P, call ?P to access a basic description and list of methods. Installed packages: \n apollonian \n base \n bqf \n visual");

\\apol.c

\\ACP=Apollonian circle packaing
	\\BASIC METHODS
		install("apol_check","iG",,"./libapol.so");
		addhelp(apol_check, "Input v, a length 4 integral vector.\n Retuns 1 if this is a Descartes quadruple, i.e. if 2(a^2+b^2+c^2+d^2)=(a+b+c+d)^2.");
		install("apol_extdepth","lG",,"./libapol.so");
		addhelp(apol_extdepth,"Input v, a Descartes quadruple.\n Returns the external depth of the quadruple v, i.e. the minimal number of swaps to reach a quadruple with negative curvature.");
		install("apol_getmatrices","",,"./libapol.so");
		addhelp(apol_getmatrices, "Returns [S1, S2, S3, S4, K], where Si generate the Apollonian group, and K*[n,A,B,C]~=theta([A, B, C]).");
		install("apol_getobstructions","",,"./libapol.so");
		addhelp(apol_getobstructions,"Returns the possible classes modulo 24 of an Apollonian circle packing.");
		install("apol_mod24","G",,"./libapol.so");
		addhelp(apol_mod24,"Input v, a Descartes quadruple.\n Returns the set of curvatures modulo 24 possible in the correponding ACP. There are 6 possible primitive sets.");
		install("apol_move","GG",,"./libapol.so");
		addhelp(apol_move, "Inputs v, ind: Descartes quadruple v, 1<=ind<=4 or a vector/vecsmall of integers between 1 and 4.\n Returns the Descartes quadruple where we replace circle ind with the other possible circle (applying this left to right if ind is a vector/vecsmall).");
		install("apol_qf","GD1,L,",,"./libapol.so");
		addhelp(apol_qf, "Inputs v, {ind=1}: Descartes quadruple v, 1<=ind<=4.\n Returns a quadratic form q where the integers primitively represented by q are a+the curvatures of the circles surrounding the circle with curvature a, a=v[ind].");
		install("apol_red","GD0,L,",,"./libapol.so");
		addhelp(apol_red,"Inputs v, {seq=0}: Descartes quadruple v.\n Reduces v and returns the reduction. If seq=1, also returns a vecsmall of the sequence of indices used to reach the reduced form.");
		install("apol_red_partial","GL",,"./libapol.so");
		addhelp(apol_red_partial,"Inputs: Descartes quadruple v and positive integer maxsteps.\n We reduce v, doing at most maxsteps. In particular, the returned quadruple may not be reduced!");

		addhelp(ap_basic,"Installed methods:\napol_check, apol_extdepth, apol_getmatrices, apol_getobstructions, apol_mod24, apol_move, apol_qf, apol_red, apol_red_partial.");

	\\CREATION OF ACPS
		install("apol_make","GD1,L,D1,L,", ,"./libapol.so");
		addhelp(apol_make,"Input q, {pos=1}, {red=1}; q a quadratic form of discriminant -4n^2.\n Returns the Descartes quadruple it corresponds to. If pos=0, the root quadruple starts with -n<0. Else, the form has a +n circle. If red=1 we reduce the form, otherwise we don't.");
		install("apol_makeall","GD1,L,p",,"./libapol.so");
		addhelp(apol_makeall,"Inputs n, {red=1}.\n Returns all Descartes quadruples containing n. We reduce the quadruples iff red=1. The output has length h^{+/-}(-4n^2), and may contain Descartes quadruples in the same packing (if n appears multiple times in ''unique ways''.");

		addhelp(ap_make,"Installed methods:\napol_make, apol_makeall.");

	\\SEARCHING FOR CURVATURES
		install("apol_circles","GG",,"./libapol.so");
		addhelp(apol_circles, "Inputs v, maxcurv: Descartes quadruple v, positive integer maxcurv.\n Computes all circles with curvature <=maxcurv in v. Returns the list, where each element is of the form [centre, radius, curvature], representing the circle centred at (x, y) with given radius/curvature. Negative radius/curvature corresponds to the outermost circle. The outer circle is centred at 0, the next largest circle is tangent at the top, and the third circle is on the left of them.");
		install("apol_curvatures","GGD0,L,",,"./libapol.so");
		addhelp(apol_curvatures,"Inputs v, bound, {countsymm=0}: Descartes quadruple v, bound>=0, countsymm=0, 1.\n Returns a sorted list of curvatures of circles in v at most bound. If countsymm=1, symmetries of the packing are counted with their multiplicity.");
		install("apol_curvatures_depth","GLD0,G,",,"./libapol.so");
		addhelp(apol_curvatures_depth,"Inputs v, depth, {bound=0}: Descartes quadruple v, depth>=1, bound>=0.\n Returns a sorted list of curvatures of circles at most depth circle replacements away from v. If bound>0, we also only save the ones of size at most bound.");
		install("apol_find","GGD0,L,",,"./libapol.so");
		addhelp(apol_find,"Inputs v, N, {countsymm=0}: Descartes quadruple v, positive integer N, countsymm=0, 1.\n Returns the Descartes quadruples in the ACP of v containing N. If countsymm=1, symmetries of the packing are counted with their multiplicity.");
		
		install("apol_orbit_layers","GLG",,"./libapol.so");
		addhelp(apol_orbit_layers,"Inputs v, maxlayers, bound.\n Returns the curvatures in the first maxlayers layers up to bound bound, with respect to the circle v[1].");
		install("apol_orbit_primes","GLG",,"./libapol.so");
		addhelp(apol_orbit_primes,"Inputs v, maxlayers, bound.\n Returns the prime curvatures in the first maxlayers layers up to bound bound, with respect to the circle v[1]. Does not include curvatures 2 or 3, and removes repeats.");

		addhelp(ap_search,"Installed methods:\napol_circles, apol_curvatures, apol_curvatures_depth, apol_find.");

	\\STRIP PACKING METHODS
		install("apol_depthelt_circle","G",,"./libapol.so");
		addhelp(apol_depthelt_circle,"Input L, an integer between 1 and 4, or a vector/vecsmall of integers between 1 and 4.\n Returns the circle/line corresponding to the depth element L. If L is an integer, this corresponds to Id_L. If L is a vecsmall/vector, this corresponds to S_L[1]*...*S_L[n].");
		install("apol_strip_qf","GD0,L,",,"./libapol.so");
		addhelp(apol_strip_qf,"Inputs L, {red=0}: L an integer between 1 and 4 or a vector/vecsmall of integers between 1 and 4, red=0 or 1.\n Returns the quadratic form corresponding to this circle in the strip packing, i.e. generating the curvatures of PSL(2, Z) times this circle.");

		addhelp(ap_strip,"Installed methods:\napol_depthelt_circle, apol_strip_qf.");

	\\VISUALIZATION
		install("printcircles_desmos","vG",,"./libapol.so");
		addhelp(printcircles_desmos,"Input c, a list of circles.\n Prints to the screen the list of equations of the circles, suitable for copying and pasting into Desmos.");
		install("printcircles_tex","GrD1,L,D0,L,D1,L,D1,L,p",,"./libapol.so");
		addhelp(printcircles_tex,"Input c, imagename, {addnumbers=1}, {modcolours=0}, {compile=1}, {open=1}: list of circles c, string imagename, addnumbers/compile/open =0, 1, modcolours>=0.\n Prints the circles in c to the tex file images/build/imagename_build.tex. If addnumbers=1, we add the curvatures to each circle. If modcolours>=1, we colour the circles based on their remainders mod modcolours. If compile=1 we compile the file and move the output to images/imagename.pdf, and if open=1 (only valid with WSL), we also open the resulting image. Returns [imagename, open].");

		addhelp(ap_visual,"Installed methods:\n.");

	\\SUPPORTING METHODS
		install("apol_words","L",,"./libapol.so");
		addhelp(apol_words,"Input d>0, the depth.\n Returns all reduced words of length d in the Apollonian group, as a Vecsmall of 1-4's (no consecutive repeats).");

		addhelp(ap_support,"Installed methods:\n.");

	\\SPECIALIZED METHODS
		install("apol_ncgp_depths","Gp",,"./libapol.so");
		addhelp(apol_ncgp_depths,"Input n.\n Computes apol_makeall(n), and returns [d0,d1,...,dk], where the number of forms at depth i is given by di.");
		install("apol_ncgp_smallcurve","GD1,L,p",,"./libapol.so");
		addhelp(apol_ncgp_smallcurve,"Input n, {red=1}: n a positive integer.\n Computes apol_makeall(n, 1, red), takes the corresponding smallest curvatures and returns the sorted list. If red=1, we reduce the quadruples, and negate the curvature (since it is <=0 always). If red=0, we do not do so, since it could be positive or negative. The return vector has length h^{+/-}(-4n^2), and each entry is from the set {0,1,...,n-1}.");
		install("apol_ncgp_smallcurve_bsteps","GLp",,"./libapol.so");
		addhelp(apol_ncgp_smallcurve_bsteps,"Input n, maxsteps.\n Does apol_ncgp_smallcurve, but does NOT reduce the forms; instead, we reduce them by at most maxsteps only. We then return the raw data of the smallest element (without negating the curvature).");
		
		addhelp(ap_special,"Installed methods:\n.");

	\\GENERAL HELP
		addhelp(apollonian,"This package is a collection of methods used to deal with integral Apollonian circle packings. Subtopics:\n Basic methods (?ap_basic)\nCreation of ACPs (?ap_make)\nSearching for curvatures (?ap_search)\nStrip packing methods (?ap_strip)\nVisualization (?ap_visual)\nSupporting methods (?ap_support)\nSpecialized methods (?ap_special)");

\\base.c

	\\INFINITY
		install("addoo","GG","addoo","./libapol.so");
		addhelp(addoo, "Inputs a,b real numbers or oo.\n Outputs a+b, where if a or b is +/- oo, returns that back. Note that this will make oo+-oo=oo and -oo+oo=-oo");
		install("divoo","GG","divoo","./libapol.so");
		addhelp(divoo, "Inputs a,b, real numbers or oo.\n Outputs a/b, where oo is output if a=oo and b>=0 or a=-oo and b<0 or b=0 and a>=0 (outputs -oo under analogous assumptions).");

	\\LINEAR ALGEBRA
		install("lin_intsolve_tc","GGG","lin_intsolve","./libapol.so");
		addhelp(lin_intsolve, "Inputs A,B,n integers.\n Outputs the general integral solutions to Ax+By=n. The format is [[s1,s2],[x0,y0]], where the general solution is x=s1*t+x0, y=s2*t+y0 for t an integer. The output is also reduced, i.e. gcd(s1,s2)=1. If A=B=0 or there are no integer solutions, returns 0.");
		install("mat3_complete_tc", "GGG", "mat3_complete", "./libapol.so");
		addhelp(mat3_complete, "Inputs A,B,C integers with gcd 1.\n Outputs a 3x3 integer matrix with determinant 1 whose top row is [A, B, C].");
	
	\\PRIMES
		install("primes_mod","GGL",,"./libapol.so");
		addhelp(primes_mod,"Input range, residues, modulus.\n Returns the primes in the given range that have the given residues modulo the modulus. Range=[a, b] where a<b are integers, and residues can either be a single integer or a vector of integers.");
	
	\\GENERAL HELP
		addhelp(base,"This package is a collection of miscellaneous methods that may be useful in a variety of settings, and not just for the programs they were originally created for \n Subtopics: \n Infinity (inf) \n Linear algebra (la) \n Primes (pr)");
		addhelp(inf,"addoo, divoo.");
		addhelp(la,"lin_intsolve, mat3_complete.");
		addhelp(pr,"primes_mod");


\\bqf.c

	\\DISCRIMINANT METHODS
		install("disclist","GGD0,L,D0,G,","disclist","./libapol.so");
		addhelp(disclist, "Inputs d1, d2, {fund=0}, {cop=0}: d1 and d2 integers with d1<=d2, fund=0,1, cop integer.\n Returns the set of proper discriminants between d1 and d2, inclusive. If fund=1, only returns fundamental discriminants. If cop!=0, only returns discriminants coprime to cop.");
		install("discprimeindex","G","discprimeindex","./libapol.so");
		addhelp(discprimeindex, "Inputs: D, a proper discriminant.\n Returns all prime divisors p of D for which D/p^2 is a proper discriminant.");
		install("discsuperorders","G","discsuperorders","./libapol.so");
		addhelp(discsuperorders,"Input D, a proper discriminant. Returns the vector of divisors D' of D such that D' is a discriminant and D/D' is a square.\n ");
		install("isdisc","iG","isdisc","./libapol.so");
		addhelp(isdisc, "Inputs: D a real number.\n Returns 1 if D is a proper discriminant, and 0 otherwise.");
		install("pell","G","pell","./libapol.so");
		addhelp(pell, "Inputs: D a positive discriminant.\n Returns [T, U], which is the smallest positive integer solution to T^2-DU^2=4 (and so (T+Usqrt(D))/2 is the fundamental unit in O_D).");
		install("posreg","Gp","posreg","./libapol.so");
		addhelp(posreg, "Inputs: D a positive discriminant.\n Returns the positive regulator of O_D, i.e. the logarithm of the fundamental totally positive unit.");
		install("quadroot","G","quadroot","./libapol.so");
		addhelp(quadroot, "Input D, a non-square integer.\n Returns sqrt(D) of type t_QUAD.");

	\\BASIC OPERATIONS ON BINARY QUADRATIC FORMS
		install("bqf_automorph_tc","G","bqf_automorph","./libapol.so");
		addhelp(bqf_automorph, "Inputs: q a BQF.\n Returns a generator of the automorph group of q in PSL(2,Z).");
		install("bqf_disc_tc","G","bqf_disc","./libapol.so");
		addhelp(bqf_disc, "Inputs: q, quadratic form.\n Returns the discriminant of q.");
		install("bqf_isequiv_tc","GGD0,L,p","bqf_isequiv","./libapol.so");
		addhelp(bqf_isequiv, "Inputs: q, S, {tmat=0}: q a BQF, S either a BQF or a set of BQFs, tmat=0,1.\n This method tests if q is PSL(2,Z) equivalent to S or any form in S. If S is a form, this returns 1 if equivalent and 0 if not (if tmat!=0, returns a possible transition matrix).\n If S is a set of forms, this returns 0 if not equivalent and an index i such that q is equivalent to S[i] otherwise. If tmat!=0, this returns [index, transition matrix].");
		install("bqf_isreduced_tc","iG","bqf_isreduced","./libapol.so");
		addhelp(bqf_isreduced,"Input q, quadratic form.\n Returns 1 if q is reduced, and 0 if not.");
		install("bqf_random","GD0,L,D1,L,","bqf_random","./libapol.so");
		addhelp(bqf_random,"Inputs maxc, {type=0}, {primitive=1}; maxc a positive integer, type=-1,0,1, and primitive=0,1.\n Returns a random BQF with coefficients bounded by maxc. If type=-1 it is positive definite, =1 is indefinite, and =0 means either. If primitive=1 the form is primitive, else it doesn't have to be.");
		install("bqf_random_D","GG","bqf_random_D","./libapol.so");
		addhelp(bqf_random_D,"Inputs maxc, D: maxc a positive integer, and D a discriminant.\n Returns a random primitive form (positive definite if D<0) of discriminant D whose B coefficient is bounded by maxc.");
		install("bqf_red_tc","GD0,L,p","bqf_red","./libapol.so");
		addhelp(bqf_red, "Inputs: q, {tmat=0}: BQF q, (tmat=0,1).\n Returns a reduced form equivalent to q, and if tmat!=0, we return [q_red, transition matrix].");
		install("bqf_roots_tc","G","bqf_roots","./libapol.so");
		addhelp(bqf_roots, "Inputs q: quadratic form q.\n Returns the roots of q with the first root first.");
		install("bqf_trans_tc","GG","bqf_trans","./libapol.so");
		addhelp(bqf_trans, "Inputs q, mtx: integral quadratic form q, 2x2 integral matrix mtx.\n Returns the form acquired by replacing (x,y)^T with m(x,y)^T.");
		install("bqf_trans_coprime_tc", "GG", "bqf_trans_coprime", "./libapol.so");
		addhelp(bqf_trans_coprime,"Inputs q, n: q a primitive integral BQF, and n an integer.\n Returns a form similar to q whose first coefficient is coprime to n.");
		install("ideal_tobqf","GG","ideal_tobqf","./libapol.so");
		addhelp(ideal_tobqf,"Inputs nf, ideal: a quadratic number field nf with ideal ideal.\n Returns the corresponding binary quadratic form.");

	\\BASIC OPERATIONS SPECIFIC TO INDEFINITE FORMS
		install("ibqf_isrecip_tc","iGp","ibqf_isrecip","./libapol.so");
		addhelp(ibqf_isrecip,"Inputs: q, a PIBQF.\n Returns 1 if q is q is reciprocal, and 0 otherwise.");
		install("ibqf_leftnbr_tc","GD0,L,p","ibqf_leftnbr","./libapol.so");
		addhelp(ibqf_leftnbr, "Inputs q, {tmat=0}: an indefinite binary quadratic form on the river and tmat=0,1. \n Returns q' or [q',mat] (if tmat=0,1 respectively), where q' is the left neighbour of q and mat is the bqf_transition matrix from q to q'. The left neighbour is the previous form along the flow of the river that is reduced (AC<0 and B>|A+C|, and occurs when the branches swap from being below to above or vice versa).");
		install("ibqf_redorbit_tc","GD0,L,D0,L,p","ibqf_redorbit","./libapol.so");
		addhelp(ibqf_redorbit, "Inputs q, (tmat), (posonly): q a PIBQF, tmat and posonly=0,1.\n Returns the reduced orbit of q. If tmat=1, also returns the corresponding transition matrices, and if posonly=1 only returns the reduced forms with A>0.");
		install("ibqf_rightnbr_tc","GD0,L,p","ibqf_rightnbr","./libapol.so");
		addhelp(ibqf_rightnbr, "Inputs q, (tmat): an indefinite binary quadratic form on the river and tmat=0,1. \n Returns q' or [q',mat] (if tmat=0,1 respectively), where q' is the right neighbour of q and mat is the bqf_transition matrix from q to q'. The right neighbour is the next form along the flow of the river that is reduced (AC<0 and B>|A+C|, and occurs when the branches swap from being below to above or vice versa).");
		install("ibqf_river_tc","Gp","ibqf_river","./libapol.so");
		addhelp(ibqf_river, "Input: q an indefinite quadratic form.\n Returns the river sequence corresponding to q, where a 1 corresponds to going right and 0 corresponds to going left.");
		install("ibqf_riverforms_tc","Gp","ibqf_riverforms","./libapol.so");
		addhelp(ibqf_riverforms, "Input q a PIBQF.\n This calculates all forms on the river of q, and returns those with A>0, in the order that they appear on the river.");
		install("ibqf_symmetricarc_tc","Gp","ibqf_symmetricarc","./libapol.so");
		addhelp(ibqf_symmetricarc,"Input q, a PIBQF.\n Returns [z,gamma_q(z)] on the root geodesic corresponding to q so that q,gamma_q are symmetric about the arc.");
		install("mat_toibqf_tc","G","mat_toibqf","./libapol.so");
		addhelp(mat_toibqf, "Inputs: mtx, a hyperbolic matrix in SL(2,Z).\n This returns the PIBQF for which it is the ibqf_automorph, namely [c,d-a,-b]/gcd(c,d-a,b) if mtx=[a,b;c,d].");

	\\CLASS GROUPS AND COMPOSITION OF FORMS
		install("bqf_allforms","GD0,L,D0,L,p","bqf_allforms","./libapol.so");
		addhelp(bqf_allforms,"Input D, {prim=0}, {GL=0}: D a proper discriminant.\n Returns the set of equivalence classes of quadratic forms of discriminant D. The returned forms are reduced. If prim=0, we include the non-primitive forms. If GL=1, we take the GL-equivalence classes, else we take SL-equivalence classes.");
		install("bqf_comp_tc","GGD1,L,p","bqf_comp","./libapol.so");
		addhelp(bqf_comp,"Inputs q1, q2, {tored=1}: BQFs q1, q2 of the same discriminant, tored=0, 1.\n Returns the composition of q1 and q2, reduced if tored=1.");
		install("bqf_identify_tc","GGp","bqf_identify","./libapol.so");
		addhelp(bqf_identify,"Inputs ncgp_lexic, q: the output of bqf_ncgp_lexic(D), and a form q of discriminant D.\n Returns [e1,...,er], where in the narrow class group we have q=g1^e1*...*gr^er.");
		install("bqf_lexicind_tobasis","GL","bqf_lexicind_tobasis","./libapol.so");
		addhelp(bqf_lexicind_tobasis,"Inputs orders (vecsmall), index ind.\n Returns [e1,...,er], where the ind element in the lexicographic ordering of the narrow class group is of the form g1^e1*...*gr^er.");
		install("bqf_ncgp","Gp","bqf_ncgp","./libapol.so");
		addhelp(bqf_ncgp, "Input D, a proper discriminant.\n Returns the narrow class group, in the format [n,[n_1,...,n_r],[g_1,...,g_r]], where it has size n, is isomorphic to c_{n_1} x ... x c_{n_r} with n_1 | n_2 | ... | n_r, and g_i is a generator of the corresponding cyclic group of order n_i.");
		install("bqf_ncgp_lexic","Gp","bqf_ncgp_lexic","./libapol.so");
		addhelp(bqf_ncgp_lexic, "Input D, a proper discriminant D.\n This returns [n,[n_1,...,n_r],[f1,f2,...,fl]], where n is the narrow class number of D, the narrow class group is c_{n_1} x ... x c_{n_r} with n_1 | n_2 | ... | n_r, and representative BQFs are the f1,f2,... written in lexicographic order: starting with the identity element, and the component with the highest order moves first.");
		install("bqf_pow_tc","GGD1,L,p","bqf_pow","./libapol.so");
		addhelp(bqf_pow,"Inputs q, n, {tored=1}: BQF q, integer n, tored=0, 1.\n Returns q^n, reduced if tored=1.");
		install("bqf_square_tc","GD1,L,p","bqf_square","./libapol.so");
		addhelp(bqf_square,"Inputs q, {tored=1}: BQF q, tored=0, 1.\n Returns q^2, reduced if tored=1.");

	\\REPRESENTATIONS OF NUMBERS BY BQFs
		install("bqf_reps_tc","GGD0,L,D1,L,p","bqf_reps","./libapol.so");
		addhelp(bqf_reps,"Inputs q, n, {proper=0}, {half=1}: BQF q, integer n, (proper=0,1), (half=0,1).\n This solves the equation q(x,y)=n over the integers. We will get a finite set of (families) of solutions, and if half=0 we only return one of the families corresponding to (x,y) and (-x,-y). If we want only coprime solutions when disc!=square, pass in proper=1. If Q has discriminant D, the return is:\n\n --------If no solutions, returns 0;\n -D>0 nonsquare and n!=0, [[1,M],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are representatives of the distinct classes of solutions and M is the invariant automorph;\n ----D>0 square and n!=0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are the solutions;\n --------------------D<0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are the solutions;\n -----------D=0 and n!=0, [[2],[[s1,s2],[x1,y1]]] where the solutions are (up to +/-) x=x1+Us1, y=y1+Us2;\n ---------n=0, D!=square, [[0],[0,0]];\n -------n=0, D=square!=0, [[2],[[s1,s2],[0,0]],[[s3,s4],[0,0]]], solutions are (x,y)=(s1k,s2k),(s3k,s4k) for integer k;\n ------------n=0 and D=0, if Q!=0 as above, if Q=0 then [[-1]] (everything is a solution).\n\n In general, -1=all, 0=finite, 1=positive, 2=linear");

	\\MORE REPRESENTATION OF NUMBERS
		install("bqf_bigreps_tc","GGp","bqf_bigreps","./libapol.so");
		addhelp(bqf_bigreps,"Inputs: Q, n, Q=[A,B,C,D,E]=Ax^2+Bxy+Cy^2+Dx+Ey and an integer n.\n Returns the solutions to Q(x,y)=n over the integers. If D=bqf_disc(Q)=B^2-4AC, then the output is: \n\n ----------If no solutions, returns 0; \n ------------D>0 nonsquare, [[1,M,[s1,s2]],[x1,y1],[x2,y2],...,[xk,yk]] where the general solution is [x;y]=M^k*[xi;yi]+[s1;s2];\n --------D>0 square EITHER, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the (xi,yi) are the solutions;\n -----------------------OR, [[2],[[s1,t1],[x1,y1]],...,[[sk,tk],[xk,yk]]] where the solutions are x=xi+Usi and y=yi+Uti for U integer;\n ----------------------D<0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the (xi,yi) are the solutions;\n ----------------------Q=0, [[-1]] if n=0 and 0 else;\n -D=0 and A=B=C=0 or D=E=0, [[2],[[s1,t1],[x1,y1]],...,[[sk,tk],[xk,yk]]] where the solutions are x=xi+Usi and y=yi+Uti for U integer (si=sj and ti=tj for all i,j in fact);\n ------------D=0 otherwise, [[-2],[[a1,b1,c1],[e1,f1,g1]],...,[[ak,bk,ck],[ek,fk,gk]]] where the solutions are x=ai*U^2+bi*U+ci, y=ei*U^2+fi*U+gi for U integer.\n\n In general, -2=quadratic, -1=all, 0=finite, 1=positive, 2=linear");
		install("bqf_linearsolve_tc","GGGGp","bqf_linearsolve","./libapol.so");
		addhelp(bqf_linearsolve,"Inputs qf, n1, lin, n2: qf a six term integer vector representing the form Ax^2+By^2+Cz^2+Dxy+Exz+Fyz, n1 an integer, lin a three term integer vector representing Ax+By+Cz, and n2 an integer.\n Solves qf(x, y, z)=n1 and lin(x, y, z)=n2 simultaneously. If there are no solutions, this returns 0. Otherwise it returns a vector v. Let v[1][1]=t, and then the format of v is:\n\n -t=-2, v=[[-2], [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]], ...], where each general solution is x=a1U^2+a2U+a3, y=b1U^2+b2U+b3, z=c1U^2+c2U+c3 for any U integral;\n -t=-1, v=[[-1], [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]]], where the solution is x=a1U+b1V+c1, y=a2U+b2V+c2, z=a3U+v3V+c3 for any U, V integral;\n --t=0, v=[[0], [a1, b1, c1], ...], where the finite set of solutions are (x,y,z)=(ai, bi, ci);\n --t=1, v=[[1, M, [s1, s2, s3]], [a1, b1, c1], ...], where the general solution is [x;y;z]=M^k[ai;bi;ci]+[s1;s2;s3] for k integral. Note that ai and si need not be integral, though M is.\n --t=2, v=[[2], [[a1, a2, a3], [b1, b2, b3]], ...], where each general solution is x=a1U+b1, y=a2U+b2, z=a3U+b3 for any U integral;\n\n In general, -2=quadratic, -1=plane, 0=finite, 1=positive, 2=linear.");

	\\GENERAL HELP
		addhelp(bqf, "This package deals with binary quadratic forms with integer coefficients. A homogeneous binary quadratic form Ax^2+Bxy+Cy^2 is stored as [A,B,C]. A proper discriminant is an integer that is equivalent to 0 or 1 modulo 4 and is not a square. \n Subtopics:\n Discriminants (disc)\n Basic operations (bqfbasic)\n Indefinite forms (ibqf)\n Class group and composition (bqfclass)\n Representation of numbers (bqfsolve)");
		addhelp(disc,"disclist, discprimeindex, discsuperorders, isdisc, pell, posreg, quadroot.");
		addhelp(bqfbasic,"bqf_automorph, bqf_disc, bqf_isequiv, bqf_isreduced, bqf_random, bqf_random_D, bqf_red, bqf_roots, bqf_trans, bqf_trans_coprime, ideal_tobqf.");
		addhelp(ibqf,"ibqf_isrecip, ibqf_leftnbr, ibqf_redorbit, ibqf_rightnbr, ibqf_river, ibqf_riverforms, ibqf_symmetricarc, mat_toibqf.");
		addhelp(bqfclass,"bqf_allforms, bqf_comp, bqf_identify, bqf_lexicind_tobasis, bqf_ncgp, bqf_ncgp_lexic, bqf_pow, bqf_square.");
		addhelp(bqfsolve,"bqf_bigreps, bqf_linearsolve, bqf_reps.");

\\farey.c
		install("fareydenom_depths","L", ,"./libapol.so");
		addhelp(fareydenom_depths,"Input n>1 an integer.\n Returns the set of fareydepth(i/n) for all 1<=i<n with gcd(n, i)=1.");
		install("fareydenom_dmode","lL", ,"./libapol.so");
		addhelp(fareydenom_dmode,"Input n>1 an integer.\n Returns the most common depth in fareydenom_depths(n).");
		install("fareydepth","lG", ,"./libapol.so");
		addhelp(fareydepth,"Input r, a rational number in (0, 1) or a Farey pair.\n Returns the depth of r in the Farey tree.");
		install("fareyup","G", ,"./libapol.so");
		addhelp(fareyup,"Input r, a rational number in (0, 1) or a Farey pair.\n Returns either the Farey pair with r as its mediant, or the (unique) Farey pair with one of the numbers as it's endpoint and the other as its mediant.");

\\visual.c
	\\DATA
		install("integerbin","GGD0,G,",,"./libapol.so");
		addhelp(integerbin,"Inputs v, binlen, {binstart=0}.\n Assumes v is a sorted list of integers, and puts them into bins of length binlen, starting with binstart (assumed to be 0). Returns [binends, counts], with binends being the last number in the bin.");
		install("integerbin_cumu","GGD0,G,",,"./libapol.so");
		addhelp(integerbin_cumu,"Inputs v, binlen, {binstart=0}.\n Assumes v is a sorted list of integers, and puts them into bins of length binlen, starting with binstart (assumed to be 0). Returns [binends, counts], with binends being the last number in the bin. This is cumulative, so counts is increasing.");
		install("veccount","G",,"./libapol.so");
		addhelp(veccount,"Input v, a vector.\n Returns [uniq, count], where uniq is the sorted vector v with repeats removed, and count is the corresponding number of times they appear in v.");
		install("vecsmallcount","G",,"./libapol.so");
		addhelp(vecsmallcount,"Input v, a vecsmall.\n Returns [uniq, count], where uniq is the sorted vecsmall v with repeats removed, and count is the corresponding number of times they appear in v.");
		install("ZV_countnonpos","lG",,"./libapol.so");
		addhelp(ZV_countnonpos,"Input v, a sorted vector of integers.\n Returns the number of entries that are nonpositive.");
		
	\\HISTOGRAMS
		install("hist_make","GrD0,L,D0,L,Drp",,"./libapol.so");
		addhelp(hist_make,"Inputs data, imagename, {compilenew=0}, {open=0}, {plotoptions=NULL}: sorted list of real numbers data, name of the tikz picture, {compilenew=0, 1}, {open=0, 1}, {plotoptions=NULL, string}.\n Automatically bins the data, and creates a pdf of the histogram using tikz and externalize. The output is in the folder /images, with the build file (named imagename_build.tex) being in the folder /images/build. If compilenew=0 assumes the LaTeX document to compile the plot is pre-made, and otherwise this method automatically writes it. If additionally, plotoptions!=NULL, this string is inserted in between \\begin{axis} and \\end{axis} in the LaTeX document (allowing one to customize how the histogram looks). If open=1, the pdf is automatically opened (only works with Linux subsystem for Windows). The returned value is used to modify the histogram, e.g. changing the bins, scaling it, and changing the range.");
		install("hist_rebin","GGGp",,"./libapol.so");
		addhelp(hist_rebin,"Inputs data, histdata, nbins: the sorted data, the length 7 vector output of a hist_ method, the number of bins.\n Rebins the data according to the new number of bins, and updates histdata.");
		install("hist_rerange","GGGGp",,"./libapol.so");
		addhelp(hist_rerange,"Inputs data, histdata, minx, maxx: the sorted data, the length 7 vector output of a hist_ method, minimum x-value, maximum x-value.\n Rebins the data according to the new minimum and maximum value. This is useful when there are outliers that skew the look of the graph. Returns the updated histdata.");
		install("hist_rescale","GGLp",,"./libapol.so");
		addhelp(hist_rescale,"Inputs data, histdata, scale: the sorted data, the length 7 vector output of a hist_ method, and scale=0, 1.\n If scale=1 scales the data so the total area is 1, and if scale=0 uses the absolute count for the y-axis. Returns the updated histdata.");

	\\REGRESSIONS & PLOTS
		install("OLS","GGD1,L,",,"./libapol.so");
		addhelp(OLS,"Inputs X, y, {retrsqr=1}:  m*n matrix X with top row being all 1's, length n column vector y, retrsqr=0, 1.\n Performs ordinary least squares regression on the data, where the n inputs are the columns of X, and the outputs are the entries of y. We must include a constant term, hence why the first row of X must be all 1's. If retrsqr=1, returns [pararms, R^2], and otherwise returns params, where params is the length m column vector of best fit parameters.");
		install("OLS_nointercept","GGD1,L,",,"./libapol.so");
		addhelp(OLS_nointercept,"Inputs X, y, {retrsqr=1}: vector X, column vector y (of same length), retrsqr=0, 1.\n Performs ordinary least squares regression on the data assuming that y[i]=c*X[i], i.e. the y-intercept is 0. Returns c if retrsqr=0, or [c, R^2] otherwise.");
		install("OLS_single","GGD1,L,",,"./libapol.so");
		addhelp(OLS_single,"Inputs x, y, {retrsqr=1}: vector x, column vector y, retrsqr=0, 1. Performs linear regression for a single variable (essentially a macro for OLS with y=mx+b.");
		install("rsquared","GGG",,"./libapol.so");
		addhelp(rsquared,"Inputs X, y, fit: X and y data supplied to OLS, and fit the proposed fit (a column vector of parameters). This returns the R^2 value for this proposal.");

	\\TEX
		install("tex_recompile","vG",,"./libapol.so");
		addhelp(tex_recompile,"Input data, the output of a tex image creation call.\n Recompiles the image, returning nothing. This is used when you edit the LaTeX document by hand.");

	\\GENERAL HELP
		addhelp(visual,"This package deals with visualizing data. Subtopics:\n Data (data)\n Histograms (hist)\n Regressions/plots (reg)\n Tex (tex)");
		addhelp(data,"integerbin, veccount, vecsmallcount.");
		addhelp(hist,"Installed methods:\n hist_make, hist_rebin, hist_rerange, hist_rescale.");
		addhelp(reg,"OLS, OLS_nointercept, OLS_single, rsquared.");
		addhelp(tex,"tex_recompile");

default(parisize, "4096M");\\Must come at the end