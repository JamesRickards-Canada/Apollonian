/*Methods for data collection and visualization.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"
#include <stdlib.h>

/*STATIC DECLARATIONS*/

/*SECTION 1: DATA*/
static GEN vecsmall_reduce(GEN v, GEN *pE);

/*SECTION 2: HISTOGRAMS*/
static void hist_autocompile(GEN minx, GEN maxx, char *imagename, char *plotoptions, int open);
static GEN hist_tobins(GEN v, GEN minx, GEN maxx, GEN nbins, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec);
static GEN hist_tobins_defaultbins(GEN data, GEN minx, GEN maxx, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec);

/*MAIN BODY*/

/*SECTION 1: DATA*/

/*If v is a sorted vector of integers, this bins the data so everything is in a bin. Returns [binends, counts]. Makes the bin length an integer.*/
GEN
integerbin(GEN v, GEN blen, GEN bstart)
{
  pari_sp av = avma;
  long lv = lg(v);
  GEN last = gel(v, lv - 1);
  long nbins = itos(divii(subis(subii(last, bstart), 1), blen)) + 1, i;
  GEN bends = cgetg(nbins + 1, t_VEC);
  gel(bends, 1) = addii(bstart, blen);
  for (i = 2; i <= nbins; i++) gel(bends, i) = addii(gel(bends, i - 1), blen);
  GEN counts = const_vecsmall(nbins, 0);
  long bind = 1;
  for (i = 1; i < lv; i++) {
	if (cmpii(gel(v, i), gel(bends, bind)) <= 0) counts[bind]++;
	else {
	  i--;//Redo this index
	  bind++;//Move to a new bin.
	}
  }
  return gerepilecopy(av, mkvec2(bends, counts));
}

/*integerbin, but cumulative.*/
GEN
integerbin_cumu(GEN v, GEN blen, GEN bstart)
{
  pari_sp av = avma;
  long lv = lg(v);
  GEN last = gel(v, lv - 1);
  long nbins = itos(divii(subis(subii(last, bstart), 1), blen)) + 1, i;
  GEN bends = cgetg(nbins + 1, t_VEC);
  gel(bends, 1) = addii(bstart, blen);
  for (i = 2; i <= nbins; i++) gel(bends, i) = addii(gel(bends, i - 1), blen);
  GEN counts = const_vecsmall(nbins, 0);
  long bind = 1;
  for (i = 1; i < lv; i++) {
	if (cmpii(gel(v, i), gel(bends, bind)) <= 0) counts[bind]++;
	else {
	  i--;/*Redo this index*/
	  bind++;/*Move to a new bin.*/
	  counts[bind] = counts[bind - 1];
	}
  }
  for (i = bind + 1; i <= nbins; i++) counts[i] = counts[i - 1];/*Filling in the rest of the bins to have the final count (shouldn't be an issue but just in case).*/
  return gerepilecopy(av, mkvec2(bends, counts));
}

/*Returns [vsort, count], where vsort is the sorted vector v with duplicates removed, and count is the Vecsmall of corresponding number of each in the original vector v.*/
GEN
vecreduce(GEN v)
{
  pari_sp av = avma;
  GEN E;
  GEN uniq = vec_reduce(v, &E);
  return gerepilecopy(av, mkvec2(uniq, E));
}

/*vecreduce for a Vecsmall*/
GEN
vecsmallreduce(GEN v)
{
  pari_sp av = avma;
  GEN E;
  GEN uniq = vecsmall_reduce(v, &E);
  return gerepilecopy(av, mkvec2(uniq, E));
}

/*Adapted from vec_reduce in bibli2.c; does the same thing, except for a Vecsmall.*/
static GEN
vecsmall_reduce(GEN v, GEN *pE)
{
  GEN E, F, P = vecsmall_indexsort(v);
  long i, m, l;
  F = cgetg_copy(v, &l);
  *pE = E = cgetg_copy(v, &l);
  for (i = m = 1; i < l;)
  {
    long u = v[P[i]];
    long k;
    for (k = i + 1; k < l; k++)
      if (v[P[k]] != u) break;
    E[m] = k - i; F[m] = u; i = k; m++;
  }
  setlg(F, m);
  setlg(E, m); return F;
}


/*SECTION 2: HISTOGRAMS*/

/*Creates LaTeX document and compiles the histogram. If plotoptions!=NULL, adds this between begin and end tikzpicture*/
static void
hist_autocompile(GEN minx, GEN maxx, char *imagename, char *plotoptions, int open)
{
  pari_sp av = avma;
  if (!pari_is_dir("images/build")) {
    int s = system("mkdir -p images/build");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY images/build");
  }
  char *autofilestr = stack_sprintf("images/build/%s_build.tex", imagename);
  FILE *f = fopen(autofilestr, "w");
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\tikzset{external/force remake}\n  \\pgfplotsset{compat=1.15}\n\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n  \\begin{axis}", imagename);
  if (!plotoptions) {
    pari_fprintf(f, "[ylabel=Frequency, xmin=%lf, xmax=%lf, ymin=0]\n  \\addplot[ycomb, blue]table[x=x, y=y, col sep=space]{%s.dat};\n", rtodbl(minx), rtodbl(maxx), imagename);
  }
  else pari_fprintf(f, "%s\n", plotoptions);
  pari_fprintf(f, "  \\end{axis}\n\\end{tikzpicture}\n\\end{document}");
  fclose(f);
  tex_compile(imagename, open);
  set_avma(av);
}

/*Makes a histogram with the inputs, defaulting to the defaults. Adjust it to your liking with the hist_re-methods. Should only be called once with the given data.*/
GEN
hist_make(GEN v, char *imagename, int compilenew, int open, char *plotoptions, long prec)
{
  return hist_tobins_defaultbins(v, gel(v, 1), gel(v, lg(v) - 1), 0, compilenew, imagename, plotoptions, open, prec);
}

/*Bins the data according to the parameters. Returns the associated parameters*/
static GEN
hist_tobins(GEN v, GEN minx, GEN maxx, GEN nbins, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec)
{
  pari_sp av = avma;
  if (!pari_is_dir("images/build")) {/*Start with the file things*/
    int s = system("mkdir -p images/build");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY images/build");
  }
  char *filetitle = stack_sprintf("images/build/%s.dat", imagename);
  FILE *f = fopen(filetitle, "w");
  minx = gtofp(minx, prec);/*Convert to t_REAL to use rtodbl*/
  maxx = gtofp(maxx, prec);
  pari_fprintf(f, "x y\n%lf ", rtodbl(minx));/*First line and start of second*/
  /*Initialize the bins*/
  long nbinslong = itos(nbins);/*long version of nbins*/
  GEN range = gsub(maxx, minx);/*Range of v*/
  GEN binlen = gdiv(range, nbins);/*bin length*/
  GEN bmax = gadd(minx, binlen);/*Upper limit of the current bin*/
  /*Outputting to the bins*/
  long ind = 1, ninbin, totinbin = 0, binno, lv = lg(v);
  if (!toscale) {
    for (binno = 0; binno < nbinslong; binno++) {
      ninbin = 0;
      while (ind < lv && gcmp(gel(v, ind), bmax) <= 0) { ind++; ninbin++; }
      pari_fprintf(f, "%d\n%lf ", ninbin, rtodbl(bmax));
      bmax = gadd(bmax, binlen);
      totinbin = totinbin + ninbin;
    }
    totinbin=lv - 1 - totinbin;
    pari_fprintf(f, "%d", totinbin);/*The last outliers*/
  }
  else {
    GEN rescalefactor = gdivsg(1, gmulgs(binlen, lv - 1));/*Rescaling factor for area 1*/
    for (binno = 0; binno < nbinslong; binno++) {
      ninbin = 0;
      while (ind < lv && gcmp(gel(v, ind), bmax) <= 0) { ind++; ninbin++; }
      pari_fprintf(f, "%lf\n%lf ", rtodbl(gmulsg(ninbin, rescalefactor)), rtodbl(bmax));
      bmax = gadd(bmax, binlen);
      totinbin = totinbin + ninbin;
    }
    totinbin=lv - 1 - totinbin;
    pari_fprintf(f, "%lf", rtodbl(gmulsg(totinbin, rescalefactor)));/*The last outliers*/
  }
  fclose(f);
  if (compilenew) hist_autocompile(minx, maxx, imagename, plotoptions, open);/*Making a new LaTeX document and compiling it*/
  else tex_compile(imagename, open);/*Just recompiling the old LaTeX file*/
  pari_printf("%d total data points\n%d total bins", lv - 1, nbinslong);
  GEN histprop = cgetg(8, t_VEC);
  gel(histprop, 1) = strtoGENstr(imagename);
  if (!open) gel(histprop, 2) = gen_0;
  else gel(histprop, 2) = gen_1;
  gel(histprop, 3) = minx; gel(histprop, 4) = maxx;
  gel(histprop, 5) = nbins; gel(histprop, 6) = stoi(toscale);
  if (!plotoptions) gel(histprop, 7) = gen_0;
  else gel(histprop, 7) = strtoGENstr(plotoptions);
  return gerepilecopy(av, histprop);
}

/*Bins the data with default options, outputting to imagename.*/
static GEN
hist_tobins_defaultbins(GEN v, GEN minx, GEN maxx, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec)
{
  pari_sp av = avma;
  long ind1 = (lg(v) - 1) >> 2, ind2 = 3 * ind1;/*The difference approximately gives the IQR of v*/
  GEN binlen = gdiv(gmulsg(2, gsub(gel(v, ind2), gel(v, ind1))), gpow(stoi(lg(v) - 1), mkfracss(1, 3), prec));/*Using the Freedman-Diaconis rule of bin width=2*IQR/(n^(1/3))*/
  GEN nbins = mpfloor(gdiv(gsub(maxx, minx), binlen));
  return gerepileupto(av, hist_tobins(v, minx, maxx, nbins, toscale, compilenew, imagename, plotoptions, open, prec));
}

/*Remake the histogram with the new bin count*/
GEN
hist_rebin(GEN v, GEN histdata, GEN nbins, long prec)
{
    if (gequal0(gel(histdata, 7))) return hist_tobins(v, gel(histdata, 3), gel(histdata, 4), nbins, itos(gel(histdata, 6)), 0, GENtostr_unquoted(gel(histdata, 1)), NULL, itos(gel(histdata, 2)), prec);
    return hist_tobins(v, gel(histdata, 3), gel(histdata, 4), nbins, itos(gel(histdata, 6)), 0, GENtostr_unquoted(gel(histdata, 1)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 2)), prec);
}

/*Remake the histogram with new minx and maxx*/
GEN
hist_rerange(GEN v, GEN histdata, GEN minx, GEN maxx, long prec)
{
    if (gequal0(gel(histdata, 7))) return hist_tobins(v, minx, maxx, gel(histdata, 5), itos(gel(histdata, 6)), 1, GENtostr_unquoted(gel(histdata, 1)), NULL, itos(gel(histdata, 2)), prec);
    return hist_tobins(v, minx, maxx, gel(histdata, 5), itos(gel(histdata, 6)), 1, GENtostr_unquoted(gel(histdata, 1)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 2)), prec);
}

/*Remake the histogram with the new scale*/
GEN
hist_rescale(GEN v, GEN histdata, int scale, long prec)
{
    if (gequal0(gel(histdata, 7))) return hist_tobins(v, gel(histdata, 3), gel(histdata, 4), gel(histdata, 5), scale, 1, GENtostr_unquoted(gel(histdata, 1)), NULL, itos(gel(histdata, 2)), prec);
    return hist_tobins(v, gel(histdata, 3), gel(histdata, 4), gel(histdata, 5), scale, 1, GENtostr_unquoted(gel(histdata, 1)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 2)), prec);
}



//REGRESSIONS


/*Perform ordinary least squares regression. X is a matrix whose columns are the parameters, and y is a column vector of results. Must include linear term as first variable of X.
The formula is B=Bhat=(X*X^T)^(-1)Xy, for the ordinary least squares regression for y=X^T*B+error (formula differs to Wikipedia due to X here being the transpose of what they define there.
Returns [best fit, R^2]*/
GEN OLS(GEN X, GEN y, int retrsqr){
  pari_sp top=avma;
  if(typ(y)!=t_COL) y=gtocol(y);
  if(lg(y)!=lg(X)) pari_err_TYPE("The inputs must have the same length.", mkvec2(X, y));
  GEN Xy=RgM_RgC_mul(X, y);
  GEN XTX=RgM_multosym(X, shallowtrans(X));//X*X^T, which is symmetric
  GEN fit=RgM_solve(XTX, Xy);//Best fit.
  if(!fit) pari_err(e_MISC, "Could not compute matrix inverse. Error with given data or with precision?");
  if(!retrsqr) return gerepileupto(top, fit);
  GEN rsqr=rsquared(X, y, fit);
  return gerepilecopy(top, mkvec2(fit, rsqr));
}

//Performs OLS where we have one independant variable and assume the intercept is 0 (so y=ax). The formua is now sum(x*y)/sum(x^2).
GEN OLS_nointercept(GEN X, GEN y, int retrsqr){
  pari_sp top=avma;
  GEN xysum=gen_0, xsqrsum=gen_0;
  if(lg(X)!=lg(y)) pari_err_TYPE("The inputs must have the same length.", mkvec2(X, y));
  for(long i=1;i<lg(X);i++){
    xysum=gadd(xysum, gmul(gel(X, i), gel(y, i)));
    xsqrsum=gadd(xsqrsum, gsqr(gel(X, i)));
  }
  GEN fit=gdiv(xysum, xsqrsum);
  if(!retrsqr) return gerepileupto(top, fit);
  long lX=lg(X);
  GEN M=cgetg(lX, t_MAT);
  for(long i=1;i<lX;i++) gel(M, i)=mkcol2(gen_1, gel(X, i));
  GEN rsqr=rsquared(M, y, mkcol2(gen_0, fit));
  return gerepilecopy(top, mkvec2(fit, rsqr));
}

//OLS, where there is only one input variable. This just puts it into a matrix form and calls OLS, and is included for convenience.
GEN OLS_single(GEN x, GEN y, int retrsqr){
  pari_sp top=avma;
  long lgx=lg(x);
  GEN xmat=cgetg(lgx, t_MAT);
  for(long i=1;i<lgx;i++) gel(xmat, i)=mkcol2(gen_1, gel(x, i));
  return gerepileupto(top, OLS(xmat, y, retrsqr));
}

//Given inputs for OLS and the proposed linear fit, this returns the R^2 value of the regression.
GEN rsquared(GEN X, GEN y, GEN fit){
  pari_sp top=avma;
  long n=lg(y)-1;//Number of observations
  GEN predicted=RgV_RgM_mul(shallowtrans(fit), X);//1xn matrix of the fitted values.
  GEN yavg=gen_0;
  for(long i=1;i<=n;i++) yavg=gadd(yavg, gel(y, i));
  yavg=gdivgs(yavg, n);//Average value of y
  GEN sstot=gen_0;
  GEN ssres=gen_0;
  for(long i=1;i<=n;i++){
    sstot=gadd(sstot, gsqr(gsub(gel(y, i), yavg)));
    ssres=gadd(ssres, gsqr(gsub(gel(y, i), gel(predicted, i))));
  }
  return gerepileupto(top, gsubsg(1, gdiv(ssres, sstot)));
}



//TEX


//Returns a ncol Latex colours. If ncol<=63, each entry of the return vector is a GEN string (HTML values), else each entry of the return vector is [a, b, c] giving the rgb values.
GEN tex_makecolours(int ncol){
  pari_sp top=avma;
  GEN rvec=cgetg(ncol+1, t_VEC);
  if(ncol<=63){
	const char *cs[] = {"00FF00", "0000FF", "FF0000", "01FFFE", "FFA6FE", "FFDB66", "006401", "010067", "95003A", "007DB5", "FF00F6", "FFEEE8", "774D00", "90FB92", "0076FF", "D5FF00", "FF937E", "6A826C", "FF029D", "FE8900", "7A4782", "7E2DD2", "85A900", "FF0056", "A42400", "00AE7E", "683D3B", "BDC6FF", "263400", "BDD393", "00B917", "9E008E", "001544", "C28C9F", "FF74A3", "01D0FF", "004754", "E56FFE", "788231", "0E4CA1", "91D0CB", "BE9970", "968AE8", "BB8800", "43002C", "DEFF74", "00FFC6", "FFE502", "620E00", "008F9C", "98FF52", "7544B1", "B500FF", "00FF78", "FF6E41", "005F39", "6B6882", "5FAD4E", "A75740", "A5FFD2", "FFB167", "009BFF", "E85EBE"};
	int first=itos(randomi(stoi(63)));
	gel(rvec, 1)=strtoGENstr(cs[first]);
	for(long i=2;i<=ncol;i++){
      first=(first+1)%63;
	  gel(rvec,i)=strtoGENstr(cs[first]);
	}
  }
  else{
    long prec=3;
    GEN shift=rdivss(1, ncol, prec);//1/ncol
    GEN a=randomr(prec), b=randomr(prec), c=randomr(prec);
    gel(rvec, 1)=mkvec3(a, b, c);
    for(long i=2;i<=ncol;i++){
	  gel(rvec, i)=mkvec3(gfrac(gadd(gmael(rvec, i-1, 1), shift)), gfrac(gadd(gmael(rvec, i-1, 2), shift)), gfrac(gadd(gmael(rvec, i-1, 3), shift)));
    }
  }
  return gerepilecopy(top, rvec);
}

//In general, an image created from a file will have a GP return value. The first entry is the image name, and the second is whether we want to open it every time or not (i.e. if we are on WSL or not). The rest of the entries will depend on the type of image.

//Compiles the latex file specified by imagename
void tex_compile(char *imagename, int open){
  char *line=pari_sprintf("(cd ./images/build && pdflatex --interaction=batchmode -shell-escape %s_build.tex)", imagename);
  int s=system(line);
  if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  pari_free(line);
  line=pari_sprintf("mv -f ./images/build/%s.pdf ./images/", imagename);//Move the file
  s=system(line);
  if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  pari_free(line);
  if(open){
    line=pari_sprintf("cmd.exe /C start images/%s.pdf", imagename);//Open the file
    s=system(line);
    if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
    pari_free(line);
  }
}

//Recompile the image. Used when we modify the TEX document explicitly, and then want to recompile it.
void tex_recompile(GEN data){
  tex_compile(GENtostr_unquoted(gel(data, 1)), itos(gel(data, 2)));
}

