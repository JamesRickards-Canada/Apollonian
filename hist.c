//IMPORTED FROM Q-QUADRATIC
//Methods for visualization and data compaliation

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "apol.h"
#endif

#ifndef STDLIB
#define STDLIB
#include <stdlib.h>
#endif


//Should also write the data to a file and offer ways to make the histogram from a file of data.



//HISTOGRAMS



//Creates LaTeX document and compiles the histogram. If plotoptions!=NULL, adds this between begin and end tikzpicture
void hist_autocompile(GEN minx, GEN maxx, char *imagename, char *autofile, char *plotoptions, int open){
  pari_sp top=avma;
  if(!pari_is_dir("images/build")){
	int s=system("mkdir -p images/build");
	if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY images/build");
  }
  char *autofilestr=pari_sprintf("images/build/%s.tex", autofile);
  FILE *f=fopen(autofilestr, "w");
  pari_free(autofilestr);//Now we have created the output file f.
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\tikzset{external/force remake}\n  \\pgfplotsset{compat=1.15}\n\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n  \\begin{axis}", imagename);
  if(plotoptions==NULL){
    pari_fprintf(f, "[ylabel=Frequency, xmin=%lf, xmax=%lf, ymin=0]\n  \\addplot[ycomb, blue]table[x=x, y=y, col sep=space]{%s.dat};\n", rtodbl(minx), rtodbl(maxx), imagename);
  }
  else pari_fprintf(f, "%s\n", plotoptions);//Printing plotoptions instead
  pari_fprintf(f, "  \\end{axis}\n\\end{tikzpicture}\n\\end{document}");
  fclose(f);
  hist_compile(imagename, autofile, open);
  avma=top;
}

//Compiles the latex file specified by filenames
void hist_compile(char *imagename, char *autoname, int open){
  char *line=pari_sprintf("(cd ./images/build && pdflatex --interaction=batchmode -shell-escape %s.tex)", autoname);
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

//Makes a histogram with the inputs, defaulting to the defaults. Adjust it to your liking with the hist_re-methods. Should only be called once with the given data.
GEN hist_make(GEN data, char *imagename, char *autofile, int compilenew, char *plotoptions, int open, long prec){
  return hist_tobins_defaultbins(data, gel(data, 1), gel(data, lg(data)-1), 0, compilenew, imagename, autofile, plotoptions, open, prec);
}

//Bins the data according to the parameters. Returns the associated parameters
GEN hist_tobins(GEN data, GEN minx, GEN maxx, GEN nbins, int toscale, int compilenew, char *imagename, char *autofile, char *plotoptions, int open, long prec){
  pari_sp top=avma;
  //Start with the file things
  if(!pari_is_dir("images/build")){
    int s=system("mkdir -p images/build");
	if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY images/build");
  }
  char *filetitle=pari_sprintf("images/build/%s.dat",imagename);
  FILE *f=fopen(filetitle, "w");
  pari_free(filetitle);//Now we have created the output file f.
  pari_fprintf(f, "x y\n%lf ", rtodbl(minx));//First line and start of second
  
  //On to initializing bin stuff
  minx=gtofp(minx, prec);//Required to make it work, otherwise this may print fractions which are sadly unreadable to latex.
  maxx=gtofp(maxx, prec);
  long nbinslong=itos(nbins);//long version of nbins
  GEN range=gsub(maxx,minx);//Range of data
  GEN binlen=gdiv(range, nbins);//bin length
  GEN bmax=gadd(minx, binlen);//Upper limit of the current bin
  
  //Outputting to the bins
  long ind=1, ninbin, totinbin=0;
  if(toscale==0){
    for(long binno=0;binno<nbinslong;binno++){
	  ninbin=0;
	  while(ind<lg(data) && gcmp(gel(data, ind), bmax)<=0){ind++;ninbin++;}
	  pari_fprintf(f, "%d\n%lf ", ninbin, rtodbl(bmax));
	  bmax=gadd(bmax, binlen);
	  totinbin=totinbin+ninbin;
    }
    totinbin=lg(data)-1-totinbin;
    pari_fprintf(f, "%d", totinbin);//The last outliers
  }
  else{
	GEN rescalefactor=gdivsg(1,gmulgs(binlen,(lg(data)-1)));//Rescaling factor for area 1
    for(long binno=0;binno<nbinslong;binno++){
	  ninbin=0;
	  while(ind<lg(data) && gcmp(gel(data, ind), bmax)<=0){ind++;ninbin++;}
	  pari_fprintf(f, "%lf\n%lf ", rtodbl(gmulsg(ninbin, rescalefactor)), rtodbl(bmax));
	  bmax=gadd(bmax, binlen);
	  totinbin=totinbin+ninbin;
    }
    totinbin=lg(data)-1-totinbin;
    pari_fprintf(f, "%lf", rtodbl(gmulsg(totinbin, rescalefactor)));//The last outliers
  }
  fclose(f);
  if(compilenew) hist_autocompile(minx, maxx, imagename, autofile, plotoptions, open);//Making new LaTeX document and compiling it
  else hist_compile(imagename, autofile, open);//Just recompiling the old LaTeX file
  pari_printf("%d total data points\n%d total bins", lg(data)-1, nbinslong);
  GEN histprop=cgetg(9, t_VEC);
  gel(histprop, 1)=gcopy(minx);gel(histprop, 2)=gcopy(maxx);
  gel(histprop, 3)=icopy(nbins);gel(histprop, 4)=stoi(toscale);
  gel(histprop, 5)=strtoGENstr(imagename);gel(histprop, 6)=strtoGENstr(autofile);
  if(plotoptions==NULL) gel(histprop, 7)=gen_0;
  else gel(histprop, 7)=strtoGENstr(plotoptions);
  if(open==0) gel(histprop, 8)=gen_0;
  else gel(histprop, 8)=gen_1;
  return gerepileupto(top, histprop);
}

//Bins the data with default options, outputting to imagename.
GEN hist_tobins_defaultbins(GEN data, GEN minx, GEN maxx, int toscale, int compilenew, char *imagename, char *autofile, char *plotoptions, int open, long prec){
  pari_sp top=avma;
  long ind1=(lg(data)-1)/4, ind2=3*ind1;//The difference approximately gives the IQR of the data
  GEN binlen=gdiv(gmulsg(2,gsub(gel(data, ind2), gel(data, ind1))),gpow(stoi(lg(data)-1), gdivgs(gen_1, 3), prec));//Using the Freedman-Diaconis rule of bin width=2*IQR/(n^(1/3))
  GEN nbins=mpfloor(gdiv(gsub(maxx, minx), binlen));
  return gerepileupto(top, hist_tobins(data, minx, maxx, nbins, toscale, compilenew, imagename, autofile, plotoptions, open, prec));
}

//Remake the histogram with the new bin count
GEN hist_rebin(GEN data, GEN histdata, GEN nbins, long prec){
	if(gequal0(gel(histdata, 7))) return hist_tobins(data, gel(histdata, 1), gel(histdata, 2), nbins, itos(gel(histdata, 4)), 0, GENtostr_unquoted(gel(histdata, 5)), GENtostr_unquoted(gel(histdata, 6)), NULL, itos(gel(histdata, 8)), prec);
	return hist_tobins(data, gel(histdata, 1), gel(histdata, 2), nbins, itos(gel(histdata, 4)), 0, GENtostr_unquoted(gel(histdata, 5)), GENtostr_unquoted(gel(histdata, 6)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 8)), prec);
}

//Recompile the histogram with the histdata. Used when we modify the TEX document explicitly, and then want to recompile it.
void hist_recompile(GEN histdata){
  hist_compile(GENtostr_unquoted(gel(histdata, 5)), GENtostr_unquoted(gel(histdata, 6)), itos(gel(histdata, 8)));
}

//Remake the histogram with new minx and maxx
GEN hist_rerange(GEN data, GEN histdata, GEN minx, GEN maxx, long prec){
	if(gequal0(gel(histdata, 7))) return hist_tobins(data, minx, maxx, gel(histdata, 3), itos(gel(histdata, 4)), 1, GENtostr_unquoted(gel(histdata, 5)), GENtostr_unquoted(gel(histdata, 6)), NULL, itos(gel(histdata, 8)), prec);
	return hist_tobins(data, minx, maxx, gel(histdata, 3), itos(gel(histdata, 4)), 1, GENtostr_unquoted(gel(histdata, 5)), GENtostr_unquoted(gel(histdata, 6)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 8)), prec);
}

//Remake the histogram with the new scale
GEN hist_rescale(GEN data, GEN histdata, int scale, long prec){
	if(gequal0(gel(histdata, 7))) return hist_tobins(data, gel(histdata, 1), gel(histdata, 2), gel(histdata, 3), scale, 1, GENtostr_unquoted(gel(histdata, 5)), GENtostr_unquoted(gel(histdata, 6)), NULL, itos(gel(histdata, 8)), prec);
	return hist_tobins(data, gel(histdata, 1), gel(histdata, 2), gel(histdata, 3), scale, 1, GENtostr_unquoted(gel(histdata, 5)), GENtostr_unquoted(gel(histdata, 6)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 8)), prec);
}



//3: REGRESSIONS



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

