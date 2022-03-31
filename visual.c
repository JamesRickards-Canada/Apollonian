//ORIGINALLY IMPORTED FROM Q-QUADRATIC
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



//DATA
//Should also write the data to a file and offer ways to make the histogram from a file of data.


//SEE ALSO vec_equiv: this does something very similar.
//Returns [vsort, count], where vsort is the sorted vector v with duplicates removed, and count is the Vecsmall of corresponding number of each in the original vector v. This is not the most efficient, but is fine.
GEN veccount(GEN v){
  pari_sp top=avma;
  GEN vsort=sort(v);//Sort it.
  long lv=lg(vsort);
  GEN uniq=vectrunc_init(lv), count=vecsmalltrunc_init(lv);
  vectrunc_append(uniq, gel(vsort, 1));
  long run=1;
  for(long i=2;i<lv;i++){
    if(gequal(gel(vsort, i), gel(vsort, i-1))){run++;continue;}//Go on
    vecsmalltrunc_append(count, run);//run is over.
    run=1;
    vectrunc_append(uniq, gel(vsort, i));//Add the new number in.
  }
  vecsmalltrunc_append(count, run);
  return gerepilecopy(top, mkvec2(uniq, count));
}

//veccount, but for a vecsmall
GEN vecsmallcount(GEN v){
  pari_sp top=avma;
  GEN vsort=gcopy(v);
  vecsmall_sort(vsort);//Sort it.
  long lv=lg(vsort);
  GEN uniq=vecsmalltrunc_init(lv), count=vecsmalltrunc_init(lv);
  vecsmalltrunc_append(uniq, vsort[1]);
  long run=1;
  for(long i=2;i<lv;i++){
    if(vsort[i]==vsort[i-1]){run++;continue;}//Go on
    vecsmalltrunc_append(count, run);//run is over.
    run=1;
    vecsmalltrunc_append(uniq, vsort[i]);//Add the new number in.
  }
  vecsmalltrunc_append(count, run);
  return gerepilecopy(top, mkvec2(uniq, count));
}



//HISTOGRAMS


//Creates LaTeX document and compiles the histogram. If plotoptions!=NULL, adds this between begin and end tikzpicture
void hist_autocompile(GEN minx, GEN maxx, char *imagename, char *plotoptions, int open){
  pari_sp top=avma;
  if(!pari_is_dir("images/build")){
    int s=system("mkdir -p images/build");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY images/build");
  }
  char *autofilestr=pari_sprintf("images/build/%s_build.tex", imagename);
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
  tex_compile(imagename, open);
  avma=top;
}

//Makes a histogram with the inputs, defaulting to the defaults. Adjust it to your liking with the hist_re-methods. Should only be called once with the given data.
GEN hist_make(GEN data, char *imagename, int compilenew, int open, char *plotoptions, long prec){
  return hist_tobins_defaultbins(data, gel(data, 1), gel(data, lg(data)-1), 0, compilenew, imagename, plotoptions, open, prec);
}

//Bins the data according to the parameters. Returns the associated parameters
GEN hist_tobins(GEN data, GEN minx, GEN maxx, GEN nbins, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec){
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
  if(compilenew) hist_autocompile(minx, maxx, imagename, plotoptions, open);//Making new LaTeX document and compiling it
  else tex_compile(imagename, open);//Just recompiling the old LaTeX file
  pari_printf("%d total data points\n%d total bins", lg(data)-1, nbinslong);
  GEN histprop=cgetg(8, t_VEC);
  gel(histprop, 1)=strtoGENstr(imagename);
  if(open==0) gel(histprop, 2)=gen_0;
  else gel(histprop, 2)=gen_1;
  gel(histprop, 3)=gcopy(minx);gel(histprop, 4)=gcopy(maxx);
  gel(histprop, 5)=icopy(nbins);gel(histprop, 6)=stoi(toscale);
  if(plotoptions==NULL) gel(histprop, 7)=gen_0;
  else gel(histprop, 7)=strtoGENstr(plotoptions);
  return gerepileupto(top, histprop);
}

//Bins the data with default options, outputting to imagename.
GEN hist_tobins_defaultbins(GEN data, GEN minx, GEN maxx, int toscale, int compilenew, char *imagename, char *plotoptions, int open, long prec){
  pari_sp top=avma;
  long ind1=(lg(data)-1)/4, ind2=3*ind1;//The difference approximately gives the IQR of the data
  GEN binlen=gdiv(gmulsg(2,gsub(gel(data, ind2), gel(data, ind1))),gpow(stoi(lg(data)-1), gdivgs(gen_1, 3), prec));//Using the Freedman-Diaconis rule of bin width=2*IQR/(n^(1/3))
  GEN nbins=mpfloor(gdiv(gsub(maxx, minx), binlen));
  return gerepileupto(top, hist_tobins(data, minx, maxx, nbins, toscale, compilenew, imagename, plotoptions, open, prec));
}

//Remake the histogram with the new bin count
GEN hist_rebin(GEN data, GEN histdata, GEN nbins, long prec){
    if(gequal0(gel(histdata, 7))) return hist_tobins(data, gel(histdata, 3), gel(histdata, 4), nbins, itos(gel(histdata, 6)), 0, GENtostr_unquoted(gel(histdata, 1)), NULL, itos(gel(histdata, 2)), prec);
    return hist_tobins(data, gel(histdata, 3), gel(histdata, 4), nbins, itos(gel(histdata, 6)), 0, GENtostr_unquoted(gel(histdata, 1)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 2)), prec);
}

//Remake the histogram with new minx and maxx
GEN hist_rerange(GEN data, GEN histdata, GEN minx, GEN maxx, long prec){
    if(gequal0(gel(histdata, 7))) return hist_tobins(data, minx, maxx, gel(histdata, 5), itos(gel(histdata, 6)), 1, GENtostr_unquoted(gel(histdata, 1)), NULL, itos(gel(histdata, 2)), prec);
    return hist_tobins(data, minx, maxx, gel(histdata, 5), itos(gel(histdata, 6)), 1, GENtostr_unquoted(gel(histdata, 1)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 2)), prec);
}

//Remake the histogram with the new scale
GEN hist_rescale(GEN data, GEN histdata, int scale, long prec){
    if(gequal0(gel(histdata, 7))) return hist_tobins(data, gel(histdata, 3), gel(histdata, 4), gel(histdata, 5), scale, 1, GENtostr_unquoted(gel(histdata, 1)), NULL, itos(gel(histdata, 2)), prec);
    return hist_tobins(data, gel(histdata, 3), gel(histdata, 4), gel(histdata, 5), scale, 1, GENtostr_unquoted(gel(histdata, 1)), GENtostr_unquoted(gel(histdata, 7)), itos(gel(histdata, 2)), prec);
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

