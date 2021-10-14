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

