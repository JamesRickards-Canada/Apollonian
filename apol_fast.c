/*Methods to quickly compute Apollonian stuff, and supporting things.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/*STATIC DECLARATIONS*/

/*SECTION 1: MISSING CURVATURES*/

/*1: C CODE*/
static void findmissing(long Bmin, long Bmax, long x[], long res[], long lenres, GEN quadfams, GEN quarfams);
static void missing_tofile(long blocks, unsigned long **rclass, GEN quadfams, GEN quarfams, long Bmin, long Bmax, long x[], long res[], long lenres, int families);


/*MAIN BODY*/

/*SECTION 1: MISSING CURVATURES*/


/*1: C CODE*/

/*Updates rclass with the new curvature.*/
inline void
missing_update(unsigned long **rclass, unsigned long *bitswap, long *maxmiss_block, int *maxmiss_bit, long Base, long *Bmax, long curv, long res[], long lenres, long Bmin)
{
  long shifted = curv - Base;
  if (shifted <= 0) return;
  long b = shifted % 24;
  long a = shifted / 24;/*shifted=24a+b. b gives the block to insert into, and we need to save "a"*/
  long v = a % 64;
  long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
  rclass[b][u] |= bitswap[v];
  if (u != maxmiss_block[b]) return;/*Not swapping in the last block.*/
  if (v != maxmiss_bit[b]) return;/*Not swapping the last bit.*/
  long i;
  for (i = v - 1; i >= 0; i--) {
	if (!(rclass[b][u] & bitswap[i])) {/*Zero, so we stop here.*/
	  maxmiss_bit[b] = i;
	  goto BMAXUPDATE;
	}
  }
  long j;/*We made it out of the block.*/
  for (j = u - 1; j >= 0; j--) {
	for (i = 63; i >= 0; i--) {
	  if (!(rclass[b][j] & bitswap[i])) {/*Zero, so we stop here.*/
		maxmiss_block[b] = j;
		maxmiss_bit[b] = i;
		goto BMAXUPDATE;
	  }
	}
  }
  maxmiss_block[b] = -1;/*Everything is gone!*/
  maxmiss_bit[b] = -1;
  BMAXUPDATE:;/*See if we update Bmax*/
  if (curv != *Bmax) return;/*Did not remove the largest exception.*/
  long worst = res[0];
  for (i = 1; i < lenres; i++) {
	if (maxmiss_block[res[i]] > maxmiss_block[worst]) {
	  worst = res[i];
	  continue;
	}
	if (maxmiss_block[res[i]] < maxmiss_block[worst]) continue;
	if (maxmiss_bit[res[i]] >= maxmiss_bit[worst]) worst = res[i];
  }
  *Bmax = Base + ((((maxmiss_block[worst] << 6) + maxmiss_bit[worst]) * 3) << 3 ) + worst;/*Update Bmax*/
  if (Bmin > *Bmax) *Bmax = 0;/*All eliminated, let's quit early!*/
}

/*Finds all missing positive curvatures in the given residue classes between B1 and B2 (inclusive), saving them to a file. Formatting of the inputs is provided by apol_missing; it is crucial that x is reduced, sorted, and res is the set of ALL residues modulo 24.*/
static void
findmissing(long Bmin, long Bmax, long x[], long res[], long lenres, GEN quadfams, GEN quarfams)
{
  unsigned long *bitswap = (unsigned long*)pari_malloc(64 * sizeof(unsigned long)), i;/*Used for swapping bits of longs.*/
  bitswap[0] = 1;
  for (i = 1; i < 64; i++) bitswap[i] = bitswap[i - 1] << 1;/*bitswap[i] = 2^i*/
  long Base = Bmin - (Bmin % 24);/*We want to start at a multiple of 24 to not ruin the mod stuff.*/
  long classmax = (Bmax - Base)/ 24 + 1;/*Maximal number of curvatures found in each class.*/
  long blocks = ((classmax - 1) / 64) + 1;/*Here is where we assume 64-bit. This is the number of 64-bit unsigned longs we need to store in each class.*/
  unsigned long **rclass = (unsigned long **)pari_malloc(24 * sizeof(unsigned long *));/*Stores pointers to the individual classes.*/
  if (!rclass) {
    printf("Insufficient memory to allocate to store the residue classes.\n");
    exit(1);
  }
  for (i = 0; i < lenres; i++) {
    rclass[res[i]] = (unsigned long *)pari_calloc(blocks * sizeof(unsigned long));/*pari_calloc the classes we want, since we want them as 0 to start.*/
    if (!rclass[res[i]]) {
      printf("Insufficient memory to allocate to store the curvatures.\n");
      exit(1);
    }
  }
  long Bmaxoriginal = Bmax;/*Save for later in case we change it.*/
  long *maxmiss_block = (long *)pari_malloc(24 * sizeof(long));/*Tracks the largest block in the residue class that still contains 0's*/
  int *maxmiss_bit = (int *)pari_malloc(24 * sizeof(int));/*Tracks the largest bit of said class that is non-zero.*/
  int Bm24 = Bmax % 24;/*We now will initialize it.*/
  long Bmaxbase = Bmax - Bm24, j;
  int foundlargest = 0;
  for (i = 0; i < lenres; i++) {
	long mc = Bmaxbase + res[i];
	if (res[i] > Bm24) {
	  if (!foundlargest) {
		foundlargest = 1;
		if (i) Bmax = Bmaxbase + res[i - 1];/*Update to the largest actually possible value.*/
		else Bmax = Bmaxbase + res[lenres - 1] - 24;
	  }
	  mc -= 24;/*The last curvature of this type.*/
	}
	if (mc < Bmin) {
	  maxmiss_bit[res[i]] = -1;
	  maxmiss_block[res[i]] = -1;
	  continue;
	}
	mc -= Base;/*Shift it back.*/
	long a = mc / 24;/*Save a in block res[i]*/
	maxmiss_bit[res[i]] = a % 64;
    maxmiss_block[res[i]] = a / 64;
  }
  if (!foundlargest) Bmax = Bmaxbase + res[lenres - 1];/*They all fit in with the +.*/
  for (i = 0; i <= lenres; i++) {
	if (i == lenres) {/*All were -1, so our range actually contains no residues.*/
	  Bmax = 0;
	}
	else if (maxmiss_bit[res[i]] != -1) break;/*We have at least one valid residue.*/
  }
  int families;
  if (quadfams) {/*Remove them right away.*/
    families = 1;
	for (i = 1; i <= lenres; i++) {
	  for (j = 1; j < lg(gel(quadfams, i)); j++) {
		long c = itos(gmael(quadfams, i, j));
		long curres = res[i - 1];/*Indices off by 1 with C array.*/
		long fa = 1, cur = c;/*cur = c * fa^2*/
		while (cur <= Bmax) {
		  missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, cur, res, lenres, Bmin);
		  do {
			fa++;
			cur = c * (fa * fa);
		  }
		  while (cur % 24 != curres);
		}
	  }
	  for (j = 1; j < lg(gel(quarfams, i)); j++) {
		long c = itos(gmael(quarfams, i, j));
		long curres = res[i - 1];/*Indices off by 1 with C array.*/
		long fa = 1, cur = c;/*cur = c * fa^4*/
		while (cur <= Bmax) {
		  missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, cur, res, lenres, Bmin);
		  do {
			fa++;
			cur = fa * fa;
			cur = c * (cur * cur);
		  }
		  while (cur % 24 != curres);
		}
	  }
	}
  }
  else families = 0;
  long maxdepth = 100;/*Maximal depth, to start.*/
  long **depthseq = (long **)pari_malloc(maxdepth * sizeof(long *));/*Tracks the sequence of Apollonian moves.*/
  if (!depthseq) {
    printf("Insufficient memory to allocate to store the depth sequence.\n");
    exit(1);
  }
  int *swaps = (int *)pari_malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps.*/
  if (!swaps) {
    printf("Insufficient memory to allocate to store the swaps.\n");
    exit(1);
  }
  for (i = 0; i < maxdepth; i++) {
    depthseq[i] = (long *)pari_malloc(sizeof(long) << 2);
    swaps[i] = -1;/*Initialize to all -1's*/
  }
  for (i = 1; i < 4; i++) {/*Do the first 3 curvatures (ignore the negative one).*/
    if (x[i] < Bmin || x[i] > Bmax) continue;
	missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, x[i], res, lenres, Bmin);
  }
  /*We adjust the starting quadruple in case of symmetries: for a+b+c=d, we put d first, and DO NOT flip it on the first iteration. If c=d, we do c, a, c, b, starting with the second one. Until one element of the depth sequence flips the third entry, we do not flip the first one, as they will be (c, c) still. There are two exceptions: [0, 0, 1, 1], and [-1, 2, 2, 3], as they have both types of symmetry. For the secon, we do 2, 3, 2, -1, and start at the third entry. For [0, 0, 1, 1], we do [1, 4, 1, 0] and also start at the third one (we do the first move, since it is forced). The variable sym keeps track of this: -1 means no symmetries, don't worry. 0 means symmetris and we have not moved beyond them, hence we cannot flip the first element. >0 means this is the index of depthseq that the first 3 occurs, so we know when we drop back into the symmetry zone.
  */
  long sym;
  if (x[0] == 0) {/*We are the [0, 0, 1, 1] packing: do [1, 4, 1, 0] and start with swapping the third.*/
    depthseq[0][0] = 1; depthseq[0][1] = 4; depthseq[0][2] = 1; depthseq[0][3] = 0;
    rclass[4][0] |= bitswap[0];/*Since we are adjusting, we did not do 4 yet.*/
    swaps[1] = 1;/*Will get incremented to 2 right away.*/
    sym = 0;/*Symmetries to avoid.*/
  }
  else if (x[0] == -1) {/*[-1, 2, 2, 3]: do [2, 3, 2, -1]*/
    depthseq[0][0] = 2; depthseq[0][1] = 3; depthseq[0][2] = 2; depthseq[0][3] = -1;
    swaps[1] = 1;
    sym = 0;/*Symmetries to avoid.*/
  }
  else if (x[1] == x[2]) {/*b=c, so do [b, a, b, d]*/
    depthseq[0][0] = x[1]; depthseq[0][1] = x[0]; depthseq[0][2] = x[2]; depthseq[0][3] = x[3];
    swaps[1] = 0;
    sym = 0;/*Symmetries to avoid.*/
  }
  else if (x[2] == x[3]) {/*c=d, so do [c, a, c, b]*/
    depthseq[0][0] = x[3]; depthseq[0][1] = x[0]; depthseq[0][2] = x[2]; depthseq[0][3] = x[1];
    swaps[1] = 0;
    sym = 0;/*Symmetries to avoid.*/
  }
  else if ((x[0] + x[1] + x[2]) == x[3]) {/*a+b+c=d, so do [d, a, b, c]*/
    depthseq[0][0] = x[3]; depthseq[0][1] = x[0]; depthseq[0][2] = x[1]; depthseq[0][3] = x[2];
    swaps[1] = 0;
    sym = -1;/*After first swap, no symmetries to worry about.*/
  }
  else {
    depthseq[0][0] = x[0]; depthseq[0][1] = x[1]; depthseq[0][2] = x[2]; depthseq[0][3] = x[3];
    sym = -1;/*No symmetries to begin with.*/
  }
  long ind = 1;/*Which depth we are working at.*/
  if (!sym) {/*Symmetries to worry about. More efficient to do this way, since we don't need to check for symmetries ever otherwise.*/
    while (ind > 0) {/*We are coming in trying to swap this circle out.*/
      int cind = ++swaps[ind];/*Increment the swapping index.*/
      if (cind == 4) {/*Overflowed, go back.*/
        swaps[ind] = -1;
        ind--;
        if (ind < sym) sym = 0;/*We moved past the first index swapping 3, so worry about symmetries again.*/
        continue;
      }
      long lastind = ind - 1;
      if (cind == swaps[lastind]) continue; /*Same thing twice, so skip it.*/
      if (!sym) {/*Worry about symmetries.*/
        if (cind == 0) continue;/*We skip the first one.*/
        else if (cind == 2) sym = ind + 1;/*First time we swap out the third one, eliminating symmetries further on in this branch.*/
      }
      long apbpc = 0;/*Now we can reasonably try a swap.*/
      for (i = 0; i < cind; i++) apbpc += depthseq[lastind][i];
      for (i = cind + 1; i < 4; i++) apbpc += depthseq[lastind][i];
      long newc = (apbpc << 1) - depthseq[lastind][cind];/*2(a+b+c)-d, the new curvature.*/
      if (newc > Bmax) {/*Too big! go back.*/
        if (ind < sym) sym = 0;/*Tried flipping out of symmetry here but it's too big.*/
        continue;
      }
      /*Do the bitswap to update the count if we are large enough.*/
	  missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, newc, res, lenres, Bmin);
      for (i = 0; i < cind; i++) depthseq[ind][i] = depthseq[lastind][i];
      depthseq[ind][cind] = newc;
      for (i = cind + 1; i < 4; i++) depthseq[ind][i] = depthseq[lastind][i];/*Add the tuple in.*/
      ind++;
      if (ind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
        long newdepth = maxdepth << 1;/*Double it.*/
        depthseq = pari_realloc(depthseq, newdepth * sizeof(long *));
        if (!depthseq) {
          printf("Insufficient memory to pari_reallocate the depth sequence.\n");
          exit(1);
        }
        swaps = pari_realloc(swaps, newdepth * sizeof(int));
        if (!swaps) {
          printf("Insufficient memory to pari_reallocate the swaps.\n");
          exit(1);
        }
        for (i = maxdepth; i < newdepth; i++) {
          depthseq[i] = (long *)pari_malloc(sizeof(long) << 2);
          swaps[i] = -1;
        }
        maxdepth = newdepth;
      }
    } 
  }
  else {/*No symmetry to worry about.*/
    while (ind > 0) {/*We are coming in trying to swap this circle out.*/
      int cind = ++swaps[ind];/*Increment the swapping index.*/
      if (cind == 4) {/*Overflowed, go back.*/
        swaps[ind] = -1;
        ind--;
        continue;
      }
      long lastind = ind - 1;
      if (cind == swaps[lastind]) continue; /*Same thing twice, so skip it.*/
      long apbpc = 0;/*Now we can reasonably try a swap.*/
      for (i = 0; i < cind; i++) apbpc += depthseq[lastind][i];
      for (i = cind + 1; i < 4; i++) apbpc += depthseq[lastind][i];
      long newc = (apbpc << 1) - depthseq[lastind][cind];/*2(a+b+c)-d, the new curvature.*/
      if (newc > Bmax) continue;/*Too big! go back.*/
      /*Do the bitswap to update the count.*/
	  missing_update(rclass, bitswap, maxmiss_block, maxmiss_bit, Base, &Bmax, newc, res, lenres, Bmin);
	  for (i = 0; i < cind; i++) depthseq[ind][i] = depthseq[lastind][i];
      depthseq[ind][cind] = newc;
      for (i = cind + 1; i < 4; i++) depthseq[ind][i] = depthseq[lastind][i];/*Add the tuple in.*/
      ind++;
      if (ind == maxdepth) {/*We are going too deep, must pari_reallocate the storage location.*/
        long newdepth = maxdepth << 1;/*Double it.*/
        depthseq = pari_realloc(depthseq, newdepth * sizeof(long *));
        if (!depthseq) {
          printf("Insufficient memory to pari_reallocate the depth sequence.\n");
          exit(1);
        }
        swaps = pari_realloc(swaps, newdepth * sizeof(int));
        if (!swaps) {
          printf("Insufficient memory to pari_reallocate the swaps.\n");
          exit(1);
        }
        for (i = maxdepth; i < newdepth; i++) {
          depthseq[i] = (long *)pari_malloc(sizeof(long) << 2);
          swaps[i] = -1;
        }
        maxdepth = newdepth;
      }
    }
  }
  missing_tofile(blocks, rclass, quadfams, quarfams, Bmin, Bmaxoriginal, x, res, lenres, families);/*Print to file.*/
  /*Time to free all of the allocated memory.*/
  pari_free(swaps);
  for (i = 0; i < maxdepth; i++) pari_free(depthseq[i]);
  pari_free(depthseq);
  pari_free(maxmiss_block);
  pari_free(maxmiss_bit);
  for (i = 0; i < lenres; i++) pari_free(rclass[res[i]]);
  pari_free(rclass);
  pari_free(bitswap);
  return;
}

/*Prints the found data to a file.*/
static void
missing_tofile(long blocks, unsigned long **rclass, GEN quadfams, GEN quarfams, long Bmin, long Bmax, long x[], long res[], long lenres, int families)
{
  long Base = Bmin - 1;
  Base = Base - (Base % 24);/*We started at a multiple of 24 to not ruin the mod stuff.*/
  char fname[200];
  int pos = 0;
  DIR* dir = opendir("missing");
  if (dir) {
    closedir(dir);/* Directory exists. */
    pos += sprintf(&fname[pos], "missing/");
  }
  else if (ENOENT == errno) {/* Directory does not exist. */
    mkdir("missing", 0777);
    pos += sprintf(&fname[pos], "missing/");
  }
  else {/* opendir() failed for some other reason. */
    printf("Could not create directory, saving to current folder");
  }
  if (x[0] < 0) pos += sprintf(&fname[pos], "m%ld_", -x[0]);
  else pos += sprintf(&fname[pos], "%ld_", x[0]);
  pos += sprintf(&fname[pos], "%ld_%ld_%ld_%ld-to-%ld", x[1], x[2], x[3], Bmin, Bmax);
  if (families) pos += sprintf(&fname[pos], "_remqq");
  pos += sprintf(&fname[pos], ".dat");
  long i;
  FILE *F;
  F = fopen(fname, "w");
  if (families) {/*Print the residue classes and the families found.*/
    fprintf(F, "[%ld", res[0]);
    for (i = 1; i < lenres; i++) fprintf(F, ", %ld", res[i]);
    fprintf(F, "]\n[");
    for (i = 1; i <= lenres; i++) {
      if (i > 1) fprintf(F, ", ");
      fprintf(F, "[");
      long j;
      for (j = 1; j < lg(gel(quadfams, i)); j++) {
        if (j > 1) fprintf(F, ", ");
        fprintf(F, "%ld", itos(gmael(quadfams, i, j)));
      }
      fprintf(F, "]");
    }
    fprintf(F, "]\n[");
    for (i = 1; i <= lenres; i++) {
      if (i > 1) fprintf(F, ", ");
      fprintf(F, "[");
      long j;
      for (j = 1; j < lg(gel(quarfams, i)); j++) {
        if (j > 1) fprintf(F, ", ");
        fprintf(F, "%ld", itos(gmael(quarfams, i, j)));
      }
      fprintf(F, "]");
    }
    fprintf(F, "]\n");
  }
  for (i = 0; i < lenres; i++) {/*The ith residue*/
    long b = res[i];
    fprintf(F, "[");
    int found = 0;
    long u;
    for (u = 0; u < blocks; u++) {/*Figure out each block.*/
      unsigned long val = rclass[b][u];
      long v;
      for (v = 0; v < 64; v++) {
        if (!(val & 1)) {/*A missing value*/
          long a = (u << 6) + v;
          long n = 24 * a + b + Base;/*The correct one!*/
          if (n <= Bmax && n >= Bmin) {/*Correct range! We only check >=Bmin due to the initial shift by up to 24.*/
            if (found) fprintf(F, ", %ld", n);/*Print it to the file.*/
            else { found = 1; fprintf(F, "%ld", n); }
          }
        }
        val >>= 1;/*Shift right one.*/
      }
    }
    fprintf(F, "]\n");
  }
  fclose(F);
}


/*1: GP ACCESS*/

/*Runs C code to find the missing curvatures up to the given bound, then returns them in a vector. If family=1, removes the known families.*/
GEN
apol_missing(GEN v, GEN B, int family, int load)
{
  pari_sp av = avma;
  long Bmin = 0, Bmax = 0, t = typ(B);
  if (t == t_INT) { Bmin = 1; Bmax = itos(B); }
  else if (t == t_VEC || t == t_COL) {
	if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	if (typ(gel(B, 1)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	if (typ(gel(B, 2)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	Bmin = itou(gel(B, 1));
	Bmax = itou(gel(B, 2));
  }
  else if (t == t_VECSMALL) {
	if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	Bmin = B[1];
	Bmax = B[2];
  }
  else pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
  GEN w = apol_red(v, 0, 0);
  w = ZV_sort(w);
  GEN modres = apol_mod24(w);
  long x[4], i;
  for (i = 1; i <= 4; i++) x[i - 1] = itos(gel(w, i));
  long lr = lg(modres), lenr = lr - 1, res[lenr];
  for (i = 1; i < lr; i++) res[i - 1] = itos(gel(modres, i));
  if (family) {
	GEN fams = apol_missingfamilies(w);
	findmissing(Bmin, Bmax, x, res, lenr, gel(fams, 1), gel(fams, 2));
  }
  else findmissing(Bmin, Bmax, x, res, lenr, NULL, NULL);
  if (!load) return gc_const(av, gen_0);/*Do not load.*/
  set_avma(av);
  return apol_missing_load(v, B, family);
}

/*Loads the saved file of curvatures.*/
GEN
apol_missing_load(GEN v, GEN B, int family)
{
  pari_sp av = avma;
  long Bmin = 0, Bmax = 0, t = typ(B);
  if (t == t_INT) { Bmin = 1; Bmax = itos(B); }
  else if (t == t_VEC || t == t_COL) {
	if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	if (typ(gel(B, 1)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	if (typ(gel(B, 2)) != t_INT) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	Bmin = itos(gel(B, 1));
	Bmax = itos(gel(B, 2));
  }
  else if (t == t_VECSMALL) {
	if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	Bmin = B[1];
	Bmax = B[2];
  }
  else pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
  v = ZV_sort(apol_red(v, 0, 0));
  char *fname;
  if (signe(gel(v, 1)) < 0) fname = pari_sprintf("m%Pd", negi(gel(v, 1)));
  else fname = pari_sprintf("%Pd", gel(v, 1));
  long i;
  for (i = 2; i <= 4; i++) fname = pari_sprintf("%s_%Pd", fname, gel(v, i));
  if (family) fname = pari_sprintf("%s_%ld-to-%ld_remqq.dat", fname, Bmin, Bmax);
  else fname = pari_sprintf("%s_%ld-to-%ld.dat", fname, Bmin, Bmax);
  if (pari_is_dir("missing")) {
    fname = pari_sprintf("missing/%s", fname);
  }
  set_avma(av);
  return gp_readvec_file(fname);/*Load them up!*/
}

/*Returns the families associated to v.*/
GEN
apol_missingfamilies(GEN v)
{
  pari_sp av = avma;
  GEN t = apol_type(v, 2), quad = NULL, quar = NULL;
  long types = (itos(gel(t, 1)) >> 1) + itos(gel(t, 2)) + itos(gel(t, 3));/*n/2+k+chi_2 where the type is (n, k): uniquely identifies it.*/
  switch (types) {
	case 5:/*6.1.1*/
	  quad = zerovec(6);
	  if (equali1(gel(t, 4))) quar = zerovec(6);
	  else quar = mkvecn(6, mkvec4s(144, 576, 1296, 5184), mkvecs(1), mkvecs(4), mkvec2s(9, 81), mkvec2s(36, 324), mkvec2s(16, 64));
	  break;
	case 3:/*6.1.-1*/
	  quad = mkvecn(6, mkvec4s(24, 48, 72, 144), mkvecs(1), mkvecs(4), mkvecs(9), mkvec2s(12, 36), mkvecs(16));
	  quar = zerovec(6);
	  break;
	case 9:/*6.5.1*/
	  quad = mkvecn(6, mkvec2s(48, 72), cgetg(1, t_VEC), mkvecs(8), mkvecs(12), cgetg(1, t_VEC), cgetg(1, t_VEC));
	  quar = zerovec(6);
	  break;
	case 7:/*6.5.-1*/
	  quad = mkvecn(6, mkvec2s(24, 144), cgetg(1, t_VEC), cgetg(1, t_VEC), mkvecs(36), cgetg(1, t_VEC), cgetg(1, t_VEC));
	  quar = zerovec(6);
	  break;
	case 17:/*6.13.1*/
	  quad = mkvecn(6, mkvec2s(24, 72), cgetg(1, t_VEC), cgetg(1, t_VEC), cgetg(1, t_VEC), cgetg(1, t_VEC), cgetg(1, t_VEC));
	  quar = zerovec(6);
	  break;
	case 15:/*6.13.-1*/
	  quad = mkvecn(6, mkvec2s(48, 144), mkvecs(4), mkvec2s(12, 36), cgetg(1, t_VEC), mkvecs(16), cgetg(1, t_VEC));
	  quar = zerovec(6);
	  break;
	case 21:/*6.17.1*/
	  quad = mkvecn(6, mkvec2s(24, 48), cgetg(1, t_VEC), cgetg(1, t_VEC), mkvecs(12), cgetg(1, t_VEC), cgetg(1, t_VEC));
	  if (equali1(gel(t, 4))) quar = mkvecn(6, mkvec2s(144, 576), cgetg(1, t_VEC), mkvecs(9), mkvecs(36), cgetg(1, t_VEC), cgetg(1, t_VEC));
	  else quar = mkvecn(6, mkvec2s(1296, 5184), cgetg(1, t_VEC), mkvecs(81), mkvecs(324), cgetg(1, t_VEC), cgetg(1, t_VEC));
	  break;
	case 19:/*6.17.-1*/
	  quad = mkvecn(6, mkvec2s(72, 144), mkvecs(8), mkvecs(9), mkvecs(36), cgetg(1, t_VEC), cgetg(1, t_VEC));
	  quar = zerovec(6);
	  break;
	case 12:/*8.7.1*/
	  quad = zerovec(8);
	  gel(quad, 1) = mkvecs(3); gel(quad, 2) = mkvecs(6);
	  quar = zerovec(8);
	  break;
	case 10:/*8.7.-1*/
	  quad = zerovec(8);
	  gel(quad, 6) = mkvecs(18);
	  quar = zerovec(8);
	  break;
	case 16:/*8.11.1*/
	  quad = zerovec(8);
	  quar = zerovec(8);
	  break;
	case 14:/*8.11.-1*/
	  quad = zerovec(8);
	  gel(quad, 1) = mkvecs(2); gel(quad, 2) = mkvecs(3);
	  gel(quad, 3) = mkvecs(6); gel(quad, 7) = mkvecs(18);
	  quar = zerovec(8);
	  break;
	default:
	  pari_err_TYPE("v did not have a recongizable type, is it a Descartes quadruple?", v);
  }
  long i;
  for (i = 1; i < lg(quad); i++) {/*Fixing the uninitialized entries.*/
	if (isintzero(gel(quad, i))) gel(quad, i) = cgetg(1, t_VEC);
	if (isintzero(gel(quar, i))) gel(quar, i) = cgetg(1, t_VEC);
  }
  return gerepilecopy(av, mkvec2(quad, quar));
}

