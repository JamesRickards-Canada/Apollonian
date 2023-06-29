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
static void findmissing(long Bmin, long Bmax, long x[], long res[], long lenres, int families);
static void missing_update_Bmin(unsigned long **rclass, long Base, unsigned long *bitswap, long curv);
static void missing_tofile(long blocks, unsigned long **rclass, unsigned long **quadfams, unsigned long **quarfams, long Bmin, long Bmax, long x[], long res[], long lenres, int families);
static void findfamilies(long Bmin, long Bmax, unsigned long **quadfams, unsigned long **quarfams, unsigned long **rclass, long res[], long lenres, unsigned long *bitswap);
static int removefamily(long d, long Bmin, long Bmax, unsigned long *curvs, unsigned long *bitswap, long u, long cop);
static int iscop(long n, long cop);

/*SECTION 1: MISSING CURVATURES*/


/*1: C CODE*/

/*Finds all missing positive curvatures in the given residue classes between B1 and B2 (inclusive), saving them to a file. Formatting of the inputs is provided by apol_missing; it is crucial that x is reduced, sorted, and res is the set of ALL residues modulo 24.*/
static void
findmissing(long Bmin, long Bmax, long x[], long res[], long lenres, int families)
{
  unsigned long *bitswap = (unsigned long*)pari_malloc(64 * sizeof(unsigned long)), i;/*Used for swapping bits of longs.*/
  bitswap[0] = 1;
  for (i = 1; i < 64; i++) bitswap[i] = bitswap[i - 1] << 1;/*bitswap[i] = 2^i*/
  long Base = Bmin - 1;
  Base = Base - (Base % 24);/*We want to start at a multiple of 24 to not ruin the mod stuff.*/
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
	missing_update_Bmin(rclass, Base, bitswap, x[i]);
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
	  missing_update_Bmin(rclass, Base, bitswap, newc);
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
	  missing_update_Bmin(rclass, Base, bitswap, newc);
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
  unsigned long **quadfams = NULL;
  unsigned long **quarfams = NULL;
  if (families) {
    quadfams = (unsigned long **)pari_malloc(lenres * sizeof(unsigned long *));/*Pointers to the families that exist. The first element of each family is the number of families in that class.*/
    quarfams = (unsigned long **)pari_malloc(lenres * sizeof(unsigned long *));
    findfamilies(Bmin, Bmax, quadfams, quarfams, rclass, res, lenres, bitswap);/*Find and remove the families.*/
  }
  missing_tofile(blocks, rclass, quadfams, quarfams, Bmin, Bmax, x, res, lenres, families);/*Print to file.*/
  /*Time to free all of the allocated memory.*/
  if (families) {
    for (i = 0; i < lenres; i++) { pari_free(quadfams[i]); pari_free(quarfams[i]); }
    pari_free(quadfams);
    pari_free(quarfams);
  }
  pari_free(swaps);
  for (i = 0; i < maxdepth; i++) pari_free(depthseq[i]);
  pari_free(depthseq);
  for (i = 0; i < lenres; i++) pari_free(rclass[res[i]]);
  pari_free(rclass);
  pari_free(bitswap);
  return;
}

/*Updates rclass with the new curvature, assuming a Bmin*/
static void
missing_update_Bmin(unsigned long **rclass, long Base, unsigned long *bitswap, long curv)
{
  long shifted = curv - Base;
  if (shifted <= 0) return;
  long b = shifted % 24;
  long a = shifted / 24;/*shifted=24a+b. b gives the block to insert into, and we need to save "a"*/
  long v = a % 64;
  long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
  rclass[b][u] |= bitswap[v];
}

/*Prints the found data to a file.*/
static void
missing_tofile(long blocks, unsigned long **rclass, unsigned long **quadfams, unsigned long **quarfams, long Bmin, long Bmax, long x[], long res[], long lenres, int families)
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
    for (i = 0; i < lenres; i++) {
      if (i) fprintf(F, ", ");
      fprintf(F, "[");
      long j;
      for (j = 1; j <= quadfams[i][0]; j++) {
        if (j > 1) fprintf(F, ", ");
        fprintf(F, "%ld", quadfams[i][j]);
      }
      fprintf(F, "]");
    }
    fprintf(F, "]\n[");
    for (i = 0; i < lenres; i++) {
      if (i) fprintf(F, ", ");
      fprintf(F, "[");
      long j;
      for (j = 1; j <= quarfams[i][0]; j++) {
        if (j > 1) fprintf(F, ", ");
        fprintf(F, "%ld", quarfams[i][j]);
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

/*Updates the families found, and updates rclass to have them removed too.*/
static void
findfamilies(long Bmin, long Bmax, unsigned long **quadfams, unsigned long **quarfams, unsigned long **rclass, long res[], long lenres, unsigned long *bitswap)
{
  long i, fampos;
  for (i = 0; i < lenres; i++) {
    switch (res[i]) {
      case 0:
        quadfams[i] = (unsigned long *)pari_malloc(5 * sizeof(unsigned long));/*Up to 4 families.*/
        fampos = 1;
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 24, 1)) {/*We have one family!*/
          quadfams[i][fampos] = 24;
          fampos++;
        }
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 48, 1)) {
          quadfams[i][fampos] = 48;
          fampos++;
        }
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 72, 1)) {
          quadfams[i][fampos] = 72;
          fampos++;
        }
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 144, 1)) {
          quadfams[i][fampos] = 144;
          fampos++;
        }
        quadfams[i][0] = fampos - 1;/*The size.*/
        quarfams[i] = (unsigned long *)pari_malloc(5 * sizeof(unsigned long));/*Up to 4 families.*/
        fampos = 1;
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 144, 1)) {
          quarfams[i][fampos] = 144;
          fampos++;
        }
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 576, 1)) {
          quarfams[i][fampos] = 576;
          fampos++;
        }
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 1296, 1)) {
          quarfams[i][fampos] = 1296;
          fampos++;
        }
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 5184, 1)) {
          quarfams[i][fampos] = 5184;
          fampos++;
        }
        quarfams[i][0] = fampos - 1;
        continue;
      case 1:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 1, 6)) {
          quadfams[i][1] = 1;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 1, 6)) {
          quarfams[i][1] = 1;
          quarfams[i][0] = 1;
        }
        else quarfams[i][0] = 0;
        continue;
      case 2:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 2, 6)) {
          quadfams[i][1] = 2;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);
        quarfams[i][0] = 0;/*None known for 2 mod 24.*/
        continue;
      case 3:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 3, 2)) {
          quadfams[i][1] = 3;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);
        quarfams[i][0] = 0;/*None known for 3 mod 24.*/
        continue;
      case 4:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 4, 6)) {
          quadfams[i][1] = 4;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 4, 6)) {
          quarfams[i][1] = 4;
          quarfams[i][0] = 1;
        }
        else quarfams[i][0] = 0;
        continue;
      case 6:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 6, 2)) {
          quadfams[i][1] = 6;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);
        quarfams[i][0] = 0;/*None known for 6 mod 24.*/
        continue;
      case 8:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 8, 3)) {
          quadfams[i][1] = 8;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);
        quarfams[i][0] = 0;/*None known for 8 mod 24.*/
        continue;
      case 9:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 9, 2)) {
          quadfams[i][1] = 9;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(3 * sizeof(unsigned long));/*Up to 2 families.*/
        fampos = 1;
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 9, 2)) {
          quarfams[i][fampos] = 9;
          fampos++;
        }
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 81, 2)) {
          quarfams[i][fampos] = 81;
          fampos++;
        }
        quarfams[i][0] = fampos - 1;
        continue;
      case 12:
        quadfams[i] = (unsigned long *)pari_malloc(3 * sizeof(unsigned long));/*Up to 2 families.*/
        fampos = 1;
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 12, 2)) {
          quadfams[i][fampos] = 12;
          fampos++;
        }
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 36, 2)) {
          quadfams[i][fampos] = 36;
          fampos++;
        }
        quadfams[i][0] = fampos - 1;
        quarfams[i] = (unsigned long *)pari_malloc(3 * sizeof(unsigned long));/*Up to 2 families.*/
        fampos = 1;
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 36, 2)) {
          quarfams[i][fampos] = 36;
          fampos++;
        }
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 324, 2)) {
          quarfams[i][fampos] = 324;
          fampos++;
        }
        quarfams[i][0] = fampos - 1;
        continue;
      case 16:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 16, 3)) {
          quadfams[i][1] = 16;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(3 * sizeof(unsigned long));/*Up to 2 families.*/
        fampos = 1;
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 16, 3)) {
          quarfams[i][fampos] = 16;
          fampos++;
        }
        if (removefamily(4, Bmin, Bmax, rclass[res[i]], bitswap, 64, 3)) {
          quarfams[i][fampos] = 64;
          fampos++;
        }
        quarfams[i][0] = fampos - 1;
        continue;
      case 18:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
        if (removefamily(2, Bmin, Bmax, rclass[res[i]], bitswap, 18, 2)) {
          quadfams[i][1] = 18;
          quadfams[i][0] = 1;
        }
        else quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);
        quarfams[i][0] = 0;/*None known for 18 mod 24.*/
        continue;
      default:
        quadfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long));/*No family*/
        quadfams[i][0] = 0;
        quarfams[i] = (unsigned long *)pari_malloc(sizeof(unsigned long) << 1);
        quarfams[i][0] = 0;
    }
  }
}

/*Looks through the family n=c*x^d, x coprime to cop, x>0, and d=2 or 4. If they truly are all missing, we return 1, and flip all the bits in curvs.*/
static int
removefamily(long d, long Bmin, long Bmax, unsigned long *curvs, unsigned long *bitswap, long c, long cop)
{
  long Base = Bmin - 1;
  Base = Base - (Base % 24);/*We started at a multiple of 24 to not ruin the mod stuff.*/
  long x = 1, n = c;
  while (n <= Bmax) {
	if (n >= Bmin) {
	  n = n - Base;
      long a = n / 24;/*n = 24a + residue*/
      long v = a % 64;
      long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
      if ((curvs[u] | bitswap[v]) == curvs[u]) return 0;/*This value was found, family not missing.*/
	}
    do {x++;} while (!iscop(x, cop));
	if (d == 4) n = x * x;
	else n = x;
	n = c * n * n;
  }
  /*If we make it here, the family exists! We now must redo the above and actually do the swaps.*/
  x = 1;
  n = c;
  while (n <= Bmax) {
	if (n >= Bmin) {
	  n = n - Base;
      long a = n / 24;/*n = 24a + residue*/
      long v = a % 64;
      long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
      curvs[u] |= bitswap[v];
    }
	do {x++;} while (!iscop(x, cop));
    if (d == 4) n = x * x;
	else n = x;
	n = c * n * n;
  }
  return 1;
}

/*Returns 1 if n is coprime to cop, where cop=1, 2, 3, or 6.*/
static int
iscop(long n, long cop)
{
  switch (cop) {
    case 1: return 1;
    case 2: return n % 2;
    case 3:
      if (n % 3) return 1;
      return 0;
    default:
      if (!(n % 2)) return 0;
      if (!(n % 3)) return 0;
      return 1;
  }
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
	Bmin = itos(gel(B, 1));
	Bmax = itos(gel(B, 2));
  }
  else if (t == t_VECSMALL) {
	if (lg(B) < 3) pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
	Bmin = B[1];
	Bmax = B[2];
  }
  else pari_err_TYPE("B must be a positive integer or a range of positive integers", B);
  GEN w = apol_red(v, 0, 0);
  GEN modres = apol_mod24(w);
  w = ZV_sort(w);
  long x[4], i;
  for (i = 1; i <= 4; i++) x[i - 1] = itos(gel(w, i));
  long lr = lg(modres), res[lr - 1];
  for (i = 1; i < lr; i++) res[i - 1] = itos(gel(modres, i));
  findmissing(Bmin, Bmax, x, res, lr - 1, family);
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
  if (family) fname = pari_sprintf("%s_%d-to-%d_remqq.dat", fname, Bmin, Bmax);
  else fname = pari_sprintf("%s_%d-to-%d.dat", fname, Bmin, Bmax);
  if (pari_is_dir("missing")) {
    fname = pari_sprintf("missing/%s", fname);
  }
  set_avma(av);
  return gp_readvec_file(fname);/*Load them up!*/
}

