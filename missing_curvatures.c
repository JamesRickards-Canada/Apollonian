/*
Compile with: 	gcc missing_curvatures.c -o missing_curvatures
Run with:		./missing_curvatures bound a b c d r1 ... rn, 
	where the reduced ACP in nondecreasing order is [a, b, c, d], we go up to bound, and r1, ..., rn are the possible residues modulo 24.
In order to do larger computations, we need a 64-bit operating system, so will assume that we have it.
*/

/*INCLUSIONS*/
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/*DECLARATIONS*/
void findmissing(long B, long x[], long res[], long lenres, int families);
unsigned long **findfamilies(long B, unsigned long **rclass, long res[], long lenres, unsigned long *bitswap);
int removefamily(long B, long *curvs, unsigned long *bitswap, long u, long cop);
static int iscop(long n, long cop);



int
main(int argc, char *argv[])
{
  if (argc < 13) return 1;/*Problem.*/
  int families = atoi(argv[1]);/*Family or not.*/
  long B = atol(argv[2]);/*Bound.*/
  long x[4], i;/*The circle packing.*/
  for (i = 3; i <= 6; i++) x[i - 3] = atol(argv[i]);
  long res[argc - 7];/*The residue classes.*/
  for (i = 7; i < argc; i++) res[i - 7] = atol(argv[i]);
  clock_t tstart = clock();
  findmissing(B, x, res, argc - 7, families);
  clock_t tend = clock();
  double time_spent = (double)(tend - tstart) / CLOCKS_PER_SEC;
  printf("Total computation time: %fs\n", time_spent);
  return 0;
}

/*Finds all missing positive curvatures in the given residue classes.*/
void
findmissing(long B, long x[], long res[], long lenres, int families)
{
  unsigned long *bitswap = (unsigned long*)malloc(64 * sizeof(unsigned long)), i;/*Used for swapping bits of longs.*/
  bitswap[0] = 1;
  for (i = 1; i < 64; i++) bitswap[i] = bitswap[i - 1] << 1;/*bitswap[i] = 2^i*/
  long classmax = B / 24 + 1;/*Maximal number of curvatures found in each class.*/
  long blocks = ((classmax - 1) / 64) + 1;/*Here is where we assume 64-bit. This is the number of 64-bit longs we need to store in each class.*/
  unsigned long **rclass = (unsigned long **)malloc(24 * sizeof(unsigned long *));/*Stores pointers to the individual classes.*/
  if (!rclass) {
	printf("Insufficient memory to allocate to store the residue classes.\n");
	exit(1);
  }
  for (i = 0; i < lenres; i++) {
	rclass[res[i]] = (unsigned long *)calloc(blocks, sizeof(unsigned long));/*calloc the classes we want, since we want them as 0 to start.*/
	if (!rclass[res[i]]) {
	  printf("Insufficient memory to allocate to store the curvatures.\n");
	  exit(1);
    }
  }
  long maxdepth = 100;/*Maximal depth, to start.*/
  long **depthseq = (long **)malloc(maxdepth * sizeof(long *));/*Tracks the sequence of Apollonian moves.*/
  if (!depthseq) {
	printf("Insufficient memory to allocate to store the depth sequence.\n");
	exit(1);
  }
  int *swaps = (int *)malloc(maxdepth * sizeof(int));/*Tracks the sequence of swaps.*/
  if (!swaps) {
	printf("Insufficient memory to allocate to store the swaps.\n");
	exit(1);
  }
  for (i = 0; i < maxdepth; i++) {
	depthseq[i] = (long *)malloc(sizeof(long) << 2);
	swaps[i] = -1;/*Initialize to all -1's*/
  }
  for (i = 1; i < 4; i++) {/*Do the first 3 curvatures (ignore the negative one).*/
	long b = x[i] % 24;
	long a = x[i] / 24;/*newc=24a+b. b gives the block to insert into, and we need to save "a"*/
	long v = a % 64;
	long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
	rclass[b][u] |= bitswap[v];
  }
  /*We adjust the starting quadruple in case of symmetries.*/
  if (x[0] == 0) {/*We are the [0, 0, 1, 1] packing: do [0, 1, 1, 0] and start with swapping the last.*/
	depthseq[0][0] = 0; depthseq[0][1] = 1; depthseq[0][2] = 1; depthseq[0][3] = 0;
	swaps[1] = 2;/*Will get incremented to 3 right away.*/
  }
  else if (x[0] == -1) {/*[-1, 2, 2, 3]: do [2, 3, -1, 2]*/
    depthseq[0][0] = 2; depthseq[0][1] = 3; depthseq[0][2] = -1; depthseq[0][3] = 2;
	swaps[1] = 1;  
  }
  else if (x[1] == x[2]) {/*b=c, so do [b, a, b, d]*/
	depthseq[0][0] = x[1]; depthseq[0][1] = x[0]; depthseq[0][2] = x[2]; depthseq[0][3] = x[3];
	swaps[1] = 0;
  }
  else if (x[2] == x[3]) {/*c=d, so do [c, a, b, c]*/
	depthseq[0][0] = x[3]; depthseq[0][1] = x[0]; depthseq[0][2] = x[1]; depthseq[0][3] = x[2];
	swaps[1] = 0;
  }
  else if ((x[0] + x[1] + x[2]) == x[3]) {/*a+b+c=d, so do [d, a, b, c]*/
	depthseq[0][0] = x[3]; depthseq[0][1] = x[0]; depthseq[0][2] = x[1]; depthseq[0][3] = x[2];
	swaps[1] = 0;
  }
  else {
	depthseq[0][0] = x[0]; depthseq[0][1] = x[1]; depthseq[0][2] = x[2]; depthseq[0][3] = x[3];
  }
  long ind = 1;/*Which depth we are working at.*/
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
	if (newc > B) continue;/*Too big! go back.*/
	/*Do the bitswap to update the count.*/
	long b = newc % 24;
	long a = newc / 24;/*newc=24a+b. b gives the block to insert into, and we need to save "a"*/
	long v = a % 64;
	long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
	rclass[b][u] |= bitswap[v];
	for (i = 0; i < cind; i++) depthseq[ind][i] = depthseq[lastind][i];
	depthseq[ind][cind] = newc;
	for (i = cind + 1; i < 4; i++) depthseq[ind][i] = depthseq[lastind][i];/*Add the tuple in.*/
	ind++;
	if (ind == maxdepth) {/*We are going too deep, must reallocate the storage location.*/
	  long newdepth = maxdepth << 1;/*Double it.*/
	  depthseq = realloc(depthseq, newdepth * sizeof(long *));
	  if (!depthseq) {
	    printf("Insufficient memory to reallocate the depth sequence.\n");
	    exit(1);
      }
	  swaps = realloc(swaps, newdepth * sizeof(int));
	  if (!swaps) {
	    printf("Insufficient memory to reallocate the swaps.\n");
	    exit(1);
      }
	  for (i = maxdepth; i < newdepth; i++) {
		depthseq[i] = (long *)malloc(sizeof(long) << 2);
		swaps[i] = -1;
	  }
	  maxdepth = newdepth;
	}
  }
  unsigned long **fams;
  if (families)	fams = findfamilies(B, rclass, res, lenres, bitswap);/*Find and remove the families.*/
  char fname[100];
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
  pos += sprintf(&fname[pos], "%ld_%ld_%ld_%ld", x[1], x[2], x[3], B);
  if (families) pos += sprintf(&fname[pos], "_fam");
  pos += sprintf(&fname[pos], ".dat");
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
	  for (j = 1; j <= fams[i][0]; j++) {
		if (j > 1) fprintf(F, ", ");
		fprintf(F, "%ld", fams[i][j]);
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
		  long n = 24 * a + b;/*The correct one!*/
		  if (n <= B && n) {/*Small enough!*/
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
  /*Time to free all of the allocated memory.*/
  if (families) {
	for (i = 0; i < lenres; i++) free(fams[i]);
	free(fams);
  }
  free(swaps);
  for (i = 0; i < maxdepth; i++) free(depthseq[i]);
  free(depthseq);
  for (i = 0; i < lenres; i++) free(rclass[res[i]]);
  free(rclass);
  free(bitswap);
  return;
}

/*Returns the families found, and updates rclass to have them removed too.*/
unsigned long **
findfamilies(long B, unsigned long **rclass, long res[], long lenres, unsigned long *bitswap)
{
  unsigned long **fams = (unsigned long **)malloc(lenres * sizeof(unsigned long *));/*Pointers to the families that exist. The first element of each family is the number of families in that class.*/
  long i, fampos;
  for (i = 0; i < lenres; i++) {
	switch (res[i]) {
	  case 0:
	    fams[i] = (unsigned long *)malloc(5 * sizeof(unsigned long));/*Up to 4 families.*/
		fampos = 1;
		if (removefamily(B, rclass[res[i]], bitswap, 24, 1)) {/*We have one family!*/
		  fams[i][fampos] = 24;
		  fampos++;
		}
		if (removefamily(B, rclass[res[i]], bitswap, 48, 1)) {
		  fams[i][fampos] = 48;
		  fampos++;
		}
		if (removefamily(B, rclass[res[i]], bitswap, 72, 1)) {
		  fams[i][fampos] = 72;
		  fampos++;
		}
		if (removefamily(B, rclass[res[i]], bitswap, 144, 1)) {
		  fams[i][fampos] = 144;
		  fampos++;
		}
		fams[i][0] = fampos - 1;/*The size.*/
		continue;
	  case 2:
	    fams[i] = (unsigned long *)malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
		if (removefamily(B, rclass[res[i]], bitswap, 2, 6)) {
		  fams[i][1] = 2;
		  fams[i][0] = 1;
		}
		else fams[i][0] = 0;
		continue;
	  case 4:
	    fams[i] = (unsigned long *)malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
		if (removefamily(B, rclass[res[i]], bitswap, 4, 6)) {
		  fams[i][1] = 4;
		  fams[i][0] = 1;
		}
		else fams[i][0] = 0;
		continue;
	  case 6:
	    fams[i] = (unsigned long *)malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
		if (removefamily(B, rclass[res[i]], bitswap, 6, 2)) {
		  fams[i][1] = 6;
		  fams[i][0] = 1;
		}
		else fams[i][0] = 0;
		continue;
	  case 8:
	    fams[i] = (unsigned long *)malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
		if (removefamily(B, rclass[res[i]], bitswap, 8, 3)) {
		  fams[i][1] = 8;
		  fams[i][0] = 1;
		}
		else fams[i][0] = 0;
		continue;
	  case 12:
	    fams[i] = (unsigned long *)malloc(3 * sizeof(unsigned long));/*Up to 2 families.*/
		fampos = 1;
		if (removefamily(B, rclass[res[i]], bitswap, 12, 2)) {
		  fams[i][fampos] = 12;
		  fampos++;
		}
		if (removefamily(B, rclass[res[i]], bitswap, 36, 2)) {
		  fams[i][fampos] = 36;
		  fampos++;
		}
		fams[i][0] = fampos - 1;
		continue;
	  case 16:
	    fams[i] = (unsigned long *)malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
		if (removefamily(B, rclass[res[i]], bitswap, 16, 3)) {
		  fams[i][1] = 16;
		  fams[i][0] = 1;
		}
		else fams[i][0] = 0;
		continue;
	  case 18:
	    fams[i] = (unsigned long *)malloc(sizeof(unsigned long) << 1);/*Up to 1 family.*/
		if (removefamily(B, rclass[res[i]], bitswap, 18, 2)) {
		  fams[i][1] = 18;
		  fams[i][0] = 1;
		}
		else fams[i][0] = 0;
		continue;
	  default:
	    fams[i] = (unsigned long *)malloc(sizeof(unsigned long));/*No family*/
		fams[i][0] = 0;
	}
  }
  return fams;
}

/*Looks through the family n=c*x^2, x coprime to cop. If they truly are all missing, we return 1, and flip all the bits in curvs.*/
int
removefamily(long B, long *curvs, unsigned long *bitswap, long c, long cop)
{
  long x = 1, n = c;
  while (n <= B) {
	long a = n / 24;/*n = 24a + residue*/
	long v = a % 64;
	long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
	if ((curvs[u] | bitswap[v]) == curvs[u]) return 0;/*This value was found, family not missing.*/
	do {x++;} while (!iscop(x, cop));
	n = c * x * x;
  }
  /*If we make it here, the family exists! We now must redo the above and actually do the swaps.*/
  x = 1;
  n = c;
  while (n <= B) {
	long a = n / 24;/*n = 24a + residue*/
	long v = a % 64;
	long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
	curvs[u] |= bitswap[v];
	do {x++;} while (!iscop(x, cop));
	n = c * x * x;
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

