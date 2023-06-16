/*
Compile with: 	gcc missing_curvatures.c -o missing_curvatures
Run with:		./missing_curvatures bound a b c d r1 ... rn, 
	where the reduced ACP in nondecreasing order is [a, b, c, d], we go up to bound, and r1, ..., rn are the possible residues modulo 24.
In order to do larger computations, we need a 64-bit operating system, so will assume that we have it.
*/

#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

void findmissing(long B, long x[], long res[], long lenres);

int
main(int argc, char *argv[])
{
  if (argc < 12) return 1;/*Problem.*/
  long B = atol(argv[1]);/*Bound.*/
  long x[4], i;/*The circle packing.*/
  for (i = 2; i <= 5; i++) x[i - 2] = atol(argv[i]);
  long res[argc - 6];/*The residue classes.*/
  for (i = 6; i < argc; i++) res[i - 6] = atol(argv[i]);
  clock_t tstart = clock();
  findmissing(B, x, res, argc - 6);
  clock_t tend = clock();
  double time_spent = (double)(tend - tstart) / CLOCKS_PER_SEC;
  printf("Total computation time: %fs\n", time_spent);
  return 0;
}

void
findmissing(long B, long x[], long res[], long lenres)
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
  for (i = 0; i < maxdepth; i++) swaps[i] = -1;/*Initialize to all -1's*/
  for (i = 1; i < 4; i++) {/*Do the first 3 curvatures (ignore the negative one).*/
	long b = x[i] % 24;
	long a = x[i] / 24;/*newc=24a+b. b gives the block to insert into, and we need to save "a"*/
	long v = a % 64;
	long u = a / 64;/*a=64u+v. u gives the entry of the array, v gives the bit to swap.*/
	rclass[b][u] |= bitswap[v];
  }
  long *xstart = (long *)malloc(4 * sizeof(long));/*We adjust the starting quadruple in case of symmetries.*/
  if (x[0] == 0) {/*We are the [0, 0, 1, 1] packing: do [0, 1, 1, 0] and start with swapping the last.*/
	xstart[0] = 0; xstart[1] = 1; xstart[2] = 1; xstart[3] = 0;
	swaps[1] = 2;/*Will get incremented to 3 right away.*/
  }
  else if (x[0] == -1) {/*[-1, 2, 2, 3]: do [2, 3, -1, 2]*/
    xstart[0] = 2; xstart[1] = 3; xstart[2] = -1; xstart[3] = 2;
	swaps[1] = 1;  
  }
  else if (x[1] == x[2]) {/*b=c, so do [b, a, b, d]*/
	xstart[0] = x[1]; xstart[1] = x[0]; xstart[2] = x[2]; xstart[3] = x[3];
	swaps[1] = 0;
  }
  else if (x[2] == x[3]) {/*c=d, so do [c, a, b, c]*/
	xstart[0] = x[3]; xstart[1] = x[0]; xstart[2] = x[1]; xstart[3] = x[2];
	swaps[1] = 0;
  }
  else if ((x[0] + x[1] + x[2]) == x[3]) {/*a+b+c=d, so do [d, a, b, c]*/
	xstart[0] = x[3]; xstart[1] = x[0]; xstart[2] = x[1]; xstart[3] = x[2];
	swaps[1] = 0;
  }
  else {
	xstart[0] = x[0]; xstart[1] = x[1]; xstart[2] = x[2]; xstart[3] = x[3];
  }
  depthseq[0] = xstart;/*Starting tuple, will always remain.*/
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
	long *newx = (long*) malloc(4 * sizeof(long));
	for (i = 0; i < cind; i++) newx[i] = depthseq[lastind][i];
	newx[cind] = newc;
	for (i = cind + 1; i < 4; i++) newx[i] = depthseq[lastind][i];
	depthseq[ind] = newx;/*Add the tuple in.*/
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
	  for (i = maxdepth; i < newdepth; i++) swaps[i] = -1;
	  maxdepth = newdepth;
	}
  }
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
  pos += sprintf(&fname[pos], "%ld_%ld_%ld_%ld.dat", x[1], x[2], x[3], B);
  FILE *F;
  F = fopen(fname, "w");
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
		  if (n <= B) {/*Small enough!*/
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

