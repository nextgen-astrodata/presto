#include "presto.h"

/* Function declarations */

int compare_powindex(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int compare_positions(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int compare_fourierprops(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int comp_bin_pow(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int comp_bin_nfftbins(void *ca, void *cb);
   /*  Used as compare function for qsort() */
int remove_dupes(position * list, int nlist);
   /*  Removes list values that are 1 unit of search away from a higher */
   /*  power candidate (dr = 0.5 or dz = 2.0).  Returns # removed.      */
int remove_dupes2(fourierprops * list, int nlist);
   /*  Removes list values that are within measurement error away from  */
   /*  a higher power candidate.  Returns # removed.                    */
int remove_dupes_bin(binaryprops * list, int nlist);
   /*  Removes list values that are within 1 Fourier bin of the PSR freq */
   /*  from a higher power candidate. Returns # removed.                 */
int remove_other(fourierprops * list, int nlist, long rlo, \
		 long rhi, double locpow);
   /*  Removes list values whose frequencies fall outside rlo and rhi */
   /*  and candidates whose local power levels are below locpow.      */ 
   /*  Returns # removed.                                             */
int remove_other_bin(binaryprops * list, int nlist);
   /*  Removes list values whose binary parameters are unrealistic.   */
   /*  For example, orbital periods under 200 sec.                    */
   /*  Returns # removed.                                             */


/* int compare_positions(const void *ca, const void *cb) */
int compare_positions(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
  position *a, *b;

  a = (position *) ca;
  b = (position *) cb;
  if ((b->pow - a->pow) < 0.0)
    return -1;
  if ((b->pow - a->pow) > 0.0)
    return 1;
  return 0;
}

/* int compare_powindex(const void *ca, const void *cb) */
int compare_powindex(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
  powindex *a, *b;

  a = (powindex *) ca;
  b = (powindex *) cb;
  if ((b->pow - a->pow) < 0.0)
    return -1;
  if ((b->pow - a->pow) > 0.0)
    return 1;
  return 0;
}

/* int comp_bin_pow(const void *ca, const void *cb) */
int comp_bin_pow(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
  binaryprops *a, *b;

  a = (binaryprops *) ca;
  b = (binaryprops *) cb;
  if ((b->pow - a->pow) < 0.0)
    return -1;
  if ((b->pow - a->pow) > 0.0)
    return 1;
  return 0;
}

/* int comp_bin_nfftbins(const void *ca, const void *cb) */
int comp_bin_nfftbins(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
  binaryprops *a, *b;

  a = (binaryprops *) ca;
  b = (binaryprops *) cb;
  if ((b->nfftbins - a->nfftbins) < 0.0)
    return -1;
  if ((b->nfftbins - a->nfftbins) > 0.0)
    return 1;
  return 0;
}

/* int compare_fourierprops(const void *ca, const void *cb) */
int compare_fourierprops(void *ca, void *cb)
/*  Used as compare function for qsort() */
{
  fourierprops *a, *b;

  a = (fourierprops *) ca;
  b = (fourierprops *) cb;
  if ((b->pow - a->pow) < 0.0)
    return -1;
  if ((b->pow - a->pow) > 0.0)
    return 1;
  return 0;
}

float percolate(position * list, int nlist)
/*  Pushes a position structure as far up a sorted list of positions */
/*  as it needs to go to keep the list sorted.  Returns the new low  */
/*  power in the list.                                               */
{
  int ct;
  position tempzz;

  for (ct = nlist - 2; ct >= 0; ct--) {
    if (list[ct].pow < list[ct + 1].pow) {
      SWAP(list[ct], list[ct + 1]);
    } else {
      break;
    }
  }
  return list[nlist - 1].pow;
}

float percolate_bin(binaryprops * list, int nlist)
/*  Pushes a binaryprops structure as far up a sorted list of structs */
/*  as it needs to go to keep the list sorted.  Returns the new low  */
/*  power in the list.                                               */
{
  int ct;
  binaryprops tempzz;

  for (ct = nlist - 2; ct >= 0; ct--) {
    if (list[ct].pow < list[ct + 1].pow) {
      SWAP(list[ct], list[ct + 1]);
    } else {
      break;
    }
  }
  return list[nlist - 1].pow;
}


int remove_dupes(position * list, int nlist)
/*  Removes list values that are 1 unit of search away from a higher */
/*  power candidate (dr = 0.5 or dz = 2.0).  Returns # removed.      */
{
  int i, j, k, ct = 0;
  position tempzz;

  for (i = 0; i < nlist - 1; i++) {
    if (list[i].pow == 0.0)
      break;
    for (j = i + 1; j < nlist; j++) {
      if (list[j].pow == 0.0)
	break;
      if ((fabs(list[j].p1 - list[i].p1) < 0.51) &&
	  (fabs(list[j].p2 - list[i].p2) < 2.01)) {
	if (j < nlist - 1) {
	  for (k = j; k < nlist - 1; k++) {
	    SWAP(list[k], list[k + 1]);
	  }
	}
	k = nlist - 1;
	list[k].pow = 0.0;
	ct++;
	j--;
      }
    }
  }
  printf("Removed %d duplicates\n", ct);
  return ct;
}

int remove_dupes_bin(binaryprops * list, int nlist)
/*  Removes list values that are within 1 Fourier bin of the PSR freq */
/*  from a higher power candidate. Returns # removed.                 */
{
  int i, j, k, ct = 0;
  binaryprops tempzz;

  for (i = 0; i < nlist - 1; i++) {
    if (list[i].pow == 0.0)
      break;
    for (j = i + 1; j < nlist; j++) {
      if (list[j].pow == 0.0)
	break;
      if ((fabs(list[j].rdetect - list[i].rdetect) < 0.6) &&
	  (fabs(list[j].rpsr - list[i].rpsr) < list[j].nfftbins / 2) &&
	  (list[j].nfftbins == list[i].nfftbins)) {
	if (j < nlist - 1) {
	  for (k = j; k < nlist - 1; k++) {
	    SWAP(list[k], list[k + 1]);
	  }
	}
	k = nlist - 1;
	list[k].pow = 0.0;
	ct++;
	j--;
      }
    }
  }
  return ct;
}

int remove_dupes2(fourierprops * list, int nlist)
/*  Removes list values that are within measurement error away from  */
/*  a higher power candidate.  Returns # removed.                    */
{
  int i, j, k, ct = 0;
  fourierprops tempzz;

  for (i = 0; i < nlist - 1; i++) {
    if (list[i].pow == 0.0)
      break;
    for (j = i + 1; j < nlist; j++) {
      if (list[j].pow == 0.0)
	break;
      if ((fabs(list[j].r - list[i].r) < list[i].rerr) &&
	  (fabs(list[j].z - list[i].z) < list[i].zerr)) {
	if (j < nlist - 1) {
	  for (k = j; k < nlist - 1; k++) {
	    SWAP(list[k], list[k + 1]);
	  }
	}
	k = nlist - 1;
	list[k].pow = 0.0;
	ct++;
	j--;
      }
    }
  }
  printf("Removed %d fine duplicates\n", ct);
  return ct;
}


int remove_other(fourierprops * list, int nlist, long rlo, \
		 long rhi, double locpow)
{
  /*  Removes list values whose frequencies fall outside rlo and rhi */
  /*  and candidates whose local power levels are below locpow.      */
  /*  Returns # removed.                                             */
  int i, j, ct = 0;
  fourierprops tempzz;

  for (i = 0; i < nlist; i++) {
    if (list[i].pow == 0.0)
      break;
    if (list[i].r < rlo || \
	list[i].r > rhi || \
	list[i].pow < locpow) {
      if (i < nlist - 1) {
	for (j = i; j < nlist - 1; j++) {
	  SWAP(list[j], list[j + 1]);
	}
      }
      j = nlist - 1;
      list[j].pow = 0.0;
      ct++;
      i--;
    }
  }
  printf("Removed %d others\n", ct);
  return ct;
}


int remove_other_bin(binaryprops * list, int nlist)
  /*  Removes list values whose binary parameters are unrealistic.   */
  /*  For example, orbital periods under 300 sec.                    */
  /*  Returns # removed.                                             */
{
  int i, j, ct = 0;
  float cutoff = 300.0;
  binaryprops tempzz;

  for (i = 0; i < nlist; i++) {
    if (list[i].pow == 0.0)
      break;
    if (list[i].pbin < cutoff) {
      if (i < nlist - 1) {
	for (j = i; j < nlist - 1; j++) {
	  SWAP(list[j], list[j + 1]);
	}
      }
      j = nlist - 1;
      list[j].pow = 0.0;
      ct++;
      i--;
    }
  }
  return ct;
}