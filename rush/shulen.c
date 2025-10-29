/***** shulen.c **************************************************************
 * Description: Compute the lengths of shortest unique substrings.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Wed Dec  6 16:12:32 2006.
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <divsufsort.h>
#include "sequenceData.h"
#include "eprintf.h"
#include "interface.h"
#include "shulen.h"

int *getSa(Sequence *seq) {
  unsigned char *t; // The text
  int *sa;  // The sa
  int n;

  // Get text
  t = (unsigned char *)seq->seq;
  // Get length of text
  n = seq->len;
  // Calculate sa.
  sa = (int *)emalloc(n * sizeof(int));
  if (divsufsort(t, sa, n) != 0) {
    printf("ERROR[getSa]: suffix sorting failed.\n");
    exit(-1);
  }
  return sa;
}

/* getLcp calculates the lcp array using Kasais' algorithm. */
int *getLcp(Sequence *seq, int *sa){
  char *t = seq->seq;
  int n = seq->len;
  int *isa = (int *)malloc(n * sizeof(int));
  int *lcp = (int *)malloc(n * sizeof(int));
  for(int i = 0; i < n; i++)
    isa[sa[i]] = i;
  lcp[0] = -1;
  int l = 0;
  for(int i = 0; i < n; i++) {
    int j = isa[i];
    if(j == 0)
      continue;
    int k = sa[j-1];
    while(t[i+l] == t[k+l])
      l++;
    lcp[j] = l;
    l--;
    if(l < 0)
      l = 0;
  }
  return lcp;
}
