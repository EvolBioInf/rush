/***** shulen.c **************************************************************
 * Description: Compute the lengths of shortest unique substrings.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Wed Dec  6 16:12:32 2006.
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "sequenceData.h"
#include "eprintf.h"
#include "interface.h"
#include "DeepShallow64/common.h"
#include "DeepShallow64/ds_ssort.h"
#include "DeepShallow64/bwt_aux.h"
#include "DeepShallow64/lcp_aux.h"
#include "shulen.h"

Int64 *getSuffixArray(Sequence *seq){
  Int64 overshoot;
  Int64 *sa;
  char *textu;

  /* init ds suffix sort routine (cf. DeepShallow64/testlcp.c) */
  overshoot = init_ds_ssort(500,2000);
  if(overshoot==0)
    eprintf("ERROR [ir]: ds initialization failed.\n");
  sa = (Int64 *)emalloc((seq->len+1)*sizeof(Int64));
  seq->seq = (char *)erealloc(seq->seq,(seq->len+overshoot)*sizeof(char));
  textu = (char *)seq->seq;
  ds_ssort((UChar *)textu, sa+1, seq->len);
  return sa;
}

Int64 *getLcp(Sequence *seq, Int64 *sa){
  Int64 occ[ALPHA_SIZE];
  UChar *textu;
  Int64 i;
  
  textu = (UChar *)seq->seq;
  for(i=0;i<ALPHA_SIZE;i++)
    occ[i] = 0;
  for(i=0;i<seq->len;i++)
    occ[textu[i]]++;
  
  return _lcp_sa2lcp_9n(textu,seq->len,sa,occ);
}

/* getShulensWithoutSentinel: return array of true shustring lengths */
Int64 *getShulensWithoutSentinel(Args *args, Sequence *seq){
  Int64 i, j, min; 
  Int64 *sl;
  Int64 *sa, *lcp;

  sa = getSuffixArray(seq);
  lcp = getLcp(seq,sa);

  sl = getRawShulens(args, seq, sa, lcp);
  free(sa);
  free(lcp);
  /* filter out shulens that cross string borders */
  min = 0;
  for(i=0;i<seq->numSeq;i++){
    if(i)
      min = seq->borders[i-1];
    for(j=seq->borders[i]-1;j>min;j--)
      if(j + sl[j] > seq->borders[i])
	sl[j] = -1;
      else
	break;
  }
  return sl;
}

Int64 *getRawShulens(Args *args, Sequence *seq, Int64 *sa, Int64 *lcp){
  Int64 i, arrayLen;
  Int64 *sl;
  
  arrayLen = seq->len;
  sl = (Int64 *)emalloc(seq->len*sizeof(Int64));
  for(i=2;i<arrayLen;i++)
    if(sa[i]<seq->len){
      if(lcp[i+1] < lcp[i])
	sl[sa[i]] = lcp[i] + 1;
      else
	sl[sa[i]] = lcp[i+1] + 1;
    }
  if(sa[i] < seq->len)
    sl[sa[i]] = lcp[i] + 1;
  if(sa[1] < seq->len)
    sl[sa[1]] = lcp[2] + 1;
  return sl;
}

/* getShulensWithSentinel: return shulens as if each sequence was terminated by a unique character. */
Int64 *getShulensWithSentinel(Args *args, Sequence *seq){
  Int64 i, j, min;
  Int64 *sl, *sa, *lcp;

  sa = getSuffixArray(seq);
  lcp = getLcp(seq,sa);
  sl = getRawShulens(args, seq, sa, lcp);
  free(sa);
  free(lcp);
  /* prune shulens that cross string borders */

  min = 0;
  for(i=0;i<seq->numSeq;i++){
    if(i)
      min = seq->borders[i-1];
    for(j=seq->borders[i]-1;j>min;j--)
      if(j + sl[j] > seq->borders[i])
	sl[j] = seq->borders[i] - j + 1;
      else
	break;
  }

  return sl;
}

void printSuffixArray(Args *args, Int64 *sa, Int64 *lcpTab, Int64 n, UChar *seq){
  Int64 i, j;

  for(i=0;i<n;i++){
    printf("%d\t%d\t%d\t",(int)i,(int)sa[i],(int)lcpTab[i]);
    for(j=sa[i];j<n;j++)
      printf("%c",seq[j]);
    printf("\n");
  }

}
