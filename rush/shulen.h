/***** shulen.h **************************************************************
 * Description: Extract shustring lengths from suffix array combined
 *   with longest common prefix array. Use in conjunction with
 *   interface.h and stringData.h
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Wed Dec  6 17:19:12 2006.
 *****************************************************************************/  
#ifndef SHULEN
#define SHULEN

#define ALPHA_SIZE _BW_ALPHA_SIZE

Int64 *getSuffixArray(Sequence *seq);
Int64 *getLcp(Sequence *seq, Int64 *sa);
Int64 *getShulensWithoutSentinel(Args *args, Sequence *seq);
Int64 *getRawShulens(Args *args, Sequence *seq, Int64 *sa, Int64 *lcp);
Int64 *getShulensWithSentinel(Args *args, Sequence *seq);
void printSuffixArray(Args *args, Int64 *sa, Int64 *lcp, Int64 n, UChar *seq);

#endif
