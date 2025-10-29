/***** shulen.h **************************************************************
 * Description: Extract shustring lengths from suffix array combined
 *   with longest common prefix array. Use in conjunction with
 *   interface.h and stringData.h
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Wed Dec  6 17:19:12 2006.
 *****************************************************************************/  
#ifndef SHULEN
#define SHULEN

int *getSa(Sequence *seq);
int *getLcp(Sequence *seq, int *sa);

#endif
