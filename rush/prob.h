/***** prob.h *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Nov 15 11:06:02 2012
 **************************************************/
#ifndef PROB
#define PROB

void initializeProb(int max, int sbjctLen, double gc);
double *shulenDist(int *sl, int n, int min, int max);
double *dshulen(double theta, int max, int sbjctLen, double *res);

#endif
