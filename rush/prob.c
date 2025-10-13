/***** prob.c *************************************
 * Description: Probability distribution of
 *   shustring lengths. Derived by Peter Pfaffel-
 *   huber as described in his memo dated
 *   August 7, 2011.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Jun 23 13:10:12 2011
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include "deepShallow64/common.h"
#include "eprintf.h"

double *pshulen(double theta, int max, int sbjctLen, double *res);
double *dshulen(double theta, int max, int sbjctLen, double *res);

double *powerX = NULL;
double *powerSbjctLen = NULL;
double *globalW = NULL;
int maxX;

void initializeProb(int max, int sbjctLen, double gc){
  int x, k;
  double p, y, logp, logp2, log_2, s, intermed;

  p = gc/2;
  logp = log(p);
  logp2 = log(0.5-p);
  log_2 = log(2);
  globalW = (double *)emalloc((max+1)*sizeof(double));
  for(x=0;x<=max;x++){
    s = 0.;
    y = x*log_2;
    for(k=0;k<=x;k++){
      intermed = y +  gsl_sf_lnchoose(x,k) + k*logp + (x-k)*logp2 +	\
	sbjctLen*log(1-exp(k*logp+(x-k)*logp2));
      s += exp(intermed);
    }
    globalW[x] = s;
    if(globalW[x] >= 1.0){
      maxX = x;
      break;
    }
  }
}

/* pshulen: Equation (2) in  version of MS dated December 21, 2011 */
double *pshulen(double pi, int max, int sbjctLen, double *res){
  int i;
  double pii, piii;

  piii = 0;
  for(i=1;i<=maxX;i++){
    pii = pi * i;
    res[i] = globalW[i] * pii / (1+pii) - 
      globalW[i-1]*piii/(1+piii);
    piii = pii;
  }
  for(i=maxX+1;i<=max;i++){
    pii = pi * i;
    res[i] = pii/(1+pii) - piii/(1+piii);
    piii = pii;
  }
  return res;
}

double *dshulen(double theta, int max, int sbjctLen, double *res){
  int i;
  double *p, sum;
    
  p = pshulen(theta, max, sbjctLen, res);
  sum = 0.;
  for(i=1;i<=max;i++)
    sum += p[i];
  for(i=1;i<=max;i++)
    p[i] /= sum;

  return p;
}

double *shulenDist(Int64 *sl, Int64 n, int min, int max){
  double *sd;
  Int64 i, c;

  sd = (double *)emalloc(sizeof(double)*(max+1));
  /* initialize frequency array */
  for(i=0;i<=max;i++)
    sd[i] = 0.;
  /* count shustrings */
  c = 0;
  for(i=0;i<n;i++)
    if(sl[i] >= min){
      c++;
      sd[sl[i]]++;
    }

  for(i=0;i<=max;i++)
    sd[i] /= c;
  
  return sd;
}

