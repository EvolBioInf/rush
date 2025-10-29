/***** varSd.c ************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Jan 10 12:39:59 2013
 **************************************************/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "rush.h"

double expVarVarSl(double l, double p){
  double vvs;
  double p5, p6, p7, p8;
  double l2, l3, l4;

  p5 = p*p*p*p*p;
  p6 = p5*p;
  p7 = p6*p;
  p8 = p7*p;

  l2 = l*l;
  l3 = l2*l;
  l4 = l3*l;

  vvs = 24./l/p5 - 264/l2/p6 + 1320/l3/p7 - 2778/l4/p8;

  return vvs;
}

double significanceVar(double meanSl, double varSl, int winLen){
  double p, x, sig, evvs;

  p = 1/meanSl;
  evvs = expVarVarSl(winLen,p);
  if(evvs > 0){
    x = (varSl-eVar(meanSl,winLen))/sqrt(evvs);
    sig = gsl_cdf_ugaussian_Q(x);
    if(sig > 1)
      return 1.0;
    else
      return sig;
  }else
    return 0.0;
}

double eVar(double meanSl, double l){
  double p, vs;
  
  p = 1./meanSl;
  vs = 1./p/p - 2./l/p/p/p + 2./l/l/p/p/p/p;

  return vs;
}
