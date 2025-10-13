/***** minimize.c *********************************
 * Description: Find the value of theta where
 *   the difference between expected and observed
 *   shulength distributions is minimized.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 31 10:03:40 2011
 **************************************************/
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "DeepShallow64/common.h"
#include "minimize.h"
#include "eprintf.h"
#include "prob.h"

Result *minimize(double *sd, Int64 mLen, Int64 sbjctLen){
  double *r;
  Result *result;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;
  int iter = 0, max_iter = 100;
  double m = 0.001, a = 0., b;
  double fm;
  int status;

  setShulenDist(sd);
  setMaxLen(mLen);
  setSbjctLen(sbjctLen);
  r = (double *)emalloc(sizeof(double)*(sbjctLen+2));
  setResult(r);

  /* Initialize method and iterate */
  F.function = &my_f;
  F.params = 0;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc(T);

  likelihood(a);
  fm = likelihood(m);
  b = m;
  while(likelihood(b) <= fm){
    b += 0.001;
    if(b>1)
      return NULL;
  }
  gsl_min_fminimizer_set(s, &F, m, a, b);
  do{
    iter++;
    status = gsl_min_fminimizer_iterate(s);
    m = gsl_min_fminimizer_x_minimum(s);
    a = gsl_min_fminimizer_x_lower(s);
    b = gsl_min_fminimizer_x_upper(s);
    status = gsl_min_test_interval(a,b,0.000001,0.0);
  }while(status == GSL_CONTINUE && iter < max_iter);

  if(status == GSL_SUCCESS){
    result = (Result *)emalloc(sizeof(Result));
    result->t = (b + a) / 2.0;
    result->l = gsl_min_fminimizer_f_minimum(s);
  }else
    result = NULL;

  gsl_min_fminimizer_free(s);
  free(r);
  return result;
}

double my_f(double t, void *params){

  if(t < 0)
    return FLT_MAX;

  return likelihood(t);
}

double likelihood(double t){
  double *ed, s;
  double *sd, lik;
  int i, c, l;

  sd = getShulenDist();
  ed = dshulen(t, getMaxLen(), getSbjctLen(), getResult());
  s = 0.;
  c = 0;
  l = getMaxLen();
  for(i=1;i<=l;i++){
    if(ed[i] > 0.0){
      c++;
      s += sd[i] * log(ed[i]);
    }
  }
  if(c){
    lik = -s;
  }else{
    lik = FLT_MAX;
  }
  return lik;
}

double *getShulenDist(){
  return globalShulenDist;
}
void setShulenDist(double *sd){
  globalShulenDist = sd;
}
double *getResult(){
  return globalResult;
}
void setResult(double *result){
  globalResult = result;
} 
int getMaxLen(){
  return globalMaxLen;
}
void setMaxLen(int maxLen){
  globalMaxLen = maxLen;
}
int getSbjctLen(){
  return globalSbjctLen;
}
void setSbjctLen(int sbjctLen){
  globalSbjctLen = sbjctLen;
}
