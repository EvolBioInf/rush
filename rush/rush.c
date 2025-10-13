/***** rush.c *************************************
 * Description: Recombination detection Using
 *   SHustrings
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Apr  4 17:17:34 2013
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include "deepShallow64/common.h"
#include "interface.h"
#include "eprintf.h"
#include "sequenceData.h"
#include "lcpTree.h"
/* #include "minimize.h" */
#include "rush.h"
#include "prob.h"

void scanFile(int fd, Args *args);
double pid(Int64 *sl, Int64 len, double gc, int min);
int  minShulen(int sbjctLen, double gc, double threshold);
double absolute(double a);

int main(int argc, char *argv[]){
  int i, sbjctDscr;
  char *version;
  Args *args;

  version = "1.4";
  setprogname2("rush");
  args = getArgs(argc, argv);
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  gsl_set_error_handler_off();
  if(args->numInputFiles == 0){
    sbjctDscr= 0;
    scanFile(sbjctDscr, args);
  }else{
    for(i=0;i<args->numInputFiles;i++){
      sbjctDscr = open(args->inputFiles[i],0);
      scanFile(sbjctDscr, args);
      close(sbjctDscr);
    }
  }
  free(args);
  free(progname());
  return 0;
}

void scanFile(int sbjctDscr, Args *args){
  Int64 *sl, i, len;
  double s, meanSl, sx, sig, ev, varSl, q, dr;
  Sequence *query, *sbjct, *seq;
  int queryDscr, r, l, winLen;

  queryDscr = open(args->q,0);
  if(queryDscr < 0)
    eprintf("ERROR: could not open query file %s\n",args->q);
  query = readFasta(queryDscr);
  sbjct = readFasta(sbjctDscr);
  close(sbjctDscr);
  prepareSeq(query);
  prepareSeq(sbjct);
  len = query->len/2 - 1;
  seq = catSeq(query,sbjct);
  freeSequence(query);
  freeSequence(sbjct);
  sl = getLcpTreeShulens(args, seq);
  s = 0.0;
  sx = 0.0;
  l = 0;
  r = 0;
  if(!args->w){
    winLen = len;
  }else{
    winLen = args->w;
  }
  /* fill first window */
  for(r=0;r<winLen;r++){
    s += sl[r];
    sx += sl[r]*sl[r];
  }
  meanSl = s/winLen;
  varSl = (sx-s*s/winLen)/(winLen-1);
  ev = eVar(meanSl,winLen);
  sig = significanceVar(meanSl, varSl, winLen);
  dr = (varSl - ev)/sqrt(24.0*pow(meanSl,5)/winLen);
  q = varSl/ev;
  if(args->w)
    printf("Pos\t%d\tQ\t%.3e\tD_r\t%.3e",(int)((r+l)/2.),q,dr);
  else
    printf("Q\t%.3e\tD_r\t%.3e",q,dr);
  if(sig != -1){
    if(dr > 0.0)
      printf("\tP\t%.3e\n",sig);
    else
      if(sig > 0.05)
	printf("\tP\t%.3e\n",sig);
      else
	printf("\tP\tfailed\n");
  }else
    printf("\tP\tfailed\n");
  /* scan remaining windows */
  while(r<len-args->s){
    for(i=0;i<args->s;i++){
      s += sl[r];
      s -= sl[l];
      sx += sl[r]*sl[r];
      sx -= sl[l]*sl[l];
      l++;
      r++;
    }
    meanSl = s/winLen;
    varSl = (sx-s*s/winLen)/(winLen-1);
    ev = eVar(meanSl,winLen);
    sig = significanceVar(meanSl, varSl, winLen);
    dr = (varSl - ev)/sqrt(24.0*pow(meanSl,5)/winLen);
    q = varSl/ev;
    printf("Pos\t%d\tQ\t%.3e\tD_r\t%3.e",(int)((r+l)/2.),q,dr);
    if(sig != -1){
      if(dr > 0.0)
	printf("\tP\t%.3e\n",sig);
      else
	if(sig > 0.05)
	  printf("\tP\t%.3e\n",sig);
	else
	  printf("\tP\tfailed\n");
    }else
      printf("\tP\tfailed\n");
  }
  freeSequence(seq);
  free(sl);
}

