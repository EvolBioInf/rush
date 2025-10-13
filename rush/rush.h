/***** recTest.h **********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 30 17:26:04 2012
 **************************************************/
#ifndef RUSH
#define RUSH

#include "interface.h"
#include "minimize.h"

void initializeProb(int max, int sbjctLen, double gc);
int isHomologous(Int64 *sl, Int64 n, int sbjctLen);
double expVarVarSl(double l, double p);
double significanceVar(double meanSl, double varSl, int winLen);
double eVar(double meanSl, double winLen);

#endif

