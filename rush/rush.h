/***** recTest.h **********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 30 17:26:04 2012
 **************************************************/
#ifndef RUSH
#define RUSH

void initializeProb(int max, int sbjctLen, double gc);
int isHomologous(int *sl, int n, int sbjctLen);
double expVarVarSl(double l, double p);
double significanceVar(double meanSl, double varSl, int winLen);
double eVar(double meanSl, double winLen);

#endif

