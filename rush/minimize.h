/***** minimize.h *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 11 12:27:58 2011
 **************************************************/
#ifndef MINIMIZE
#define MINIMIZE

double *globalShulenDist;
double *globalResult;
int globalMaxLen;
int globalSbjctLen;

typedef struct result{
  double t; /* theta */
  double l; /* likelihood */
} Result;

Result *minimize(double *sd, Int64 mLen, Int64 sbjctLen);

double *getShulenDist();
void setShulenDist(double *shulenDist);
double *getResult();
void setResult(double *result);
int getMaxLen();
void setMaxLen(int maxLen);
int getSbjctLen();
void setSbjctLen(int sbjctLen);
double my_f(double t, void *params);
double likelihood(double t);
#endif
