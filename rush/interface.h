#ifndef INTERFACE
#define INTERFACE

#include <float.h>

#define DEFAULT_D 1000000

/* define argument container */
typedef struct args{
  double p; /* pi */
  float t;  /* threshold */
  float g;  /* deviation in GC-content */
  int m;    /* minimum shustring length */
  char *q;  /* query sequence */
  int D;    /* maximum depth of suffix tree */
  int w;    /* window length */
  int s;    /* step length */
  char d;   /* print debug messages */
  char h;   /* help message? */
  char v;   /* version message? */
  char e;   /* error message? */
  char **inputFiles;
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);

#endif
