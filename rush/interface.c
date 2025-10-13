/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "interface.h"
#include "eprintf.h"

Args *args;

Args *getArgs(int argc, char *argv[]){
  char c;
  char *optString = "hvq:w:s:";

  args = (Args *)emalloc(sizeof(Args));
  args->w = 0;
  args->s = 0;
  args->D = DEFAULT_D;
  args->d = 0;
  args->q = NULL;
  args->h = 0;
  args->v = 0;
  args->e = 1;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'w':                           /* window length */
      args->w = atoi(optarg);
      break;
    case 's':                           /* step length */
      args->s = atoi(optarg);
      break;
    case 'q':                           /* query file */
      args->q = optarg;
      args->e = 0;
      break;
    case '?':                           /* fall-through is intentional */
    case 'h':                           /* print help */
      args->h = 1;
      break;
    case 'v':                           /* print version */
      args->v = 1;
      break;
    default:
      printf("# unknown argument: %c\n",c);
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  args->s = (int)(args->w/10+0.5);
  if(!args->s)
    args->s = 1;
  if(!args->q && !(args->e || args->h || args->v))
    printf("ERROR: Please provide name of query file.\n");
  return args;
}


void printUsage(char *version){
  printf("Usage: %s [options] [inputFiles]\n",progname());
  printf("Recombination detection Using SHustrings\n");
  printf("rush -q query.fasta sbjct.fasta\n");
  printf("Options:\n");
  printf("\t-q <FILE> file containing query sequence\n");
  printf("\t[-w window length; default: whole sequence]\n");
  printf("\t[-s step length; default: window_len/10]\n");
  printf("\t[-h print this help message and exit]\n");
  printf("\t[-v print program information and exit]\n");
  exit(0);
}

void printSplash(char *version){
  printf("%s %s\n",progname(),version);
  printf("Written by Bernhard Haubold, Linda Krause, Thomas Horn, and Peter Pfaffelhuber.\n");
  printf("Distributed under the GNU General Public License.\n");
  printf("Please send bug reports to haubold@evolbio.mpg.de\n");
  exit(0);
}
