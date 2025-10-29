/***** lcpTree.h *************************************************************
 * Description: Header file for lcp-interval tree processing.
 * Reference: Abouelhoda, M. I., Kurtz, S., and Ohlebusch, E. (2002).
 *   The enhanced suffix array and its applications to genome analysis.
 *   Proceedings of the Second Workshop on Algorithms in Bioinformatics,
 *   Springer-Verlag, Lecture Notes in Computer Science.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Mon Aug  6 15:37:55 2007.
 *****************************************************************************/
#ifndef LCPTREE
#define LCPTREE

typedef struct interval{
  int lcp;                        /* longest common prefix of members of this interval */
  int lb;                         /* left border */
  int rb;                         /* right border */
  struct interval **children;     /* children of this lcp-interval */
  int numChildren;                /* number of children */
  unsigned int isQuery : 1;       /* query flag */
  unsigned int isSubject : 1;     /* subject flag */
  unsigned int isNull : 1;
  unsigned int isFree : 1;
  int numQueryLeaves;
  int maxNumChildren;
  int maxNumLeaves;
  int *queryLeaves;
  int id;
}Interval;

int *getLcpTreeShulens(Args *args, Sequence *sequence);

void process(Args *args, Interval *interval, int *shulens);

int *traverseLcpTree(Args *args, int *lcpTab, int *sa, Sequence *seq);

Interval *getInterval(int lcp, int lb, int rb, Interval *interval);

void addChild(Interval *parent, Interval *child);
void checkLeaves(Interval *interval);
void freeInterval(Interval *interval);

#endif
