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
  Int64 lcp;                        /* longest common prefix of members of this interval */
  Int64 lb;                         /* left border */
  Int64 rb;                         /* right border */
  struct interval **children;     /* children of this lcp-interval */
  int numChildren;                /* number of children */
  unsigned int isQuery : 1;       /* query flag */
  unsigned int isSubject : 1;     /* subject flag */
  unsigned int isNull : 1;
  unsigned int isFree : 1;
  Int64 numQueryLeaves;
  Int64 maxNumChildren;
  Int64 maxNumLeaves;
  Int64 *queryLeaves;
  Int64 id;
}Interval;

Int64 *getLcpTreeShulens(Args *args, Sequence *sequence);

void process(Interval *interval, Int64 *shulens);

Int64 *traverseLcpTree(Int64 *lcpTab, Int64 *sa, Sequence *seq);

Interval *getInterval(Int64 lcp, Int64 lb, Int64 rb, Interval *interval);

void addChild(Interval *parent, Interval *child);
void checkLeaves(Interval *interval);
void freeInterval(Interval *interval);

#endif
