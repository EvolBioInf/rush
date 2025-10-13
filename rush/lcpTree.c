/***** lcpTree.c *************************************************************
 * Description: Functions for lcp-interval tree processing.
 * Reference: Abouelhoda, M. I., Kurtz, S., and Ohlebusch, E. (2002).
 *   The enhanced suffix array and its applications to genome analysis.
 *   Proceedings of the Second Workshop on Algorithms in Bioinformatics,
 *   Springer-Verlag, Lecture Notes in Computer Science.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Mon Aug 6 15:39:58 2007.
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eprintf.h"
#include "interface.h"
#include "sequenceData.h"
#include "deepShallow64/common.h"
#include "deepShallow64/ds_ssort.h"
#include "deepShallow64/bwt_aux.h"
#include "deepShallow64/lcp_aux.h"
#include "shulen.h"
#include "lcpTree.h"
#include "intervalStack.h"

Int64 *suffixArray;
Int64 queryStart, queryEnd;
Int64 numInterval = 0;
Int64 nextId = 0;
Int64 maxDepth;

/* getLcpTreeShulens: compute shulens from lcp tree traversal.
 * This is the only entry point to the functions in this file.
 */
Int64 *getLcpTreeShulens(Args *args, Sequence *seq){

  Int64 *sa, *lcpTab, *sl;
  Int64 min, i, j;

  maxDepth = args->D;

  sa = getSuffixArray(seq);
  lcpTab = getLcp(seq, sa);
  sl = traverseLcpTree(args, lcpTab, sa, seq);

  /* filter out shulens that cross string borders */
  min = 0;
  for(i=0;i<seq->numSeq;i++){
    if(i)
      min = seq->borders[i-1] + 1;
    for(j=seq->borders[i]-1;j>=min;j--){
      if(j + sl[j] > seq->borders[i])
	sl[j] = seq->borders[i] - j + 1;
      else
	break;
    }
  }
  if(args->d){
    for(i=0;i<seq->len;i++)
      printf("sl[%d] = %d\n",(int)i,(int)sl[i]);
  }
  free(sa);
  free(lcpTab);
  return sl;
}

/* traverseLcpTree: bottom-up traversal of lcp-interval tree */
Int64 *traverseLcpTree(Args *args, Int64 *lcpTab, Int64 *sa, Sequence *seq){

  Interval *lastInterval, *interval;
  Int64 i, j, lb, rightEnd, lastIsNull, *shulens;
  Stack *treeStack, *reserveStack;

  sa++;
  lcpTab++;

  suffixArray = sa;
  queryStart = seq->queryStart;
  queryEnd = seq->queryEnd;
  rightEnd = seq->len-1;
  lcpTab[0] = 0;
  /* allocate space for shulens of forwad & reverse strand */
  shulens = (Int64 *)emalloc((2*seq->len+1)*sizeof(Int64));
  
  treeStack = createStack();
  reserveStack = createStack();
  lastIsNull = 1; 
  lastInterval = NULL;
  interval = NULL;
  push(treeStack,getInterval(0,0,rightEnd,NULL));
  for(i=1; i<seq->len; i++){
    lb = i - 1;
    while(lcpTab[i] < treeStack->top->lcp){
      treeStack->top->rb = i - 1;
      lastInterval = pop(treeStack);
      /* save child invervals of popped interval */
      for(j=0;j<lastInterval->numChildren;j++){
	lastInterval->children[j]->numChildren = 0;
	push(reserveStack,lastInterval->children[j]);
      }
      lastIsNull = 0;
      process(args, lastInterval, shulens);
      lb = lastInterval->lb;
      if(lcpTab[i] <= treeStack->top->lcp){
	addChild(treeStack->top, lastInterval);
	lastIsNull = 1;
      }
    }
    if(lcpTab[i] > treeStack->top->lcp){
      if(isEmpty(reserveStack)){
	interval = getInterval(lcpTab[i],lb,rightEnd,NULL);
      }else{
	interval = pop(reserveStack);
	interval->lcp = lcpTab[i];
	interval->lb = lb;
	interval->rb = rightEnd;
	interval->numChildren = 0;
      }
      if(!lastIsNull){
	addChild(interval,lastInterval);
	lastIsNull = 1;
      }
      push(treeStack,interval);
    }
  }
  while(!isEmpty(treeStack)){
    interval = pop(treeStack);
    process(args, interval, shulens);
    for(i=0;i<interval->numChildren;i++)
      freeInterval(interval->children[i]);
    freeInterval(interval);
  }

  while(!isEmpty(reserveStack)){
    interval = pop(reserveStack);
    for(i=0;i<interval->numChildren;i++)
      freeInterval(interval->children[i]);
    freeInterval(interval);
  }
  freeStack(treeStack);
  freeStack(reserveStack);
  return shulens;
}

Interval *getInterval(Int64 lcp, Int64 lb, Int64 rb, Interval *child){
  Interval *interval;
  Int64 i;

  if(++numInterval > maxDepth){
    printf("WARNING [getInterval]: program terminated with Warning Code 1;\n");
    printf("see documentation for further details\n");
    exit (0);
  }

  interval = (Interval *)emalloc(sizeof(Interval));
  interval->lcp = lcp;  /* longest common prefix */
  interval->lb = lb;    /* left border */
  interval->rb = rb;    /* right border */
  interval->id = ++nextId;
  interval->maxNumChildren = 0;
  interval->children = (Interval **)emalloc(interval->maxNumChildren*sizeof(Interval *));
  for(i=0;i<interval->maxNumChildren;i++)
    interval->children[i] = (Interval *)emalloc(sizeof(Interval *)); 
  if(child == NULL)
    interval->numChildren = 0;
  else{
    interval->children[0] = child;
    interval->numChildren = 1;
  }
  interval->numQueryLeaves = 0;
  interval->maxNumLeaves = 0;
  interval->queryLeaves = (Int64 *)emalloc(interval->maxNumLeaves*sizeof(Int64));
  interval->isNull = 0;
  interval->isFree = 0;
  return interval;
}

void addChild(Interval *parent, Interval *child){

  child->isFree = 0;
  if(parent->numChildren >= parent->maxNumChildren){
    parent->maxNumChildren++;
    parent->children = (Interval **)erealloc(parent->children,parent->maxNumChildren * sizeof(Interval *));
  }
  parent->children[parent->numChildren] = child;
  parent->numChildren++;
}

/* process: this function does two things:
 * 1) it labels an interval according to whether it has suffix tree leaves from query (isQuery)
 *    or from subject (isSubject) or both (isQuery && isSubject).
 * 2) if(isQuery && isSubject) it determines the corresponding query shustring lengths (if any)
 * Some rough notes on this procedure can be found in the black notebook, p. 79, August 8, 2007.
 */
void process(Args *args, Interval *interval, Int64 *shulens){
  Int64 i, len;
  Int64 ind;

  interval->isQuery = 0;
  interval->isSubject = 0;
  interval->numQueryLeaves = 0;
  if(interval->numChildren == 0){                       /* leaf of lcp tree */
    for(i=interval->lb;i<=interval->rb;i++){
      if(suffixArray[i] >= queryStart && suffixArray[i] <= queryEnd)
	interval->isQuery = 1;
      else 
	interval->isSubject = 1;
    } 
  }else{
    for(i=0;i<interval->numChildren;i++){                /* internal node of lcp tree */
      if(interval->children[i]->isQuery)
	interval->isQuery = 1;
      if(interval->children[i]->isSubject)
	interval->isSubject = 1;
    }
    checkLeaves(interval);
  }


  if(interval->isQuery && interval->isSubject){
    if(interval->numChildren == 0){
      for(i=interval->lb;i<=interval->rb;i++){
	if(suffixArray[i] >= queryStart && suffixArray[i] < queryEnd){
	  ind = suffixArray[i] - queryStart;
	  len = interval->lcp + 1;
	  shulens[ind] = ((ind + len - 1) <= queryEnd) ? len : queryEnd - ind + 1;
	}
      }
    }else if(interval->numQueryLeaves > 0){
      for(i=0;i<interval->numQueryLeaves;i++){
	ind = interval->queryLeaves[i]-queryStart;
	len = interval->lcp + 1;
	shulens[ind] = ((ind + len - 1) <= queryEnd) ? len : queryEnd - ind + 1;
      }
    }
  }
  if(args->d){
    if(interval->isQuery && interval->isSubject)
      printf("%d-[%d..%d]s\n",(int)interval->lcp,(int)interval->lb,(int)interval->rb);
    else if(interval->isQuery)
      printf("%d-[%d..%d]c\n",(int)interval->lcp,(int)interval->lb,(int)interval->rb);
    else if(interval->isSubject)
      printf("%d-[%d..%d]t\n",(int)interval->lcp,(int)interval->lb,(int)interval->rb);
  }
}

/* checkLeaves: note the suffix tree leaves from query attached to interval.
 */
void checkLeaves(Interval *interval){
  Int64 i, j;
  
  if(interval->numChildren < 1)
    return;

  /* do the child intervals cover the parent interval? */
  /* left border */
  if(!interval->children[0]->isSubject){
    for(i=interval->children[0]->lb;i<=interval->children[0]->rb;i++){
      if(interval->numQueryLeaves >= interval->maxNumLeaves){
	interval->maxNumLeaves++;
	interval->queryLeaves = (Int64 *)erealloc(interval->queryLeaves, interval->maxNumLeaves*sizeof(Int64));
      }
      interval->queryLeaves[interval->numQueryLeaves++] = suffixArray[i];
    }
  }
  if(interval->lb < interval->children[0]->lb){
    for(i=interval->lb;i<interval->children[0]->lb;i++)
      if(suffixArray[i] >= queryStart && suffixArray[i] < queryEnd){
	interval->isQuery = 1;
	if(interval->numQueryLeaves >= interval->maxNumLeaves){
	  interval->maxNumLeaves++;
	  interval->queryLeaves = (Int64 *)erealloc(interval->queryLeaves, interval->maxNumLeaves*sizeof(Int64));
	}
	interval->queryLeaves[interval->numQueryLeaves++] = suffixArray[i];
      }else
	interval->isSubject = 1;
  }
  /* intervals covered by subject children */
  for(i=1;i<interval->numChildren;i++){
    if(!interval->children[i]->isSubject){
      for(j=interval->children[i]->lb;j<=interval->children[i]->rb;j++){
	if(interval->numQueryLeaves >= interval->maxNumLeaves){
	  interval->maxNumLeaves++;
	  interval->queryLeaves = (Int64 *)erealloc(interval->queryLeaves, interval->maxNumLeaves*sizeof(Int64));
	}
	interval->queryLeaves[interval->numQueryLeaves++] = suffixArray[j];
      }
    }
    if(interval->children[i-1]->rb+1 < interval->children[i]->lb){
      for(j=interval->children[i-1]->rb+1;j<interval->children[i]->lb;j++){
	if(suffixArray[j] >= queryStart && suffixArray[j] < queryEnd){
	  interval->isQuery = 1;
	  if(interval->numQueryLeaves >= interval->maxNumLeaves){
	    interval->maxNumLeaves++;
	    interval->queryLeaves = (Int64 *)erealloc(interval->queryLeaves, interval->maxNumLeaves*sizeof(Int64));
	  }
	  interval->queryLeaves[interval->numQueryLeaves++] = suffixArray[j];
	}else
	  interval->isSubject = 1;
      }
    }
  }
  /* right border */
  if(interval->children[interval->numChildren-1]->rb < interval->rb){
    for(i=interval->children[interval->numChildren-1]->rb+1;i<=interval->rb;i++)
      if(suffixArray[i] >= queryStart && suffixArray[i] < queryEnd){
	interval->isQuery = 1;
	if(interval->numQueryLeaves >= interval->maxNumLeaves){
	  interval->maxNumLeaves++;
	  interval->queryLeaves = (Int64 *)erealloc(interval->queryLeaves, interval->maxNumLeaves*sizeof(Int64));
	}
	interval->queryLeaves[interval->numQueryLeaves++] = suffixArray[i];
      }else
	interval->isSubject = 1;
  }
}

void freeInterval(Interval *interval){
  free(interval->queryLeaves);
  free(interval->children);
  free(interval);
  numInterval--;
}

