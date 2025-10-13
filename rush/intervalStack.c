/***** intervalStack.c *******************************************************
 * Description: Functions from manipulating the interval stack.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Mon Aug  6 15:47:30 2007.
 *****************************************************************************/
#include <stdlib.h>
#include "eprintf.h"
#include "interface.h"
#include "sequenceData.h"
#include "deepShallow64/common.h"
#include "lcpTree.h"
#include "intervalStack.h"

Stack *createStack(){
  Stack *stack;

  stack = (Stack *)emalloc(sizeof(Stack));
  stack->curSize = 0;
  stack->intervals = (Interval **)emalloc(sizeof(Interval *));
  stack->maxSize = 1;
  stack->top = NULL;
  return stack;
}

/* freeStack: free memory occupied by a stack */
void freeStack(Stack *stack){
  Interval *interval;

  while(!isEmpty(stack)){
    interval = pop(stack);
    freeInterval(interval);
  }

  free(stack->intervals);
  free(stack);
}

void push(Stack *stack, Interval *interval){
  if(stack->curSize == stack->maxSize){
    stack->maxSize *= 2;
    stack->intervals = (Interval **)erealloc(stack->intervals,sizeof(Interval *)*stack->maxSize);
  }
  stack->intervals[stack->curSize] = interval;
  stack->top = stack->intervals[stack->curSize];
  stack->curSize++;
}

Interval *pop(Stack *stack){ 
  Interval *interval;

  interval = stack->intervals[stack->curSize-1];
  stack->curSize--;
  if(stack->curSize > 0)
    stack->top = stack->intervals[stack->curSize-1];
  else
    stack->top = NULL;
  return interval;
}

int isEmpty(Stack *stack){
  if(stack->curSize > 0)
    return 0;
  else
    return 1;
}
