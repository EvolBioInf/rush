/***** intervalStack.h *******************************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Fri Feb  8 08:43:20 2008.
 *****************************************************************************/
#ifndef STACK
#define STACK

typedef struct stack{
  Interval **intervals;
  int curSize;
  int maxSize;
  Interval *top; 
}Stack;

Interval *pop(Stack *stack);
void push(Stack *stack, Interval *interval);
Stack *createStack(void);
int isEmpty(Stack *stack);
void freeStack(Stack *stack);

#endif
