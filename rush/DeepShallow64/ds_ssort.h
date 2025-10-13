#ifndef DS_SORT_H
#define DS_SORT_H

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Prototypes for the Deep Shallow Suffix Sort routines
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */ 
/* void ds_ssort(unsigned char *t, int *sa, int n); */
void ds_ssort(UChar *t, Int64 *sa, Int64 n); /* 64 bit version */
/* int init_ds_ssort(int adist, int bs_ratio); */
Int64 init_ds_ssort(Int64 adist, Int64 bs_ratio);

#endif /* DS_SORT_H */

