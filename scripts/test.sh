seq=$(mktemp)
# The following command breaks down as
#   ms 2 1: run Hudson's ms to generate 1 set of 2 haplotypes
#   -s 1000: the sequences are separated by 1000 mutations
#   -r 100 100000: the sequences undergo an expected 100
#      recombination events and are 100000 bp long
ms 2 1 -s 1000 -r 100 100000 |
    ms2dna > $seq; # Convert the sequences to DNA
# Extract query & subject sequences
query=$(mktemp)
sbjct=$(mktemp)
getSeq S1 $seq > $query
getSeq S2 $seq > $sbjct
# Run rush
rush -q $query $sbjct
# Clean up
rm $seq $query $sbjct
