./rush -q ../data/testQ.fasta ../data/testS.fasta > tmp.txt
d=$(diff r.txt tmp.txt)
if [[ $d != "" ]]; then
    echo "failed"
else
    echo "passed"
fi
rm tmp.txt
