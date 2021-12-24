# HLA gene graphs built with pggb

Uniformly built using `pggb -s 1000 -p 70 -n $(cat $f.fai | wc -l) -G 2000,2000,2000,2000 -P 1,7,11,2,33,1 -k 19`

```
ls ~/HLA-zoo/seqs/*.fa | while read f; do echo $f; time pggb -i $f -t 16 -s 1000 -p 70 -n $(cat $f.fai | wc -l) -G 2000,2000,2000,2000 -P 1,7,11,2,33,1 -k 19 -o $(basename $f .fa) -Z; done
rm */*.smooth.{1,2,3}.gfa.gz
```

The graphs and other redundant files are gzipped.
