# the HLA zoo

_HLA variation graphs_

This is a collection of different renderings of a set of variation graphs made from different versions of the GRCh38 alts.

These are meant for testing variation graph methods.

## input sequences

These are ALT sequences of HLA genes from GRCh38, collected in the Human Genome Variation / Map project.

```txt
gene            alts  size
A-3105          11    147718
B-3106          9     30751
C-3107          10    33810
DMA-3108        11    49624
DMB-3109        10    64515
DOA-3111        10    54307
DOB-3112        10    42837
DPA1-3113       11    178303
DPB1-3115       11    151390
DQA1-3117       10    73280
DQB1-3119       10    73913
DRA-3122        11    57216
DRB1-3123       12    163416
DRB3-3125       4     52271
DRB4-3126       5     64760
DRB5-3127       3     38568
E-3133          9     43200
F-3134          10    49951
G-3135          11    45653
H-3136          10    33111
J-3137          10    39844
K-3138          9     28660
L-3139          8     59036
MICA-100507436  8     116989
MICB-4277       11    168093
TAP1-6890       11    96391
TAP2-6891       11    185580
V-352962        10    9865
```

## output graphs

A variation graph represents a set of alignments between sequences and the sequences themselves in a graph.

A number of methods are applied to build variation graphs from these sequences.

### progressive partial order alignment with `spoa`

The classic partial order alignment algorithm, sped up with SIMD instructions and written into GFA format rather than MSA, with some light post-processing

```
ls seqs/*fa | sort | while read f
do
    echo $f
    in=$f
    out=graphs/spoa/$(basename $f .fa).gfa
    spoa -GR $in >$out.1 && odgi build -g $out.1 -o - | odgi unchop -i - -o - | odgi sort -i - -o - | odgi view -i - -g >$out
    rm -f $out.1
done
```

(Some of these require a lot of memory to compute, so not all were created.)

### progressive variation graph construction with `vg msga`

A progressive graph construction using the `vg map` algorithm, `vg msga` is similar to `spoa`, but it allows for structural variation in the graph by "chunking" the alignment problem, threading the results back together using dynamic programming, and further cleaning up the alignment result by locally aligning any remaining unaligned fragments of the sequence.

This heuristic approach works well enough to build these graphs, but it does not scale efficiently to large problems.

```
ls seqs/*fa | sort | while read f
do
    echo $f
    in=$f
    out=graphs/vg-msga/$(basename $f .fa).gfa
    vg msga -f $in | vg view - >$out
done
```

### `seqwish` variation graph induction

[`seqwish`](https://github.com/ekg/seqwish) reads a set of alignments and sequences and renders the variation graph that they imply.
It is the WYSIWYG of variation graph construction methods, in that the resulting graph perfectly reflects the input alignment set.
To change the graph, we simply change the input alignments.

`seqwish` depends on an indepednent alignment process to produce a set of alignments (in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md)).
Two alignment approaches have been validated.

#### `minimap2` → `seqwish`

[`minimap2`](https://github.com/lh3/minimap2) is a popular and efficient aligner for long sequences based on minimizer seeding/chaining and [adaptive banded alignment](https://github.com/ocxtal/libgaba) to derive base pair alignments.

```
ls seqs/*fa | sort | while read f
do
    echo $f
    in=$f
    out=graphs/seqwish/minimap2/$(basename $f .fa)
    minimap2 -c -x asm20 $in $in >$out.paf \
        && seqwish -s $in -p <(fpa drop -l 3000 <$out.paf) -g $out.gfa -t 16 \
        && odgi build -g $out.gfa -o - \
            | odgi sort -i - -o - -p sYgYs -k 1000 -G 1 -A -t 16 -P \
            | odgi view -i - -g >$out.sort.gfa
done
```

This builds the graph and also does a sort of it.
`seqwish` graphs are not sorted by default, and the nodes in them occur in the order of first appearance in the input sequence set.


#### `mashmap` → `seqwish`

[`mashmap`](https://github.com/marbl/MashMap) uses [`mash`](http://mash.readthedocs.org/) kmer similarity estimates to find matching blocks of long sequences.
These are then chained to produce alignments.
The resulting alignments can be recomputed with edlib to obtain base-exact descriptions of them, which allows the method to be used as an input to seqwish.
In contrast to `minimap2`, `mashmap` does not suffer from problems with repetitive minimizer seeds or low sequence complexity, and in general runs extremely fast, particularly in its approximate mapping mode.

Here, we'll use [a fork `mashmap` that supports output in PAF format and multithreaded alignment](https://github.com/ekg/MashMap).

```
ls seqs/*fa | sort | while read f
do
    echo $f
    in=$f
    out=graphs/seqwish/minimap2/$(basename $f .fa)
    mashmap -r $in -q $in -o $out.mashmap.paf --pi 70 -k 11 -s 500 -t 16 -n 10 \
        && mashmap-align -s $in -q $in --mappingFile $out.mashmap.paf --pi 0 -t 16 -o $out.paf
        && seqwish -s $in -p <(fpa drop -l 3000 <$out.paf) -g $out.gfa -t 16 \
        && odgi build -g $out.gfa -o - \
            | odgi sort -i - -o - -p sYgYs -k 1000 -G 1 -A -t 16 -P \
            | odgi view -i - -g >$out.sort.gfa
done
```

### smoothing tangles with `smoothxg`

Graphs produced by these methods (aside from spoa) can have local complex structures that frustrate many kinds of downstream analyses.
We can flatten these tangles using [`smoothxg`](https://github.com/ekg/smoothxg).
These often result from inducing a graph over a region of low-complexity sequence where many alignment descriptions are equivalent and differentially preferred based on aspects of the aligned sequences and their relative orientation.
To mitigate these issues, we can apply "smoothing" to the graph, which basically collects bins of nodes up to a given amount of embedded path sequence and then subjects these to spoa.
The resulting graphs are then "laced" together by walking back through the original paths.

This method is somewhat experimental, and we need to modify the input graph in a number of ways (sorting it and chopping the nodes to a short length) to ensure that it works well (sorted) and without memory exhaustion (in spoa, due to the quadratic costs associated with aligning very long sequences).

```
find graphs/seqwish/mashmap/*sort.gfa | while read f
do
    echo $f
    o=graphs/smoothxg/seqwish/mashmap/$(basename $f .sort.gfa)
    odgi build -g $f -o - \
        | odgi chop -c 100 -i - -o - \
        | odgi sort -i - -o - -O \
        | odgi sort -i - -o - -p sYgYs -k 1000 -G 1 -A -t 16 \
        | odgi view -i - -g >$o.pre.gfa \
    && time smoothxg -g $o.pre.gfa -w 1000 -j 100 >$o.smooth.gfa.1 \
    && odgi build -g $o.smooth.gfa.1 -o - \
        | odgi chop -c 100 -i - -o - \
        | odgi sort -i - -o - -O \
        | odgi sort -i - -o - -p sYgYs -k 1000 -G 1 -A -t 16 \
        | odgi view -i - -g >$o.smooth.gfa
    rm -f $o.smooth.gfa.1
done
```

The same can be repeated for `minimap2`-based `seqwish` graphs and the `vg msga` graphs.
