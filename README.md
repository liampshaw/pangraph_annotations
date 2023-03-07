
[Pangraph](https://github.com/neherlab/pangraph) is an annotation-free approach to building a pangenome graph data structure from sets of genomes. This repository contains a helper script to connect existing genome annotations (in [GFF3](https://www.ensembl.org/info/website/upload/gff3.html) format) to the 'pancontigs' of the graph.

It can add pancontig information to the attributes of gff features or produce a new gff with pancontigs as sequence regions.  

## Dependencies

Written in python, requires `pandas`:

```
conda create -n pangraph_annotations pandas
conda activate pangraph_annotations
```

## Usage

```
# Add pancontigs as attributes to original gff
python add_pancontigs_to_gff.py --pangraph {pangraph.json} \
    --input_gff {genome.gff} \
    --mode keep_original \
    --output_gff {pancontigs_as_attributes.gff}
# Add pancontigs as regions
# WARNING: abuse of the GFF3 format 
python add_pancontigs_to_gff.py --pangraph {pangraph.json}\
    --input_gff {genome.gff} \
    --mode make_new \
    --output_gff {pancontigs_as_regions.gff}
```

Output files will have the original header with an additional header-string e.g.

```
#!pancontig information relative to data/pangraph.json added on 03/07/2023, 10:33:55
```

## Example

Example test data are two *E. coli* genomes: [GCA_000597845.1](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/597/845/GCA_000597845.1_ASM59784v1) and [GCA_000599625.1](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/599/625/GCA_000599625.1_ASM59962v1/). 

Their genomes are in `data/input_genomes.fa`. Also there are GFF3 annotations for `GCA_000597845.1`. 

First, we build a pangraph with

```
pangraph build --circular data/input_genomes.fa > data/pangraph.json
```

Then, we can then add the pancontig information 

```
python add_pancontigs_to_gff.py --pangraph data/pangraph.json \
    --input_gff data/GCA_000597845.1.gff \
    --mode keep_original \
    --output_gff pancontigs_as_attributes.gff
```

This adds pancontig information as an attribute to each annotated feature. Here is a (simplified attributes) example from that output file of a gene that stretches across multiple pancontigs:

```
CP007265.1  Genbank gene    4699050 4700590 .   -   .   ID=gene-BU34_30355;pancontigs=TVFWZRNOTJ-_6,DCKHVIHAKN-_6,LGQMDYQNWO-_7,WFWWUGUICI-_7,TZQFPNNGZQ+_7,FUBHNOVGNG+_7
```

This 1540bp gene (`gene-BU34_30355`) has been fragmented across six pancontigs in the genome. The pancontig string in the attributes (last column of gff) gives a comma-separated list of these pancontigs with the ID, strand, and occurrence of each pancontig e.g. `TVFWZRNOTJ-_6`: pancontig `TVFWZRNOTJ` on the negative (`-`) strand in its sixth (`_6`) occurrence.

If alternatively we wish to know the positions of the annotations on top of the pancontigs, we can use

```
python add_pancontigs_to_gff.py --pangraph data/pangraph.json\
    --input_gff data/GCA_000597845.1.gff \
    --mode make_new \
    --output_gff pancontigs_as_regions.gff
```

Our 1540bp gene (`gene-BU34_30355`) stretches across six pancontigs and so has been split into six lines. **This is an abuse of the GFF3 format so use with caution.**

```
TVFWZRNOTJ  CP007265.1  gene    1   90  .   -   .   ID=gene-BU34_30355-fragment6;parent=gene-BU34_30355;pancontigID=TVFWZRNOTJ;pancontigStrand=-;pancontigN=6
DCKHVIHAKN  CP007265.1  gene    1   636 .   -   .   ID=gene-BU34_30355-fragment5;parent=gene-BU34_30355;pancontigID=DCKHVIHAKN;pancontigStrand=-;pancontigN=6
LGQMDYQNWO  CP007265.1  gene    1   224 .   -   .   ID=gene-BU34_30355-fragment4;parent=gene-BU34_30355;pancontigID=LGQMDYQNWO;pancontigStrand=-;pancontigN=7
WFWWUGUICI  CP007265.1  gene    1   117 .   -   .   ID=gene-BU34_30355-fragment3;parent=gene-BU34_30355;pancontigID=WFWWUGUICI;pancontigStrand=-;pancontigN=7
TZQFPNNGZQ  CP007265.1  gene    1   202 .   -   .   ID=gene-BU34_30355-fragment2;parent=gene-BU34_30355;pancontigID=TZQFPNNGZQ;pancontigStrand=+;pancontigN=7
FUBHNOVGNG  CP007265.1  gene    1   272 .   -   .   ID=gene-BU34_30355-fragment1;parent=gene-BU34_30355;NpancontigID=FUBHNOVGNG;pancontigStrand=+;pancontigN=7
```

Coordinates are now in terms of the pancontigs i.e. `1` is first base of pancontig. The attributes include pancontig `ID`, `Strand` (`-/+`)and `N` (occurrence). The `ID` in the attributes has ID of the 'parent' annotation in the original GFF, with `-fragmentX` added to indicate that it is a fragment.

For a gene that is not fragmented, here is what that looks like (simplifying attributes)

In `pancontigs_as_attributes.gff`:

```
# In pancontigs_as_attributes.gff
CP007265.1  Genbank gene    47  1450    .   +   .   ID=gene-BU34_00005;pancontigs=QSHBRKARJW-_1
```

In `pancontigs_as_regions.gff`:

```
QSHBRKARJW      CP007265.1      gene    21756   23159   .       +       .       ID=gene-BU34_00005;pancontigID=QSHBRKARJW;pancontigStrand=-;pancontigN=1

```

We can use this output file to work out things like how many of the annotated genes are fragmented, and across how many pancontigs:

```
awk '$3=="gene"' pancontigs_as_attributes.gff | sed -e 's/.*pancontigs=//g' | awk -F "," ' { print NF } ' | sort -n | uniq -c 
# 4549 1
#   72 2
#   18 3
#    2 4
#    8 6
```

The majority of the gene features in our chosen genome (`4549`) were on a single pancontig in the pangenome we built with `pangraph`. (This uses the fact that the `pancontigs` attribute is a comma-separated list and counts up how many items are in the list.)

```
# What proportion of genes were split across more than one pancontig?
fragmented=$(awk '$3=="gene"' pancontigs_as_attributes.gff | sed -e 's/.*pancontigs=//g' | grep "," | wc -l)
total=$(awk '$3=="gene"' pancontigs_as_attributes.gff | wc -l)
echo $fragmented / $total | bc -l
# 0.0215 
# i.e. ~2%
```





The majority of genes 
