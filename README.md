
Test data are two *E. coli* from Antwerp:
```
GCA_000597845.1
GCA_000599625.1
```

These are in `data/input_genomes.fa`. 

We build a pangraph with

```
pangraph build --circular data/input_genomes.fa > data/pangraph.json
```

```
conda create -n pangraph_annotations pandas
conda activate pangraph_annotations
```

We can add the pancontig information for an example set of gff annotations like 

```
python add_pancontigs_to_gff.py --pangraph data/pangraph.json --input_gff data/GCA_000597845.1.gff --output_gff test_output.gff
```

This adds pancontig information as an attribute to each annotated feature e.g.

```
CP007265.1	Genbank	gene	1455	2555	.	+	.	ID=gene-BU34_00010;Name=BU34_00010;gbkey=Gene;gene_biotype=protein_coding;locus_tag=BU34_00010;pancontigs=QSHBRKARJW-_1
```

The pancontig string includes the ID, the strand, and the occurrence of the pancontig in the specified genome (here, the first occurrence of QSHBRKARJW and the block is on the negative strand).

(NOT DONE) Go the other way, and create a gff from a fasta file of pancontigs and a set of annotations i.e. a gff relative to pancontigs in a particular strain

```
PANCONTIG	CP007265.1	gene    10	1110   .       +       .             ID=gene-BU34_00010
```

