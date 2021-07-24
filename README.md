# FragMap
Map and count sequencing reads to genomic regions

It requires 6 positional arguments to run:

1. Path to flowcell
2. A path to a bedfile with the reference genomic regions
3. A FASTA file with the adapters to be trimmed
4. Path to the reference genome
5. MAPQ value
6. The SAM flag `-f` to filter the BAM file

It returns a tab-delimted table with counts per genomic region for each sample
