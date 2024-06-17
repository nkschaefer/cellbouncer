# utils/bam_split_bcs
The program `utils/bam_split_bcs` is designed for special cases where you want to divide up the data from a single-cell sequencing experiment by individual of origin. It takes as input a BAM file and outputs a set of BAM files, one per individual. It requires an [`assignments` file](../README.md#output-files) and will omit any cells with barcodes missing from the `.assignments` file. 

## How to run

```
utils/bam_split_bcs -b [BAM file] \
  -a [assignments file] \
  -o [output_prefix] \
  (-d) \
  (-l LLR)
```

The optional argument `-d` keeps doublets and will create a new BAM file for each type of doublet in the `assignments` file.

The optional argument `-l` sets a minimum log likelihood ratio of assignment (a measure of confidence) for a cell to be included in the BAM file created for its assigned individual.

## One more thing
Don't forget to index the resulting BAM files after they're created, using [`samtools index`](https://github.com/samtools/samtools).

[Back to main README](../README.md)
