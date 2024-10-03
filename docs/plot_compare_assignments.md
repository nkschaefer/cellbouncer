# plot/compare_assignments.R

<img src="../img/vcf_vs_mito.png" width=400 alt="compare_assignments" />

`plot/compare_assignments.R` compares two sets of assignments on the same cells. A typical use case for this program might be to compare the results of `demux_vcf` and `demux_mt` run on the same data set, to see whether the programs worked, and potentially to assign an individual of origin (from the VCF) to each mitochondrial haplotype. To generate this plot, run 

```
plot/compare_assignments.R [output_prefix1.assignments] [output_prefix2.assignments] [plot_out] (S)
```

where `[output_prefix1.assignments]` is the assignments file from the first program run, `[output_prefix2.assignments]` is the assignments file from the other program run, `[plot_out]` is the base output name for plots (both `.pdf` and `.png` format plots will be created), and `S` is an optional final argument that, if present, will limit the plot to show singlet identifications only.

All possible arguments for the last argument are:

```
S = only use singlets from both files
S1 = only use singlets from the first file
S2 = only use singlets from the second file
D = treat doublets as unique categories (do not consider as combinations of singlets)
D1 = treat doublets as unique categories only for the first file
D2 = treat doublets as unique categories only for the second file
S1D2, D2S1, S2D1, D1S2 = combinations of above
```

In addition to creating a plot, `plot/compare_assignments.R` prints a table showing how many cells receive each combination of labels from file 1 and file 2, limited to only singlets.

If labels with different names are passed (i.e. VCF-based sample names in one label set, and numeric mitochondrial haplotype IDs in the second label set), the program will attempt to identify corresponding labels in the two data sets. To do this, it measures the [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) between each pair of labels in the first and second file. Any pair of labels that both produce the highest Jaccard index for the label in the other file are chosen as a corresponding pair. 

Given these corresponding labels, `plot/compare_assignments.R` computes the sensitivity, specificity, and F statistic conditional on file 1 representing the truth, and conditional on file 2 representing the truth.

It also computes the above values for singlets and doublets separately.

It will also report cases where a label from one file has a highest-Jaccard partner in the other file for which the first label is not the reciprocal highest-Jaccard partner for the second label. This can be indicative of problems -- i.e. one mitochondrial cluster might correspond to two genotype-based identities, possibly because two individuals share a recent maternal ancestor. 

