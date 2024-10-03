# utils/compare_props.md

If you have run a program (like [`bulkprops`](bulkprops.md)) that models the RNA as a mixture of individuals, you may wish to compare the output pool composition to another set of proportions (i.e. the same cells at a later time point or under another treatment). A difficulty here is how to compute significance: two sets of proportions can be visually compared, but it is less straightforward to determine whether any differences between them are statistically significant.

To this end, the CellBouncer programs that infer pool composition ([`bulkprops`](bulkprops.md) and the `.contam_prof` file output by [`quant_contam`](quant_contam.md)) can perform bootstrap resampling (and do by default). This allows these programs to compute many estimates of the mixture proportions using different sets of input allele frequencies. From these estimates, they then fit a Dirichlet distribution, and the concentration parameters of this distribution indicate how sensitive the maximum likelihood mixture proportions are to fluctuations in the input data. If bootstrapping was done, an output file will have a third column appended containing the maximum likelihood Dirichlet concentration parameters.

You can then assess the significance of the differences in mixture proportions between two files using `utils/compare_props.R`, as follows:

```
utils/compare_props.R [reference] [test]
```
`[reference]` is a file of mixture proportions with bootstrap resampling performed. It serves as the null hypothesis against which the `[test]` file is tested for differences.

`[test]` can be a file of mixture proportions or an `.assignments` file that contains the same labels (i.e. if `[reference]` was produced by running `bulkprops` with a VCF file of named samples and `[test]` was produced by running `demux_vcf` with the same VCF, this will work. On the other hand, if `[test]` is an `.assignments` file containing numeric indices for inferred mitochondrial haplotypes, this will not work unless you determine which individual corresponds to which mitochondrial haplotype and then rename the labels in the `.assignments` file to match the individuals). 

By the aggregation property of the Dirichlet distribution, we can then test for over/underrepresentation of each individual component using one of the marginal beta distributions. `utils/compare_props.R` reports the two-sided p-value for the null hypothesis that each individual's proportion in file 2 is the same as in file 1, corrected for multiple hypothesis testing using the Benjamini-Hochberg (FDR) method.


[Back to main README](../README.md)
