## bulkprops
<p>
<img src="../img/bulkprops.png" width=200, alt="bulkprops" />
</p>

If you have a VCF of genotype data for a set of individuals and a BAM file of aligned sequencing data from those individuals, which can be either bulk or single-cell (in which case cell barcodes will be ignored), `bulkprops` can infer the proportion of reads originating from each of the individuals in the pool. It works the same way as the ambient RNA profiling step in [`quant_contam`](quant_contam.md), but it can only consider SNPs for which there are no missing genotypes.

To run:
```
bulkprops -b [mapped_data.bam] -v [variants.vcf.gz] -o [output_prefix]
```
Where `[mapped_data.bam]` is your sequencing data and `[variants.vcf.gz]` is a VCF file (gzipped or not) containing genotype data for the individuals in the pool, containing at least some sites at which there are no missing genotypes (`./.`). `bulkprops` will write a file called `[output_prefix].bgprops`, which will contain these fields (tab-separated) on each line:
* An individual's name in the VCF
* The inferred proportion of reads originating from that individual

#### Optional arguments
```
--qual/q Minimum variant quality to consider
```
Set this to a high number to more strictly filter variants in the VCF.
```
--ids/-i Individuals to include
```
If your VCF file contains individuals not in the pool, make a text file listing the individuals to include, one per line. Specify the file with this option (`--ids/-i`).
```
--n_trials/-N Number of random starting proportions
```
To avoid getting stuck in a local (but not global) likelihood maximum, this option allows `bulkprops` to randomly shuffle initial guesses of mixture proportions a fixed number of times, computing the maximum likelihood solution each time, and reporting the likeliest of all these solutions. When used in conjunction with bootstrap replicates (next option), this can cost lots of execution time, so it is disabled by default. Default initial guess for starting proportions is an even mixture of all individuals.

```
--bootstrap/-B Number of bootstrap samples
```
This option controls the number of bootstrap samples performed in order to compute Dirichlet concentration parameters on mixture proportions, which provide a measure of how variable they are when different sets of SNPs are supplied. If this is done, the resulting file can be used to assess significance when compared to another file of proportions (see [`utils/compare_props.R`](utils_compare_props.md)). If this option is disabled (`-B 0`), then the third column that would normally list Dirichlet concentration parameters will be omitted from the output file. 

```
--error_rate/-e Sequencing/alignment error rate
```
How often you expect reads to mismatch their expected genotypes, due to sequencing error and misalignment. If you do not provide this value, it will be estimated from the data. You might want to provide this value if comparing pool composition across different conditions, to eliminate differences in results that come from differences in this estimated value.

```
--threads/-T Number of processing threads
```
Number of threads to use for parallel processing in likelihood calculations. Reading the BAM file will still be done by a single thread.

### Alternative run mode

It might be the case that you have already computed maximum likelihood mixture proportions and want to infer how variable the pool composition is at different SNPs or genes. In this case, you can provide a file of mixture proportions (or .assignments, which will be converted into mixture proportions by the program) along with a BAM and VCF. If the individuals in the proportions file match those in the VCF, then this program will output the log likelihood of the observed data at various genomic loci, given the mixture of individuals.

To run in this mode, provide:
```
--props/-p The previous file of proportions (i.e. from bulkprops) or assignments (i.e. from demux_vcf)
```
Optionally also provide:
```
--genes/-g Report mean likelihood per gene
```
Without `-g` enabled, `bulkprops` will report the log likelihood of each SNP, in BED format. If `-g` is enabled, then it will instead summarize the log likelihood across genes (reporting the mean value). This option requires that `GN` and `GX` tags are present in the BAM, which are inserted by programs like CellRanger and STARsolo to indicate which gene a read mapped to.

[Back to main README](../README.md)
