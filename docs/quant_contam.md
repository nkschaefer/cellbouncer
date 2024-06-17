# quant_contam
<p>
<img src="../img/quant_contam.png" width=225 alt="quant_contam" />
</p>

After you have run [`demux_vcf`](demux_vcf.md) on single-cell RNA-seq data for which you have pre-existing variant call data for the individuals of origin, you can run `quant_contam` to profile and quantify the ambient RNA contamination in your data set.

## How this differs from other methods
Other methods for profiling ambient RNA often rely on the idea that certain genes will be more highly expressed in cells affected by ambient RNA contamination. This may be the case, but it requires ambient RNA contamination to be relatively low and heterogeneous in your data set. `quant_contam`, however, quantifies ambient RNA contamination by seeing how often cells mismatch alleles that they should possess, according to a preexisting set of variant calls (external data). This means that `quant_contam` can quantify ambient RNA contamination, even when it is extensive and uniform across all cells in your data set.

In addition to quantifying ambient RNA contamination, `quant_contam` models the ambient RNA as a mixture of individuals in your sequencing pool. This means that you can see whether or not specific individuals contribute an outsize proportion of the ambient RNA. If this is the case, and contamination is extensive, this could bias results in differential expression analyses downstream.

Rather than attempt to correct your count matrices, `quant_contam` reports the percent of each cell's RNA likely to have originated from ambient RNA contamination. With this information, you can filter individual cells, include the values as covariates in a single-cell level differential expression analysis, and/or test which genes have their expression most positively and negatively correlated with ambient RNA contamination rates in your data set.

## Running the program
