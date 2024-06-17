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
You can run `quant_contam` in the most basic way easily:
```
quant_contam -o [output_prefix]
```
Where `[output_prefix]` is the `--output_prefix/-o` argument you passed to `demux_vcf`.
### Optional arguments
Some useful optional arguments you can provide are:
```
--ids -i If you limited demux_vcf to specific IDs using the --ids/-i argument, provide
  the same value here.

--ids_doublet -I If you limited demux_vcf to specific IDs using the --ids_doublet/-I argument,
  provide the same value here.

--dump_freqs -d Dump the inferred frequencies of each type of allele in ambient RNA to a file
  called [output_prefix].contam.dat. You probably don't need to do this.

--llr -l Filter the .assignments file to this minimum log likelihood ratio before inferring
  ambient RNA contamination

--no_weights -w Do not weight assignments by log likelihood ratio (see below)

--disable_profile -p If the program is taking too long to run, you can skip the (expensive)
  step of modeling the ambient RNA as a mixture of individuals, and only output estimates of
  percent ambient RNA per cell.

--max_cells -c The (expensive) step of modeling ambient RNA as a mixture of individuals
  subsamples cells for calculating the expected allele frequencies per individual. This
  parameter sets how many cells are sampled (default 50).

--n_mixprop_trials -n When modeling ambient RNA as a mixture of individuals, several trials
  of maximizing the likelihood function are done with different random starting values. This
  sets the number of random trials (default 10)

--other_species -s If your pool contained members of another species that didn't make it into
  the data set analyzed here, attempt to model other species' contribution to ambient RNA
  (see below)

--error_ref -e The underlying true error rate for misreading ref alleles as alt (default = 0.001,
  should be set to sequencer error rate)

--error_alt -E The underlying true error rate for misreading alt alleles as ref (default = 0.001,
  should be set to sequencer error rate)
```
* The `--no_weights/-w` option disables weighting cell-individual assignments by log likelihood ratio (confidence of assignment). This default setting will allow more confident assignments to contribute more to inference of the ambient RNA profile. You might want to disable this if, for example, you have very different numbers of cells per different individuals in your assignments and you are worried that some individuals might mostly be noise. This would allow all cells to contribute equally to the solution.
* The `--other_species/-s` argument is designed for cases in which you have pooled together cells from multiple species, but then separated those cells by species (for example, using [`demux_species`](demux_species.md)) and mapped each group to a different reference genome. In this case, ambient RNA in each data set could have originated from another species that is not present in the data being examined. In this case, this option can include a "fake" individual in the pool consisting of only reference alleles (an approximation for another species, which should generally match the most common allele, wherever applicable).

### Output files
This program (by default) outputs two files:
* `[output_prefix].contam_prof` lists individuals (from the VCF) and the fraction of the ambient RNA pool made up of their RNA. Each line is one individual name followed by a decimal between 0 and 1 indicating their contribution to the pool, tab separated.
* `[output_prefix].contam_rate` lists cell barcodes and the fraction of their RNA likely to have originated from ambient RNA contamination. Each line is one cell barcode followed by a decimal between 0 and 1 indicating the percent ambient RNA contamination in that cell, tab separated.
### Plotting
The program `plot/contam.R` can create plots showing information about ambient RNA contamination in a data set. To run, just pass the `[output_prefix]`:
```
plot/contam.R [output_prefix]
```
Two plots will be created: `[output_prefix].contam.pdf` (vector) and `[output_prefix].contam.png` (rasterized).

#### Example plots
Low contamination data set (10X 40k NSCLC)| High contamination data set (Tetraploid composite iPSCs)|
:--------------------------:|:---------------------------------------:|
![](../img/nsclc.contam.png)  |  ![](../img/tet_ipsc.contam.png) | 

#### Information in the plots
* The top left corner of the plot lists the mean and standard deviation of contamination rate per cell in the data set. 
* The bar plot in the top left shows the fraction of ambient RNA inferred to have originated from each individual, in the form of stacked bars. If you have very many individuals, this information may overflow the space allocated in the plot; in that case, you will need to load `[output_prefix].contam_prof` and plot in your favorite plotting program.
* The top right corner is a kernel density curve showing the contamination rate per cell across the data set. `quant_contam` uses an Empirical Bayes prior to shrink per-cell estimates toward the mean; provided per-cell estimates are maximum *a posteriori* rather than maximum likelihood estimates.
* The panel below the density curve shows the contamination rate in each cell (Y-axis) against the log likelihood ratio of the cell's identity (X-axis). Error bars show uncertainty in contamination rate estimates (using Fisher information). Low confidence assignments may have higher inferred contamination rates, reflecting some incorrectly-assigned cells.
* The bottom left panel shows each cell's contamination rate as a horizontal bar, colored by the individual ID of each cell. Error bars show uncertainty in contamination rate estimates (using Fisher information).
* The bottom right panel shows the distribution of contamination rate estimates for cells asssigned to each individual ID as box plots, colored the same way as the bottom left panel and top left panel.

[Back to main README](../README.md)
