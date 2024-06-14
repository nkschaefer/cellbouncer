# demux_vcf
<p>
<img src="../img/demux_vcf.png" width=250 alt="demux_vcf" />
</p>
A common task in pooled culture experiments is to identify every cell's individual of origin, given genotype data for each individual in the pool that was collected beforehand. Unfortunately, some methods do not perform well when too many or too few variant sites are included, or when the depth of sequencing of the single cell data falls outside an ideal range. 

`demux_vcf` is designed to identify cells' individuals of origin quickly and accurately and is also capable of identifying specific combinations of two individuals (as appears in 2-cell "doublet" droplets). Additionally, the output of `demux_vcf` can be provided to another `cellbouncer` program, `quant_contam`, which quantifies the ambient RNA contamination per cell and models the contribution of each individual in the pool to the ambient RNA.

## Input data
`demux_vcf` requires variant data to be in the [VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf). 

### Converting from PLINK
Sometimes, genotyping data is provided in [PLINK](https://zzz.bwh.harvard.edu/plink/) format (`.map` and `.ped` files). To convert these files to VCF, first download and install `plink`. 

Then, run 
```
plink --file [base] --recode vcf --out [out_base] --real-ref-alleles
```
where `[base]` is the beginning of your `.map` and `.ped` file names, and `[out_base]` will be the beginning of the output file (before the `.vcf` extension).

The output file may contain some weird/phony variants (sequence name and position 0). To get around that, do this: 
```
cat [out_base].vcf | awk '{if (substr($1, 1, 1) == "#"){ print $line; } else if ( $2 != 0 && $3 != 0){ print $line; }}' > [out_base].filt.vcf
```

Next, if you plan on filtering with a tool like [`bcftools`](https://samtools.github.io/bcftools/) (see below), compress with [`bgzip`](http://www.htslib.org/doc/bgzip.html) and index with [`tabix`](http://www.htslib.org/doc/tabix.html): 
```
bgzip [out_base].filt.vcf
tabix -p vcf filt.vcf.gz
```

### Filtering variants
Sometimes, and especially if you have called variants yourself from whole-genome sequencing or exome sequencing data, you will need to filter variants to eliminate low-quality ones that result from sequencing error. We recommend [`bcftools`](https://samtools.github.io/bcftools/). A basic filtering step would be as follows:
```
bcftools view -O z -m 2 -M 2 -v snps -i "(QUAL >= 100) & (F_MISSING < 0.5)" [input.vcf.gz] > [output.vcf.gz]
```
This will limit the variant calls to biallelic SNPs with variant quality of at least 100 and no more than 50% of individuals missing genotypes.

Another important step is to exclude all individuals in the panel that were not used in your single cell sequencing experiment. This can be accomplished by running `bcftools view` with the argument `-S [file.txt]` where `file.txt` is a text file listing each sample you want to include, one per line. If you include extra individuals in the VCF, you increase execution time and memory usage as well as risking identifying a few of your cells as individuals you know not to be present.

### Is imputation required?
Sometimes (especially when working with SNP chip genotyping data), relatively few variants will make it into the VCF. Some tools recommend using a reference panel like the [high-coverage genomes from the HGDP panel](https://pubmed.ncbi.nlm.nih.gov/32193295/) or [1000 Genomes Project data](https://www.internationalgenome.org/) to impute genotypes at other sites in the genome. In our experience, this is unnecessary; `demux_vcf` seems to work well even with only about 100,000 variants (although more is better). 

In the case of non-human data, imputation is generally impossible due to the lack of phased genomes. 

## Running the program
Simply run as follows:
```
demux_vcf -b [input.bam] -v [input.vcf] -B [filtered_barcodes.tsv] -o [output_base]
```
Where `[input.bam]` is the aligned BAM file of single cell sequence data, `[input.vcf]` is the VCF format data of variants segregating among the individuals used in the experiment, `[output_base]` is the base name to use for output files, and `[filtered_barcodes.tsv]` is an optional but recommended argument: the filtered list of barcodes determined by the alignment program to represent true cells (i.e. usually located in the `filtered_feature_bc_matrix` directory within [CellRanger](https://www.10xgenomics.com/support/software/cell-ranger/latest) output. This file can be gzipped.

This process can take several hours, depending on how deeply sequenced the single-cell data is and how many variants are in the VCF.

This will create the following output files:
* `[output_base].assignments` contains the most likely identity assigned to each cell.
* `[output_base].counts` is gzipped data containing the counts of each type of allele in each cell used to make the assignments. If you run `demux_vcf` again with the same `[output_base]`, this file will be loaded instead of going through the expensive step of computing these counts again.
* `[output_base].samples` contains the names of the samples in the VCF, in the same order, for use by the `quant_contam` program, if you choose to run it.
* `[output_base].summary` contains some summary information about the run:
  * Where the second column is `param`, the third and fourth columns list parameters used by `demux_vcf` on this run.
  * Where the second column is `result`, the third and fourth columns list inferred parameters/results computed in this run.
    * `error_ref_posterior` gives the inferred rate at which reference alleles were misread as alt alleles in this data set. Rates much higher than 1-5% are indicative of high ambient RNA contamination.
    * `error_alt_posterior` gives the inferred rate at which alternate alleles were misread as ref alleles in this data set. Rates much higher than 1-5% are indicative of high ambient RNA contamination.
    * `tot_cells` gives the total number of cells assigned identites.
    * `frac_doublets` tells how many cell barcodes were inferred to be droplets containing two different individuals (note that this is lower than the actual doublet rate, which includes self-self doublets.
    * `doublet.chisq.p` gives the p value of a chi-squared test: given the prevalence of each individual in the pool (computed from singlet identifications) and the doublet rate (approximated as the number of inter-individual doublets divided by the number of singlets), it computes the expected number of doublets of each type and compares to the identified number of doublets of each type. A high p-value means that the results agreed well with expectation.
  * Where the second column is `perc_doublets_include`, the fourth column gives the percent of inter-individual doublets that include the individual in the third column. 
  * Where the second column is `p_fewer_cells`, the third and fourth columns give the p-value that the ID in the third column has a significantly lower cell count than the average ID across the data set. This is computed using the Poisson CDF.
  * Where the second column is `p_lower_LLRs`, the third and fourth columns give the p-value that the ID in the third column has a significantly lower distribution of log likelihood ratios of assignments (a measure of confidence) than the rest of the data set as a whole. This is computed using a Mann-Whitney U-test of the distribution of LLRs for the individual against the distribution of all other LLRs.
  * Where the second column is `num_cells`, the fourth column gives the total number of cells assigned to the individual in the third column.
    
