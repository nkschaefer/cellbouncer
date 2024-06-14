<p>
<img src="img/logo.png", width=300, alt="CellBouncer" />
</p>

Tools for demultiplexing and keeping the riffraff out of pooled single cell sequencing data sets. 

|I want to...|I have...|Tool to use|
|------------|---------|-----------|
|Demultiplex cells by **species** before mapping to a reference genome|A transcriptome (FASTA) or annotation (GTF) and genome (FASTA) per species|[`demux_species`](#demux_species)|
|Demultiplex cells by **individual of origin**, and I hope individuals are unrelated enough to have different mitochondrial haplotypes|Nothing|[`demux_mt`](#demux_mt)|
|Demultiplex cells by **individual of origin**|VCF of known variants|[`demux_vcf`](#demux_vcf)|
|Demultiplex individuals by **custom label** or **treatment**|MULTIseq/HTO/CITE-seq data|`demux_tags`|
|Assign **sgRNAs** to cells|sgRNA capture data|`demux_tags`|
|Quantify **ambient RNA** per cell and infer its origins|VCF of known variants|`quant_contam`|
|Infer global **doublet rate**|Output from a `CellBouncer` program|`doublet_dragon`|

# Installation
The non-plotting programs in `cellbouncer` require only [HTSLib](https://github.com/samtools/htslib) and [zlib](https://www.zlib.net/). You can install these dependencies yourself, or to make it easier, get these and plotting-related dependencies by creating a [conda](https://github.com/conda-forge/miniforge/releases) environment from the included `cellbouncer.yml` file. Then all that's left to do is `make`:

```
git clone --recurse-submodules https://github.com/nkschaefer/cellbouncer.git
cd cellbouncer
conda env create --file=cellbouncer.yml
conda activate cellbouncer
make
```
You've now got all the programs compiled, and you can run them as long as you remember to `conda activate cellbouncer` first.

# Overview
The programs in `cellbouncer` are standalone command line tools. If you run one of them with no arguments or with `-h`, it will give you detailed information about how to run it. In general, though, each program writes output files that start with a name you give via `--output_prefix` or `-o`. 

## Output files
Demultiplexing tools all write a file called `[output_prefix].assignments`, which tells you information about each cell's identity. These files are 4 columns, tab separated: cell barcode, most likely identity, droplet type: S (for singlet), D (for doublet), or in some cases M (for multiplet, 3+ individuals, but not all programs look for these), and a ratio of the log likelihood of the best to the second best assignment. 

In output files, cell barcodes will by default be printed without any additional text (they will consist only of DNA sequences). When there are multiple data sets to be analyzed together, additional text must be appended to barcodes from each data set to prevent barcode collisions. Different programs have different conventions for handling this, such as [CellRanger](https://www.10xgenomics.com/support/software/cell-ranger/latest) appending "-" and a numeric ID (starting from 1) to cell barcodes. In [scanpy](https://scanpy.readthedocs.io/en/stable/), the [anndata.concatenate](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.concatenate.html) command also follows this convention, unless the `batch_categories` argument is used. Some `CellBouncer` programs have a `--batch_id` argument that allows users to append unique identifiers in the same format (separated from the barcode sequence with `-`).

# Programs

## [demux_species](docs/demux_species.md)
<p>
<img src="img/demux_species.png", width=200, alt="demux_species" />
</p>
Before mapping data, infer the species of origin of each cell barcode by counting k-mers unique to each species' transcriptome. Separate FASTQ files by species and optionally plot species abundances.

[[more]](docs/demux_species.md)


## [demux_mt](docs/demux_mt.md)
<p>
<img src="img/mito.png", width=200, alt="demux_mt" />
</p>
Using a BAM file of aligned scATAC-seq (ideally) or whole-cell scRNA-seq data containing cells originating from multiple individuals, infer the set of mitochondrial haplotypes in the mixture, as well as the number of individuals. Assign each cell an identity based on its likeliest mitochondrial haplotype. These assignments can then be used to label individuals of origin in the BAM, and a variant caller can then identify genomic SNPs and their genotypes in the inferred individuals.

[[more]](docs/demux_mt.md)


## [demux_vcf](docs/demux_vcf.md)
<p>
<img src="img/demux_vcf.png" width=250, alt="demux_vcf" />
</p>
Given genotype data for the individuals in a pool and a BAM file of aligned single cell sequencing data, quickly infer the individual (or doublet) of origin of each cell in the pool. Confidently identifies specific doublets of origin where they occur and has been shown to be accurate even in identifying the correct contributor cell lines in the case of composite cell lines formed through inter-species cell fusions.

[[more]](docs/demux_vcf.md)

## [demux_tags](docs/demux_tags.md)

[more](docs/demux_tags.md)

## Plotting
In the `plot` directory, there are R scripts to plot output from some of the programs. If you run one with no arguments, it will tell you how to run it. Plotting programs are described in more detail on the README pages for specific tools.

`plot/assignment_llr.R`| `plot/compare_assignments.R` | 
:--------------------------:|:---------------------------------------:
![](img/assn_llr.png)  |  ![](img/vcf_vs_mito.png)

Two plotting programs are useful across multiple `cellbouncer` programs, however.

`plot/assignment_llr.R` plots a set of cell-to-identity assignments as stacked bars, to show the proportion of the pool made up of cells of each assigned identity. Along the X-axis, the heights of these bars change to show how many cells of each category would remain if the user filtered the data using the log likelihood ratio cutoff on the X-axis. This is useful for visualizing the total number of cells, pool composition, and feasibility of different filtering thresholds. To generate this plot, simply run 
```
plot/assignment_llr.R [output_prefix]
```
where `[output_prefix]` is the output prefix given to the `cellbouncer` program you just ran.

`plot/compare_assignments.R` compares two sets of assignments on the same cells. A typical use case for this program might be to compare the results of `demux_vcf` and `demux_mt` run on the same data set, to see whether the programs worked, and potentially to assign an individual of origin (from the VCF) to each mitochondrial haplotype. To generate this plot, run 
```
plot/compare_assignments.R [output_prefix1.assignments] [output_prefix2.assignments] [plot_out] (S)
```
where `[output_prefix1.assignments]` is the assignments file from the first program run, `[output_prefix2.assignments]` is the assignments file from the other program run, `[plot_out]` is the base output name for plots (both `.pdf` and `.png` format plots will be created), and `S` is an optional final argument that, if present, will limit the plot to show singlet identifications only.

## Other tools
The `utils` directory is the junk drawer of `cellbouncer`. It contains several programs meant for specific tasks, which can aid the above programs.
|Name|Purpose|
|----|-------|
|[`utils/get_unique_kmers`](docs/demux_species.md#preparing-data)|For use in demultiplexing reads using species-specific k-mers. Takes lists of k-mers in different species' transcriptomes from [FASTK](https://github.com/thegenemyers/FASTK) and outputs lists of high-complexity k-mers unique to each species.|
|[`utils/split_read_files`](docs/demux_species.md#running-in-parallel)|Splits FASTQ files (paired or not) into a set number of approximately evenly-sized chunks|
|[`utils/combine_species_counts`](docs/demux_species.md#running-in-parallel)|Combines the output of multiple runs of [`demux_species`](docs/demux_species.md) on chunks of data into a single file that can then be used to demultiplex reads|
|[`utils/bam_indiv_rg`](docs/demux_mt.md#utilsbam_indiv_rg)|Add read groups to a BAM file to mark cells' individuals of origin, so a variant caller can find variant sites that segregate among identified individuals|
|`utils/bam_split_bcs`|Takes BAM file of aligned single cell sequencing data and an [`.assignments` file](#output-files) from a `cellbouncer` program and outputs one BAM file per identity in the `.assignments` file.|
|`utils/atac_fq_preprocess`|Finds valid cell barcodes in scATAC-seq data (including from [10X multiome](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression) experiments) and outputs new FASTQ files with cell barcodes inserted as sequence comments, where they can be transformed into BAM tags by some aligners.|
