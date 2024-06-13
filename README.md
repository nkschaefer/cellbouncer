<p>
<img src="img/logo.png", width=300, alt="CellBouncer" />
</p>

Tools for demultiplexing and keeping the riffraff out of pooled single cell sequencing experiments. 

|I want to...|I have...|Tool to use|
|------------|---------|-----------|
|Demultiplex cells by **species**|A transcriptome (FASTA) or annotation (GTF) and genome (FASTA) per species|[`demux_species`](#demux_species)|
|Demultiplex cells by **individual of origin**, where individuals have different mitochondrial haplotypes|Nothing|`demux_mt`|
|Demultiplex cells by **individual of origin**|VCF of known variants|`demux_vcf`|
|Demultiplex individuals by **individual of origin** or **treatment**|MULTIseq/HTO/CITE-seq data|`demux_tags`|
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

# The main I/O idea
The programs in `cellbouncer` are standalone command line tools. If you run one of them with no arguments or with `-h`, it will give you detailed information about how to run it. In general, though, each program writes output files that start with a name you give via `--output_prefix` or `-o`. Demultiplexing tools all write a file called `[output_prefix].assignments`, which tells you information about each cell's identity. These files are 4 columns, tab separated: cell barcode, most likely identity, droplet type: S (for singlet), D (for doublet), or in some cases M (for multiplet, 3+ individuals, but not all programs look for these), and a ratio of the log likelihood of the best to the second best assignment. 

In the `plot` directory, there are R scripts to plot output from some of the programs. If you run one with no arguments, it will tell you how to run it.

# Programs
<p>
<img src="img/demux_species.png", width=200, alt="demux_species" />
</p>

## demux_species
This program can be used before aligning your data, to assign cells to species of origin and separate reads by species, so you can then map data from each species to its own reference genome and annotation. To accomplish this, it counts species-specific k-mers. 

To run, you must first build a set of reference k-mers. To accomplish this,
* Obtain a transcriptome (in FASTA format) for each species in the pool.
  * If you do not have one, you can create one from a [GTF](http://genome.ucsc.edu/FAQ/FAQformat#format4) or [GFF](http://genome.ucsc.edu/FAQ/FAQformat#format3)-format annotation and a FASTA reference genome indexed with [`samtools faidx`](https://github.com/samtools/samtools). Install [gffread](https://github.com/gpertea/gffread) and run it like this: `gffread -F -w [output_tx.fa] -g [genome_fasta.fa] [unzipped_gtf.gtf]` where `[output_tx.fa]` is the output file of transcripts to be created, `genome_fasta.fa` is an unzipped genome in FASTA format, and `unzipped_gtf.gtf` is the unzipped gene annotation.
* Run [FASTK](https://github.com/thegenemyers/FASTK) to count k-mers in each species' transcriptome. You should use a value of k that is large enough to be unique but small enough to allow for multiple k-mers per read and limit the amount of RAM required; a value between 25-45 should work well. To run FASTK, do this:
  * `FastK -N[output_prefix] -k[kmer_size] -t1 [output_tx.fa]` where `[output_prefix]` will be the beginning of the name of the output files and `[output_tx.fa]` is the file you created in the last step. Note that there should be no space between argument names and argument values.
* Compare each species' k-mer sets to find k-mers unique to each species' transcriptome, and create the files needed by `demux_species`
  * Run `utils/get_unique_kmers` in this package, with `-k [name.ktab]` specified once for each species you ran `FastK` on and `-n [name]` once again for each species, in the same order. For example, if you have the files `human.ktab`, `chicken.ktab`, and `mouse.ktab` representing human, chicken, and mouse k-mers, you could run `utils/get_unique_kmers -k human.ktab -k chicken.ktab -k mouse.ktab -n Human -n Chicken -n Mouse -o hcm_kmers`. This will create a set of "unique k-mers" files beginning with `hcm_kmers`.
* Run `demux_species` using the unique k-mer lists you have just created as the `-k` argument: in the example above, `-k hcm_kmers`. If you run `demux_species` with no arguments, there will be a detailed usage screen. At minimum, you will need to provide a file of forward and reverse single cell RNA-seq reads, an output prefix, a [cell barcode whitelist](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist), and a unique k-mer set. If you have multiple types of data from cells that share the same barcodes (i.e. [10X Genomics multiome data](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression), or [feature barcoding](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-what-is-feature-bc) data for the same cells), you can provide these files as well. The RNA-seq alone will be checked against the species-specific k-mers, but once cells are assigned to species, the other types of reads will be split into separate files by species as well.
  * This program will output the following:
    * `species_counts.txt` A table of cell barcodes and counts of unique k-mers from each species
    * `species_names.txt` A file listing species names in the same order as k-mer counts are listed in `species_counts.txt`
    * `species.assignments` A file mapping cell barcodes to species (or inter-species doublets), along with log likelihood ratios/confidence of assignments
    * `species.filt.assignments` The same as above, but filtered to (hopefully) exclude many non-cell barcodes, in order to improve plotting and give a better idea of the proportions of each species in the pool
    * Subdirectories for each species, containing the input read files (with the same names), subset to only reads from cells assigned to that species
    * Species-specific [10X Genomics-format library files](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-libraries-csv), to aid in running [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) on data from each species separately



