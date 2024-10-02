<p>
<img src="../img/demux_species.png", width=200, alt="demux_species" />
</p>

# demux_species
This program can be used before aligning your data, to assign cells to species of origin and separate reads by species, so you can then map data from each species to its own reference genome and annotation. This approach uses counts of species-specific k-mers to determine each cell's species of origin. Because there are so many potential genomic k-mers but a more limited number of transcriptomic k-mers, we recommend only using this approach on transcriptomic data (if you have multiome, feature barcoding, or other data where cell barcodes are shared across modalities, you can use the transcriptomic reads to demultiplex the other data types as well). 

Alternatively, if you map reads to a composite reference genome, you can use [`utils/composite_bam2counts`](#counting-reads-on-a-composite-reference-genome) to count reads per cell mapped to each species' reference genome, and then run `demux_species` to assign cells to species of origin and optionally demultiplex reads.

### Preparing input data (k-mer counting method)

To run, you must first build a set of reference k-mers. To accomplish this,
* Obtain a transcriptome (in FASTA format) for each species in the pool.
  * If you do not have one, you can create one from a [GTF](http://genome.ucsc.edu/FAQ/FAQformat#format4) or [GFF](http://genome.ucsc.edu/FAQ/FAQformat#format3)-format annotation and a FASTA reference genome indexed with [`samtools faidx`](https://github.com/samtools/samtools). Install [gffread](https://github.com/gpertea/gffread) and run it like this: `gffread -F -w [output_tx.fa] -g [genome_fasta.fa] [unzipped_gtf.gtf]` where `[output_tx.fa]` is the output file of transcripts to be created, `genome_fasta.fa` is an unzipped genome in FASTA format, and `unzipped_gtf.gtf` is the unzipped gene annotation.
* Run [FASTK](https://github.com/thegenemyers/FASTK) to count k-mers in each species' transcriptome. You should use a value of k that is large enough to be unique but small enough to allow for multiple k-mers per read and limit the amount of RAM required; a value between 25-45 should work well. To run FASTK, do this:
  * `FastK -N[output_prefix] -k[kmer_size] -t1 [output_tx.fa]` where `[output_prefix]` will be the beginning of the name of the output files and `[output_tx.fa]` is the file you created in the last step. Note that there should be no space between argument names and argument values.
* Compare each species' k-mer sets to find k-mers unique to each species' transcriptome, and create the files needed by `demux_species`
  * Run `utils/get_unique_kmers` in this package, with `-k [name.ktab]` specified once for each species you ran `FastK` on and `-n [name]` once again for each species, in the same order. For example, if you have the files `human.ktab`, `chicken.ktab`, and `mouse.ktab` representing human, chicken, and mouse k-mers, you could run:
    
    ```
    utils/get_unique_kmers -k human.ktab -k chicken.ktab -k mouse.ktab -n Human -n Chicken -n Mouse -o hcm_kmers
    ```
    
    This will create a set of "unique k-mers" files beginning with `hcm_kmers`.

### Running the program (k-mer counting method)
* Run `demux_species` using the unique k-mer lists you have just created as the `-k` argument: in the example above, `-k hcm_kmers`. If you run `demux_species` with no arguments, there will be a detailed usage screen. At minimum, you will need to provide a file of forward and reverse single cell RNA-seq reads, an output prefix, a [cell barcode whitelist](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist), and a unique k-mer set. If you have multiple types of data from cells that share the same barcodes (i.e. [10X Genomics multiome data](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression), or [feature barcoding](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-what-is-feature-bc) data for the same cells), you can provide these files as well. The RNA-seq alone will be checked against the species-specific k-mers, but once cells are assigned to species, the other types of reads will be split into separate files by species as well.

Some additional, optional arguments are:
```
--doublet_rate -D The prior expectation of a cell being an inter-species doublet (default = 0.5)
--disable_umis -u Do not collapse UMIs (can slightly improve speed at the risk of duplicate reads inflating species counts)
--exact -e Require exact matches to cell barcode list (will improve speed at the risk of missing some cells)
--rna_r1/-r, --rna_r2/-R = RNA-seq read pairs to demultiplex (will be used to count k-mers)
--atac_r1/-1, --atac_r2/-2, --atac_r3/-3 = ATAC-seq read triplets to demultiplex (will not be used to count k-mers)
--custom_r1/-x, --custom_r2/-X, --names_custom/-N = Custom-format read pairs, plus names of the data types they contain
--whitelist_rna/-w The allowed cell barcode list for RNA-seq reads (also used for custom reads)
--whitelist_atac/-W The allowed cell barcode list for ATAC-seq (lines in -w and -W should correspond if multiome data)
--atac_preproc/-A If demultiplexing ATAC-seq, rather than just writing out the read triplets as they came in, outputs forward and reverse read pairs of genomic sequence (with barcodes removed), with error-corrected cell barcodes inserted in sequence comments as SAM tags.
```
Cell barcodes will be error corrected in output, demultiplexed reads, to speed up downstream processing.

  * If you set the `-A` option above, ATAC-seq reads will be output as genomic read pairs with barcodes removed, error corrected, and appended to sequence comments in SAM tag form (i.e. `CB:Z:AGTGGACGTGAAGTGA`). These can then be aligned to a reference genome using a paired read aligner capable of inserting sequence comments as BAM tags: 
    * `minimap2 -a -x -sr -y`
    * `bwa mem -C`
    * `bowtie2 --sam-append-comment`
  * This program will output the following:
    * `species_counts.txt` A table of cell barcodes and counts of unique k-mers from each species
    * `species_names.txt` A file listing species names in the same order as k-mer counts are listed in `species_counts.txt`
    * `species.assignments` A file mapping cell barcodes to species (or inter-species doublets), along with log likelihood ratios/confidence of assignments
    * `species.filt.assignments` The same as above, but filtered to (hopefully) exclude many non-cell barcodes, in order to improve plotting and give a better idea of the proportions of each species in the pool
    * `dists.txt` contains parameters from the fit multinomial mixture model used to assign cells to species.
    * Subdirectories for each species, containing the input read files (with the same names), subset to only reads from cells assigned to that species
    * Species-specific [10X Genomics-format library files](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-libraries-csv), to aid in running [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) on data from each species separately
#### Extra options

### Note about cell barcode lists
`demux_species` needs to check for valid cell barcodes in reads. To do this, it requires one or more [cell barcode whitelists](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist) for the FASTQ files you provide. If using RNA-seq data only, only one whitelist is required (`-w` option). If using 10X Genomics multiome data, however, two whitelists are required, one for RNA-seq (`-w` option) and one for ATAC-seq (`-W` option). This is because these kits use separate barcodes for ATAC-seq and RNA-seq, where each ATAC-seq barcode corresponds to an RNA-seq barcode and is converted to the matching RNA-seq barcode in output data. For help finding the multiome whitelist files, see [here](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist).
### Running in parallel
Unfortunately, matching cell barcodes to large whitelists and counting k-mers can be slow. To help speed things up, this program is designed to be multi-threaded (both the k-mer counting and read demultiplexing steps are designed to use multiple threads). You can specify the number of threads with the `-T` argument to `demux_species`.

Alternatively, if you have access to a compute cluster, `cellbouncer` has a way to split up data, run in parallel on a cluster, and join results.
* Chop up the input read files using `utils/split_read_files`
  * Usage:

    ```
    utils/split_read_files -1 MyLibrary_S1_L001_R1_001.fastq.gz \
        -2 MyLibrary_S1_L001_R2_001.fastq.gz \
        -o [output_directory] \
        -n [number of chunks]
    ```
    
  * This will create files in `[output_directory]` with the same names as the input read files, but with a 1-based numeric index appended to the end.
* Run `demux_species` on each chunk in batch mode, using the same output directory for all runs
  * Pass one forward/reverse read pair file chunk in: i.e. `-r MyLibrary_S1_L001_R1_001.1.fastq.gz -R MyLibrary_S1_L001_R2_001.1.fastq.gz` and add the chunk number, so it can be appended to output files: i.e. `--batch_num 1`
* Ensure all runs have completed successfully: there should be a `species_counts.[batch_num].txt` and `species_names.[batch_num].txt` for each run in the output directory
* Join all runs together using `utils/combine_species_counts` with `-o` set to the output directory you used and `-n` set to the number of chunks/batch numbers
* Re-run `demux_species` with the same output directory you used for all runs, but now provide all reads you would like to separate by species.

**NOTE**: k-mers are stored in a suffix tree, which can become very large. Because of this, only one species' k-mers are counted at a time. Small k-mer sizes can help curtail memory usage (but will also become more likely to become non-species-specific as sizes decrease), but you should expect k-mer counting runs to take a lot of RAM (potentially 60-80 GB or more). In our experience, a k of 20 is sufficient for distinguishing between distantly related species (i.e. human and mouse), but a higher number (30 or 35) is preferable when separating more closely-related species.

You should choose the lowest value of k that gives good results. To check how well a run worked, you can [plot](#plotting) the results and see how the cells cluster.

### Counting reads on a composite reference genome
If you would rather use a composite reference genome mapping (i.e. if you only have scATAC-seq data and cannot use the transcriptomic k-mer counting method), you can use the program `utils/composite_bam2counts` to create a counts table in the format expected by `demux_species`. Run it like this:
```
utils/composite_bam2counts -b [bamfile] -o [output_directory] (-s [separator] -e [end])
```
Where `[bamfile]` is your data aligned to a composite reference genome and `-o` is the directory where the counts table and species names will be written. This assumes that the composite reference genome has the name of each species appended or prepended to the beginning or end of each chromosome/scaffold name. By default, it's assumed that species names are prepended to sequence names, separated by an underscore (`_`), but you can change the delimiter with the `-s` option and specify that species names are appended to the end, rather than the beginning, of sequence names, using the `-e` option.

After this runs, you can run `demux_species` with `-o` set to the `--output_directory` you gave above. It will then load the counts and proceed as usual. You can either use it to just assign species (it will create `species.assignments` and `species.filt.assignments` files in the output directory), or to assign species and demultiplex reads (if you give it read files to demultiplex).

One thing to be cautious of when using data mapped to a composite reference genome is that misassemblies specific to one genome can attract repetitive reads from all species -- i.e. if there is a long run of G bases in one but not other genomes, then all low-quality reads containing runs of Gs risk aligning to that region and inflating the counts for that species. To avoid this, you should try masking low-complexity regions in your composite reference by using a tool like [GenMap](https://github.com/cpockrandt/genmap). 

### Plotting
<p>
<img src="../img/species.png" width=600 alt="demux_species plot" />
</p>

Run `plot/species.R [output_directory]` after a `demux_species` run to create a summary plot.

The left panel contains a heatmap of species-specific k-mer counts in cells, clustered by cell-cell similarity and with the inferred identity of each cell shown as color bars on the left side. If the color bars do not line up well with groups of cells, something may have gone wrong with assigning cells. The right panel shows the proportion of cells assigned each identity, both in the full set and the filtered set of cell barcodes. The unfiltered/left bar plot represents what will go into the separated FASTQ files (filtering will be done later, during QC after alignment). The filtered/right bar plot is likely a better reflection of the true proportions of species in the pool.

[Back to main README](../README.md)
