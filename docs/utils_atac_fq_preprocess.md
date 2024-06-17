# utils/atac_fq_preprocess
`utils/atac_fq_preprocess` is a helper program for users who want to align scATAC data without using a complicated pipeline. It takes as input raw FASTQs of scATAC-seq data, along with a cell barcode whitelist (or two whitelists, in the case of 10X Genomics multiome data). It outputs copies of the input FASTQs with validated cell barcodes appended as sequence comments, which can be mapped to a reference genome using an aligner capable of interpreting sequence comments as SAM/BAM tags.

## File name and order
Single cell ATAC-seq FASTQs come in groups of three: files for a given library might look like `Lib1_S1_L001_R1_001.fastq.gz`, `Lib1_S1_L001_R2_001.fastq.gz`, and `Lib1_S1_L001_R3_001.fastq.gz`. The first file (`R1`) contains the forward read (which contains genomic sequence), the second file (`R2`) is the barcode read (and contains no genomic sequence), and the third file (`R3`) is the reverse read (which contains genomic sequence). When mapping to a reference genome, however, it is typical to align two files: forward and reverse reads of genomic sequence. For this reason, `utils/atac_fq_preprocess` takes three files as input but outputs two: the data in the output `R1` file matches input `R1`, but the data in output `R2` matches input `R3`. Input `R2` is used only to extract cell barcodes, which are inserted into sequence comments in output `R1` and `R2`. 

## Cell barcode lists
`utils/atac_fq_preprocess` requires a [list of valid cell barcodes](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist) against which to check barcode sequences in reads. It allows up to one mismatch, as is the case with other popular single cell programs like `cellranger` and `STARsolo`. In the case of 10X Genomics [multiome](https://www.10xgenomics.com/support/single-cell-multiome-atac-plus-gene-expression) data, however, two barcode lists are required. This is because different barcodes are used for ATAC-seq and RNA-seq in these kits, but every ATAC-seq barcode has a corresponding RNA-seq barcode (identified by being on the same number line in the barcode list file) which is reported in output data. 

## How to run
To run, first create an output directory in which to store the processed reads. Do not use the directory where input reads are stored as the output directory, because `utils/atac_fq_preprocess` creates output FASTQ files with the same names as input FASTQ files. Once you have an output directory, you can run as follows:

```
utils/atac_fq_preprocess \
  -1 [input]_R1_001.fastq.gz \
  -2 [input]_R2_001.fastq.gz \
  -3 [input]_R3_001.fastq.gz \
  -o [output_directory] \
  -w [path_to_barcode_whitelist \
  (-W [path_to_ATAC_whitelist])
```
Where the `-W` argument is only required in the case of multiome data; in this case `-w` should be the multiome RNA-seq barcode list.

## What to do next
With reads processed, you can now align your ATAC data to a reference genome using an aligner that can insert sequence comments into SAM-format output as tags. We recommend [`minimap2`](https://github.com/lh3/minimap2), which can output SAM format with the `-a` option and insert sequence comments as tags with the `-y` option enabled. Remember to pipe the output to [`samtools`](https://github.com/samtools/samtools) to sort and compress to [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf):

```
minimap2 -a -y -x sr [genome.fa|genome.idx] [output]_R1_001.fastq.gz [output]_R2_001.fastq.gz | samtools sort -o [output2].bam
samtools index [output2].bam
```
where `[genome.fa|genome.idx]` is either a FASTA-format reference genome or a `minimap2` index you built for one, `[output]` signifies the output of your run of `utils/atac_fq_preprocess`, and `[output2]` is the name you choose for your sorted BAM file.

### Calling peaks and generating fragment files
If you need to obtain a `fragments` file from your BAM (to use with single cell ATAC analysis software like [ArchR](https://www.archrproject.com/)), we recommend the program [`sinto`](https://timoast.github.io/sinto/basic_usage.html).

If you plan on calling peaks from your BAM file, we recommend [Genrich](https://github.com/jsh58/Genrich), although this will require you to re-sort the BAM by name instead of coordinate, again using `samtools`.

[Back to main README](../README.md)
