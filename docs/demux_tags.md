# demux_tags
<p>
<img src="../img/demux_tags.png" width=150, alt="demux_tags" />
</p>

Some experiments involve labeling cells with barcode sequences that have a custom meaning. In these data sets, a set of reads is collected that can be used to map cell barcodes to these "tag" barcodes, and then counts of each "tag" barcode per cell can be used to assign the cell an identity from these tags. Examples of these types of experiments include [MULTIseq](https://www.nature.com/articles/s41592-019-0433-8), [cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1), [CITE-seq](https://emea.illumina.com/techniques/sequencing/rna-sequencing/cite-seq.html), and [10X Feature Barcoding](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-what-is-feature-bc). In single cell CRISPR experiments, [sgRNA capture](https://www.nature.com/articles/s41587-020-0470-y) data can be processed in a similar way. For these purposes, the program `demux_tags` can count occurrences of barcodes in reads and then assign cell barcodes to identities based on these counts.

## Basic usage
`demux_tags` can be run several different ways.

### Starting from FASTQs
If you have FASTQs from one of the above data types, `demux_tags` can count the occurrences of each feature (tag or sgRNA) in each cell barcode and then assign cells to identities. To run in this mode, you must provide the following arguments:
```
--read1 -1 Forward read FASTQ file(s) (can specify multiple times)
--read2 -2 Reverse read FASTQ file(s) (can specify multiple times)
--whitelist -w List of allowed cell barcode sequences
--seqs -s File mapping sequences to identities
```
**Notes**
* The `--read1` and `--read2` argument names must be provided before each file. For example, to provide two read file pairs `Lib1_S1_L001_R1_001.fastq.gz, Lib1_S1_L001_R2_001.fastq.gz` and `Lib2_S2_L001_R1_001.fastq.gz, Lib2_S2_L001_R2_001.fastq.gz`, pass them like this:
  ```
  --read1 Lib1_S1_L001_R1_001.fastq.gz --read2 Lib1_S1_L001_R2_001.fastq.gz --read1 Lib2_S2_L001_R1_001.fastq.gz --read2 Lib2_S2_L001_R2_001.fastq.gz
  ```
* The `--whitelist` argument is specific to the droplet-based sequencing technology used; for more information on the 10X Genomics barcode lists, see [here](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist)
* The `--seqs` file should be tab-separated, with one sequence and one identity per line, i.e.:
  ```
  ACGTGGAGCTTG  Name1
  GTGGACGTGAGT  Name2
  TTGAGGTGAGGT  Name3
  ```

