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
--names -N (OPTIONAL) convert second column in the --seqs file to a final identity (see below)
```
**Notes**
* The `--read1` and `--read2` argument names must be provided before each file. For example, to provide two read file pairs `Lib1_S1_L001_R1_001.fastq.gz, Lib1_S1_L001_R2_001.fastq.gz` and `Lib2_S2_L001_R1_001.fastq.gz, Lib2_S2_L001_R2_001.fastq.gz`, pass them like this:
  ```
  --read1 Lib1_S1_L001_R1_001.fastq.gz \
  --read2 Lib1_S1_L001_R2_001.fastq.gz \
  --read1 Lib2_S2_L001_R1_001.fastq.gz \
  --read2 Lib2_S2_L001_R2_001.fastq.gz
  ```
* The `--whitelist` argument is specific to the droplet-based sequencing technology used; for more information on the 10X Genomics barcode lists, see [here](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist)
* The `--seqs` file should be tab-separated, with one sequence and one identity per line, i.e.:
  ```
  ACGTGGAGCTTG  Name1
  GTGGACGTGAGT  Name2
  TTGAGGTGAGGT  Name3
  ```
* You can include an optional file via the `--names/-N` parameter for use cases, such as [MULTIseq](https://pubmed.ncbi.nlm.nih.gov/31209384/), where there is a consistent DNA barcode sequence -> plate well mapping across all experiments, but in each experiment, the meaning attached to each plate well changes. In these types of cases, you can use the same `--seqs` file in every run (which maps MULTIseq barcode sequences to well names), but use a unique `--names/-N` file for every individual run, to assign a unique final meaning to each intermediate name (i.e. well ID).
  * This file should be tab-separated, with two values per line: the intermediate ID, i.e. well name, followed by the final name (user-supplied). For example, the `--seqs` file might contain:
    ```
    TGGTAGCT	A1
    CCGTAGGT	A2
    CGGGAAGA	A3
    GATACGAT  A4
    ```
    and the `--names` file might contain:
    ```
    A1  Treatment1
    A2  Treatment2
    A3  Treatment3
    ```
    * In this example, the same `--seqs` file could be provided for every experiment of this type, but a unique `--names` file would be provided for every individual experiment, where meanings of tags change. One advantage of 
