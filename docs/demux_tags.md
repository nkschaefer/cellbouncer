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
--mismatches -m The number of allowable mismatches to a tag sequence to count it (default 2, -1 = best match with no limit)
--umi_len -u Length of UMI sequences in read 1 (assumed to immediately follow cell barcode), in bp (default 12, cannot exceed 16)
--cell_barcodes -B A filtered list of cell barcodes (i.e. from an RNA-seq alignment tool like STARsolo or CellRanger)
--exact -e Do not allow mismatches in cell barcode sequences (default: allow one mismatch)
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
  * In this example, the same `--seqs` file could be provided for every experiment of this type, but a unique `--names` file would be provided for every individual experiment, where meanings of tags change.
  * In addition to the convenience of reusing one `--seqs` file for multiple experiments, this setup allows the user to detect mixups. In the above example, one well (`A4`) is not used in the experiment and has no identity in the `--names` file. When this is the case, `demux_tags` will automatically check to see whether the sequence tied to `A4` appears more often in the data than any of the named sequences and will alert the user when this is the case.
  * If the `--names` file is omitted, output files from `demux_tags` will use the second column of the `--seqs` file (in this case, well IDs) in output files.
  * In some cases, including sgRNA capture data, there is no intermediate identity: each sgRNA sequence has a unique name specific to the experiment. In those cases, users should omit the `--names/-N` argument.
* If you have already processed other data (i.e. scRNA-seq) from this experiment and have a filtered list of barcodes likely representing true cells (i.e. from CellRanger or STARsolo), passing that list to the `--cell_barcodes/-B` argument will greatly speed up processing by reducing the number of barcodes each read is checked against.
* If you experience slowness counting data, the `--exact` argument can significantly speed up processing, especially if the list of barcodes in the `--whitelist` file is very long. This risks missing some legitimate counts, though. If you have provided a pre-filtered barcode list with `--cell_barcodes/-B`, this argument is likely unnecessary.
* `--umi_len/-u` is set to 12 by default. Older versions of kits often used 10 bp UMIs, however. Be careful when setting this parameter.
     
### Starting from Market exchange format (MEX) output
If users prefer to have another program, such as [10X Genomics CellRanger](https://www.10xgenomics.com/support/software/cell-ranger/latest) or [kallisto/bustools kite](https://github.com/pachterlab/kite) count tags in reads, such programs offer output in [MEX format](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices), among other options. This format usually involves a specific output directory that contains (at least) three files:
* Cell barcodes file (file name should include `barcodes.tsv` or `barcodes.txt` and can end with/without `.gz`)
* Genes or features file (file name should include `features.txt/tsv` or `genes.txt/tsv` and can end with/without `.gz`)
* Matrix file (file name should include `.mtx` and can end with/without `.gz`)

If a user wishes to use `demux_tags` to assign identities from counts produced by another such program, `demux_tags` provides three input options to load these files (described above):
```
--cell_barcodes -B Cell barcodes, optionally gzipped
--features -F Features (sometimes called genes), optionally gzipped
--mtx/-M Matrix, or .mtx file, optionally gzipped
--feature_type -f (OPTIONAL), see below
```
#### Multiple data types together
If a data set includes multiple data types (i.e. a CellRanger run including both gene expression and feature barcoding, or gene expression and sgRNA capture), the MEX-format data will contain all these types, and `demux_tags` must be told which data type to load. For this, the `--feature_type/-f` argument must be used.
* For 10X cell hashing, for example, specify `Multiplexing\ Capture`
* For 10X antibody capture specify `Antibody\ Capture`
* For 10X sgRNA capture, specify `CRISPR\ Guide\ Capture`

### Starting from the output of a previous run
If you already ran `demux_tags` to count tags in reads, you can simply re-run with the same `--output_prefix` argument with no input options (see above), and it will automatically load the data written to `[output_prefix].counts`. This is useful for when you want to change parameters than can affect final assignments, but you have already obtained the count data.

### Assigning cells to identities
A few arguments affect how `demux_tags` will run after counts are computed or loaded.

```
--cell_barcodes -B Filters the list of cell barcodes under consideration. If you counted reads without this argument, it will limit the number of cells assigned identities.
--filt -f Applies an optional filtering step at the end (see below) which can improve precision at the cost of decreasing recall.
--comma -c By default, in all CellBouncer programs, cells identified as doublets receive two names, separated by +. This changes the separator to ,. If you specify --sgRNA, this is activated by default
--sgRNA -s Assume input is sgRNA counts instead of cell hashing/MULTIseq type data (see below)
```

* The `--filt/-f` argument performs an extra filtering step at the end of cell-identity assignment.
  * After all cells have been assigned a most likely identity, `demux_tags` infers (via maximum likelihood estimation) the percent of tag counts for each cell that likely originated from background/ambient tags (these values are output in the `[output_prefix].bg` file).
  * With the `--filt` option set, `demux_tags` then fits a two-way mixture model to the total ambient tag counts per cell. If there is reasonable evidence that this distribution is truly bimodal (the mean of the upper distribution must be at or above the point where the lower distribution's CDF reaches 0.9), then cells likelier to have their total ambient counts fall in the upper distribution are removed from the data set and do not receive assignments.
  * The reasonining behind this is that incorrect assignments should produce unreasonably high numbers of ambient/off-target tag counts compared to correct assignments.

## Output files
`demux_tags` creates several output files, with names defined by `[output_prefix]`, followed by an extension.
* `[output_prefix].assignments` contains the cell-identity assignments, in a similar format to other `cellbouncer` programs. The fourth column, as in other files, gives the log likelihood ratio of the best to second best assignment; this is a measure of confidence and can be used to filter results.
  * Because `demux_tags` considers not only singlet and doublet assignments (as other programs do) but also considers combinations of three or more assignments per cell, the third column by default can take on the values of `S` (singlet), `D` (doublet), or `M` (multiplet, meaning three or more identities), whereas other programs only assign `S` or `D`.
  * If you specified that the data is from sgRNA capture (`--sgRNA` option), there are some slight differences.
    * First, multiple identities in the second column are separated by `,` instead of `+`. This is also the case if you did not set the `--sgRNA` option but set the `--comma` option. This is to avoid confusion in cases where guide names contain the `+` character.
    * Second, to make things easier in cases where multiple guides are expected (i.e. high-MOI CRISPR screens), the third column does not contain `S` (for singlet) or `D` (for doublet) as in other data types, but instead contains the number of guides assigned. Note that cells with zero guides assigned will not appear in the `.assignments` file.
* If you set the `--sgRNA` argument, another file is created called `[output_prefix].table`. This file can become very large, but its purpose is to be helpful in high-MOI CRISPR screens, where every cell can receive multiple guides. This file is a table where the first column is the cell barcode, and all other columns are sgRNAs (named in the first line of the file). Every row below the header row is then a cell barcode, followed by `0` (for not assigned) or `1` (for assigned) in every column, signifying whether or not the cell received that guide RNA.
* The `[output_prefix].bg` file contains maximum likelihood estimates of the percent of tag counts consisting of ambient/background tag counts for every cell. This includes cells that were filtered out of the data set and do not appear in `[output_prefix].assignments`.
  * It may be useful to filter using this file in cases (such as sgRNA capture data) where one wishes to be conservative in assigning cells to identities. This can be used together with the log likelihood ratio of assignment (the fourth column in the `.assignments` file).

[Back to main README](../README.md)
