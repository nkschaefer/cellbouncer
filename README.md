<p>
<img src="img/logo.png", width=350, alt="CellBouncer" />
</p>

Tools for checking cell identities and keeping the riffraff out of pooled single cell sequencing data sets. 

| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |I want to...|I have...|Tool to use|
|-|------------|---------|-----------|
| <img src="img/demux_species_mini.png" alt="demux_species" /> | Demultiplex cells by **species** |Raw reads, plus a transcriptome (FASTA) or annotation (GTF) and genome (FASTA) per species <p align="center">*OR*</p> A BAM file of reads mapped to a composite reference genome|[`demux_species`](#demux_species)|
| <img src="img/demux_mt_mini.png" alt="demux_mt" /> | Demultiplex cells by **individual of origin**, and I hope individuals are unrelated enough to have different mitochondrial haplotypes|A BAM file of aligned scATAC-seq or whole cell scRNA-seq data|[`demux_mt`](#demux_mt)|
| <img src="img/demux_vcf_mini.png" alt="demux_vcf" /> | Demultiplex cells by **individual of origin**|VCF of known variants, plus a BAM file of aligned single cell sequencing data|[`demux_vcf`](#demux_vcf)|
| <img src="img/demux_tags_mini.png" alt="demux_tags" /> | Demultiplex individuals by **custom label** or **treatment**|FASTQs containing MULTIseq/HTO/CITE-seq data, or a table of pre-computed counts, optionally in [MEX format](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices)|[`demux_tags`](#demux_tags)|
| <img src="img/demux_tags_mini.png" alt="demux_tags" /> | Assign **sgRNAs** to cells|FASTQs containing sgRNA capture data, or a table of pre-computed counts, optionally in [MEX format](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices)|[`demux_tags`](#demux_tags)|
| <img src="img/quant_contam_mini.png" alt="quant_contam" /> | Quantify **ambient RNA** per cell, infer its origins, and optionally adjust gene counts|Output from `demux_vcf` (plus optional single-cell expression data to adjust, in [MEX format](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices))|[`quant_contam`](#quant_contam)|
| <img src="img/doublet_dragon_mini.png" alt="doublet_dragon" /> | Infer global **doublet rate** and proportions of individuals|Output from one or more `CellBouncer` programs run on the same cells|[`doublet_dragon`](#doublet_dragon)|
| <img src="img/bulkprops_mini.png" alt="bulkprops" /> | Determine **proportion of individuals** in a pool|A VCF of known variants, plus a BAM of aligned sequence data (can be bulk)|[`bulkprops`](#bulkprops)|

<br>
<p >
 &nbsp;&nbsp;<img src="img/viz_compare3.png" width=50 />
</p>

### Visualizing and comparing results

|I want to...|I have...|Tool to use|
|------------|---------|-----------|
| **Visualize** a set of labels and the pool compositions they produce at different confidence cutoffs | An `.assignments` file from a CellBouncer program |[`plot/assignment_llr.R`](docs/plot_assignment_llr.md)|
| **Compare** two sets of labels on the same cells | Two `.assignments` files from CellBouncer programs run on the same data | [`plot/compare_assignments.R`](docs/plot_compare_assignments.md)|
| **Merge** two sets of labels on the same cells into one set of labels | Two `.assignments` files from CollBouncer programs run on the same data | [`utils/merge_assignments.R`](docs/utils_merge_assignments.md)
| **Compare** two sets of pool proportions and assess significance if possible | Two files describing pool composition (i.e. from [`bulkprops`](#bulkprops) or contamination profile from [`quant_contam`](#quant_contam)), or one file describing pool composition and an `.assignments` file describing cell labels | [`utils/compare_props.R`](docs/utils_compare_props.md) |
| **Refine** genotype calls to better match cell-individual labels | A preexisting set of genotypes in VCF format, a BAM file of aligned single-cell data, and an `.assignments` file mapping cells to individuals of origin | [`utils/refine_vcf`](docs/utils_refine_vcf.md) |

<br>
<p>
<img src="img/manipulate.png" width=50 />
</p>

### Manipulating input files

|I want to...|I have...|Tool to use|
|------------|---------|-----------|
|**Split** a BAM file into one file per cell identity | A BAM file of aligned single-cell sequencing data and a CellBouncer-format `.assignments` file | [`utils/bam_split_bcs`](docs/utils_bam_split_bcs.md) |
|**Tag** reads in a BAM file to mark individual of origin | A BAM file of aligned single-cell sequencing data and a CellBouncer-format `.assignments` file | [`utils/bam_indiv_rg`](docs/demux_mt.md#utilsbam_indiv_rg) |
|**Convert** 10X or Scanpy (AnnData) data from `.h5` to [MEX](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices) format | A CellRanger-format [`.h5`](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices) or Scanpy-format [`.h5ad`](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) file | [`utils/h5tomex.py`](docs/mex_format.md#h5-formats) |
|**Split** MEX-format data into one data set per library/run | Single-cell expression data in [MEX format](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices) | [`utils/split_mex_libs.py`](docs/mex_format.md#splitting-libraries)  |
|**Subset** MEX-format data to specific cell barcodes | Single-cell expression data in [MEX format](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices) | [`utils/subs_mex_bc.py`](docs/mex_format.md#subsetting-to-a-data-type-or-barcode-list)  |
|**Subset** MEX-format data to a specific feature type | Single-cell expression data in [MEX format](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices) | [`utils/subs_mex_featuretype.py`](docs/mex_format.md#subsetting-to-a-data-type-or-barcode-list) |

# Installation
To install ([see below](#get-the-repository)):

* Clone the repository (and its submodules)
* Choose a `conda` environment file to install
  * All necessary dependencies: `cellbouncer_minimum.yml`
  * All necessary dependencies plus extra helper programs mentioned in documentation:
    * Mac OS X: `cellbouncer_extra_osx.yml`
    * Linux: `cellbouncer_extra.yml`
* Run `make`

For more information about installing or updating CellBouncer, see [here](docs/installation_notes.md).

### Get the repository
```
git clone --recurse-submodules git@github.com:nkschaefer/cellbouncer.git
cd cellbouncer
```
### Create conda environment
#### Linux
```
[conda/mamba] env create --file=[cellbouncer_minimum/cellbouncer_extra].yml
conda activate cellbouncer
conda env config vars set LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}" -n cellbouncer
```
#### Mac OS X (M1)
```
CONDA_SUBDIR=osx-arm64 [conda/mamba] env create --file=[cellbouncer_minimum/cellbouncer_extra_osx].yml
conda activate cellbouncer
conda env config vars set DYLD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${DYLD_LIBRARY_PATH}" -n cellbouncer
```
#### Mac OS X (Intel)
```
[conda/mamba] env create --file=[cellbouncer_minimum/cellbouncer_extra_osx].yml
conda activate cellbouncer
conda env config vars set DYLD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${DYLD_LIBRARY_PATH}" -n cellbouncer
```
### Compile
```
make
```
You've now got all the programs compiled, and you can run them as long as you remember to `conda activate cellbouncer` first.

# Test data set
You can get a test data set from [this link](https://ucsf.box.com/s/wvrzl1tvdrozojlj0z4sp5fhmvqe05mu). It contains an example `.bam` file, `vcf` file, and cell hashing `.counts` file. The `README` in the linked directory will explain everything, but these give you the opportunity to test-run `demux_species`, `demux_mt`, `demux_tags`, `demux_vcf`, `quant_contam`, and `doublet_dragon`.

# Overview
The programs in `cellbouncer` are standalone command line tools. If you run one of them with no arguments or with `-h`, it will give you detailed information about how to run it. Each program uses the concept of an `--output_prefix/-o`, which is a base name that will be used for all output files.

## Output files
Demultiplexing tools all write a file called `[output_prefix].assignments`, which tells you information about each cell's identity. These files are 4 columns, tab separated: 
* cell barcode (optionally with unique ID appended; see below)
* most likely identity (doublets are two names in alphabetical order separated by `+`)
* droplet type: `S` (for singlet), `D` (for doublet), or in some cases `M` (for multiplet, 3+ individuals, so far only considered by `demux_tags`)
* ratio of the log likelihood of the best to the second best assignment (a measure of confidence in the assignment)

## Cell barcode format and merging with other data
To load data from CellBouncer into a single cell analysis tool like [Seurat](https://satijalab.org/seurat/) or [scanpy](https://scanpy.readthedocs.io/en/stable/), you will need to load CellBouncer's (text format) output files and merge with your single cell data set. This requires ensuring that cell barcodes are formatted the same way by CellBouncer as in your data set. Read more [here](docs/merging.md). 

## Note about multi-library data sets
`CellBouncer` programs take input from single libraries. If you have concatenated multiple single cell sequencing data sets, `CellBouncer` will interpret all cells with the same barcode sequence as the same cell, ignoring any unique IDs you have added to barcodes. If you need to load data in [MEX format](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices) (i.e. for `demux_tags` or `quant_contam`) and that data comes from multiple libraries that were concatenated together, you can separate the data by library using the program [`utils/split_mex_libs.py`](docs/mex_format.md#splitting-libraries).

## Acknowledgments
We thank [Helena Pinheiro](https://www.hpinheiro.com/) for creating the icons for `demux_species`, `demux_mt`, `demux_vcf`, `demux_tags`, `quant_contam`, and `bulkprops`, and the cell drawing used throughout.
