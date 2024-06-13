# CellBouncer
Tools for demultiplexing and keeping the riffraff out of pooled single cell sequencing experiments. 

|I want to...|I have...|Tool to use|
|------------|---------|-----------|
|Demultiplex species|A transcriptome (FASTA) or annotation (GTF) and genome (FASTA) per species|`demux_species`|
|Demultiplex individuals with different mitochondrial haplotypes|Nothing|`demux_mt`|
|Demultiplex individuals|VCF of known variants|`demux_vcf`|
|Demultiplex individuals|MULTIseq/HTO/CITE-seq data|`demux_tags`|
|Assign sgRNAs|sgRNA capture data|`demux_tags`|
|Quantify ambient RNA per cell and infer its origins|VCF of known variants|`quant_contam`|
|Infer global doublet rate|Output from a `CellBouncer` program|`doublet_dragon`|

<p align="center">
<img src="logo.png", width=300, alt="CellBouncer" />
</p>

## Installation
The non-plotting programs in `CellBouncer` require only [HTSLib](https://github.com/samtools/htslib) and [zlib](https://www.zlib.net/). You can install these dependencies yourself, or to make it easier, get these and plotting-related dependencies by creating a [conda](https://github.com/conda-forge/miniforge/releases) environment from the included `cellbouncer.yml` file. Then all that's left to do is `make`:

```
git clone --recurse-submodules https://github.com/nkschaefer/cellbouncer.git
cd cellbouncer
conda env create --file=cellbouncer.yml
conda activate cellbouncer
make
```
You've now got all the programs compiled, and you can run them as long as you remember to `conda activate cellbouncer` first.

