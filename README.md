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
