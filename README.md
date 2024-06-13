# CellBouncer
Tools for keeping the riffraff out of pooled single cell sequencing experiments. 

|I want to...|I have...|Tool to use|
|------------|---------|-----------|
|Demultiplex species|A transcriptome (FASTA) or annotation (GTF) and genome (FASTA) per species|`demux_species`|
|Demultiplex individuals|Nothing|`demux_mt`|
|Demultiplex individuals|VCF of known variants|`demux_vcf`|
|Demultiplex individuals|MULTIseq/HTO/CITE-seq data|`demux_tags`|
|Assign sgRNAs|sgRNA capture data|`demux_tags`|
|Quantify ambient RNA|VCF of known variants|`quant_contam`|
|Infer global doublet rate|Output from a `CellBouncer` program|`doublet_dragon`|

CellBouncer can:
* Demultiplex cells 
  * by species
  * by mitochondrial haplotype
  * by individual of origin
  * by cell hashing identity
  * using sgRNA capture data
* Quantify ambient RNA
  * Infers proportion per cell
  * Infers the fraction of ambient RNA originating from each individual
* Infer a global doublet rate
  * Use results from one or more demultiplexing types
<p>
<img src="logo.png", width=300, alt="CellBouncer" />
</p>
