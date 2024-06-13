# demux_mt
<p>
<img src="../img/mito.png" width="200px" alt="demux_mt" />
</p>

## Overview
This program attempts to identify the individual of origin of each cell in a pool by inferring the mixture of mitochondrial haplotypes. It does not require any prior knowledge of variant sites or genotypes. Here are some potential use cases:

* Individuals were pooled together and sequenced with no prior knowledge of genotypes, and you want to quickly learn whether they all made it into the data set and whether they are fairly evenly represented
* Individuals that have already been sequenced were pooled together, and `demux_vcf` was used to infer the individual of origin of each cell. To ensure that this worked correctly, you want to have a second, independent identification of each cell to compare to.
* You would like to know something about the mitochondrial haplotypes present in the cells
* You want to obtain genomic variants for the individuals in your pool, without prior knowledge of genotypes of the individuals. You can use this program to get an idea of which cells come from which individuals. With this first round of cell-to-individual assignments, you can use the program `utils/bam_indiv_rg` to add read groups to the BAM that mark which individual a cell comes from, and then run a variant caller to find genomic SNPs that segregate among the individuals. If desired, you can then run `demux_vcf` using these variants, to get more confident assignments for more cells. 

## Pitfalls

* [NUMTs](https://en.wikipedia.org/wiki/Nuclear_mitochondrial_DNA_segment) are very common in some species and will attract reads away from the mitochondrial sequence, making it hard for this program to work. NUMTs are more of a problem in [non-human primates](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7671390/) than in humans. To deal with them, we recommend an approach similar to that described [here](https://github.com/caleblareau/mitoblacklist/): simulate mitochondrial reads, map to the reference genome with the mitochondrial sequence removed, run a peak caller like [MACS2](https://github.com/macs3-project/MACS) to find candidate NUMTs, and use the resulting BED file to mask the reference genome using a tool like [bedtools maskfasta](https://gensoft.pasteur.fr/docs/bedtools/2.29.2/content/tools/maskfasta.html). This must be done before aligning your single-cell sequencing data.
* To confidently identify the mitochondrial haplotype of a cell, there must be a large number of reads in that cell that align to the mitochondrial sequence. Unfortunately, high numbers of mitochondrial reads are also often a [red flag](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html) during quality control in single cell sequencing analysis. This means that the cells most confidently assigned to an individual of origin using `demux_mt` are also often low-quality cells that will later be filtered out. Because of this, if you need to identify individuals of origin and do not have variant call data, we recommend using `demux_mt`, inserting read groups to mark these inferred individuals in the BAM file, running a variant caller to genotype these individuals at genomic variants, and then running `demux_vcf` to re-identify your cells using these genomic variants.

## Running the program

* Run `demux_mt` 
  
[Back to main README](../README.md)
