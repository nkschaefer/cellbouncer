# utils/refine_vcf

If you have an `.assignments` file mapping cell barcodes to identities, as well as a BAM file of single-cell sequencing data and a VCF file of genotype data per individual (the individual names in the VCF must match those in the `.assignments` file), you can run this tool to refine the genotype calls in the VCF to maximize the likelihood of the allele counts in the BAM, given the cell identities. This is sort of the reverse thing that [`demux_vcf`](demux_vcf.md) does. 

To run, just do this:

```
utils/refine_vcf -b [bamfile] -a [assignments] -v [vcf] (-T [num_threads])
```
Where `[bamfile]` is your aligned single-cell sequencing data, `[assignments]` is the `.assignments` file containing cell identities, and `[vcf]` is the VCF file of genotype data (this can optionally be bgzipped or compressed to bcf).

The optional `-T` argument controls the number of threads used for parallel processing (default = 1).

This program outputs the refined variants to `stdout` as it goes. Variants are written in gzipped VCF format, but there may be a future option to write BCF instead.

[Back to main README](../README.md)
