# utils/refine_vcf

If you have an `.assignments` file mapping cell barcodes to identities, as well as a BAM file of single-cell sequencing data and a VCF file of genotype data per individual (the individual names in the VCF must match those in the `.assignments` file), you can run this tool to refine the genotype calls in the VCF to maximize the likelihood of the allele counts in the BAM, given the cell identities. This is sort of the reverse thing that [`demux_vcf`](demux_vcf.md) does. 

To run, just do this:

```
utils/refine_vcf -b [bamfile] -a [assignments] -v [vcf] (-j)
```
Where `[bamfile]` is your aligned single-cell sequencing data, `[assignments]` is the `.assignments` file containing cell identities, and `[vcf]` is the VCF file of genotype data (this can optionally be bgzipped or compressed to bcf).

The `-j` argument is only useful when there are very few variants: it index-jumps through the BAM file instead of streaming it. This can speed things up when there aren't many variants to look up, but it will slow things down when there are many.

This program outputs the refined variants to `stdout` as it goes. 

[Back to main README](../README.md)
