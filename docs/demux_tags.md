# demux_tags
<p>
<img src="../img/demux_tags.png" width=150, alt="demux_tags" />
</p>

Some experiments involve labeling cells with barcode sequences that have a custom meaning. In these data sets, a set of reads is collected that can be used to map cell barcodes to these "tag" barcodes, and then counts of each "tag" barcode per cell can be used to assign the cell an identity from these tags. Examples of these types of experiments include [MULTIseq](https://www.nature.com/articles/s41592-019-0433-8), [cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1), [CITE-seq](https://emea.illumina.com/techniques/sequencing/rna-sequencing/cite-seq.html), and [10X Feature Barcoding](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-what-is-feature-bc). In single cell CRISPR experiments, [sgRNA capture](https://www.nature.com/articles/s41587-020-0470-y) data can be processed in a similar way. For these purposes, the program `demux_tags` can count occurrences of barcodes in reads and then assign cell barcodes to identities based on these counts.

## Basic usage
To run `demux_tags`, 
