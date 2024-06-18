# utils/merge_assignments.R
This program takes as input two `.assignments` files describing the same labels on the same data. For example, if you have assigned the same cells to the same individual identities using both a VCF (`demux_vcf`) and mitochondrial haplotypes (`demux_mt`), and the individuals have matching names, then this program outputs consensus assignments for all cells.

## How it makes decisions
* If a cell is missing from one but not the other file, the identity is taken from the file where it is not missing.
* If a cell is present in both files, and it has a higher log likelihood ratio of assignment in one file, then the assignment is taken from the file with the higher LLR.
* If a cell is present in both files and has different assignments in both, and the LLR is the same in both files, then the cell is *not* output in the consensus file.

## How to run
Just run 
```
utils/merge_assignments.R [file1.assignments] [file2.assignments] > [consensus.assignments]
```
Where `[consensus.assignments]` is the new file to create (the program writes to `stdout`).

[Back to main README](../README.md)
