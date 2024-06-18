# utils/split_mex_libs
CellBouncer programs expect to be run on single libraries. One reason for this is that cell barcodes are represented in memory by a numeric hash of their sequences, so multiple unique cells with the same barcode sequence will be treated as the same cell.

If you want to run a program like `demux_tags` that loads [market exchange format](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices) data, but the matrices you have come from multiple libraries (where library names are appended to the cell barcode sequences), this program can split the data into one matrix per library.

To run:
```
utils/split_mex_libs.py -m [input].mtx(.gz) \
  -f [input features/genes].tsv(.gz) \
  -b [input barcodes].tsv(.gz) \
  -o [output_prefix]
```

## Output files

This will detect library names tacked onto cell barcodes and create a new set of files for each library, where the names are:
* `[output_prefix]_[libname].barcodes.tsv.gz`
* `[output_prefix]_[libname].featues.tsv.gz`
* `[output_prefix]_[libname].matrix.mtx.gz`

[Back to main README](../README.md)
