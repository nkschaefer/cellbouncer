# MEX format, converting to it, and manipulating it

The program [`demux_tags`](demux_tags.md) can load tag or sgRNA counts per cell computed by an outside program. These data must be in [market exchange format](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices), however. 

## h5 formats
Some programs, such as [`kite`](https://github.com/pachterlab/kite) from `kallisto/bustools` and [`cellranger`](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-feature-bc-analysis) can also output data in the `.h5` file format. Additionally, `AnnData` objects used by `scanpy` can be stored in [a specific `h5` format](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.write_h5ad.html) as well.

CellBouncer provides a program, `utils/h5tomex.py`, that can convert `.h5` format output files to market exchange format. To do this, simply run

```
utils/h5tomex.py -H [file.h5] -o [output_directory] (--scanpy)
```

where `[file.h5]` is the input file, `[output_directory]` is a new directory that will contain `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` after the program runs, and `--scanpy` is an optional flag to set if the data is an `.h5ad` file saved by `scanpy/Anndata` rather than an `.h5` file from a data processing program.

### Note
`demux_tags` needs to know which pieces of the `MEX`-format data represent tag or sgRNA capture data. If you converted an `h5` file from CellRanger, it likely contains gene expression as well as other data types. If you're not sure what data types exist in the file, you can try this:

```
zcat [output_directory]/features.tsv.gz | cut -f3 | sort | uniq
```

If the data came from 10X and there are multiple data types, you should see `CRISPR Guide Capture` in the list (for sgRNA capture data), `Multiplexing Capture` (for cell hashing data), `Antibody Capture` (for antibody capture data), or `Custom` (for general feature barcoding data). 

It is important to specify the name of your custom data type of interest when loading the MEX format data with `demux_tags` by passing the name to `--feature_type/-t` with quotes; for example:

```
demux_tags ... -M [output_dir]/matrix.mtx.gz \
  -B [output_dir]/barcodes.tsv.gz \
  -F [output_dir]/features.tsv.gz \
  -t "CRISPR guide capture"
```

## Splitting libraries

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
