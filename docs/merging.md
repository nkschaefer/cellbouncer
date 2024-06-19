# Cell barcode format and merging with other data sets

## Overview
After running one or more `CellBouncer` programs, you will probably want to load the results into your single cell analysis tool of choice. Such tools have slightly different conventions for representing data, but at their core, they store (sparse) matrices of single cell count or expression data, along with one table of metadata for features (or genes) and another table of metadata for cells. 

`CellBouncer` output should be merged with cell metadata in a way that does not discard cells that are missing from `CellBouncer` output files. The first thing to do is ensure that cell barcodes written by `CellBouncer` are the same format as those in your other data set.

## Cell barcode format
In output files, cell barcodes will by default be printed without any additional text (they will consist only of DNA sequences). When you plan to combine data from multiple single-cell libraries, however, additional text must be appended to barcodes from each data set to prevent barcode collisions. 

Different programs have different conventions for handling this:
* [CellRanger](https://www.10xgenomics.com/support/software/cell-ranger/latest) appends `-` and a numeric ID (starting from 1) to cell barcodes.
  * If you are not using `cellranger aggr` to combine libraries, then this will result in `-1` being at the end of all barcodes.
* In [scanpy](https://scanpy.readthedocs.io/en/stable/), the [anndata.concatenate](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.concatenate.html) command also follows this convention, unless the `batch_categories` argument is used. In this case, names passed to `batch_categories` are appended to barcode sequences, separated by `-`. The `index_unique` argument can also change `-` to another character.
* [Seurat](https://satijalab.org/seurat/archive/v4.3/merge) appends unique IDs given by the `add_cell_ids` argument to the `merge` function, but these are prepended to barcode sequences, separated by `_`. If loading CellRanger data, this may result in barcodes that begin with unique IDs and end with `-1`.

`CellBouncer` programs that write files containing cell barcodes have a `--libname/-n` argument that allows users to append unique identifiers to barcode sequences, by default at the end and with the separator `-`.
  * To mimic multi-library `cellranger aggr` output, or `anndata.concatenate` output without the `batch_categories` argument, set `--libname` to the 1-based numeric index of this library among the set of all being concatenated.
  * To mimic `anndata.concatenate` output with `batch_categories` set, set `--libname` to the unique name for this library in `batch_categories`.
  * To mimic `Seurat`, use the `--seurat/-S` argument to instead append `--libname` to the beginning of each cell barcode, with the separator `_`.
  * To use `_` as a separator (at the end of barcode sequences) instead of `-`, use the `--underscore/-U` argument.
  * To mimic single-library `cellranger` output, use the `--cellranger/-C` argument to append `-1` to the end of each barcode sequence. This can be used in conjunction with other options: with `--libname`, output will be `[barcode seq]-1-[libname]`, and with `--seurat`, output will be `[libname]_[barcode seq]-1`.

`CellBouncer` programs each expect data from a single library: extra characters like library identifiers are removed from cell barcode sequences when loading them. This means that if you do provide data from multiple libraries together, any cells that share the same barcode sequence will be treated as the same cell. 

## Merging data
Once you have your file with properly-formatted cell barcodes, you can merge with your single cell data set. Here is some example code for merging an `.assignments` file from a `CellBouncer` program with an exsting single cell data set.

### In scanpy
```
import pandas as pd
import scanpy as sc
import anndata

# Load scanpy project
adata = sc.read_h5ad('[scanpy_saved].h5ad')

# Load CellBouncer assignments file
assn = pd.read_csv('[output_prefix].assignments', \
  sep='\t', index_col=0, header=['individual', 'droplet_type', 'individual_llr'])

# Merge CellBouncer data into cell metadata, without dropping any cells
adata.obs = adata.obs.merge(assn, how='left', left_index=True, right_index=True)
```
`adata.obs` should now contain the columns `individual`, `droplet_type`, and `individual_llr`, with values set to `NA` for all cells missing from the `CellBouncer` output file (which could not be identified). This will not work if your barcode format does not match that in `scanpy` (see above section).

### In Seurat


[Back to main README](../README.md)
