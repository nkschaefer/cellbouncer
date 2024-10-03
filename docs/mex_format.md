# Working with MEX-format input files

The program [`demux_tags`](demux_tags.md) can load tag or sgRNA counts per cell computed by an outside program. 

The program [`quant_contam`](quant_contam.md) can also optionally attempt to remove contamination from single-cell gene expression data.

These data must be in [market exchange format](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices). This format consists of three files: one listing cell barcodes, one listing genes (and potentially other features, like antibody or cell hashing barcodes), and a third listing counts of features per cell barcodes (a matrix file). 

## h5 formats
Some programs, such as [`kite`](https://github.com/pachterlab/kite) from `kallisto/bustools` and [`cellranger`](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-feature-bc-analysis) can also output data in the `.h5` file format. Additionally, `AnnData` objects used by `scanpy` can be stored in [a specific `h5` format](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.write_h5ad.html) as well.

CellBouncer provides a program, `utils/h5tomex.py`, that can convert `.h5` format output files to market exchange format. To do this, simply run

```
utils/h5tomex.py -H [file.h5] -o [output_directory] (--list) (--feature_type) (--scanpy)
```

where:
* `[file.h5]` is the input file
* `[output_directory]` is a new directory that will contain `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` after the program runs
* `--scanpy/-s` is an optional flag to set if the data is an `.h5ad` file saved by `scanpy/Anndata` rather than an `.h5` file from a data processing program
* For CellRanger `.h5` files, the matrix will often contain multiple data types. If you have processed feature barcoding/sgRNA capture together with gene expression data, for example, both will be output together in this file. Two other options help with this scenario:
  * `--list/-l` lists all data types in the `.h5` file and quits
  * `--feature_type/-t` lets you subset the data to only the type you specify here. The value should be in double quotes and should match one of the options output when you run with the `--list` option. For example, if your data contains both RNA-seq and sgRNA capture, running with `--list` will show:

    ```
    CRISPR Guide Capture
    Gene Expression
    ```
    And you can run again with `--feature_type "CRISPR Guide Capture"` to include only sgRNA capture data in the output files, which will be loaded by `demux_tags`.

### Multiple data types
`demux_tags` needs to know which pieces of the `MEX`-format data represent tag or sgRNA capture data. The same is true for `quant_contam`. If you converted an `h5` file from CellRanger and did not filter it to a specific data type (see `--feature_type` option above), it likely contains gene expression as well as other data types. If you're not sure what data types exist in the file, you can try this:

```
zcat [output_directory]/features.tsv.gz | cut -f3 | sort | uniq
```

If the data came from 10X and there are multiple data types, you should see `CRISPR Guide Capture` in the list (for sgRNA capture data), `Multiplexing Capture` (for cell hashing data), `Antibody Capture` (for antibody capture data), or `Custom` (for general feature barcoding data). For 10X Genomics data sets containing multiple data types, RNA-seq is usually called "Gene Expression."

It is important to specify the name of your custom data type of interest when loading the MEX format data with `demux_tags` by passing the name to `--feature_type/-t` with quotes; for example:

```
demux_tags ... -M [output_dir]/matrix.mtx.gz \
  -B [output_dir]/barcodes.tsv.gz \
  -F [output_dir]/features.tsv.gz \
  -t "CRISPR guide capture"
```
The same `-t` argument also exists for `quant_contam` when loading gene expression data (with `-B`, `-F`, and `-M` arguments).

If you already filtered your data, or the `.h5` file only contains the data type of interest, you can omit the `--feature_type/-t` argument to `demux_tags` and `quant_contam`.

## Subsetting to a data type or barcode list

CellBouncer has two programs to subset a MEX-format data set:

`utils/subs_mex_bc.py` will subset a MEX data set given a list of barcodes to include. To run:
```
utils/subs_mex_bc.py -m [matrix] -f [features] -b [barcodes] -B [barcodes_filter] -o [output]
```
where `-m`, `-f`, and `-b` are input files, `-B` is a file listing barcodes to retain (one per line), and `-o` is the desired output directory (it will create `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` within this directory).

`utils/subs_mex_featuretype.py` will subset a MEX data set to the data set of interest. To run:
```
utils/subs_mex_featuretype.py -m [matrix] -f [features] -b [barcodes] -t [feature_type] -o [output]
```
where `-m`, `-f`, and `-b` are input files, `-t` is the feature type to retain (see above, *e.g.* "Gene Expression"), and `-o` is the output directory.

If you subset your data using `utils_subs_mex_featuretype.py` so that it only contains gene expression data, then you can omit the `-t` argument to `demux_tags` and `quant_contam` when loading it.

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
