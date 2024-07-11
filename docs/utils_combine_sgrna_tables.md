# utils/combine_sgrna_tables.py
If you have run [`demux_tags`](demux_tags.md) with the `--sgRNA` option to assign sgRNAs to cells, for multiple libraries from the same experiment, you can use this program to combine the resulting `.table` files from each run into one table. For this to work properly, you must have kept cell barcode names unique by supplying a different `--libname` argument to each `demux_tags` run.

### Usage
```
utils/combine_sgrna_tables.py [output_prefix1].table [output_prefix2].table... > [output_combined].table
```
Where `[output_prefix1]`, `[output_prefix2]`, and so on are `--output_prefix` options given to `demux_tags` for each run. You can provide as many input `.table` files as you like. Output will be printed to `stdout`, so to save a file you must use the `>` character and the name of the output file to create.

[Back to main README](../README.md)
