# Installation notes

## For show-offs
The non-plotting programs in `cellbouncer` require only [HTSLib](https://github.com/samtools/htslib) and [zlib](https://www.zlib.net/). If you want to be cocky about it, you can just install those two dependencies and then clone and compile this repository:
```
git clone --recurse-submodules git@github.com:nkschaefer/cellbouncer.git
cd cellbouncer
make
```
You will not be able to run the plotting scripts unless you also install all needed R packages.

## With conda
[`conda`](https://github.com/conda-forge/miniforge) is a useful package manager that can be used to install CellBouncer dependencies. We also recommend installing [`mamba`](https://anaconda.org/conda-forge/mamba), which is an add-on for `conda` that is able to locate packages more quickly.

We provide several configuration files that can be used to create `conda` environments for use with CellBouncer. The file `cellbouncer_minimum.yml` contains all dependencies you need to run `cellbouncer` programs and make plots. The file `cellbouncer_extra.yml` also includes other programs useful for processing/analyzing data outside of `cellbouncer` that are mentioned in these help pages. 

**Note for Mac users**: The file `cellbouncer_extra_osx.yml` is a version of `cellbouncer_extra.yml` that omits the program [`FastK`](https://github.com/thegenemyers/FASTK), which is not available through `conda` for Mac OSX. If you want to use that program on a Mac, you will need to install it through other means. Additionally, M1 Macs can run programs compiled for both x86_64 and arm64 architectures (and conda doesn't care which by default), but attempting to link against programs compiled for one architecture while compiling for the other will cause problems. Therefore, you need to set a special variable before creating the conda environment if using an M1 Mac (see below).

Installing a conda environment
```
(CONDA_SUBDIR=osx-arm64) [conda/mamba] env create --file=[cellbouncer_minimum/cellbouncer_extra/cellbouncer_extra_osx].yml
conda activate cellbouncer
```
The first part (`CONDA_SUBDIR=osx-arm64`) should only be included if you are installing on a silicon M1+ Mac.

The environment you choose should either be `cellbouncer_minimum.yml` if you only want necessary dependencies, `cellbouncer_extra.yml` if you want to include extra programs, or `cellbouncer_extra_osx.yml` if you want to include extra programs and are on a Mac.

## Compiling
The only thing you need to do to compile is run `make`, and all compiled programs will be in the root-level directory. In the root (`cellbouncer`) directory:
```
make
```
or, if you have compiled old code and are updating:

```
make clean
make
```

## Running
If you have installed a `conda` environment, be sure to run `conda activate cellbouncer` to activate the environment before you run any CellBouncer programs.

## Note about git submodules

`CellBouncer` depends on several other repositories included as git submodules. If you forget the `--recurse-submodules` option in your `git pull` command above, you can get all the submodules with
```
git submodule update --init --recursive
```
If you later need to update a local `CellBouncer` to the latest version, you can do so like this:
```
git pull --recurse-submodules
git submodule update --remote
make clean
make
```

[Back to main README](../README.md)
