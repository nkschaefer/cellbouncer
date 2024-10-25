## Building a unique k-mer list for demux_species
If you would like to build a reference k-mer list for demux_species without using the helper tool `utils/demux_species_kmers.py`, here are the steps involved:

* Obtain a transcriptome (in FASTA format) for each species in the pool.
  * If you do not have one, you can create one from a [GTF](http://genome.ucsc.edu/FAQ/FAQformat#format4) or [GFF](http://genome.ucsc.edu/FAQ/FAQformat#format3)-format annotation and a FASTA reference genome indexed with [`samtools faidx`](https://github.com/samtools/samtools). Install [gffread](https://github.com/gpertea/gffread) and run it like this: `gffread -F -w [output_tx.fa] -g [genome_fasta.fa] [unzipped_gtf.gtf]` where `[output_tx.fa]` is the output file of transcripts to be created, `genome_fasta.fa` is an unzipped genome in FASTA format, and `unzipped_gtf.gtf` is the unzipped gene annotation.
* Run [FASTK](https://github.com/thegenemyers/FASTK) to count k-mers in each species' transcriptome. You should use a value of k that is large enough to be unique but small enough to allow for multiple k-mers per read and limit the amount of RAM required; a value between 25-45 should work well. To run FASTK, do this:
  * `FastK -N[output_prefix] -k[kmer_size] -t1 [output_tx.fa]` where `[output_prefix]` will be the beginning of the name of the output files and `[output_tx.fa]` is the file you created in the last step. Note that there should be no space between argument names and argument values.
* Compare each species' k-mer sets to find k-mers unique to each species' transcriptome, and create the files needed by `demux_species`
  * Run `utils/get_unique_kmers` in this package, with `-k [name.ktab]` specified once for each species you ran `FastK` on and `-n [name]` once again for each species, in the same order. For example, if you have the files `human.ktab`, `chicken.ktab`, and `mouse.ktab` representing human, chicken, and mouse k-mers, you could run:
    
    ```
    utils/get_unique_kmers -k human.ktab -k chicken.ktab -k mouse.ktab -n Human -n Chicken -n Mouse -o hcm_kmers
    ```
    
    This will create a set of "unique k-mers" files beginning with `hcm_kmers`.
    * Note that you can also sample a set number of k-mers instead of choosing all of them. This can save memory and processing time at the cost of potentially identifying fewer cells. In our hands, we recommend sampling at least 10 million (`-N 10000000`) k-mers per species when doing this. With this option set, k-mers will be sorted in decreasing order of their frequency in each transcriptome before sampling. k-mers with equal frequency above the threshold for selection will be randomly chosen.
   
[Back to demux_species](demux_species.md)
