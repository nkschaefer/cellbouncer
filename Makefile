SHELL=bash
COMP=g++
CCOMP=gcc
PREFIX ?=/usr/local
FLAGS=-std=c++11 --std=gnu++11 -fPIC -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2)
CFLAGS=-Wall
IFLAGS=-I$(PREFIX)/include -Iinclude
LFLAGS=-L$(PREFIX)/lib -Llib
MAX_SITES ?= 2000
MAKE=make
PROJROOT=$(shell pwd)
BC_LENX2=32
KX2=16
DEPS=lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a

all: demux_vcf demux_mt demux_tags demux_species quant_contam doublet_dragon utils/bam_indiv_rg utils/bam_split_bcs utils/get_unique_kmers utils/fastq_cell_bcs utils/split_read_files utils/atac_fq_preprocess

demux_vcf: src/demux_vcf.cpp build/common.o build/demux_vcf_io.o build/demux_vcf_hts.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) -g src/demux_vcf.cpp -o demux_vcf build/common.o build/demux_vcf_io.o build/demux_vcf_hts.o -lz -lhts lib/libmixturedist.a lib/liboptimml.a lib/libhtswrapper.a

demux_mt: src/demux_mt.cpp src/common.h build/common.o $(DEPS)
	$(COMP) -D MAX_SITES=$(MAX_SITES) $(IFLAGS) $(LFLAGS) $(FLAGS) src/demux_mt.cpp -o demux_mt build/common.o lib/libhtswrapper.a lib/libmixturedist.a lib/liboptimml.a -lhts -lz

demux_species: src/demux_species.cpp src/kmsuftree.h src/common.h build/kmsuftree.o build/common.o build/demux_species_io.o build/species_kmers.o build/reads_demux.o
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) -pthread src/demux_species.cpp -o demux_species build/kmsuftree.o build/common.o build/demux_species_io.o build/species_kmers.o build/reads_demux.o lib/libhtswrapper.a lib/libmixturedist.a -lhts -lz

demux_tags: src/demux_tags.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) -D PROJ_ROOT=$(PROJROOT) src/demux_tags.cpp -o demux_tags build/common.o lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a -lz

quant_contam: src/common.h src/quant_contam.cpp src/ambient_rna.h build/common.o build/demux_vcf_io.o build/ambient_rna.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) src/quant_contam.cpp -o quant_contam build/common.o build/demux_vcf_io.o build/ambient_rna.o lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a -lz

doublet_dragon: src/doublet_dragon.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) -g src/doublet_dragon.cpp -o doublet_dragon build/common.o lib/liboptimml.a lib/libmixturedist.a lib/libhtswrapper.a -lz

utils/bam_indiv_rg: src/bam_indiv_rg.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) src/bam_indiv_rg.cpp -o utils/bam_indiv_rg build/common.o $(DEPS) -lz -lhts

utils/bam_split_bcs: src/bam_split_bcs.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) src/bam_split_bcs.cpp -o utils/bam_split_bcs build/common.o $(DEPS) -lz -lhts

utils/get_unique_kmers: src/get_unique_kmers.c src/FASTK/libfastk.c build/libfastk.o
	$(CCOMP) -fPIC -g src/get_unique_kmers.c -o utils/get_unique_kmers build/libfastk.o -lz

utils/kmsuftree_test: src/kmsuftree_test.cpp src/kmsuftree.h build/kmsuftree.o
	$(COMP) $(FLAGS) src/kmsuftree_test.cpp -o utils/kmsuftree_test build/kmsuftree.o

utils/atac_fq_preprocess: src/atac_fq_preprocess.cpp src/common.h build/common.o $(DEPS)
	$(COMP) $(FLAGS) $(IFLAGS) $(LFLAGS) src/atac_fq_preprocess.cpp -o utils/atac_fq_preprocess build/common.o $(DEPS) -lz -lhts

utils/split_read_files: src/split_read_files.cpp src/common.h build/common.o
	$(COMP) $(FLAGS) $(IFLAGS) $(LFLAGS) src/split_read_files.cpp -o utils/split_read_files build/common.o $(DEPS) -lz -lhts

utils/combine_species_counts: src/combine_species_counts.cpp src/common.h build/common.o
	$(COMP) $(FLAGS) $(IFLAGS) $(LFLAGS) src/combine_species_counts.cpp -o utils/combine_species_counts build/common.o $(DEPS) -lz -lhts

build/common.o: src/common.cpp src/common.h lib/libhtswrapper.a lib/libmixturedist.a
	$(COMP) $(IFLAGS) $(FLAGS) src/common.cpp -c -o build/common.o

build/demux_vcf_io.o: src/demux_vcf_io.cpp src/demux_vcf_io.h src/common.h
	$(COMP) $(IFLAGS) $(FLAGS) src/demux_vcf_io.cpp -c -o build/demux_vcf_io.o

build/demux_vcf_hts.o: src/demux_vcf_hts.cpp src/demux_vcf_hts.h src/common.h lib/libhtswrapper.a
	$(COMP) $(IFLAGS) $(FLAGS) src/demux_vcf_hts.cpp -c -o build/demux_vcf_hts.o

build/ambient_rna.o: src/ambient_rna.cpp src/ambient_rna.h src/common.h $(DEPS)
	$(COMP) $(IFLAGS) $(FLAGS) src/ambient_rna.cpp -c -o build/ambient_rna.o

build/species_kmers.o: src/species_kmers.cpp src/species_kmers.h src/common.h src/kmsuftree.h
	$(COMP) $(IFLAGS) $(FLAGS) src/species_kmers.cpp -c -o build/species_kmers.o 

build/reads_demux.o: src/reads_demux.cpp src/reads_demux.h src/common.h lib/libhtswrapper.a
	$(COMP) $(IFLAGS) $(FLAGS) src/reads_demux.cpp -c -o build/reads_demux.o

build/demux_species_io.o: src/demux_species_io.cpp src/demux_species_io.h src/common.h lib/libhtswrapper.a
	$(COMP) $(IFLAGS) $(FLAGS) src/demux_species_io.cpp -c -o build/demux_species_io.o

build/kmsuftree.o: src/kmsuftree.c src/kmsuftree.h
	$(CCOMP) -fPIC -c src/kmsuftree.c -o build/kmsuftree.o

build/libfastk.o: src/FASTK/libfastk.c src/FASTK/libfastk.h
	$(CCOMP) -fPIC -c src/FASTK/libfastk.c -o build/libfastk.o

build/gene_core.o: src/FASTK/gene_core.c src/FASTK/gene_core.h
	$(CCOMP) -fPIC -c src/FASTK/gene_core.c -o build/gene_core.o

lib/libhtswrapper.a:
	#cd dependencies/htswrapper && $(MAKE) clean
	cd dependencies/htswrapper && $(MAKE) PREFIX=../.. BC_LENX2=$(BC_LENX2) KX2=$(KX2)
	cd dependencies/htswrapper && $(MAKE) install PREFIX=../..

lib/libmixturedist.a:
	#cd dependencies/mixtureDist && $(MAKE) clean
	cd dependencies/mixtureDist && $(MAKE) PREFIX=../..
	cd dependencies/mixtureDist && $(MAKE) install PREFIX=../..

lib/liboptimml.a:
	#cd dependencies/optimML && $(MAKE) clean
	cd dependencies/optimML && $(MAKE) PREFIX=../..
	cd dependencies/optimML && $(MAKE) install PREFIX=../..

clean:
	rm build/*.o
	rm lib/*.a
	rm utils/bam_indiv_rg
	rm utils/bam_split_bcs
	rm utils/get_unique_kmers
	rm utils/kmsuftree_test
	rm utils/atac_fq_preprocess
	rm utils/split_read_files
	rm utils/combine_species_counts
	rm demux_vcf
	rm demux_mt
	rm demux_species
	rm demux_tags
	rm quant_contam
	rm doublet_dragon
