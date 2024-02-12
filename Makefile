SHELL=bash
COMP=g++
CCOMP=gcc
PREFIX ?=/usr/local
FLAGS=-std=c++11 --std=gnu++11 -fPIC
CFLAGS=-Wall
IFLAGS=-I$(PREFIX)/include -Iinclude
LFLAGS=-L$(PREFIX)/lib -Llib
MAX_SITES ?= 2000
MAKE=make
PROJROOT=$($(SHELL) pwd)
DEPS=lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a

all: demux/demux_vcf demux/demux_mt analysis/quant_contam utils/bam_indiv_rg utils/bam_split_bcs

demux/demux_vcf: src/demux_vcf.cpp src/robin_hood.h build/common.o build/demux_vcf_io.o build/demux_vcf_hts.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) -std=c++11 --std=gnu++11 -fPIC src/demux_vcf.cpp -o demux/demux_vcf build/common.o build/demux_vcf_io.o build/demux_vcf_hts.o -lz -lhts lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a

demux/demux_mt: src/demux_mt.cpp src/robin_hood.h src/common.h build/common.o $(DEPS)
	$(COMP) -D MAX_SITES=$(MAX_SITES) $(IFLAGS) $(LFLAGS) $(FLAGS) src/demux_mt.cpp -o demux/demux_mt build/common.o -lz -lhts lib/libmixturedist.a lib/libhtswrapper.a

analysis/quant_contam: src/quant_contam.cpp src/ambient_rna.h build/common.o build/demux_vcf_io.o build/ambient_rna.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) -std=c++11 --std=gnu++11 -fPIC src/quant_contam.cpp -o analysis/quant_contam build/common.o build/demux_vcf_io.o build/ambient_rna.o lib/libmixturedist.a lib/libhtswrapper.a lib/liboptimml.a -lz

utils/bam_indiv_rg: src/bam_indiv_rg.cpp build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) src/bam_indiv_rg.cpp -o utils/bam_indiv_rg build/common.o -lz -lhts $(DEPS)

utils/bam_split_bcs: src/bam_split_bcs.cpp build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) src/bam_split_bcs.cpp -o utils/bam_split_bcs build/common.o -lz -lhts $(DEPS)

build/common.o: src/common.cpp src/common.h lib/libhtswrapper.a lib/libmixturedist.a
	$(COMP) $(IFLAGS) $(FLAGS) src/common.cpp -c -o build/common.o

build/demux_vcf_io.o: src/demux_vcf_io.cpp src/demux_vcf_io.h src/common.h
	$(COMP) $(IFLAGS) $(FLAGS) src/demux_vcf_io.cpp -c -o build/demux_vcf_io.o

build/demux_vcf_hts.o: src/demux_vcf_hts.cpp src/demux_vcf_hts.h src/common.h lib/libhtswrapper.a
	$(COMP) $(IFLAGS) $(FLAGS) src/demux_vcf_hts.cpp -c -o build/demux_vcf_hts.o

build/ambient_rna.o: src/ambient_rna.cpp src/ambient_rna.h src/robin_hood.h $(DEPS)
	$(COMP) $(IFLAGS) $(FLAGS) src/ambient_rna.cpp -c -o build/ambient_rna.o

lib/libhtswrapper.a:
	cd dependencies/htswrapper && $(MAKE) PREFIX=../..
	cd dependencies/htswrapper && $(MAKE) install PREFIX=../..

lib/libmixturedist.a:
	cd dependencies/mixtureDist && $(MAKE) PREFIX=../..
	cd dependencies/mixtureDist && $(MAKE) install PREFIX=../..

lib/liboptimml.a:
	cd dependencies/optimML && $(MAKE) PREFIX=../..
	cd dependencies/optimML && $(MAKE) install PREFIX=../..

clean:
	rm build/*.o
	rm lib/*.a
	rm utils/bam_indiv_rg
	rm utils/bam_split_bcs
	rm demux/demux_mt
	rm demux/demux_vcf
