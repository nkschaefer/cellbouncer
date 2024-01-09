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
DEPS=lib/libmixturedist.a lib/libhtswrapper.a

all: demux/demux_vcf demux/demux_mt utils/bam_indiv_rg utils/bam_split_bcs

demux/demux_vcf: src/demux_vcf.cpp src/robin_hood.h src/nnls.h build/nnls.o build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) src/demux_vcf.cpp -o demux/demux_vcf build/nnls.o build/common.o -llapack -lz -lhts lib/libmixturedist.a lib/libhtswrapper.a

demux/demux_mt: src/demux_mt.cpp src/robin_hood.h src/common.h build/common.o $(DEPS)
	$(COMP) -D MAX_SITES=$(MAX_SITES) $(IFLAGS) $(LFLAGS) $(FLAGS) src/demux_mt.cpp -o demux/demux_mt build/nnls.o build/common.o -lz -lhts lib/libmixturedist.a lib/libhtswrapper.a

utils/bam_indiv_rg: src/bam_indiv_rg.cpp build/nnls.o build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) src/bam_indiv_rg.cpp -o utils/bam_indiv_rg build/nnls.o build/common.o -lz -lhts $(DEPS)

utils/bam_split_bcs: src/bam_split_bcs.cpp build/nnls.o build/common.o $(DEPS)
	$(COMP) $(IFLAGS) $(LFLAGS) $(FLAGS) src/bam_split_bcs.cpp -o utils/bam_split_bcs build/nnls.o build/common.o -lz -lhts $(DEPS)

build/nnls.o: src/nnls.c src/nnls.h
	$(CCOMP) $(CFLAGS) -c src/nnls.c -o build/nnls.o

build/common.o: src/common.cpp src/common.h lib/libhtswrapper.a lib/libmixturedist.a
	$(COMP) $(IFLAGS) $(FLAGS) src/common.cpp -c -o build/common.o

lib/libhtswrapper.a: dependencies/htswrapper/lib/libhtswrapper.a
	cd dependencies/htswrapper && $(MAKE) PREFIX=../..
	cd dependencies/htswrapper && $(MAKE) install PREFIX=../..

lib/libmixturedist.a: dependencies/mixtureDist/lib/libmixturedist.a
	cd dependencies/mixtureDist && $(MAKE) PREFIX=../..
	cd dependencies/mixtureDist && $(MAKE) install PREFIX=../..

clean:
	rm build/*.o
	rm lib/*.a
	rm utils/bam_indiv_rg
	rm utils/bam_split_bcs
	rm demux/demux_mt
	rm demux/demux_vcf
