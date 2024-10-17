#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// Define programs
split_reads = baseDir + "/../utils/split_read_files"
demux_species = baseDir + "/../demux_species"
combine_species_counts = baseDir + "/../utils/combine_species_counts"
plot = baseDir + "/../plot/species.R"

// Error check parameters
if (!params.kmers){
    error("No kmers specified. Please provide base name given to get_unique_kmers")
}
if (!params.whitelist){
    error("No barcode whitelist specified. Provide whitelist = RNA-seq barcode whitelist")
}
if (!params.output_directory){
    error("No output_directory specified.")
}
if (!params.rna_dir){
    error("No rna_dir specified. Please provide a directory containing RNA-seq reads")
}
if (!params.libs){
    error("No libs specified. Provide libs = a file listing library names, one per line.")
}
if (params.atac_dir){
    if (!params.whitelist_atac){
        error("ATAC data provided, but no whitelist_atac was given.")
    }
    if (!params.atac_map){
        println "No atac_map given. Assuming ATAC-seq library names match RNA-seq library names."
    }
}
if (params.custom_dir){
    if (!params.custom_map){
       error("Custom reads dir given without custom_map. Provide custom_map = file mapping custom library names to data types (tab-separated), with optional third column listing corresponding RNA-seq library name, if library names do not match.")
    }
}
if (params.num_chunks <= 1){
    params.demux_pieces = false
    params.num_chunks = 0
}

// Get absolute file paths
def abs_path(filename){
    File f = new File(filename)
    return f.getAbsoluteFile().toString()
}

def abs_path_kmers(filename){
    names = filename + ".names"
    File f = new File(names)
    abspath = f.getAbsoluteFile().toString()
    return abspath.substring(0, abspath.length()-6)
}

abs_wl = abs_path(params.whitelist)
abs_kmers = abs_path_kmers(params.kmers)
abs_out = abs_path(params.output_directory)
abs_rna = abs_path(params.rna_dir)
abs_atac = ""
if (params.atac_dir){
    abs_atac = abs_path(params.atac_dir)
}
abs_custom = ""
if (params.custom_dir){
    abs_custom = abs_path(params.custom_dir)
}

// Get maps of ATAC/custom lib names -> RNA lib names
def atac_map = [:]
def custom_map = [:]

// Get a map of custom lib names -> custom lib types
def custom_names = [:]

if (params.atac_map){
    new File(params.atac_map).eachLine { line ->
        def (atac, gex) = line.split('\t')
        atac_map[gex] = atac
    }
    def countmatch = 0
    new File(params.libs).eachLine { line -> 
        if ( line.trim() in atac_map ){
            countmatch++
        }
    }
    if (countmatch == 0){
        error("""
No RNA-seq libraries found in ${params.libs} matching those listed in ${params.atac_map}.
Please check ${params.atac_map} and ensure you have named libraries correctly.
""")
    }
}
else{
    new File(params.libs).eachLine { line ->
        def linetrim = line.trim()
        atac_map[linetrim] = linetrim
}
}

if (params.custom_map){
    new File(params.custom_map).eachLine { line -> 
        def fields = line.split('\t')
        if (fields[1].split(' ')[0] == "CRISPR"){
            custom_names[fields[0]] = "CRISPR"
        }
        else if (fields[1].split(' ')[0] == "Antibody"){
            custom_names[fields[0]] = "Ab"
        }
        else{
            custom_names[fields[0]] = fields[1]
        }
        if (fields.size() == 3){
            if (fields[2] in custom_map){
                custom_map[fields[2]].add(fields[0])
            }
            else{
                custom_map[fields[2]] = [fields[0]]
            }
        }
    }
    if (custom_map.size() == 0){
        new File(params.lib).eachLine{ line -> 
            def linetrim = line.trim()
            custom_map[linetrim] = linetrim
        }
    }
    else{
        def countmatch = 0
        new File(params.libs).eachLine { line -> 
            if ( line.trim() in custom_map ){
                countmatch++
            }
        }
        if (countmatch == 0){
            error("""
    No RNA-seq libraries found in ${params.libs} matching those listed in ${params.custom_map}.
    Please check ${params.custom_map} and ensure you have named libraries correctly.
    """)
        }
    }
}

def get_lib_idx(lib){
    def lib_idx = [:]
    new File(params.libs).eachWithIndex { line, index ->
        lib_idx[line.trim()] = index
    }
    return lib_idx[lib]
}

// Count the number of species in the index
def n_species = file(params.kmers + ".names").readLines().size()

/*
 * Split RNA-seq read pairs into chunks.
 */
process split_rna_reads{
    input:
    tuple val(file_idx), 
        val(lib_id), 
        val(basename_R1), 
        val(basename_R2), 
        file(R1), 
        file(R2)
    
    output:
    tuple val(file_idx), 
        val(lib_id), 
        file("${basename_R1}.*.fastq.gz"), 
        file("${basename_R2}.*.fastq.gz")
    
    script:
    """
    ${split_reads} -1 ${R1} -2 ${R2} -n ${params.num_chunks} -o .
    """
}

/*
 * Split ATAC read pairs into chunks.
 */
process split_atac_reads{
    input:
    tuple val(file_idx), 
        val(name), 
        val(lib_id), 
        val(R1base), 
        val(R2base), 
        val(R3base),
        file(R1), 
        file(R2),
        file(R3)
    
    output:
    tuple val(file_idx), 
        val(name), 
        val(lib_id), 
        file("${R1base}.*.fastq.gz"), 
        file("${R2base}.*.fastq.gz"),
        file("${R3base}.*.fastq.gz") 
    script:
    """
    ${split_reads} -1 ${R1} -2 ${R2} -3 ${R3} -n ${params.num_chunks} -o .
    """
}

/*
 * Split custom-type read pairs into chunks.
 */
process split_custom_reads{
    input:
    tuple val(file_idx), 
        val(name), 
        val(lib_id), 
        val(type), 
        val(R1base), 
        val(R2base), 
        file(R1), 
        file(R2)
    
    output:
    tuple val(file_idx), 
        val(name), 
        val(lib_id), 
        val(type), 
        file("${R1base}.*.fastq.gz"), 
        file("${R2base}.*.fastq.gz")
    
    script:
    """
    ${split_reads} -1 ${R1} -2 ${R2} -n ${params.num_chunks} -o .
    """
}

/*
 * Count species specific k-mers across an entire file pair.
 */
process count_kmers{
    cpus params.threads
    memory params.memgb + ' GB'
    
    input:
    tuple val(file_idx), 
        val(lib_id),
        file(R1), 
        file(R2),
        file(wl),
        val(kmerbase),
        file(kmerfiles)
            
    output:
    tuple val(lib_id), 
        file("species_counts.*.txt"), 
        file("species_names.*.txt") 
    
    script:
    
    batch_idx = file_idx + 1
    
    """
    ${demux_species} -o . -d --batch_num ${batch_idx} -k ${kmerbase}\
         -w ${wl} -T ${params.threads} -r ${R1} -R ${R2} 
    """

}
/*
 * Count species specific k-mers in one input file chunk.
 */
process count_kmers_chunk{
    cpus params.threads
    memory params.memgb + ' GB'

    input:
    tuple val(file_idx), 
        val(lib_id), 
        file(R1), 
        file(R2),
        file(wl),
        val(kmerbase),
        file(kmerfiles)
            
    output:
    tuple val(lib_id), 
        file("species_counts.*.txt"), 
        file("species_names.*.txt") 
    
    script:
    
    // Extract batch ID/index from file name
    def (whole, begin, extra, idx) = \
        ( R1 =~ /(.*)\_R1(_\d+)?\.(\d+)\.fastq\.gz/ )[0]
    
    batch_idx = (file_idx * params.num_chunks) + idx.toInteger()
    
    """
    ${demux_species} -o . -d --batch_num ${batch_idx} -k ${kmerbase}\
         -w ${wl} -T ${params.threads} -r ${R1} -R ${R2} 
    """

}

def boolean isCollectionOrArray(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

/*
 * Combine together species-specific k-mer counts from all chunks.
 */
process join_counts{
    input:
    tuple val(libname), file(countsfiles), file(namefiles)
    
    output:
    tuple val(libname), file("species_counts.txt"), file("species_names.txt")

    script:
    def n_files = countsfiles.size()
    if (!isCollectionOrArray(countsfiles)){
        n_files = 1
    }
    """
    ${combine_species_counts} -o . -n ${n_files} 
    """
}

def get_fit_model_script_extra(libname){
    def extra=""
    if (params.append_libname){
        if (params.cellranger){
            libidx = get_lib_idx(libname) + 1
            extra += " -n ${libidx}"
        }
        else{
            extra += " -n ${libname}"
        }
        if (params.seurat){
            extra += " -S"
        }
        if (params.underscore){
            extra += " -U"
        }
    }
    return extra
}

def get_fn_arg_str(fns, arg){
    str = ""
    first = true
    for (fn in fns){
        if (!first){
            str += " "
        }
        first = false
        str += arg + ' ' + fn.toString()
    }
    return str
}

/*
 * After combining all k-mer counts, fit the model and write species assignments
 * to disk.
 * 
 * For use when there are no reads to separate.
 */
process fit_model_dump{
    input:
    tuple val(libname),
        file(countsfile),
        file(namefile)
    
    publishDir "${params.output_directory}/${libname}", mode: 'copy'

    output:
    tuple val(libname),
        file("dists.txt"),
        file("species.assignments"),
        file("species.filt.assignments"),
        file("species.pdf"),
        file("species.png"),
        file(countsfile),
        file(namefile)
    
    script:
    def extra = get_fit_model_script_extra(libname)
    def r1str = get_fn_arg_str(gexR1, '-r')
    def r2str = get_fn_arg_str(gexR2, '-R')
    """
    ${demux_species} -d -o . ${r1str} ${r2str}${extra}
    if [ \$( cat species.filt.assignments | wc -l ) -gt 0 ]; then
        ${plot} . 
    else
        echo "" > species.pdf
        echo "" > species.png
    fi
    """    
}

/*
 * After combining all k-mer counts, fit the model and write species assignments
 * to disk.
 * 
 * For use when there are only regular RNA-seq reads to separate.
 */
process fit_model_gex{
    input: 
    tuple val(libname), 
        file(gexR1), 
        file(gexR2), 
        file(countsfile), 
        file(namefile)
    
    publishDir "${params.output_directory}/${libname}", mode: 'copy'
    
    output:
    tuple val(libname), 
        file("dists.txt"), 
        file("species.assignments"), 
        file("species.filt.assignments"), 
        file("species.pdf"), 
        file("species.png"), 
        file(countsfile), 
        file(namefile), 
        file("*library")
    
    script:
    def extra = get_fit_model_script_extra(libname)
    def r1str = get_fn_arg_str(gexR1, '-r')
    def r2str = get_fn_arg_str(gexR2, '-R')
    """
    ${demux_species} -d -o . ${r1str} ${r2str}${extra}
    if [ \$( cat species.filt.assignments | wc -l ) -gt 0 ]; then
        ${plot} . 
    else
        echo "" > species.pdf
        echo "" > species.png
    fi
    """
}

/*
 * After combining all k-mer counts, fit the model and write species assignments
 * to disk.
 * 
 * For use when there are RNA-seq and ATAC-seq reads to separate.
 */
process fit_model_gex_atac{
    input: 
    tuple val(libname), 
        file(gexR1), 
        file(gexR2),
        val(atacname),
        file(atacR1),
        file(atacR2),
        file(atacR3),
        file(countsfile), 
        file(namefile)
    
    publishDir "${params.output_directory}/${libname}", mode: 'copy'
    
    output:
    tuple val(libname), 
        file("dists.txt"), 
        file("species.assignments"), 
        file("species.filt.assignments"), 
        file("species.pdf"), 
        file("species.png"), 
        file(countsfile), 
        file(namefile), 
        file("*library")
    
    script:
    def extra = get_fit_model_script_extra(libname)
    def r1str = get_fn_arg_str(gexR1, '-r')
    def r2str = get_fn_arg_str(gexR2, '-R')
    def ar1str = get_fn_arg_str(atacR1, '-1')
    def ar2str = get_fn_arg_str(atacR2, '-2')
    def ar3str = get_fn_arg_str(atacR3, '-3')
    """
    ${demux_species} -d -o . ${r1str} ${r2str} ${ar1str} ${ar2str} ${ar3str}${extra}
    if [ \$( cat species.filt.assignments | wc -l ) -gt 0 ]; then
        ${plot} . 
    else
        echo "" > species.pdf
        echo "" > species.png
    fi
    """    
}

/*
 * After combining all k-mer counts, fit the model and write species assignments
 * to disk.
 * 
 * For use when there are RNA-seq and custom-type reads to separate.
 */
process fit_model_gex_custom{
    input: 
    tuple val(libname), 
        file(gexR1), 
        file(gexR2), 
        val(custom_name), 
        val(custom_type), 
        file(customR1), 
        file(customR2), 
        file(countsfile), 
        file(namefile)
    
    publishDir "${params.output_directory}/${libname}", mode: 'copy'
    
    output:
    tuple val(libname), 
        file("dists.txt"), 
        file("species.assignments"), 
        file("species.filt.assignments"), 
        file("species.pdf"), 
        file("species.png"), 
        file(countsfile), 
        file(namefile), 
        file("*library")
    
    script:
    def extra = get_fit_model_script_extra(libname)
    def r1str = get_fn_arg_str(gexR1, '-r')
    def r2str = get_fn_arg_str(gexR2, '-R')
    def c1str = get_fn_arg_str(customR1, '-x')
    def c2str = get_fn_arg_str(customR2, '-X')
    def nstr = get_fn_arg_str(custom_type, '-N')
    """
    ${demux_species} -d -o . ${r1str} ${r2str} ${c1str} ${c2str} ${nstr}${extra}
    if [ \$( cat species.filt.assignments | wc -l ) -gt 0 ]; then
        ${plot} . 
    else
        echo "" > species.pdf
        echo "" > species.png
    fi
    """    
}

/*
 * After combining all k-mer counts, fit the model and write species assignments
 * to disk.
 * 
 * For use when there are RNA-seq, ATAC-seq, and custom-type reads to separate.
 */
process fit_model_gex_atac_custom{
    input: 
    tuple val(libname), 
        file(gexR1), 
        file(gexR2),
        val(atac_name),
        file(atacR1),
        file(atacR2),
        file(atacR3),
        val(custom_name), 
        val(custom_type), 
        file(customR1), 
        file(customR2), 
        file(countsfile), 
        file(namefile)
    
    publishDir "${params.output_directory}/${libname}", mode: 'copy'
    
    output:
    tuple val(libname), 
        file("dists.txt"), 
        file("species.assignments"), 
        file("species.filt.assignments"), 
        file("species.pdf"), 
        file("species.png"), 
        file(countsfile), 
        file(namefile), 
        file("*library")
    
    script:
    def extra = get_fit_model_script_extra(libname)
    def r1str = get_fn_arg_str(gexR1, '-r')
    def r2str = get_fn_arg_str(gexR2, '-R')
    def a1str = get_fn_arg_str(atacR1, '-1')
    def a2str = get_fn_arg_str(atacR2, '-2')
    def a3str = get_fn_arg_str(atacR3, '-3')
    def c1str = get_fn_arg_str(customR1, '-x')
    def c2str = get_fn_arg_str(customR2, '-X')
    def nstr = get_fn_arg_str(custom_type, '-N')
    """
    ${demux_species} -d -o . ${r1str} ${r2str} ${a1str} ${a2str} ${a3str} ${c1str} ${c2str} ${nstr}${extra}
    if [ \$( cat species.filt.assignments | wc -l ) -gt 0 ]; then
        ${plot} . 
    else
        echo "" > species.pdf
        echo "" > species.png
    fi
    """    
}
/* 
 * Splits an entire (non-split) file of RNA-seq reads by species
 */
process demux_rna_reads{
    cpus params.threads

    input:
    tuple val(libname), 
        file(R1), 
        file(R2), 
        file(assn), 
        file(counts), 
        file(names),
        file(wl)
 
    publishDir "${params.output_directory}/${libname}", mode: "copy"
    
    output:
    tuple file("*/GEX_${R1}"), file("*/GEX_${R2}")

    script:
    """
    ${demux_species} -o . -w ${wl} -T ${params.threads} -r ${R1} -R ${R2}
    """
}

/*
 * Splits an RNA-seq read pair chunk by species
 */
process demux_rna_reads_chunk{
    cpus params.threads

    input:
    tuple val(uid), 
        val(libname), 
        file(R1), 
        file(R2), 
        file(assn), 
        file(counts), 
        file(names),
        file(wl)
    
    output:
    tuple val(uid), 
        val(libname), 
        file("*/GEX_${R1}"), 
        file("*/GEX_${R2}")

    script:
    """
    ${demux_species} -o . -w ${wl} -T ${params.threads} -r ${R1} -R ${R2}
    """
}

/* 
 * Splits an entire (non-split) file of ATAC-seq reads by species
 */
process demux_atac_reads{
    cpus params.threads

    input:
    tuple val(libname), 
        file(R1), 
        file(R2),
        file(R3),
        file(assn), 
        file(counts), 
        file(names),
        file(wl),
        file(wl_atac)
 
    publishDir "${params.output_directory}/${libname}", mode: "copy"
    
    output:
    tuple file("*/ATAC_${R1}"), file("*/ATAC_${R2}"), file("*/ATAC_${R3}")

    script:
    """
    ${demux_species} -o . -w ${wl} -W "${wl_atac}" -T ${params.threads} -1 ${R1} -2 ${R2} -3 ${R3}
    """
}

/*
 * Splits an ATAC_seq read triplet chunk by species
 */
process demux_atac_reads_chunk{
    cpus params.threads

    input:
    tuple val(uid), 
        val(libname_atac), 
        val(libname_gex), 
        file(R1), 
        file(R2), 
        file(R3),
        file(assn), 
        file(counts), 
        file(names),
        file(wl),
        file(wl_atac)
 
    output:
    tuple val(uid), 
        val(libname_gex),
        file("*/ATAC_${R1}"),
        file("*/ATAC_${R2}"),
        file("*/ATAC_${R3}")
 
    script:
    """
    ${demux_species} -o . -w ${wl} -W ${wl_atac} -T ${params.threads} -1 ${R1} -2 ${R2} -3 ${R3}
    """
}

/* 
 * Splits an entire (non-split) file of custom type reads by species
 */
process demux_custom_reads{
    cpus params.threads

    input:
    tuple val(libname),
        val(libtype), 
        file(R1), 
        file(R2),
        file(assn), 
        file(counts), 
        file(names),
        file(wl)
 
    publishDir "${params.output_directory}/${libname}", mode: "copy"
    
    output:
    tuple file("*/${libtype}_${R1}"), file("*/${libtype}_${R2}")

    script:
    """
    ${demux_species} -o . -w ${wl} -T ${params.threads} -x ${R1} -X ${R2} -N ${libtype}
    """
}

/*
 * Splits a custom-type read pair chunk by species
 */
process demux_custom_reads_chunk{
    cpus params.threads

    input:
    tuple val(uid), 
        val(libname_custom), 
        val(type_custom), 
        val(libname_gex), 
        file(R1), 
        file(R2), 
        file(assn), 
        file(counts), 
        file(names),
        file(wl)
 
    output:
    tuple val(uid), 
        val(libname_gex),
        file("*/${type_custom}_${R1}"),
        file("*/${type_custom}_${R2}")    
    
    script:
    """
    ${demux_species} -o . -w ${wl} -T ${params.threads} -x ${R1} -X ${R2} -N ${type_custom}
    """
}

/*
 * Combines a set of RNA-seq read pair chunks, after splitting into chunks
 * and separating chunks by species
 */
process cat_read_chunks_rna{
    input:
    tuple val(uid), 
        val(libname), 
        val(species), 
        file(R1s), 
        file(R2s), 
        val(R1out), 
        val(R2out)
    
    publishDir "${params.output_directory}/${libname}/${species}", mode: "copy"
    
    output:
    tuple val(uid), 
        val(libname), 
        val(species), 
        file(R1out), 
        file(R2out)

    script:
    def r1s_str = R1s.join(' ')
    def r2s_str = R2s.join(' ')
    """
    zcat ${r1s_str} | gzip -c - > ${R1out}
    zcat ${r2s_str} | gzip -c - > ${R2out}
    """ 
}

/*
 * Combines a set of ATAC read triple chunks, after splitting into chunks
 * and separating chunks by species
 */
process cat_read_chunks_atac{
    input:
    tuple val(uid), 
        val(libname), 
        val(species), 
        file(R1s), 
        file(R2s),
        file(R3s),
        val(R1out), 
        val(R2out),
        val(R3out)
    
    publishDir "${params.output_directory}/${libname}/${species}", mode: "copy"
    
    output:
    tuple val(uid), val(libname), val(species), file(R1out), file(R2out), file(R3out)

    script:
    def r1s_str = R1s.join(' ')
    def r2s_str = R2s.join(' ')
    def r3s_str = R3s.join(' ')
    """
    zcat ${r1s_str} | gzip -c - > ${R1out}
    zcat ${r2s_str} | gzip -c - > ${R2out}
    zcat ${r3s_str} | gzip -c - > ${R3out}
    """ 
}
/*
 * Combines a set of custom-type read pair chunks, after splitting into chunks
 * and separating chunks by species
 */
process cat_read_chunks_custom{
    input:
    tuple val(uid), 
        val(libname), 
        val(species), 
        file(R1s), 
        file(R2s), 
        val(R1out), 
        val(R2out)
    
    publishDir "${params.output_directory}/${libname}/${species}", mode: "copy"
    
    output:
    tuple val(uid), val(libname), val(species), file(R1out), file(R2out)

    script:
    def r1s_str = R1s.join(' ')
    def r2s_str = R2s.join(' ')
    """
    zcat ${r1s_str} | gzip -c - > ${R1out}
    zcat ${r2s_str} | gzip -c - > ${R2out}
    """ 
}

def cartesian_product(A, B) {
    A.collectMany{ a -> B.collect { b -> [a, b] } }
}

def extensions = [
    '.fastq.gz',
]

def suffixes = [
    '*_R{1,2}_001',
    '*_R{1,2}',
    '*_{1,2}',
]

def suffixes_triple = [
    '*_R{1,2,3}_001',
    '*_R{1,2,3}',
    '*_{1,2,3}',
]

/**
 * Returned channel:
 * include_index == false: grouped
 *   [lib, [r1s], [r2s]]
 * include_index == true:
 *   [file_index, lib, [r1, r2]]
 */
def get_rna_channel(extensions, suffixes, include_index){
    // Read library IDs from params.libs
    def library_ids = Channel.fromPath(params.libs).splitText().map { id ->
        id = id.trim()
        return id
    }
    // Get all RNA-seq read pairs 
    def patterns = cartesian_product(suffixes, extensions).collect { 
        "${params.rna_dir}/${it.join()}"
    }     
    def rna_pairs = Channel.fromFilePairs(patterns).map{ id, reads -> 
        [ id.replaceFirst(/_S(\d+)_L(\d+)/, ''), reads ]
    }
    
    if (include_index){
        // Return only RNA-seq read pairs with library IDs in the expected list
        // Return format:
        //   0-based file index (in library)
        //   library name
        //   base name R1 (when files are merged again)
        //   base name R2 (when files are merged again)
        //   R1 file
        //   R2 file
        def file_idx = -1 
        return library_ids.cross(rna_pairs).map{ x ->
            file_idx++
            def match = (x[1][1][0].toString() =~ /(.*)\/(.*)\_R1(_\d+)?\.(fastq|fq)(\.gz)?/)[0]
            def end = ""
            if (match[3] != null){
                end = match[3]
            }
            def r1base = match[2] + '_R1' + end
            def r2base = match[2] + '_R2' + end
            return [file_idx, x[0], r1base, r2base, x[1][1][0], x[1][1][1]]
        }
    }
    else{
        // Transform lib ID, [[r1A, r2A], [r1B, r2B]] into
        // lib ID, [r1A, r1B], [r2A, r2B]
        return rna_pairs.groupTuple().map{ lib, reads ->
            def r1s = []
            def r2s = []
            for (elt in reads){
                r1s.add(elt[0])
                r2s.add(elt[1])    
            }
            return [lib, r1s, r2s]
        }
    }
}

/**
 * Returned channel:
 * include_index == false:
 *   [libname, name, [r1, r2, r3]]
 * include_index == true:
 *   [file_index, libname, name, [r1, r2, r3]]
 */
def get_atac_channel(atac_map, include_index){
    def atac_library_ids = Channel.fromPath(params.libs).splitText().map { id ->
        def idtrim = id.trim()
        if (idtrim in atac_map){
            return [atac_map[idtrim], idtrim]
        }
    }
    
    def atac_triples = Channel.fromPath("${params.atac_dir}/*_R1*.fastq.gz").map{ fn ->
        def r1 = fn.toString().trim()
        def match = (r1 =~ /(.*)\/(.*)\_R1(_\d+)?\.(fastq|fq)(\.gz)?/)[0]
        def libname_atac = match[2].replaceFirst(/_S(\d+)_L(\d+)/, '')
        def dirn = ""
        if (match[1] != null && match[1] != ""){
            dirn += match[1] + '/'
        }
        def end = ""
        if (match[3] != null){
            end = match[3]
        }
        def gz = ""
        if (match[5] != null){
            gz = match[5]
        }
        def r2 = dirn + match[2] + "_R2" + end + '.' + match[4] + gz
        def r3 = dirn + match[2] + "_R3" + end + '.' + match[4] + gz
        def r1base = match[2] + "_R1" + end
        def r2base = match[2] + "_R2" + end
        def r3base = match[2] + "_R3" + match[3]
        return [libname_atac, r1base, r2base, r3base, file(r1), file(r2), file(r3)]
    }
    
    if (include_index){
        // Return only ATAC read triples with library IDs in the expected list
        // Return format:
        //   0-based file index (in library)
        //   library name (ATAC-specific)
        //   library name (corresponding RNA-seq library name)
        //   base name R1 (when files are merged again)
        //   base name R2 (when files are merged again)
        //   base name R3 (when files are merged again)
        //   R1 file
        //   R2 file
        //   R3 file
        def file_idx = -1 
        return atac_library_ids.cross(atac_triples).map{ x, y -> 
            file_idx++
            return [file_idx, x[0], x[1], y[1], y[2], y[3], y[4], y[5], y[6]]
        }
    }
    else{
        // Transform lib ID, [[r1A, r2A, r3A], [r1B, r2B, r3B]] into
        // lib ID, [r1A, r1B, r1C], [r2A, r2B, r2C], [r3A, r3B, r3C]
        return atac_library_ids.join(atac_triples.groupTuple())\
            .map{ lib_atac, lib, base1, base2, base3, r1, r2, r3 -> 
                return [lib, lib_atac, r1, r2, r3]
            }
    }
}

/*
 * Returned channel:
 * include_index == false: grouped
 *   [lib, [names], [types], [r1s], [r2s]]
 * include_index == true:
 *   [file_index, lib, name, type, [r1, r2]]
 */
def get_custom_channel(extensions, suffixes, custom_map, custom_names, include_index){
    def custom_library_ids = Channel.fromPath(params.libs).splitText().flatMap { id ->
        def ret = []
        id = id.trim();
        if (id in custom_map){
            for (elt in custom_map[id]){
                ret.add([elt, id, custom_names[elt]])
            }
        }
        return ret;
    }

    // Get all RNA-seq read pairs 
    def custom_patterns = cartesian_product(suffixes, extensions).collect { 
        "${params.custom_dir}/${it.join()}"
    }     
    def custom_pairs = Channel.fromFilePairs(custom_patterns).map{ id, reads -> 
        [ id.replaceFirst(/_S(\d+)_L(\d+)/, ''), reads ]
    }
    
    if (include_index){
        // Return only custom read pairs with library IDs in the expected list
        // Return format:
        //   0-based file index (in library)
        //   library name (custom library-specific)
        //   library name (corresponding RNA-seq library name)
        //   custom library type
        //   base name R1 (when files are merged again)
        //   base name R2 (when files are merged again)
        //   R1 file
        //   R2 file
        def file_idx = -1 
        def indexed_custom_pairs = custom_library_ids.cross(custom_pairs).map{ x, y ->
            file_idx++
            def match = (y[1][0].toString() =~ /(.*)\/(.*)\_R1(_\d+)?\.(fastq|fq)(\.gz)?/)[0]
            def end = ""
            if (match[3] != ""){
                end = match[3]
            }
            def r1base = match[2] + '_R1' + end
            def r2base = match[2] + '_R2' + end
            return [file_idx, x[0], x[1], x[2], r1base, r2base, y[1][0], y[1][1]]
        }
        return indexed_custom_pairs
    }
    else{  
        return custom_library_ids.cross(custom_pairs).map{ x, y ->
            [x[1], x[0], x[2], y[1]]
        }.groupTuple().map{ x -> 
            def r1s = []
            def r2s = []
            for (readpair in (x[3])){
                r1s.add(readpair[0])
                r2s.add(readpair[1])
            }
            return [x[0], x[1], x[2], r1s, r2s]
        }
    }
}

/*
 * Creates a channel of files to concatenate, after separating by species
 * reads that were originally split into chunks.
 */
def cat_job(demuxed, has_R3){
    
    to_cat = demuxed.flatMap{ x ->
        def uid = x[0]
        def lib = x[1]
        def r1s = x[2]
        def r2s = x[3]
        def species_R1 = [:]
        def species_R2 = [:]
        def species_R3 = [:]
        def species_R1_out = ""
        def species_R2_out = ""
        def species_R3_out = ""
        for (r1 in r1s){
            def r1split = r1.toString().split('/')
            def r1spec = r1split[r1split.size()-2]
            species_R1[r1spec] = r1
            if (species_R1_out == ""){
                species_R1_out = r1split[r1split.size()-1]\
                    .replaceFirst(/_R1(_\d+)?\.(\d+).fastq.gz/, '_R1.fastq.gz')
            }
        }
        for (r2 in r2s){    
            def r2split = r2.toString().split('/')
            def r2spec = r2split[r2split.size()-2]
            species_R2[r2spec] = r2 
            if (species_R2_out == ""){
                species_R2_out = r2split[r2split.size()-1]\
                    .replaceFirst(/_R2(_\d+)?\.(\d+).fastq.gz/, '_R2.fastq.gz')
            }
        }
        if (has_R3){
            for (r3 in x[4]){
                def r3split = r3.toString().split('/')
                def r3spec = r3split[r3split.size()-2]
                species_R3[r3spec] = r3
                if (species_R3_out == ""){
                    species_R3_out = r3split[r3split.size()-1]\
                        .replaceFirst(/_R3(_\d+)?\.(\d+).fastq.gz/, '_R3.fastq.gz')
                }
            }
        }
        
        def to_return = []
        for (species in species_R1.keySet()){
            def row = [uid + "_" + species, lib, species, \
                species_R1[species], species_R2[species]]
            if (has_R3){
                row.add(species_R3[species])
            }
            row.add(species_R1_out)
            row.add(species_R2_out)
            if (has_R3){
                row.add(species_R3_out)
            }
            to_return.add(row)
        }
        return to_return;
    }.groupTuple().map{x -> 
        uid = x[0]
        libnames = x[1]
        species = x[2]
        R1s = x[3]
        R2s = x[4]
        R1outs = x[5]
        R2outs = x[6]
        if (has_R3){
            R3s = x[5]
            R1outs = x[6]
            R2outs = x[7]
            R3outs = x[8]
            return [uid, libnames[0], species[0], R1s, R2s, R3s, R1outs[0], \
                R2outs[0], R3outs[0]]
        }
        else{
            return [uid, libnames[0], species[0], R1s, R2s, R1outs[0], R2outs[0]]
        }
    }
    return to_cat
}

workflow{
    
    kmerfiles = Channel.fromPath(params.kmers + ".*{names,kmers}").collect().map{ x ->
        [x]
    }
     
    indexed_rna_pairs = get_rna_channel(extensions, suffixes, true) 
    
    joined = null
    
    def kmerSplit = params.kmers.toString().split("/")
    def kmerBase = kmerSplit[kmerSplit.size()-1]
    
    if (params.num_chunks <= 1){
        to_count = indexed_rna_pairs.map{ x -> 
            return [x[0], x[1], x[4], x[5], file(params.whitelist), kmerBase]
        }.combine(kmerfiles)

        filestojoin = count_kmers(to_count).groupTuple()
        
        // Combine counts
        joined = join_counts(filestojoin)        
    }
    else{
        rsplit = split_rna_reads(indexed_rna_pairs).transpose().map({ tup ->
            def tup2 = tup
            tup2.add(file(params.whitelist))
            return tup2
        })
        
        // Add in k-mer files before counting k-mers    
        to_count = rsplit.map{ rsplitelts ->
            def to_return = []
            for (elt in rsplitelts){
                to_return.add(elt)
            }
            to_return.add(kmerBase)
            return to_return
        }.combine(kmerfiles)
        
        // Group files by library name
        filestojoin = count_kmers_chunk(to_count).groupTuple()

        // Combine counts
        joined = join_counts(filestojoin)
    }

    reads_rna = null
    reads_atac = null
    reads_custom = null
    
    if (params.demux_reads){
        
        reads_rna = get_rna_channel(extensions, suffixes, false) 
        
        if (params.atac_dir){
            reads_atac = get_atac_channel(atac_map, false)
            atac_plus_rna = reads_atac.cross(reads_rna).map{ x, y ->
                return [y[0], y[1], y[2], x[1], x[2], x[3], x[4]]
            }
        }
        if (params.custom_dir){
            reads_custom = get_custom_channel(extensions, suffixes, custom_map,
                custom_names, false)
            custom_plus_rna = reads_custom.cross(reads_rna).map{ x, y ->
                return [x[0], y[1], y[2], x[1], x[2], x[3], x[4]]
            }
        }
         
        // Fit model, assign species, and create library files
        
        if (!params.atac_dir && !params.custom_dir){
            // Just RNA.
            rna_joined = reads_rna.cross(joined).map{ x -> 
                def ret = x[0]
                ret.add(x[1][1])
                ret.add(x[1][2])
                return ret
            }
            models = fit_model_gex(rna_joined)
        }
        else if (params.custom_dir && !params.atac_dir){
            
            custom_plus_rna_joined = custom_plus_rna.cross(joined.collect()).collect().map{ x ->
                def ret = x[0]
                ret.add(x[1][1])
                ret.add(x[1][2])
                return ret
            }
            models = fit_model_gex_custom(custom_plus_rna_joined)
        }
        else if (params.atac_dir && !params.custom_dir){
            atac_plus_rna_joined = atac_plus_rna.cross(joined.collect()).collect().map{ x ->
                def ret = x[0]
                ret.add(x[1][1])
                ret.add(x[1][2])
                return ret
            }
            models = fit_model_gex_atac(atac_plus_rna_joined)    
        }
        else if (params.atac_dir && params.custom_dir){ 
            all_together = atac_plus_rna.cross(custom_plus_rna).map{ x, y ->
                return [x[0], x[1], x[2], x[3], x[4], x[5], x[6], y[3], y[4], y[5], y[6]]
            }
            all_together_joined = all_together.cross(joined.collect()).collect().map{ x-> 
                def ret = x[0]
                ret.add(x[1][1])
                ret.add(x[1][2])
                return ret
            }

            models = fit_model_gex_atac_custom(all_together_joined)
        }
        
        if (params.demux_pieces){
            // Demux each chunk of RNA files
            rsplit2 = rsplit.map{ idx, lib, r1, r2, wl -> 
                [lib, lib + "_" + idx, r1, r2, wl]
            } 
            
            rna_chunk_demuxed = models.cross(rsplit2).map{ x, y ->
                [y[1], y[0], y[2], y[3], x[2], x[6], x[7], y[4]]
            } | demux_rna_reads_chunk
            
            cat_job(rna_chunk_demuxed, false) | cat_read_chunks_rna
            
            if (params.atac_dir){
                // Split, demux, and join ATAC-seq read files
                readtrios_atac = get_atac_channel(atac_map, true)
                asplit = split_atac_reads(readtrios_atac).transpose()\
                    .map{ idx, libname_atac, libname_gex, r1, r2, r3 ->
                        [libname_gex, libname_atac + "_" + idx, libname_atac, r1, r2, r3]
                    }
                atac_reads_chunk_demuxed = models.cross(asplit).map{x, y ->
                    [y[1], y[2], y[0], y[3], y[4], y[5], x[2], x[6], x[7], \
                    file(params.whitelist), file(params.whitelist_atac)]
                } | demux_atac_reads_chunk
                
                cat_job(atac_reads_chunk_demuxed, true) | cat_read_chunks_atac 
            }
            if (params.custom_dir){
                // Split, demux, and join custom read files
                readpairs_custom = get_custom_channel(extensions, suffixes,
                    custom_map, custom_names, true)
                csplit = split_custom_reads(readpairs_custom).transpose()\
                    .map{ idx, libname_custom, libname_gex, type, r1, r2 ->
                        [libname_gex, libname_custom + "_" + idx, libname_custom, type, r1, r2]
                    }
                
                custom_reads_chunk_demuxed = models.cross(csplit).map{x, y ->
                    [y[1], y[2], y[3], y[0], y[4], y[5], x[2], x[6], x[7], file(params.whitelist)]
                } | demux_custom_reads_chunk
                
                cat_job(custom_reads_chunk_demuxed, false) | cat_read_chunks_custom
            }
        }
        else{
            // Demux read files intact.
            non_indexed_rna_pairs = indexed_rna_pairs.map{ x -> 
                def ret = []
                for (int i = 1; i < x.size(); ++i){
                    ret.add(x[i])
                }
                return ret
            }
            models.cross(non_indexed_rna_pairs).map{ x, y -> 
                [x[0], y[3], y[4], x[2], x[6], x[7], file(params.whitelist)]
            } | demux_rna_reads
            
            if (params.atac_dir){
                readtrios_atac = get_atac_channel(atac_map, true)
                non_indexed_atac_trios = readtrios_atac.map{ x ->
                    [x[2], x[6], x[7], x[8]] 
                }
                models.cross(non_indexed_atac_trios).map{ x, y ->
                    [x[0], y[1], y[2], y[3], x[2], x[6], x[7], file(params.whitelist), file(params.whitelist_atac)]
                } | demux_atac_reads
            }
            if (params.custom_dir){
                readpairs_custom = get_custom_channel(extensions, suffixes,
                    custom_map, custom_names, true)
                non_indexed_readpairs_custom = readpairs_custom.map{ x -> 
                    [x[2], x[3], x[6], x[7]]
                }
                models.cross(non_indexed_readpairs_custom).map{ x, y ->
                    [x[0], y[1], y[2], y[3], x[2], x[6], x[7], file(params.whitelist)]
                } | demux_custom_reads
            }
        }
    }
    else{
        // Just need to join model and dump.
        fit_model_dump(joined)     
    }
}

