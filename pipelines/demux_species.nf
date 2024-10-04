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
}

if (params.custom_map){
    new File(params.custom_map).eachLine { line -> 
        def fields = line.split('\t')
        custom_names[fields[0]] = fields[1]
        if (fields.size() == 3){
            if (fields[2] in custom_map){
                custom_map[fields[2]].add(fields[0])
            }
            else{
                custom_map[fields[2]] = [fields[0]]
            }
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

// Function to count files starting with a specific prefix
def count_lib_readfiles(String libname) {
    def dir = new File(params.rna_dir) 
    if (!dir.exists() || !dir.isDirectory()) {
        throw new IllegalArgumentException("The specified directory does not exist or is not a directory.")
    }

    // Count files that start with the given prefix
    def count = dir.listFiles().count { file ->
        file.name.startsWith(libname) && file.isFile() && file.name.contains("_R1") && 
            file.name.endsWith(".fastq.gz") 
    }
    return count
}

def copy_kmers(kmer_base){
    def idx = 0
    
    fn = file(kmer_base + ".names")

    // Create a File object
    def fileToCheck = new File(params.kmers + '.' + idx + '.kmers')
    while (fileToCheck.exists()){
        fn = file(params.kmers + '.' + idx + '.kmers')
        println "Exists " + idx
        idx += 1
        fileToCheck = new File(params.kmers + '.' + idx + '.kmers')
    }
}

// Count the number of species in the index
def n_species = file(params.kmers + ".names").readLines().size()

process split_rna_reads{
    input:
    tuple val(file_idx), val(lib_id), file(filepair)
    
    output:
    tuple val(file_idx), val(lib_id), file("${lib_id}*_R1.*.fastq.gz"), file("${lib_id}*_R2.*.fastq.gz")
    
    script:
    """
    ${split_reads} -1 ${filepair[0]} -2 ${filepair[1]} -n ${params.num_chunks} -o .
    """
}

process count_kmers{

    input:
    tuple val(file_idx), val(lib_id), file(R1), file(R2)
    
    //publishDir "${workDir}/${lib_id}", mode: 'copy'

    output:
    tuple val(lib_id), file("species_counts.*.txt"), file("species_names.*.txt") 
    
    script:

    // Extract batch ID/index from file name
    def (whole, begin, idx) = ( R1 =~ /(.*)\_R1\.(\d+)\.fastq\.gz/ )[0]
    
    batch_idx = (file_idx * params.num_chunks) + idx.toInteger()
    
    """
    ${demux_species} -o . -d --batch_num ${batch_idx} -k ${abs_kmers} -w ${abs_wl} -T ${params.threads} -r ${R1} -R ${R2} 
    """

}

process join_counts{
    input:
    tuple val(libname), file(countsfiles), file(namefiles)
    
    //publishDir "${params.output_directory}/${libname}", mode: 'copy'

    output:
    tuple val(libname), file("species_counts.txt"), file("species_names.txt")

    script:
    def n_files = countsfiles.size()
    
    """
    ${combine_species_counts} -o . -n ${n_files} 
    """
}

process fit_model{
    input:
    tuple val(libname), file(countsfile), file(namefile)
    
    publishDir "${params.output_directory}/${libname}", mode: 'copy'

    output:
    tuple val(libname), file("dists.txt"), file("species.assignments"), file("species.filt.assignments"), file("species.pdf"), file("species.png"), file(countsfile), file(namefile), file("*.library")

    script:
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
    def libname_atac=""
    if (params.atac_map){
        libname_atac = atac_map[libname]
    }
    def libnames_custom=[]
    def types_custom=[]
    if (params.custom_map){
        libnames_custom = custom_map[libname]
        for (n in libnames_custom){
            types_custom.add(custom_names[n])
        }
    }
    def libnames_custom_txt=libnames_custom.join(' ')
    def types_custom_txt=types_custom.join(' ')
    """
    rna_files1=""
    rna_files2=""
    atac_files1=""
    atac_files2=""
    atac_files3=""
    custom_files1=""
    custom_files2=""
    custom_types=""
    if [ ${params.demux_reads} ]; then
        rna_files1=\$(ls ${abs_rna}/${libname}*_R1*.fastq.gz | sort -V | sed 's/^/-r /' | tr '\\n' ' ')
        rna_files2=\$(ls ${abs_rna}/${libname}*_R2*.fastq.gz | sort -V | sed 's/^/-R /' | tr '\\n' ' ')
        rna_files1=" \${rna_files1}"
        rna_files2=" \${rna_files2}"  
        
        if [ '${abs_atac}' != "" ]; then
            atac_files1=\$(ls ${abs_atac}/${libname_atac}*_R1*.fastq.gz | sort -V | sed 's/^/-1 /' | tr '\\n' ' ' )
            atac_files2=\$(ls ${abs_atac}/${libname_atac}*_R2*.fastq.gz | sort -V | sed 's/^/-2 /' | tr '\\n' ' ' )
            atac_files3=\$(ls ${abs_atac}/${libname_atac}*_R3*.fastq.gz | sort -V | sed 's/^/-3 /' | tr '\\n' ' ' )
            atac_files1=" \${atac_files1}"
            atac_files2=" \${atac_files2}"
            atac_files3=" \${atac_files3}"
        fi
        
        if [ '${abs_custom}' != "" ]; then
            idx=0
            for lnc in ${libnames_custom_txt}; do
                custom_files1=\$(ls ${abs_custom}/\${lnc}*_R1*.fastq.gz | sort -V | sed 's/^/-x /' | tr '\\n' ' ' )
                custom_files2=\$(ls ${abs_custom}/\${lnc}*_R2*.fastq.gz | sort -V | sed 's/^/-X /' | tr '\\n' ' ' )
                custom_types=\$(ls ${abs_custom}/\${lnc}*_R1*.fastq.gz | sed 's/^/-N \${types_custom_txt[\$idx]}\\t/' | cut -f1 | tr '\\n' ' ' )
                idx=\$(( \$idx + 1 ))
            done
            custom_files1=" \${custom_files1}"
            custom_files2=" \${custom_files2}"
            custom_types=" \${custom_types}"
        fi
    fi
    ${demux_species} -d -o .${extra}\${rna_files1}\${rna_files2}\${atac_files1}\${atac_files2}\${atac_files3}\${custom_files1}\${custom_files2}\${custom_types}
    ${plot} .
    """
}

// For demultiplexing an entire file of F/R reads
process demux_rna_reads{
    input:
    tuple val(libname), file(R1), file(R2), val(libname2), file(dists), file(assn), file(assn_filt), file(pdf), file(png), file(counts), file(names)
    
    publishDir "${params.output_directory}/${libname}", mode: "copy"
    output:
    tuple file("*/GEX_${R1}"), file("*/GEX_${R2}")

    script:
    """
    ${demux_species} -o . -w ${abs_wl} -T ${params.threads} -r ${R1} -R ${R2}
    """
}

// For demultiplexing a F/R file chunk
process demux_rna_reads_chunk{
    input:
    tuple val(libname), val(file_idx), file(R1), file(R2), val(libname2), file(dists), file(assn), file(assn_filt), file(pdf), file(png), file(counts), file(names)
    
    output:
    tuple val(libname), val(file_idx), file("*/GEX_${R1}"), file("*/GEX_${R2}")

    script:
    """
    ${demux_species} -o . -w ${abs_wl} -T ${params.threads} -r ${R1} -R ${R2}
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

workflow{
    
    // Read library IDs from params.libs
    library_ids = Channel.fromPath(params.libs).splitText().map { it.trim() } | view()  
    
    // Get all RNA-seq read pairs 
    def patterns = cartesian_product(suffixes, extensions).collect { 
        "${params.rna_dir}/${it.join()}"
    }     
    rna_pairs = Channel.fromFilePairs(patterns)
    
    // Return only RNA-seq read pairs with library IDs in the expected list
    // Also append a 0-based file index (which gives the number of the file
    // for that specific library), to help set batch IDs
    def file_idx = -1 
    indexed_rna_pairs = rna_pairs.join(library_ids).map { tup  ->
        file_idx++
        return [file_idx, tup[0], tup[1]]
    }
    
    rsplit = split_rna_reads(indexed_rna_pairs).transpose()
    
    // Group files by library name
    filestojoin = count_kmers(rsplit).groupTuple()
    //filestojoin = (split_rna_reads(indexed_rna_pairs).transpose() | count_kmers).groupTuple()
    
    // Combine counts
    joined = join_counts(filestojoin)

    // Fit the model, assign species, and plot
    models = fit_model(joined) 
    
    if (params.demux_reads){
        
        if (params.demux_pieces){
            //modMerge = models.collect().flatten().filter{ !it.toString().endsWith(".library") }.collate(12)
            modMerge = models.flatten().filter{ !it.toString().endsWith(".library") }.collate(12)
            
            rsplitMerge = rsplit.collect().flatten().collate(4).map{ tup -> [tup[1], tup[0], tup[2], tup[3]]}
            modMerge.cross(rsplitMerge).map{ tup ->
                [tup[1][0] + "_" + tup[1][1], tup[1][2], tup[1][3], tup[0][1], tup[0][5], tup[0][6]]
            }.groupTuple().transpose().view()

            //rsplit.collect().collate(4).map{ tup -> [tup[1], tup[0], tup[2], tup[3]] } | view
            //rsplit.collect().set().map{ tup -> [tup[1], tup[0], tup[2], tup[3]]}.view()
            /*
            rsplit.map{ tup -> [tup[1], tup[0], tup[2], tup[3]]}\
                .cross(models)\
                .view()
            */
            /*
            rsplit.map{ tup -> [tup[1], tup[0], tup[2], tup[3]]}\
                .cross(models)\
                .flatten()\
                .filter{ !it.toString().endsWith(".library") }\
                .collate(12)\
                .collect()\
                .map{ tup -> [ tup[0] + '_' + tup[1], tup[2][0], tup[3][0] ] }\
                .groupTuple()\
                .view()
            */
            //concat_rna(demuxed_chunks)
        }
        else{
            // Also need to re-obtain input read files
            def patterns2 = cartesian_product(suffixes, extensions).collect { 
                "${params.rna_dir}/${it.join()}"
            }     
            rna_pairs2 = Channel.fromFilePairs(patterns)

            demux_rna_reads(rna_pairs2.cross(models).flatten().filter{ !it.toString().endsWith(".library") }.collate(11))
            
            if (params.atac_dir){

            }
            if (params.custom_dir){

            }
        }
    }
     

    //Channel.fromPath("${params.output_directory}/*").filter { f -> 
    //    f.isDirectory()
    //}.set { subdirectories }
    
    //join_counts(library_ids)
}

