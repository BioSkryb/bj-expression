nextflow.enable.dsl=2
params.timestamp = ""

process PETASUITE_DECOMPRESS {
    tag "${compressed_file}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    path compressed_file
    val(publish_dir)
    val(enable_publish)
    
    output:
    path("*.fastq.gz"), emit: fastq, optional: true
    path("*.bam"), emit: bam, optional: true
    path("petasuite_version.yml"), emit: version
    
    script:
    """


    petasuite --decompress --numthreads $task.cpus --md5match *
        
    export PETASUITE_VER=\$(petasuite --version 2>&1 | head -n1 | sed 's/^.*version //;s/ .*\$//')
    echo PETASUITE DECOMPRESS: \$PETASUITE_VER > petasuite_version.yml
    """
}

workflow PETASUITE_DECOMPRESS_WF {
    take:
        ch_files
        ch_publish_dir
        ch_enable_publish
    main:
        PETASUITE_DECOMPRESS (
                                ch_files,
                                ch_publish_dir,
                                ch_enable_publish
                            )
    emit:
        version = PETASUITE_DECOMPRESS.out.version
        fastq = PETASUITE_DECOMPRESS.out.fastq
        bam = PETASUITE_DECOMPRESS.out.bam
}

workflow{
    ch_pgcompressed = Channel.fromPath(params.reads)
    ch_pgcompressed.view()
    
    PETASUITE_DECOMPRESS_WF(
                            ch_pgcompressed,
                            params.publish_dir,
                            params.enable_publish
                        )
                        
    ch_fastq = PETASUITE_DECOMPRESS_WF.out.fastq
                    .map{ tag -> subtags = ( tag =~ /.*\/(.*)_(S\d+)_(L+\d+)/)[0]; [ subtags[1], subtags[2], subtags[3], tag ] }
                    .groupTuple()
                    .map{ tag, sample, lane, pair -> [ tag, sample[1], lane[1], pair ] }

    ch_fastq.view()
    
}