nextflow.enable.dsl=2
params.timestamp = ""


process CUSTOM_FASTQ_MERGE {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/output/${sample_name}/", enabled:"$enable_publish"
    


    input:
    tuple val(sample_name), val(sample_number), val(lane_number), path(reads)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), path("${sample_name}_R1.fastq.gz"), path("${sample_name}_R2.fastq.gz"), emit: reads
    path("custom_fastq_merge_version.yml"), emit: version
    //path("*.fastq.gz"), emit: merged_fastqs_to_compress
    
    script:
    
    """
    # echo ${sample_number.join('-')}
    # echo ${lane_number.join('-')}
    
    cat *R1_001* > ${sample_name}_R1.fastq.gz
    cat *R2_001* > ${sample_name}_R2.fastq.gz
    
    echo custom_fastq_merge: v0.0.1 > custom_fastq_merge_version.yml
    """
}

workflow CUSTOM_FASTQ_MERGE_WF{
    take:
        ch_publish_dir
        ch_enable_publish
        ch_reads
        
    main:
        CUSTOM_FASTQ_MERGE ( 
                            ch_reads,
                            ch_publish_dir,
                            ch_enable_publish
                           )
                           
    emit:
        reads = CUSTOM_FASTQ_MERGE.out.reads
        version = CUSTOM_FASTQ_MERGE.out.version
        //merged_fastqs_to_compress = CUSTOM_FASTQ_MERGE.out.merged_fastqs_to_compress
        
    
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], subtags[2], subtags[3], pair] }
                            
    CUSTOM_FASTQ_MERGE_WF(
                            params.publish_dir,
                            params.enable_publish,
                            ch_reads
                        )
}
