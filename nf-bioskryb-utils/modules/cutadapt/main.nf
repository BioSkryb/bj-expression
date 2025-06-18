nextflow.enable.dsl=2
params.timestamp = ""

process CUTADAPT {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/primary_analyses/metrics/${sample_name}/cutadapt/", enabled:"$enable_publish"
    
    input:
    tuple val(sample_name), path(reads)
    val(adapter_r1)
    val(adapter_r2)
    val(publish_dir)
    val(enable_publish)

    
    output:
    tuple val(sample_name), path("*cutadapt.fastq.gz"), emit: reads
                                              
    """

    cutadapt -b ${adapter_r1} -B ${adapter_r2} -a ${adapter_r1} -o ${sample_name}_R1_cutadapt.fastq.gz -p ${sample_name}_R2_cutadapt.fastq.gz ${reads[0]} ${reads[1]}

    
    """
}

