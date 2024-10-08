nextflow.enable.dsl=2
params.timestamp = ""


process CUSTOM_AWK_DEMUX_CBC_FASTQ {
    tag "CUSTOM_AWK_DEMUX_CBC_FASTQ_${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"



    input:
    tuple val(sample_name), path(fastq_cbc)
    val(publish_dir)
    val(enable_publish)


    output:
    path("*.demux.fastq"), emit: fastqs
    
    script:
    
    """
    
    cat ${fastq_cbc} | paste - - - - | awk -v FS="\\t" -v OFS="\\n" 'BEGIN { reader=0; } { reader++; if(reader % 50000 == 0){print "Reader " reader }; split(\$1,a,"_");barcode=a[2];print \$1,\$2,\$3,\$4>> barcode".demux.fastq" ; close(barcode".demux.fastq")  }'

    
    
    """
}

