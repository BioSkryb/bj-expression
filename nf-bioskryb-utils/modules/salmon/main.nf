nextflow.enable.dsl=2
params.timestamp = ""

process SALMONQUANT {
  tag "${sample_name}"
  publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
  
  input:
  tuple val(sample_name), path(reads)
  path salmonindex
  val(publish_dir)
  val(enable_publish)
  

  output:
  
 // tuple val(sample_name), val(sample_number), val(lane_number), path("*.salmon.sam"), emit: sam
  path("salmon_outdir_${sample_name}"), emit: outdir
  path("salmon_version.yml"), emit: version
  
  script:
 
 
  reads_paired_cmd = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${reads[0]}"
 
  """
  
    salmon quant \\
    -i ${salmonindex} \\
    -l A \\
    --minAssignedFrags 1 \\
     ${reads_paired_cmd} \\
    -o salmon_outdir_${sample_name} \\
    -p $task.cpus \\
    -z ${sample_name}.salmon.sam 
    
    export SALMON_VER=\$(salmon --version | sed -e "s/salmon //g")
    echo Salmon: \$SALMON_VER > salmon_version.yml   
    
  """

}

workflow SALMON_WF{
    take:
        ch_reads
        ch_salmonindex
        ch_publish_dir
        ch_enable_publish
        
    main:
        SALMONQUANT (
                       ch_reads,
                       ch_salmonindex,
                       ch_publish_dir,
                       ch_enable_publish
                     )
                     
    emit:
        version = SALMONQUANT.out.version
        sam = SALMONQUANT.out.sam
        outdir = SALMONQUANT.out.outdir
    
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], subtags[2], subtags[3], pair] }
    SALMON_WF (
                        ch_reads,
                        params.salmon_index,
                        params.publish_dir,
                        params.enable_publish
                     )
}