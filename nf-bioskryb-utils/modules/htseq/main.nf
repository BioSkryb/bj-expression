nextflow.enable.dsl=2
params.timestamp = ""

process HTSEQ_COUNTS {
  tag "${sample_name}"
  publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$disable_publish"

  input:
  tuple val(sample_name), path(bam), path(bai)
  path(gtf_file)
  val(publish_dir)
  val(enable_publish)

  output:
  path("${sample_name}.htseq_counts.tsv"), emit: htseq_counts
  path("htseq_version.yml"), emit: version

  script:
  """
  htseq-count \\
  -f bam -r pos -s no -t exon -i gene_id --additional-attr=gene_name -m union \\
  --secondary-alignments=ignore \\
  ${bam} ${gtf_file} > ${sample_name}.htseq_counts.tsv
  
  export HTSEQ_VER=\$(htseq-count . . --version)
  echo HTSeq: \$HTSEQ_VER > htseq_version.yml  
  """
}

workflow HTSEQ_COUNTS_WF{
    take:
        ch_bam
        ch_gtf
        ch_publish_dir
        ch_enable_publish
        
    main:
        HTSEQ_COUNTS (
                       ch_bam,
                       ch_gtf,
                       ch_publish_dir,
                       ch_enable_publish
                     )
                     
    emit:
        version = HTSEQ_COUNTS.out.version
        htseq_counts = HTSEQ_COUNTS.out.htseq_counts
        
    
}

workflow{
   Channel
            .fromFilePairs( params.bam_dir + '/*/*.{bam,bai}', size: -1 )
                          { file -> file.name.replaceAll(/.bam|.bai$/,'') }
            .mix ( Channel
                        .fromFilePairs( params.bam_dir + '/*.{bam,bai}', size: -1 )
                        { file -> file.name.replaceAll(/.bam|.bai$/,'') } )
            .set { ch_bam }
            
    HTSEQ_COUNTS_WF (
                        ch_bam,
                        params.gtf_file,
                        params.publish_dir,
                        params.enable_publish
                     )
}