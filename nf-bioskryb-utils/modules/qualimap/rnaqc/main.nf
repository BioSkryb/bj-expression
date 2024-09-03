nextflow.enable.dsl=2
params.publish_dir = ""

process QUALIMAP_BAMRNA {
  tag "${sample_name}"
  publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$disable_publish"

    
  input:
  tuple val (sample_name), path(bam)
  path(gtf_file)
  val(publish_dir)
  val(enable_publish)
    
  output:
  path("qualimap_outdir_${sample_name}"), emit: outdir
  path("qualimap_outdir_${sample_name}/rnaseq_qc_results.txt"), emit: rnaseq_qc
  path("qualimap_version.yml"), emit: version
    
  script:
    
  """
  unset DISPLAY
  mkdir tmp
  export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
  #JAVA_MEM_SIZE=${task.memory.toGiga()}G
  #java_options="-Djava.awt.headless=true -Xmx\$JAVA_MEM_SIZE -XX:MaxPermSize=1024m"
  #java_options="-Djava.awt.headless=true -XX:MaxPermSize=1024m"
  #qualimap --java-mem-size=${task.memory.toGiga()}G
  
  qualimap \\
  --java-mem-size=${task.memory.toGiga()}G \\
  rnaseq \\
  -gtf ${gtf_file} \\
  -bam ${bam} \\
  --sorted \\
  -outfile ${sample_name}.report.pdf \\
  -oc ${sample_name}.oc \\
  -outdir qualimap_outdir_${sample_name}
  
 
  echo QualiMap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//') > qualimap_version.yml
  
  """
	
}
workflow QUALIMAP_BAMRNA_WF{
    take:
        ch_bam
        ch_gtf
        ch_publish_dir
        ch_enable_publish
        
    main:
        QUALIMAP_BAMRNA (
                       ch_bam,
                       ch_gtf,
                       ch_publish_dir,
                       ch_enable_publish
                     )
                     
     emit:
        results = QUALIMAP_BAMRNA.out.outdir
        version = QUALIMAP_BAMRNA.out.version
        
    
}

workflow{
   Channel
            .fromFilePairs( params.bam_dir + '/*/*.{bam,bai}', size: -1 )
                          { file -> file.name.replaceAll(/.bam|.bai$/,'') }
            .mix ( Channel
                        .fromFilePairs( params.bam_dir + '/*.{bam,bai}', size: -1 )
                        { file -> file.name.replaceAll(/.bam|.bai$/,'') } )
            .set { ch_bam }
            
    QUALIMAP_BAMRNA (
                        ch_bam,
                        params.gtf_file,
                        params.publish_dir,
                        params.enable_publish
                     )
}