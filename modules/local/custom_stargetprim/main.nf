nextflow.enable.dsl=2
params.timestamp = ""


process STARGETPRIM {
  tag "${sample_name}"
  publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/alignment_htseq/", enabled:"$enable_publish", pattern: "${sample_name}.bam*"
  publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/alignment_htseq/", enabled:"$enable_publish", pattern: "${sample_name}_Chimeric.out.junction"
  
  
  input:
  tuple val (sample_name), path(bam), path(junction)
  val(publish_dir)
  val(enable_publish)

  output:
  tuple val(sample_name), path("${sample_name}.raw.bam"), emit: raw_bam
  tuple val (sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"),  emit: primary_bundle
  tuple val(sample_name), path("${sample_name}_Chimeric.out.junction"), emit: junction

  path("samtools_version.yml"), emit: version
  
  script:
  """
 
  samtools index ${bam}
  samtools view -h ${bam}  | sed 's/^\\(@RG\\tID:.*\\)/\\1\\tSM:${sample_name}\\tPL:Illumina/' | samtools view -h -b - | samtools sort -O BAM -o ${sample_name}.raw.bam
  samtools view -F 256 -h -o ${sample_name}.bam ${sample_name}.raw.bam
  samtools index ${sample_name}.bam
  
  export SAMTOOLS_VER=\$(samtools --version 2>&1 |  sed -n -e '1p' | grep -Eo [0-9][.]*[0-9]*)
  echo Samtools: \$SAMTOOLS_VER > samtools_version.yml

  touch ${junction}
  """
	
}

workflow  STARGETPRIM_WF{
    take:
        ch_bam
        ch_publish_dir
        ch_enable_publish

    main:
        STARGETPRIM (
                       ch_bam,
                       ch_publish_dir,
                       ch_enable_publish
                     )

    emit:
        version = STARGETPRIM.out.version
        primary_bam= STARGETPRIM.out.primary_bundle
        bam_raw = STARGETPRIM.out.raw_bam
        junction = STARGETPRIM.out.junction


}

workflow{
   
   Channel
            .fromFilePairs( params.bam_dir + '/*/*.{bam,bai}', size: -1 )
                          { file -> file.name.replaceAll(/.bam|.bai$/,'') }
            .mix ( Channel
                        .fromFilePairs( params.bam_dir + '/*.{bam,bai}', size: -1 )
                        { file -> file.name.replaceAll(/.bam|.bai$/,'') } )
            .set { ch_bam }
   
    STARGETPRIM_WF(
                       ch_bam,
                       params.publish_dir,
                       params.enable_publish
                     )
}
