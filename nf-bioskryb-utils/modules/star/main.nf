nextflow.enable.dsl=2
params.timestamp = ""

process STARALIGN {
  tag "${sample_name}"
  publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$disable_publish"

  input:
  tuple val(sample_name), path(reads)
  path starindex
  val(publish_dir)
  val(enable_publish)

  output:
  path("star_outdir_${sample_name}"), emit: outdir
  tuple val(sample_name), path("star_outdir_${sample_name}/*Aligned.sortedByCoord.out.bam"), path("*Chimeric.out.junction"),  emit: bam
  tuple val(sample_name), path("star_outdir_${sample_name}/*Chimeric.out.junction"), emit: junction
  path("star_version.yml"), emit: version    

  script:
  
  reads_paired_cmd = reads.size() == 2 ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
  
  """
 
  mkdir "star_outdir_${sample_name}";
    
    STAR \\
    --genomeDir ${starindex} \\
    --runThreadN $task.cpus \\
    --outSAMunmapped Within KeepPairs \\
    --twopassMode Basic \\
    --outReadsUnmapped None \\
    --outSAMtype BAM SortedByCoordinate \\
    --outFileNamePrefix "star_outdir_${sample_name}/${sample_name}_" \\
    --readFilesCommand zcat \\
    --readFilesIn ${reads_paired_cmd}\\
    --outSAMstrandField intronMotif  \\
    --outSAMattributes All \\
    --outFilterMultimapNmax 10 \\
    --outSAMprimaryFlag OneBestScore \\
    --chimSegmentMin 12im \\
    --chimJunctionOverhangMin 8 \\
    --chimOutJunctionFormat 1 \\
    --outSAMattrRGline "ID:${sample_name}" 
    
    export STAR_VER=\$(STAR --version | sed -e "s/STAR_//g")
    echo STAR: \$STAR_VER > star_version.yml
    
   rm -r star_outdir_${sample_name}/${sample_name}__STARgenome
   rm -r star_outdir_${sample_name}/${sample_name}__STARpass1

  cp star_outdir_${sample_name}/*Chimeric.out.junction .
 
  """
	
}
workflow STAR_WF{
    take:
        ch_reads
        ch_starindex
        ch_publish_dir
        ch_enable_publish

    main:
        STARALIGN (
                       ch_reads,
                       ch_starindex,
                       ch_publish_dir,
                       ch_enable_publish
                     )

    emit:
        version = STARALIGN.out.version
        bam = STARALIGN.out.bam
        junction = STARALIGN.out.junction
        outdir = STARALIGN.out.outdir

}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], subtags[2], subtags[3], pair] }
    STAR_WF (
                        ch_reads,
                        params.star_index,
                        params.publish_dir,
                        params.enable_publish
                     )
}
