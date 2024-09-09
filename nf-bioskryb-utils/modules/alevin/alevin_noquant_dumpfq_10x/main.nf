nextflow.enable.dsl=2
params.timestamp = ""

process ALEVIN_NOQUANT_DUMPFQ {
  tag "alevin_noquant_dumpfastq_${sample_name}"
  publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
  
  input:
  tuple val(sample_name), path(reads)
  path (salmonindex)
  path (tx2gene)
  val (lib_type)
  val (lib_protocol)
  val (publish_dir)
  val (enable_publish)
  

  output:
  path("alevin_dumpfastq_outdir"), emit: outdir
  tuple val(sample_name), path("*.cbc.fastq"), emit:fastq_cbc
  path("salmon_version.yml"), emit: version

  
  script:
 
 
  """
  
  cat ${tx2gene} | cut -f1,3 | sort -u > df_tx2gene.tsv
  
  mlib_type=`echo -e "--${lib_protocol}"`;
  
  salmon alevin -l ${lib_type} -1  ${reads[0]} -2 ${reads[1]}  --noQuant --dumpfq -i index  -p ${task.cpus} -o alevin_dumpfastq_outdir --tgMap df_tx2gene.tsv \${mlib_type}  > ${sample_name}_L001_R1_001.cbc.fastq

    
  export SALMON_VER=\$(salmon --version | sed -e "s/salmon //g")
  echo salmon: \$SALMON_VER > salmon_version.yml   
    
  """

}

workflow ALEVIN_NOQUANT_DUMPFQ_WF{
    take:
        ch_reads
        ch_salmonindex
        ch_tx2gene
        ch_lib_type
        ch_lib_protocol
        ch_publish_dir
        ch_enable_publish
        

  
        
        
    main:
        ALEVIN_NOQUANT_DUMPFQ (
                       ch_reads,
                       ch_salmonindex,
                       ch_tx2gene,
                       ch_lib_type,
                       ch_lib_protocol,
                       ch_publish_dir,
                       ch_enable_publish
                     )
                     
    emit:
        version = ALEVIN_NOQUANT_DUMPFQ.out.version
        fastq_cbc = ALEVIN_NOQUANT_DUMPFQ.out.fastq_cbc
        outdir = ALEVIN_NOQUANT_DUMPFQ.out.outdir
    
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], subtags[2], subtags[3], pair] }
    ALEVIN_NOQUANT_DUMPFQ_WF  (
                        ch_reads,
                        params.salmon_index,
                        params.tx2gene,
                        params.lib_protocol_10x,
                        params.lib_type_10x,
                        params.publish_dir,
                        params.enable_publish
                     )
}