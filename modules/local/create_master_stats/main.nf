nextflow.enable.dsl=2

process CREATE_MASTER_STATS {
    tag "create_master_stats"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$disable_publish"

    input:
    path(df_files)
    path(input_csv)
    val(publish_dir)
    val(disable_publish)

    output:
    path("overall_stats_mqc.csv"), emit: stats
    
 
    script:
    """
    
    
    Rscript  /usr/local/bin/create_master_stats_scrnaseqwf.R ${input_csv}

    
    """
}
  
