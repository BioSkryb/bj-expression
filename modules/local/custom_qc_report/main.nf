nextflow.enable.dsl=2
params.timestamp = ""

process CREATE_QC_REPORT {
    tag "QC Report"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/secondary_metrics/", enabled:"$enable_publish"
    
    input:
    path(qualimap_outdir)
    val(publish_dir)
    val(enable_publish)

    output:
    path("pipeline_metrics_summary_percents.csv"), emit: pstats
    path("qualimap_stats_mqc.csv"), emit: stats
   

    script:
    """
    echo -e "Working on ${qualimap_outdir}"
    Rscript /usr/local/bin/parse_qualimap.R 
    
    """
}

workflow CREATE_QC_REPORT_WF{
    take:
        ch_qualimap_outdir
        ch_publish_dir
        ch_enable_publish
        
    main:
        CREATE_QC_REPORT (
                       ch_qualimap_outdir,
                       ch_publish_dir,
                       ch_enable_publish
                     )
    
}

workflow{
   ch_qualimap_outdir = Channel.from(params.qualimap_outdir)
            
   CREATE_QC_REPORT (ch_qualimap_outdir,
                    params.publish_dir, 
                    params.enable_publish)
                     
}
