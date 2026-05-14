nextflow.enable.dsl=2

process CREATE_MASTER_STATS {
    tag "create_master_stats"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$disable_publish"

    input:
    path(df_files)
    path(reads_csv)
    path(insert_size_summary)
    val(large_files)
    val(publish_dir)
    val(disable_publish)

    output:
    path("SelectedMetrics_mqc.csv"), emit: stats
    
    script:
    if (large_files) {
        """
        Rscript /usr/local/bin/create_master_stats_scrnaseqwf.R ${reads_csv} ${insert_size_summary}
        """
    } else {
        """
        touch SelectedMetrics_mqc.csv
        echo "Skipping process as large_files is set to false."
        """
    }
}