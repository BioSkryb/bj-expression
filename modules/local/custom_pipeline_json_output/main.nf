params.publish_dir = ""
custom_json = "597246834581.dkr.ecr.us-east-1.amazonaws.com/miscellaneous:nfrnaseq_json_output_0.2"
params.timestamp = ""

process CUSTOM_PIPELINE_JSON_OUTPUT {
    echo true
    tag 'pipeline_output'
    label 'process_low'
    container custom_json
    publishDir "${publish_dir}_${params.timestamp}/EXECUTION_INFO/", enabled:"$enable_publish"

    
    input:
    path(files)
    val(publish_dir)
    val(enable_publish)


    output:
    path "pipeline_json_output*", emit: file
    path("custom_pipeline_json_output_version.yml"), emit: version

    
    
    script:
    """
   
    python3 /scripts/custom_json_output_rnaseq.py -b '*.bam*' \\
                                           -f '*.fastq.gz' \\
                                           -csh 'df_gene_star*htseq_counts.tsv' \\
                                           -cstt 'df_tx_salmon_outdir_*tsv' \\
                                           -cstg 'df_gene_salmon_outdir_*tsv' \\
                                           -msh 'df_gene_counts_starhtseq.tsv' \\
                                           -mstt 'df_transcript_counts_salmon.tsv' \\
                                           -mstg 'df_gene_counts_salmon.tsv' \\
                                           -mshmt 'df_mt_gene_counts_starhtseq.tsv' \\
                                           -msmxt 'df_mt_gene_counts_salmon.tsv' \\
                                           --multiqc_report "multiqc_report.html" \\
                                           -o 'pipeline_json_output.json' \\
                                           -p "${publish_dir}" \\
                                           --run_name ${workflow.runName} \\
                                           --session_id ${workflow.sessionId} \\
                                           --pipeline_name ${workflow.manifest.name} \\
                                           --pipeline_version ${workflow.manifest.version} 
                                           
    echo custom_pipeline_json: v0.0.1 > custom_pipeline_json_output_version.yml
   
    
    """
    
}