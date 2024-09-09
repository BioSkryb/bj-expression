nextflow.enable.dsl=2
params.timestamp = ""


process PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY {
    tag "plotter_pcaheatmap"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    input:
    path(tsv_files)
    val(publish_dir)
    val(enable_publish)

    output:
    path("*.png"), emit: pcaheatmap_plot
   
   
    script:
    """
    Rscript /usr/local/bin/pca_heatmap_qc_salmon_htseq_nometadata.R
    """
}

workflow PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY_WF {
    take:
        ch_tsv_files
        ch_publish_dir
        ch_enable_publish
        
    main:
        PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY (
                       ch_tsv_files,
                       ch_publish_dir,
                       ch_enable_publish
                     )
     emit:
        pcaheatmap_plot = PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY.out.pcaheatmap_plot
   
        
    
}

workflow{
   ch_tsv_files = Channel.from( params.tsv_files)
            
   PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY_WF  ( ch_tsv_files,
                                    params.publish_dir, 
                                    params.enable_publish)
                     
}
