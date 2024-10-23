nextflow.enable.dsl=2

process MERGE_HTSEQ_SUMMARY {
    tag "merge_htseq_summary"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/quantification_htseq/", enabled:"$enable_publish"
    
    input:
    path(df_files)
    val(publish_dir)
    val(enable_publish)

    output:
    path("df*starhtseq.tsv"), emit: merge_tsv
    
 
    script:
    """

    
    echo -e "File\\tgene_id\\tcountHTSeq\\tgene_biotype\\tgene_symbol\\tgene_symbol_gene_id"> df_gene_counts_starhtseq.tsv
    cat df_gene_star_htseq* | grep -v "^File" >> df_gene_counts_starhtseq.tsv
    
    echo -e "File\\tgene_biotype\\tNumFeatures\\tPropFeatures\\tcountHTSeq\\tPropcountHTSeq" >  df_gene_types_detected_summary_starhtseq.tsv
    cat df_sum_detected_gene_* | grep -v "^File" >>  df_gene_types_detected_summary_starhtseq.tsv
    
    
    echo -e "File\\tTotalFeatures\\tMT_NumFeatures\\tMT_Counts\\tTotal_Counts\\tPropMT" >  df_mt_gene_counts_starhtseq.tsv
    cat df_mtcounts_star_htseq_* | grep -v "^File" >>  df_mt_gene_counts_starhtseq.tsv
    
    
    """
}
  


workflow MERGE_HTSEQ_SUMMARY_WF{
    take:
        ch_df_files
        ch_publish_dir
        ch_enable_publish
        
    main:
        MERGE_HTSEQ_SUMMARY (
                       ch_df_files,
                       ch_publish_dir,
                       ch_enable_publish
                     )
    emit:
        
       
        merge_tsv = MERGE_HTSEQ_SUMMARY.out.merge_tsv
}

workflow{
   ch_df_files = Channel.from( params.df_files)
            
   MERGE_HTSEQ_SUMMARY_WF (  ch_df_files,
                              params.publish_dir, 
                              params.enable_publish)
                     
}
