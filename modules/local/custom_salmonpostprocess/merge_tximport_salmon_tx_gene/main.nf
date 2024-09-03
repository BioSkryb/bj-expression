nextflow.enable.dsl=2
params.timestamp = ""

process MERGE_TXIMPORT_SALMON_TX_GENE {
   tag "merge_tximport_salmon_tx_gene"
   publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/quantification_salmon/", enabled:"$enable_publish"
    
    input:
    path(df_files)
    val(publish_dir)
    val(enable_publish)
  

    output:
    path("df*_salmon.tsv"), emit: quant_merge_tsv
  

    script:
    """
   

    
    echo -e "File\\ttranscript_id\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\tcountsFromAbundanceNo\\tcountsFromAbundancescaledTPM\\tcountsFromAbundancelengthScaledTPM\\tcountsFromAbundancedtuScaledTPM\\ttranscript_biotype\\tgene_id\\tgene_biotype\\tgene_symbol" > merged_df_tx_salmon.txt
    cat df_tx* | grep -v "^File" >> merged_df_tx_salmon.txt
    
    echo -e "File\\tgene_id\\tcountsFromAbundanceNo\\tcountsFromAbundancescaledTPM\\tcountsFromAbundancelengthScaledTPM\\tgene_biotype\\tgene_symbol\\tgene_symbol_gene_id"> merged_df_gene_salmon.txt
    cat df_gene* | grep -v "^File" >> merged_df_gene_salmon.txt
    
    echo -e "File\\tgene_biotype\\ttranscript_biotype\\tNumFeatures\\tPropFeatures\\tcountsFromAbundanceNo\\tPropcountsFromAbundanceNo" > merged_df_num_detected_tx.txt
    cat df_sum_detected_tx_* | grep -v "^File" >> merged_df_num_detected_tx.txt
    
    echo -e "File\\tgene_biotype\\tNumFeatures\\tPropFeatures\\tcountsFromAbundanceNo\\tPropcountsFromAbundanceNo" > merged_df_num_detected_gene.txt
    cat df_sum_detected_gene_* | grep -v "^File" >> merged_df_num_detected_gene.txt
    
    echo -e "File\\tTotalFeatures\\tMT_NumFeatures\\tMT_Counts\\tTotal_Counts\\tPropMT" >  df_mt_gene_counts_salmon.tsv
    cat df_mtcounts_* | grep -v "^File" >>  df_mt_gene_counts_salmon.tsv
    
    mv merged_df_tx_salmon.txt df_transcript_counts_salmon.tsv
    mv merged_df_gene_salmon.txt df_gene_counts_salmon.tsv
    mv merged_df_num_detected_tx.txt  df_transcript_types_detected_summary_salmon.tsv
    mv merged_df_num_detected_gene.txt  df_gene_types_detected_summary_salmon.tsv
   
   
    """
}


workflow MERGE_TXIMPORT_SALMON_TX_GENE_WF{
    take:
        ch_df_files
        ch_publish_dir
        ch_enable_publish
        
    main:
        MERGE_TXIMPORT_SALMON_TX_GENE (
                       ch_df_files,
                       ch_publish_dir,
                       ch_enable_publish
                     )
                
    
    emit:
        
        quant_merge_tsv= MERGE_TXIMPORT_SALMON_TX_GENE.out.quant_merge_tsv
   
}

workflow{
   ch_df_files = Channel.from( params.df_files)
            
   MERGE_TXIMPORT_SALMON_TX_GENE (  ch_df_files,
                                    params.tx2gene, 
                                    params.publish_dir, 
                                    params.enable_publish)
                     
}
