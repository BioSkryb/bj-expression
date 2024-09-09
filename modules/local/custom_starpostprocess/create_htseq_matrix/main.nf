nextflow.enable.dsl=2
params.timestamp = ""

process CREATE_HTSEQ_MATRIX {
    tag "create_htseq_matrix"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/quantification_htseq/", enabled:"$enable_publish"


    input:
    path(merged_df)
    val(publish_dir)
    val(enable_publish)

    output:
    path("matrix*"), emit: htseq_matrix
    
    script:
    """
    Rscript /usr/local/bin/htseq_summarydf_to_matrix.R
    
    """
}
workflow CREATE_HTSEQ_MATRIX_WF{
    take:
        ch_htseq_file
        ch_publish_dir
        ch_enable_publish
        
    main:
        CREATE_HTSEQ_MATRIX (
                       ch_htseq_file,
                       ch_publish_dir,
                       ch_enable_publish
                     )
    emit:
        htseq_matrix = CREATE_HTSEQ_MATRIX.out.htseq_matrix
        
    
}

workflow{
   ch_htseq_file = Channel.from( params.htseq_file)
            
   CREATE_HTSEQ_MATRIX_WF  ( ch_htseq_file,
                                    params.publish_dir, 
                                    params.enable_publish)
                     
}

