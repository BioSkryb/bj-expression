nextflow.enable.dsl=2
params.timestamp = ""

process CREATE_MATRIX_SALMON_TX_GENE {
    tag "create_matrix_salmon_tx_gene"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/quantification_salmon/", enabled:"$enable_publish"


    input:
    path(merged_df)
    val(publish_dir)
    val(enable_publish)

    output:
    path("*"), emit: salmon_matrix
    

    script:
    """
    Rscript /usr/local/bin/salmon_summarydf_to_matrix.R
    
    """
}
workflow CREATE_MATRIX_SALMON_TX_GENE_WF{
    take:
        ch_salmon_files
        ch_publish_dir
        ch_enable_publish
        
    main:
        CREATE_MATRIX_SALMON_TX_GENE (
                       ch_salmon_files,
                       ch_publish_dir,
                       ch_enable_publish
                     )
    emit:
        salmon_matrix = CREATE_MATRIX_SALMON_TX_GENE.out.salmon_matrix
        
    
}

workflow{
   ch_salmon_files = Channel.from( params.salmon_files)
            
   CREATE_MATRIX_SALMON_TX_GENE_WF  ( ch_salmon_files,
                                    params.publish_dir, 
                                    params.enable_publish)
                     
}

