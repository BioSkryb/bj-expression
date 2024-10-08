nextflow.enable.dsl=2
params.timestamp = ""

process CREATE_TXIMPORT_SALMON_TX_GENE {
    tag "CREATE_TXIMPORT_SALMON_TX_GENE"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    input:
    path(salmon_outdir)
    path(tx2gene)
    val(publish_dir)
    val(enable_publish)

    output:
    path("*tsv")

    script:
    """
    echo -e "Working on ${salmon_outdir}"
    Rscript /usr/local/bin/salmon_to_tx_gene_df.R ${salmon_outdir} ${tx2gene} 
    
    """
}

workflow CREATE_TXIMPORT_SALMON_TX_GENE_WF{
    take:
        ch_salmon_outdir
        ch_tx2gene
        ch_publish_dir
        ch_enable_publish
        
    main:
        CREATE_TXIMPORT_SALMON_TX_GENE (
                       ch_salmon_outdir,
                       ch_tx2gene,
                       ch_publish_dir,
                       ch_enable_publish
                     )
    
}

workflow{
   ch_salmon_outdir = Channel.from( params.salmon_outdir)
            
   CREATE_TXIMPORT_SALMON_TX_GENE ( ch_salmon_outdir,
                                    params.tx2gene, 
                                    params.publish_dir, 
                                    params.enable_publish)
                     
}
