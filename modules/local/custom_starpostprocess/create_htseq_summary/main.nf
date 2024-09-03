nextflow.enable.dsl=2
r_nf_rnaseq_cntr = "597246834581.dkr.ecr.us-east-1.amazonaws.com/miscellaneous:custom_r_nf_rnaseq_0.8"
params.timestamp = ""

process CREATE_HTSEQ_SUMMARY {
    tag "create_htseq_summary"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    container r_nf_rnaseq_cntr

    input:
    path(htseq_file)
    path(tx2gene)
    val(publish_dir)
    val(enable_publish)

    output:
    path("*tsv"), emit: htseq_summary
    

    script:
    """
    echo -e "Working on ${htseq_file}"
    
    Rscript /usr/local/bin/htseq_to_gene_df.R ${htseq_file} ${tx2gene}
    """
}
workflow CREATE_HTSEQ_SUMMARY_WF{
    take:
        ch_htseq_file
        ch_tx2gene
        ch_publish_dir
        ch_enable_publish
        
    main:
        CREATE_HTSEQ_SUMMARY (
                       ch_htseq_file,
                       ch_tx2gene,
                       ch_publish_dir,
                       ch_enable_publish
                     )
    emit:
        htseq_summary = CREATE_HTSEQ_SUMMARY.out.htseq_summary
        
    
}

workflow{
   ch_htseq_file = Channel.from( params.htseq_file)
            
   CREATE_HTSEQ_SUMMARY_WF  ( ch_htseq_file,
                                    params.tx2gene, 
                                    params.publish_dir, 
                                    params.enable_publish)
                     
}

