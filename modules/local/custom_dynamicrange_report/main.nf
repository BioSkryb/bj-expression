nextflow.enable.dsl=2
params.timestamp = ""

process CALC_DYNAMICRANGE {
    tag "CALC_DYNAMICRANGE"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/secondary_metrics/", enabled:"$enable_publish"
    
    input:
    path(count_matrix)
    val(publish_dir)
    val(enable_publish)

    output:
    path ("df_dynamicrange_expression.tsv")

    script:
    """
    
    
    Rscript /usr/local/bin/scrnaseq_dynmicrange.R ${count_matrix}
    
    """
}
