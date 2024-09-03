nextflow.enable.dsl=2
params.timestamp = ""

process CELL_TYPING {
    tag "CELL_TYPING"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/classification_cell_typing/", enabled:"$enable_publish"
    
    input:
    path(tsv_files)
    path(ref_files)
    val(publish_dir)
    val(enable_publish)

    output:
    path("df*.tsv"), emit: annotations
//  path("seurat.rds"), emit: seurat_object

    script:
    """
    
    Rscript /usr/local/bin/scrnaseq_celltype_refactored_pipeline.R
    
    """
}



