nextflow.enable.dsl=2

process GENE_BODY_COVERAGE_RNA_PLOT {
  tag "GENE_BODY_COVERAGE_RNA_PLOT"
  publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$disable_publish"

  input:
  path (df_files)
  val (disable_publish)


  output:
  path("*.png"), emit: gene_body_png

  script:
  """
  
  cat df_sum* | grep Percentile | head -n1  > header
  cat df_sum* | grep -v Percentile | cat header -  > df_all.tsv
  geneBodyCoveragePlot.py  
 
  
  """
}
