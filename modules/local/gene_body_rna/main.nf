nextflow.enable.dsl=2

process GENE_BODY_COVERAGE_RNA {
  tag "${sample_name}_GENE_BODY_COVERAGE_RNA"
  publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$disable_publish"


  input:
  tuple val(sample_name), path (bam), path(bai)
  path (reference_bed)
  val (disable_publish)


  output:
  path("df_*.tsv"), emit: df
  path("*.png"), emit: gene_body_png

  script:
  """
  
  geneBodyCoverageIsai.py -r ${reference_bed} -i .
  
  mv df_sum.geneBodyCoverage.tsv df_sum_${sample_name}.geneBodyCoverage.tsv
  mv df_all.genebodypercentile.tsv df_all_${sample_name}.genebodypercentile.tsv
  mv df_skewness.genebodypercentile.tsv df_skewness_${sample_name}.genebodypercentile.tsv 
  mv transcript_body_coverage_mqc.png transcript_body_coverage_mqc_${sample_name}.png
  
  """
}
