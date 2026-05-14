/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//

include { PUBLISH_INPUT_DATASET_WF } from '../nf-bioskryb-utils/modules/bioskryb/publish_input_dataset/main.nf'
include { CUSTOM_FASTQ_MERGE_WF } from '../modules/local/custom_fastq_merge/main.nf'
// include { PETASUITE_DECOMPRESS_WF } from '../nf-bioskryb-utils/modules/petasuite/decompress/main.nf'
include { SEQTK_WF } from '../nf-bioskryb-utils/modules/seqtk/sample/main.nf'
include { FastpFull_WF } from '../nf-bioskryb-utils/modules/fastp/main.nf'
include { SALMONQUANT } from '../nf-bioskryb-utils/modules/salmon/main.nf'
include { CREATE_TXIMPORT_SALMON_TX_GENE } from '../modules/local/custom_salmonpostprocess/create_tximport_salmon_tx_gene/main.nf'
include { MERGE_TXIMPORT_SALMON_TX_GENE_WF } from '../modules/local/custom_salmonpostprocess/merge_tximport_salmon_tx_gene/main.nf'
include { STARALIGN } from '../nf-bioskryb-utils/modules/star/main.nf'
include { QUALIMAP_BAMRNA } from '../nf-bioskryb-utils/modules/qualimap/rnaqc/main.nf'
include { STARGETPRIM_WF } from '../modules/local/custom_stargetprim/main.nf'
include { HTSEQ_COUNTS } from '../nf-bioskryb-utils/modules/htseq/main.nf'
include { CREATE_HTSEQ_SUMMARY_WF } from '../modules/local/custom_starpostprocess/create_htseq_summary/main.nf'
include { MERGE_HTSEQ_SUMMARY_WF } from '../modules/local/custom_starpostprocess/merge_htseq_summary/main.nf'
include { PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY_WF } from '../modules/local/custom_starpostprocess/plotter_pcaheatmapqc_htseq_summary/main.nf'
include { CELL_TYPING } from '../modules/local/custom_cell_typing/main.nf'
include { CREATE_QC_REPORT } from '../modules/local/custom_qc_report/main.nf'
include { MULTIQC_WF } from '../nf-bioskryb-utils/modules/multiqc/main.nf'
include { REPORT_VERSIONS_WF } from '../nf-bioskryb-utils/modules/bioskryb/report_tool_versions/main.nf'
// include { CUSTOM_PIPELINE_JSON_OUTPUT } from '../modules/local/custom_pipeline_json_output/main.nf'
include { CREATE_HTSEQ_MATRIX } from '../modules/local/custom_starpostprocess/create_htseq_matrix/main.nf'
include { CREATE_MATRIX_SALMON_TX_GENE } from '../modules/local/custom_salmonpostprocess/create_matrix_salmon_tx_gene/main.nf'
include { GENE_BODY_COVERAGE_RNA } from '../modules/local/gene_body_rna/main.nf'
include { GENE_BODY_COVERAGE_RNA_PLOT } from '../modules/local/gene_body_rna_plot_all/main.nf'
include { CREATE_MASTER_STATS } from '../modules/local/create_master_stats/main.nf'
include { ALEVIN_NOQUANT_DUMPFQ } from  '../nf-bioskryb-utils/modules/alevin/alevin_noquant_dumpfq_10x/main.nf'
include { CUSTOM_AWK_DEMUX_CBC_FASTQ } from  '../modules/local/custom_10x_demux/awk_demux_cbc_fastq/main.nf'
include { CALC_DYNAMICRANGE } from '../modules/local/custom_dynamicrange_report/main.nf'
include { COUNT_READS_FASTQ_WF } from '../nf-bioskryb-utils/modules/bioskryb/custom_read_counts/main.nf'
include { CUTADAPT } from '../nf-bioskryb-utils/modules/cutadapt/main.nf'
include { RNA_QC_PLOTS } from '../nf-bioskryb-utils/subworkflows/qc_plots/main.nf'
include { INSERT_SIZE_ANALYSIS_WF } from '../modules/local/insert_size_analysis/main.nf'

//include { STARFUSION_WF } from  '../nf-bioskryb-utils/modules/starFusion/main.nf'
//include { MERGE_STAR_FUSION } from '../modules/local/custom_starpostprocess/merge_starfusion/main.nf'




workflow RNASEQ_WF {
    
    take:
        ch_publish_dir
        ch_enable_publish
        ch_disable_publish
        ch_reads
        ch_input_csv 
        ch_adapter_sequence
        ch_adapter_sequence_r2
        ch_salmon_index
        ch_star_index
        ch_gtf
        ch_tx2gene
        ch_reference_celltype
        ch_multiqc_config
        ch_genebody_ref
        ch_min_reads
        ch_skip_10X
        ch_lib_type_10x
        ch_lib_protocol_10x
        ch_skip_subsampling
        ch_n_reads
        ch_read_length
        ch_seqtk_sample_seed
        ch_two_color_chemistry
        ch_celltype_ref
        ch_instrument
        ch_genome
        ch_tmp_dir
        ch_timestamp
        ch_skip_qc_plots
        ch_bed_files_directory

    main:
    
    if( ! ch_skip_10X ){
    
      ALEVIN_NOQUANT_DUMPFQ ( ch_reads, 
                                ch_salmon_index,
                                ch_tx2gene,
                                ch_lib_type_10x,
                                ch_lib_protocol_10x,
                                ch_publish_dir,
                                ch_disable_publish
       )
       
        CUSTOM_AWK_DEMUX_CBC_FASTQ ( ALEVIN_NOQUANT_DUMPFQ.out.fastq_cbc,
                                ch_publish_dir,
                                ch_disable_publish
        
                              )
                              
        ch_fastqs = CUSTOM_AWK_DEMUX_CBC_FASTQ.out.fastqs
        
        
    }else{
    
        if ( ch_input_csv  ) {
            
            PUBLISH_INPUT_DATASET_WF (
                                        ch_input_csv,
                                        ch_publish_dir,
                                        ch_enable_publish
                                    )

        } else {
            ch_input_csv = Channel.of(file("NO_FILE"))
        }
        COUNT_READS_FASTQ_WF (
                                    ch_reads,
                                    ch_publish_dir,
                                    ch_enable_publish
                                )
        COUNT_READS_FASTQ_WF.out.read_counts
            .map { sample_id, files, read_count_file -> 
                def read_count = read_count_file.text.trim().toLong()
                [sample_id, files, read_count]
            }
            .branch { read ->
                small: read[2] < ch_min_reads
                large: read[2] >= ch_min_reads
            }
            .set { branched_reads }
        
        ch_fastqs = ch_reads
    }
    
    ch_seqtk_version = Channel.empty()
    ch_fastp_report = Channel.empty()
    ch_fastp_version = Channel.empty()
    
    if ( !ch_skip_subsampling ) {
        // Run sub sampling followed by fastp only when subsampling is enabled 
        // and when there are reads with min threshold
        if (branched_reads.large.ifEmpty { false }) {
            ch_fastqs_with_nreads = ch_fastqs.map { sample_name, reads ->
                tuple(sample_name, reads, ch_n_reads)
            }
            SEQTK_WF (ch_fastqs_with_nreads, false, ch_read_length, ch_seqtk_sample_seed, ch_publish_dir, ch_disable_publish)
            ch_seqtk_version = SEQTK_WF.out.version

            if (!params.skip_cutadapt) {
                CUTADAPT ( SEQTK_WF.out.reads, ch_adapter_sequence, ch_adapter_sequence_r2, ch_publish_dir, ch_disable_publish)
                FastpFull_WF( CUTADAPT.out.reads, ch_two_color_chemistry, ch_adapter_sequence, ch_adapter_sequence_r2, ch_publish_dir, ch_disable_publish)
            } else {
                FastpFull_WF( SEQTK_WF.out.reads, ch_two_color_chemistry, ch_adapter_sequence, ch_adapter_sequence_r2, ch_publish_dir, ch_disable_publish)
            }
            ch_fastp_report = FastpFull_WF.out.report
            ch_fastp_version = FastpFull_WF.out.version
        }
    } else {

        if (!params.skip_cutadapt) {
            CUTADAPT ( ch_fastqs, ch_adapter_sequence, ch_adapter_sequence_r2, ch_publish_dir, ch_disable_publish)
            FastpFull_WF(CUTADAPT.out.reads, ch_two_color_chemistry, ch_adapter_sequence, ch_adapter_sequence_r2, ch_publish_dir, ch_disable_publish)
        } else {
            FastpFull_WF(ch_fastqs, ch_two_color_chemistry, ch_adapter_sequence, ch_adapter_sequence_r2, ch_publish_dir, ch_disable_publish)
        }
        
        ch_fastp_report = FastpFull_WF.out.report
        ch_fastp_version = FastpFull_WF.out.version
        
     }

    // Process fastp output and branch based on filtered read count
    FastpFull_WF.out.reads
        .join(FastpFull_WF.out.json)
        .map { sample_id, reads, json -> 
            def fastp_data = new groovy.json.JsonSlurper().parseText(json.text)
            def read_count = fastp_data.summary.after_filtering.total_reads
            [sample_id, reads, read_count]
        }
        .branch { read ->
            small: read[2] < ch_min_reads
            large: read[2] >= ch_min_reads
        }
        .set { branched_reads_after_filter }
    ch_fastqs_filtered_poor_quality_samples = branched_reads_after_filter.large.map { sample_id, reads, _read_count ->
            tuple(sample_id, reads)
    }
              
    SALMONQUANT (ch_fastqs_filtered_poor_quality_samples, ch_salmon_index, ch_publish_dir, ch_disable_publish)
    ch_salmon_report = SALMONQUANT.out.outdir
    ch_salmon_version = SALMONQUANT.out.version
    
    CREATE_TXIMPORT_SALMON_TX_GENE ( SALMONQUANT.out.outdir, ch_tx2gene, ch_publish_dir, ch_disable_publish)
    df_tximport = CREATE_TXIMPORT_SALMON_TX_GENE.out.collect()
    MERGE_TXIMPORT_SALMON_TX_GENE_WF (df_tximport, ch_publish_dir, ch_enable_publish)
    ch_qunt_merge_tsv = MERGE_TXIMPORT_SALMON_TX_GENE_WF.out.quant_merge_tsv
    CREATE_MATRIX_SALMON_TX_GENE ( ch_qunt_merge_tsv, ch_publish_dir, ch_enable_publish )

    STARALIGN ( ch_fastqs_filtered_poor_quality_samples, ch_star_index, ch_publish_dir,ch_disable_publish)
    ch_star_report = STARALIGN.out.outdir
    ch_star_version = STARALIGN.out.version
    // ch_junction_file = STARALIGN.out.junction
    
    STARGETPRIM_WF ( STARALIGN.out.bam, ch_publish_dir, ch_enable_publish)
    ch_stargetprim_version = STARGETPRIM_WF.out.version
    // ch_raw_bam = STARGETPRIM_WF.out.primary_bam.collect()

    // STARGETPRIM_WF.out.junction
    //     .collectFile( name: "junction_files.txt", newLine: true, sort: { it[0] }, storeDir: "${ch_tmp_dir}" )
    //     { it[0].replaceFirst(/_.*/,"") + "\t" + "${ch_publish_dir}_${ch_timestamp}/secondary_analyses/alignment_htseq/" + it[1].getName() }

    STARGETPRIM_WF.out.primary_bam
        .collectFile( name: "bam_files.txt", newLine: true, sort: { item -> item[0] }, storeDir: "${ch_tmp_dir}" )
            { item -> item[0] + "\t" + "${ch_publish_dir}_${ch_timestamp}/secondary_analyses/alignment_htseq/" + item[1].getName() }

    GENE_BODY_COVERAGE_RNA (STARGETPRIM_WF.out.primary_bam, ch_genebody_ref, ch_disable_publish )
    ch_gene_body_df = GENE_BODY_COVERAGE_RNA.out.df.collect()
    GENE_BODY_COVERAGE_RNA_PLOT (ch_gene_body_df, ch_disable_publish)
    ch_gene_body_report = GENE_BODY_COVERAGE_RNA_PLOT.out.gene_body_png

    // Run INSERT_SIZE_ANALYSIS_WF only for GRCh38 genome
    ch_insert_size_summary = Channel.empty()
    if (ch_genome == "GRCh38") {
        INSERT_SIZE_ANALYSIS_WF(STARGETPRIM_WF.out.primary_bam,
                ch_bed_files_directory,
                ch_publish_dir,
                ch_enable_publish
        )
        ch_insert_size_summary = INSERT_SIZE_ANALYSIS_WF.out.aggregated_summary
    }
    

    QUALIMAP_BAMRNA ( STARGETPRIM_WF.out.bam_raw , ch_gtf, ch_publish_dir, ch_disable_publish)
    ch_qualimap_report = QUALIMAP_BAMRNA.out.outdir
    ch_qualimap_version = QUALIMAP_BAMRNA.out.version
                        
                    
    HTSEQ_COUNTS ( STARGETPRIM_WF.out.primary_bam, ch_gtf, ch_publish_dir, ch_disable_publish)
    ch_htseq_counts_version = HTSEQ_COUNTS.out.version
    
    CREATE_HTSEQ_SUMMARY_WF ( HTSEQ_COUNTS.out.htseq_counts, ch_tx2gene, ch_publish_dir, ch_disable_publish)
    df_htseqsum = CREATE_HTSEQ_SUMMARY_WF.out.htseq_summary.collect()
    MERGE_HTSEQ_SUMMARY_WF ( df_htseqsum, ch_publish_dir, ch_enable_publish)
    ch_merge_tsv = MERGE_HTSEQ_SUMMARY_WF.out.merge_tsv.collect()
    CREATE_HTSEQ_MATRIX ( ch_merge_tsv, ch_publish_dir, ch_enable_publish)
    PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY_WF (ch_merge_tsv, ch_publish_dir, ch_disable_publish)
    ch_pcaheatmap_plot=  PLOTTER_PCAHEATMAP_HTSEQ_SUMMARY_WF.out.pcaheatmap_plot
    ch_reference_celltype =  Channel.fromPath(ch_celltype_ref)
    CELL_TYPING (ch_merge_tsv, ch_reference_celltype.collect(), ch_publish_dir, ch_enable_publish)
     
    //STARFUSION_WF (ch_junction_file, params.ctat_resource_lib, ch_publish_dir, ch_disable_publish)
    //ch_starfusion_all = STARFUSION_WF.out.fusion_file.collect()
    //MERGE_STAR_FUSION ( ch_starfusion_all, ch_publish_dir, ch_enable_publish )
      

    /*
    ========================================================================================
        REPORTING
    ========================================================================================
    */

    
    CREATE_QC_REPORT (ch_qualimap_report.collect(), ch_publish_dir, ch_enable_publish)
    ch_qc_report_stats= CREATE_QC_REPORT.out.stats
    ch_qc_report_pstats = CREATE_QC_REPORT.out.pstats
    
    CALC_DYNAMICRANGE (CREATE_HTSEQ_MATRIX.out.htseq_matrix ,ch_publish_dir, ch_enable_publish)

    if ( ch_skip_subsampling || branched_reads.large.ifEmpty { true } ) {
    	collect_master_stats = ch_merge_tsv.collect().ifEmpty([])
	    .combine( ch_fastp_report.collect().ifEmpty([]))
	    .combine( ch_qunt_merge_tsv.collect().ifEmpty([]) )
	    .combine( ch_qc_report_stats.collect().ifEmpty([]) )
	    .combine( CELL_TYPING.out.annotations.collect().ifEmpty([]) )
	    .combine( ch_gene_body_df.collect().ifEmpty([]) )
	    .combine( CALC_DYNAMICRANGE.out.collect().ifEmpty([]) )
	    
    }else{

       collect_master_stats = ch_merge_tsv.collect().ifEmpty([])
            .combine( ch_fastp_report.collect().ifEmpty([]))
            .combine( SEQTK_WF.out.metadata.collect().ifEmpty([]))
            .combine( ch_qunt_merge_tsv.collect().ifEmpty([]) )
            .combine( ch_qc_report_stats.collect().ifEmpty([]) )
            .combine( CELL_TYPING.out.annotations.collect().ifEmpty([]) )
            .combine( ch_gene_body_df.collect().ifEmpty([]) )
            .combine( CALC_DYNAMICRANGE.out.collect().ifEmpty([]) )


    }

    CREATE_MASTER_STATS (collect_master_stats, COUNT_READS_FASTQ_WF.out.combined_read_counts, ch_insert_size_summary.ifEmpty([]), branched_reads_after_filter.large.ifEmpty{}, ch_publish_dir, ch_disable_publish)

    ch_tool_versions = ch_fastp_version.take(1)
                                    .combine(ch_seqtk_version.take(1).ifEmpty([]))
                                    .combine(ch_salmon_version.take(1).ifEmpty([]))
                                    .combine(ch_star_version.take(1).ifEmpty([]))
                                    .combine(ch_qualimap_version.take(1).ifEmpty([]))
                                    .combine(ch_stargetprim_version.take(1).ifEmpty([]))
                                    .combine(ch_htseq_counts_version.take(1).ifEmpty([]))


    REPORT_VERSIONS_WF(
                            ch_tool_versions,
                            ch_publish_dir,
                            ch_enable_publish)

    qc_plots_composition_jpg = Channel.empty()
    if ( !ch_skip_qc_plots ) {
        RNA_QC_PLOTS ( 
            CREATE_HTSEQ_MATRIX.out.htseq_matrix,
            CREATE_MASTER_STATS.out.stats,
            ch_input_csv,
            ch_publish_dir,
            ch_enable_publish
        )
        qc_plots_composition_jpg = RNA_QC_PLOTS.out.composition_jpg
    }

    collect_mqc = CREATE_MASTER_STATS.out.stats.collect().ifEmpty([])
                    .combine( ch_fastp_report.collect().ifEmpty([]))
                    .combine( ch_star_report.collect().ifEmpty([]))
                    .combine( ch_qualimap_report.collect().ifEmpty([]))
                    .combine( ch_salmon_report.collect().ifEmpty([]))
                    .combine( ch_qc_report_stats.collect().ifEmpty([]))
                    .combine( ch_qc_report_pstats.collect().ifEmpty([]))
                    .combine( ch_qunt_merge_tsv.collect().ifEmpty([]))
                    .combine( ch_merge_tsv.collect().ifEmpty([]))
                    .combine( ch_pcaheatmap_plot.collect().ifEmpty([]))
                    .combine( ch_gene_body_report.collect().ifEmpty([]))
                    .combine( REPORT_VERSIONS_WF.out.versions.collect().ifEmpty([]))
                    .combine( qc_plots_composition_jpg.collect().ifEmpty([]))
                    .combine( CREATE_HTSEQ_MATRIX.out.housekeeping_genes_CV.collect().ifEmpty([]))
                    .combine( CREATE_HTSEQ_MATRIX.out.housekeeping_genes_counts.collect().ifEmpty([]))
                    .combine( CREATE_HTSEQ_MATRIX.out.housekeeping_genes_clustergram.collect().ifEmpty([]))

    params_meta = [
            session_id: workflow.sessionId,
            instrument: ch_instrument,
            genome: ch_genome,
            read_length: ch_read_length,
            adapter_sequence: ch_adapter_sequence,
            adapter_sequence_r2: ch_adapter_sequence_r2,
            skip_subsampling: ch_skip_subsampling,
            skip_10X: ch_skip_10X
    ]

    if (!ch_skip_subsampling) {
        params_meta['subsample_reads'] = ch_n_reads
    }

    if (!ch_skip_10X) {
        params_meta['lib_protocol_10x'] = ch_lib_protocol_10x
        params_meta['lib_type_10x'] = ch_lib_type_10x
    }

    MULTIQC_WF ( collect_mqc,
                    params_meta,
                    ch_multiqc_config.collect(),
                    ch_publish_dir,
                    ch_enable_publish
                  )
        
                            
    //ch_deliverable_files = CUSTOM_FASTQ_MERGE_WF.out.merged_fastqs_to_compress.collect().ifEmpty([])
    //                                .combine(ch_raw_bam.ifEmpty([]))
    //                                .combine(df_htseqsum.ifEmpty([]))
    //                                .combine(df_tximport.ifEmpty([]))
    //                                .combine(ch_merge_tsv.ifEmpty([]))
    //                               .combine(ch_qunt_merge_tsv.ifEmpty([]))

    //CUSTOM_PIPELINE_JSON_OUTPUT (   ch_deliverable_files,
    //                              ch_publish_dir,
    //                              ch_enable_publish
    //                   )
                                       
                                       
    
   
  
}

workflow.onComplete {
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["bam"]         = [:]
    // output["output"]["junction"]    = [:]


    bam_outfile = file("$params.tmp_dir/bam_files.txt")
    bam_outfile_lines = bam_outfile.exists() ? bam_outfile.readLines() : []
    for ( bam_line : bam_outfile_lines ) {
        def (sample_name, bam_path) = bam_line.split('\t')
        output["output"]["bam"][sample_name] = [:]
        output["output"]["bam"][sample_name]["bam"] = bam_path
    }
    
    // junction_outfile = file("$params.tmp_dir/junction_files.txt")
    // junction_outfile_lines = junction_outfile.readLines()
    // for ( junction_line : junction_outfile_lines ) {
    //     def (sample_name, junction_path) = junction_line.split('\t')
    //     output["output"]["junction"][sample_name] = [:]
    //     output["output"]["junction"][sample_name]["junction"] = junction_path
    // }
    
    def output_json = groovy.json.JsonOutput.toJson(output)
    def output_json_pretty = groovy.json.JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)
}

workflow.onError {
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["bam"]         = [:]
    // output["output"]["junction"]    = [:]
    
    def output_json = groovy.json.JsonOutput.toJson(output)
    def output_json_pretty = groovy.json.JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)
    

    def subject = """\
        [nf-rnaseq-pipeline] FAILED: ${workflow.runName}
        """

    def msg = """\
        Pipeline execution summary 
        --------------------------------
        Script name       : ${workflow.scriptName ?: '-'}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Workflow repo     : ${workflow.repository ?: '-' }
        Workflow revision : ${workflow.repository ? "$workflow.revision ($workflow.commitId)" : '-'}
        Workflow profile  : ${workflow.profile ?: '-'}
        Workflow cmdline  : ${workflow.commandLine ?: '-'}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
        Error Report      : ${workflow.errorReport}
        """
        .stripIndent()
    
    log.info ( msg )
    
    if ( "${params.email_on_fail}" && workflow.exitStatus != 0 ) {
        sendMail(to: "${params.email_on_fail}", subject: subject, body: msg)
    }
}
