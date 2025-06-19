include { printHeader; helpMessage } from './help'
include { RNASEQ_WF } from './workflow/rnaseq.nf'

workflow {

    if ( params.help ) {
        helpMessage()
        workflow.onComplete { 
            System.exit(0) 
        }
    }
    
    printHeader()
    if ( params.reads ) {
                
        ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
        ch_reads.ifEmpty{ error "ERROR: cannot find any fastq files matching the pattern: ${params.reads}\nMake sure that the input file exists!" }

    } else if ( params.input_csv  ) {
                
        ch_reads = Channel.fromPath( params.input_csv  ).splitCsv( header:true )
                        .map { row -> [ row.biosampleName, [ row.read1, row.read2 ] ] }
        ch_reads.ifEmpty{ error "ERROR: Input csv file is empty." }
    }


                        
     ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)
     ch_reference_celltype = Channel.fromPath(params.celltype_ref)    

 
       
    RNASEQ_WF( 
                params.publish_dir,
                params.enable_publish,
                params.disable_publish,
                ch_reads, 
                params.input_csv,
                params.adapter_sequence,
                params.adapter_sequence_r2,
                params.salmon_index,
                params.star_index, 
                params.gtf_file, 
                params.tx2gene,
                ch_reference_celltype,
                ch_multiqc_config,
                params.genebody_ref,
                params.min_reads,
                params.skip_10X,
                params.lib_type_10x,
                params.lib_protocol_10x,
                params.skip_subsampling,
                params.n_reads,
                params.read_length,
                params.seqtk_sample_seed,
                params.two_color_chemistry,
                params.celltype_ref,
                params.instrument,
                params.genome,
                params.tmp_dir,
                params.timestamp


             )
}
