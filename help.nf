nextflow.enable.dsl=2

def printHeader() {
  
  log.info """\
  NF-RNASEQ   P I P E L I N E
  ===================================
  fastq files    : ${ params.reads }
  fastq files    : ${ params.input_csv }
  time stamp     : ${ params.timestamp }
  \n
  """

}

def helpMessage() {

  yellow = "\033[0;37m"
  blue = "\033[0;35m"
  white = "\033[0m"
  red = "\033[0;31m"

  log.info """\
    nf-rnaseq-pipeline
    Usage:
        nextflow run main.nf [options]
    Script Options: see nextflow.config
        ${yellow}
        [required]
        --reads             FILE    Path to fastq files specified as a glob pattern
        OR
        --reads_csv         FILE    Path to input csv file
        
        --organization      STR     BioSkryb organization/client id
                                    DEFAULT: ${params.organization}
    
        --workspace         STR     Organization/client's workspace
                                    DEFAULT: ${params.workspace}
    
        --project           STR     Organization/client's project
                                    DEFAULT: ${params.project}
                                    
        --salmon_index      DIR     Path to a folder containg a reference salmon indexes 
                                    DEFAULT: ${params.salmon_index}
                                    
        --star_index        DIR     Path to a folder containg a reference star indexes
                                    DEFAULT: ${params.star_index}   
                                    
        --gtf_file          FILE    GTF file containg gene locations
                                    DEFAULT: ${params.gtf_file}  
        
        --tx2gene           FILE    kdasda
                                    
                                    
        [optional]
        
        --publish_dir       DIR     Path to run output directory
                                    DEFAULT: ${params.publish_dir}
       
        --n_reads           VAL     Number of reads to sample for analysis eg. 2.5M == 5M paired reads
                                    DEFAULT: ${params.n_reads}
        --read_length       VAL     Desired read length for analysis and excess to be trimmed
                                    DEFAULT: ${params.read_length}
                                    
        --email_on_fail     STR     Email to receive upon failure
                                    DEFAULT: ${params.email_on_fail}
                                    
        --skip_subsampling  STR     Skip Qualimap module
                                    DEFAULT: ${params.skip_qualimap}
                                    
        --instrument        STR     Specify instrument. If 'NextSeq' or 'NovaSeq' set two_color_chemistry param true
                                    DEFAULT: ${params.instrument}
                                    
        --help              BOOL    Display help message
        
    ${yellow}
    """.stripIndent()


}

workflow{
  printHeader()
  helpMessage()
}