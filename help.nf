nextflow.enable.dsl=2

def printHeader() {
  
  log.info """\
  BJ-EXPRESSION   PIPELINE
  ===================================
  fastq files    : ${ params.input_csv }
  time stamp     : ${ params.timestamp }
  \n
  """

}

def helpMessage() {

  def yellow = "\033[0;33m"
  def blue = "\033[0;34m"
  def white = "\033[0m"
  def red = "\033[0;31m"

  log.info """\
  ${blue}bj-expression pipeline
  Usage:
      nextflow run main.nf [options]
  Script Options: see nextflow.config
  ${red}
    [required]
    --reads_csv         FILE    Path to input csv file

    --genome            STR     Reference genome to use. Available options - GRCh38, GRCm39
                                DEFAULT: ${params.genome}
                                
                                
  ${yellow}                              
    [optional]
    
    --publish_dir       DIR     Path to run output directory
                                DEFAULT: ${params.publish_dir}
   
    --n_reads           VAL     Number of reads to sample for analysis
                                DEFAULT: ${params.n_reads}

    --read_length       VAL     Desired read length for analysis and excess to be trimmed
                                DEFAULT: ${params.read_length}

    --min_reads         VAL       Minimum number of reads required for analysis.
                                DEFAULT: ${params.min_reads}
                                
    --skip_subsampling  STR     Skip Qualimap module
                                DEFAULT: ${params.skip_subsampling}

    --skip_cutadapt     BOOL    Skip cutadapt module
                                DEFAULT: ${params.skip_cutadapt}
                                
    --help              BOOL    Display help message
    
${white}
""".stripIndent()

}

workflow{
  printHeader()
  helpMessage()
}