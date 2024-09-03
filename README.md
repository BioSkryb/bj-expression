# BJ-Expression

Pipeline processes RNASeq data and assess transcript-level and gene-level quantification


**Input Options**

The input for the pipeline is fastq files. The input can be passed either directly as path to the input file directory or via a input.csv with a meta data.

- **Reads Directory Input**: Use the `--reads` parameter to specify the path to a directory containing the input files. By default, this parameter is set to `null`. For example, to use fastq files from a specific directory, you would use: 
`--reads 's3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/rnaseq/*R{1,2}_001.fastq.gz'`

- **CSV Metadata Input**: Alternatively, you can use the `--input_csv` parameter to specify a CSV file containing metadata. This parameter is also `null` by default. The CSV file should have 6 columns: `biosampleName`, `sampleId`, `reads`, `readLength`, `read1` and `read2`. 
The `biosampleName` column contains the name of the biosample, the `sampleId` column contains the sample name in Illumina specified name format, `reads` have the number of reads, `readLength` with length of the reads and `read1` and `read2` has the path to the input reads. For example:

```
biosampleName,sampleId,reads,readLength,read1,read2
Expression-test1,Expression-test1_S1_L001,1000000,76,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/rnaseq/Expression-test1_S1_L001_R1_001.fastq.gz,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/rnaseq/Expression-test1_S1_L001_R2_001.fastq.gz
```


**Optional Modules**


This pipeline includes several optional modules. You can choose to include or exclude these modules by adjusting the following parameters:

- `--skip_subsampling`: Set this to `true` to exclude the subsampling module. By default, it is set to `true`.
- `--skip_fastq_merge`: Set this to `true` to exclude the fastq merge module. By default, it is set to `false`.

**Outputs**


The pipeline saves its output files in the designated "publish_dir" directory. The bam files after htseq alignment are stored in the "secondary_analyses/alignment_htseq/" subdirectory and the metrics files are saved in the "secondary_analyses/secondary_metrics/" subdirectory.

**nf-test**


The BioSkryb BJ-Expression nextflow pipeline run is tested using the nf-test framework.

Installation:

nf-test has the same requirements as Nextflow and can be used on POSIX compatible systems like Linux or OS X. You can install nf-test using the following command:
```
wget -qO- https://code.askimed.com/install/nf-test | bash
```
It will create the nf-test executable file in the current directory. Optionally, move the nf-test file to a directory accessible by your $PATH variable.

Usage:

```
nf-test test
```

The nf-test for this repository is saved at tests/ folder.

```
    test("scrna-seq test") {

        when {
            params {
                // define parameters here. Example: 
                publish_dir = "${outputDir}/results"
                timestamp = "test"
            }
        }

        then {
            assertAll(
                // Check if the workflow was successful
                { assert workflow.success },

                // Verify existence of the multiqc report HTML file
                {assert new File("${outputDir}/results_test/multiqc/multiqc_report.html").exists()},

                // Check for a match in the pipeline_metrics_summary csv file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/secondary_metrics/pipeline_metrics_summary.csv")).match("pipeline_metrics_summary")},

                // Check for a match in the pipeline_metrics_summary_percents csv file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/secondary_metrics/pipeline_metrics_summary_percents.csv")).match("pipeline_metrics_summary_percents")},

                // Check for a match in the df_dynamicrange_expression tsv file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/secondary_metrics/df_dynamicrange_expression.tsv")).match("df_dynamicrange_expression")},

                // Verify existence of the df_gene_counts_salmon file
                {assert new File("${outputDir}/results_test/secondary_analyses/quantification_salmon/df_gene_counts_salmon.tsv").exists()},

                // Verify existence of the df_gene_counts_starhtseq file
                {assert new File("${outputDir}/results_test/secondary_analyses/quantification_htseq/df_gene_counts_starhtseq.tsv").exists()},

                // Verify existence of the bam file
                {assert new File("${outputDir}/results_test/secondary_analyses/alignment_htseq/Expression-test1.bam.bai").exists()}

            )
        }

    }
```
