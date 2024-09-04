# BJ-Expression

Pipeline processes RNASeq data and assess transcript-level and gene-level quantification

# Running locally

Following are instructions for running BJ-Expression in a local Ubuntu 18.04 server

## Install Java 11

```
sudo apt-get install default-jdk

java -version
#openjdk version "11.0.18" 2023-01-17
#OpenJDK Runtime Environment (build 11.0.18+10-post-Ubuntu-0ubuntu118.04.1)
#OpenJDK 64-Bit Server VM (build 11.0.18+10-post-Ubuntu-0ubuntu118.04.1, mixed mode, sharing)
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

```

## Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

```

## Download and Unzip Pipeline Repository
```
cd bj-expression
```

## Sentieon License Setup

The Sentieon license is based on a lightweight floating license server process running on one node, and serving licenses though TCP to all other nodes. Normally this license server is running in a special non-computing node on the cluster periphery that has unrestricted access to the outside world through HTTPS, and serves the licenses to the rest of the nodes in the cluster by listening to a specific TCP port that needs to be open within the cluster. The license server needs to have external https access to validate the license, while the computing nodes do not need to have access to the internet.

Client will need to provide FQDN (hostname) or IP address of the machine that they would like to use to host the license server along with the port the license server will listen at to create a license file by Sentieon.

Sentieon license server supports connection through the proxy server. Set the standard `http_proxy` environment before starting the license server.

In order to run Sentieon software will need to start a Sentieon license server on eg. `b06x-pbs01.inet.xxxxxxxxx`; running the following command will setup the license server as a running daemon in your system:

```
export http_proxy=<proxy_server_name_and_port>
<SENTIEON_DIR>/bin/sentieon licsrvr --start --log <LOCATION_OF_LOG_FILE> <LICENSE_FILE> 
```
To run Sentieon software in the computing nodes, will need to set an environmental variable to tell Sentieon software the location of the license server and port. This can be added to bash profile or to the scripts that will drive the pipelines:
```
export SENTIEON_LICENSE=b06x-pbs01.inet.xxxxxxxxxxxx:xxxx
```

## Test Pipeline Execution

**Input Options**

The input for the pipeline is passed via a input.csv with a meta data.

- **CSV Metadata Input**: The CSV file should have 4 columns: `biosampleName`, `reads`, `read1` and `read2`. 
The `biosampleName` column contains the name of the biosample, `reads` have the number of reads and `read1` and `read2` has the path to the input reads. For example:

```
biosampleName,reads,read1,read2
Expression-test1,1000000,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/rnaseq/Expression-test1_S1_L001_R1_001.fastq.gz,s3://bioskryb-data-share/BioSkryb-Testing-Data/genomics/homo_sapiens/GRCh38/illumina/fastq/small/rnaseq/Expression-test1_S1_L001_R2_001.fastq.gz
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
