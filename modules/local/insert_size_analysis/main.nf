nextflow.enable.dsl=2

process INSERT_SIZE_ANALYSIS {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/insert_size/", enabled:"$enable_publish", pattern: "${sample_name}_insert_size_summary.tsv"
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    path(bed_files_dir)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}_insert_size_summary.tsv"), emit: summary_table
    path("samtools_version.yml"), emit: version

    script:
    """
    # Create summary table header
    echo -e "Sample\\tBiotype\\tCount\\tMean\\tMedian" > ${sample_name}_insert_size_summary.tsv
    
    # Process each BED file
    for bed_file in ${bed_files_dir}/*_transcripts.bed; do
        if [ -s "\$bed_file" ]; then
            # Extract biotype name from filename
            biotype=\$(basename "\$bed_file" _transcripts.bed)
            
            # Step 1: Subset BAM file using BED file
            samtools view -b -h -L "\$bed_file" ${bam} -o ${sample_name}_\${biotype}.bam
            
            # Step 2: Index the subset BAM
            samtools index ${sample_name}_\${biotype}.bam
            
            # Step 3: Calculate insert sizes for properly paired reads
            samtools view -f 3 -F 4 -F 8 ${sample_name}_\${biotype}.bam | \\
                awk '{if(\$9>0) print \$9}' > ${sample_name}_\${biotype}_insert_sizes.txt
            
            # Step 4: Calculate statistics
            if [ -s "${sample_name}_\${biotype}_insert_sizes.txt" ]; then
                stats=\$(sort -n ${sample_name}_\${biotype}_insert_sizes.txt | \\
                    awk '{
                        arr[NR]=\$1; 
                        sum+=\$1; 
                    } END {
                        n=NR;
                        if(n>0) {
                            mean=sum/n;
                            if(n%2==1) {
                                median=arr[(n+1)/2]
                            } else {
                                median=(arr[n/2]+arr[n/2+1])/2
                            }
                            print n "\\t" mean "\\t" median;
                        } else {
                            print "0\\t0\\t0";
                        }
                    }')
                
                # Add to summary table
                echo -e "${sample_name}\\t\${biotype}\\t\${stats}" >> ${sample_name}_insert_size_summary.tsv
            else
                # No reads found for this biotype
                echo -e "${sample_name}\\t\${biotype}\\t0\\t0\\t0" >> ${sample_name}_insert_size_summary.tsv
            fi
            
            # Clean up intermediate BAM files to save space
            rm -f ${sample_name}_\${biotype}.bam ${sample_name}_\${biotype}.bam.bai
        else
            # Empty BED file
            biotype=\$(basename "\$bed_file" _transcripts.bed)
            echo -e "${sample_name}\\t\${biotype}\\t0\\t0\\t0" >> ${sample_name}_insert_size_summary.tsv
        fi
    done
    
    # Get samtools version
    export SAMTOOLS_VER=\$(samtools --version 2>&1 | sed -n -e '1p' | grep -Eo [0-9][.]*[0-9]*)
    echo "Samtools: \$SAMTOOLS_VER" > samtools_version.yml
    """
}

process COMBINE_INSERT_SIZE_SUMMARIES {
    tag "Combining insert size summaries"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/insert_size/", enabled:"$enable_publish", pattern: "aggregated_insert_size_summary.tsv"
    
    input:
    path(summary_files)
    val(publish_dir)
    val(enable_publish)

    output:
    path("combined_insert_size_summary.tsv"), emit: combined_summary
    path("aggregated_insert_size_summary.tsv"), emit: aggregated_summary

    script:
    """
    # Combine all individual summary files
    echo -e "Sample\\tBiotype\\tCount\\tMean\\tMedian" > combined_insert_size_summary.tsv
    
    # Skip headers and combine all files
    for file in ${summary_files}; do
        if [ -s "\$file" ]; then
            tail -n +2 "\$file" >> combined_insert_size_summary.tsv
        fi
    done
    
    # Create aggregated RNA types summary
    python3 - << 'EOF'
    import pandas as pd
    import numpy as np

    # Read the combined data
    df = pd.read_csv('combined_insert_size_summary.tsv', sep='\\t')

    # Define RNA type aggregations
    rna_groups = {
        'Protein_Coding': ['protein_coding'],
        'Ribosomal_RNA': ['rRNA', 'Mt_rRNA', 'rRNA_pseudogene'],
        'Mitochondrial_RNA': ['Mt_rRNA', 'Mt_tRNA'], 
        'Small_RNAs': ['snRNA', 'snoRNA', 'scRNA', 'scRNAs', 'scaRNA', 'ribozyme', 'vault_RNA', 'sRNA'],
        'MicroRNAs': ['miRNA']
    }

    # Get unique samples
    samples = df['Sample'].unique()

    # Initialize results
    results = []

    for sample in samples:
        sample_data = df[df['Sample'] == sample]
        result_row = {'Sample': sample}
        
        for group_name, biotypes in rna_groups.items():
            # Filter data for this RNA group
            group_data = sample_data[sample_data['Biotype'].isin(biotypes)]
            
            if len(group_data) > 0 and group_data['Count'].sum() > 0:
                # Calculate weighted means using counts as weights
                total_count = group_data['Count'].sum()
                
                if total_count > 0:
                    # Weighted mean calculation
                    weighted_mean = np.average(group_data['Mean'], weights=group_data['Count'])
                    
                    # For median, we need to be more careful - use weighted approach
                    # Simple approach: weight by counts for existing medians
                    non_zero_data = group_data[group_data['Count'] > 0]
                    if len(non_zero_data) > 0:
                        weighted_median = np.average(non_zero_data['Median'], weights=non_zero_data['Count'])
                    else:
                        weighted_median = 0
                else:
                    weighted_mean = 0
                    weighted_median = 0
            else:
                weighted_mean = 0
                weighted_median = 0
            
            result_row[f'{group_name}_Mean'] = round(weighted_mean, 2)
            result_row[f'{group_name}_Median'] = round(weighted_median, 2)
        
        results.append(result_row)

    # Create DataFrame and save
    result_df = pd.DataFrame(results)

    # Reorder columns
    column_order = ['Sample', 
                    'Protein_Coding_Mean', 'Protein_Coding_Median',
                    'Ribosomal_RNA_Mean', 'Ribosomal_RNA_Median',
                    'Mitochondrial_RNA_Mean', 'Mitochondrial_RNA_Median', 
                    'Small_RNAs_Mean', 'Small_RNAs_Median',
                    'MicroRNAs_Mean', 'MicroRNAs_Median']

    result_df = result_df[column_order]

    # Save to file
    result_df.to_csv('aggregated_insert_size_summary.tsv', sep='\\t', index=False)

    print("Aggregated summary created successfully!")
    print(f"Processed {len(samples)} samples")
    print("RNA type aggregations:")
    for group, types in rna_groups.items():
        print(f"  {group}: {', '.join(types)}")
    EOF
    """
}

workflow INSERT_SIZE_ANALYSIS_WF {
    take:
        ch_bam_bundle  // tuple val(sample_name), path(bam), path(bai)
        ch_bed_files_dir
        ch_publish_dir
        ch_enable_publish

    main:
        INSERT_SIZE_ANALYSIS(
            ch_bam_bundle,
            ch_bed_files_dir,
            ch_publish_dir,
            ch_enable_publish
        )

        COMBINE_INSERT_SIZE_SUMMARIES(
            INSERT_SIZE_ANALYSIS.out.summary_table.map { _sample_name, summary_file -> summary_file }.collect(),
            ch_publish_dir,
            ch_enable_publish
        )

    emit:
        version = INSERT_SIZE_ANALYSIS.out.version
        summary_table = INSERT_SIZE_ANALYSIS.out.summary_table
        combined_summary = COMBINE_INSERT_SIZE_SUMMARIES.out.combined_summary
        aggregated_summary = COMBINE_INSERT_SIZE_SUMMARIES.out.aggregated_summary
}