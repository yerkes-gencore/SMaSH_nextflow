process SUBSAMPLE_BAM {
    label 'samtools'

    input:
    tuple val(sample_id), path(full_bam)
    
    output:
    tuple val(sample_id), path('sub.bam')
    
    script:
    """
    filesize=\$(samtools view $full_bam --no-header -c)
    if [[ \$filesize -lt 5000000 ]]
    then
        samtools view $full_bam --with-header -b > sub.bam
    else
        frac=\$(echo "scale=3;$params.n_reads_smash / \$filesize" | bc)
        samtools view $full_bam --with-header --subsample \$frac -b > sub.bam
    fi
    """
}

process SORT_BAM {
    label 'samtools'
    publishDir "$params.outdir/subset_bams", mode: 'copy', overwrite: false
    
    input:
    tuple val(sample_id), path(subsampled_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_smash.bam")
    
    script:
    """
    samtools sort $subsampled_bam -@ $params.cpus -o ${sample_id}_smash.bam
    """
}

process INDEX_BAM {
    label 'samtools'
    publishDir "$params.outdir/subset_bams", mode: 'copy', overwrite: false
    
    input:
    tuple val(sample_id), path("${sample_id}_smash.bam")
    
    output:
    path("${sample_id}_smash.*", includeInputs: true)
    
    script:
    """
    samtools index ${sample_id}_smash.bam -@ $params.cpus
    """
}

process SMASH {
    publishDir "$params.outdir/smash_out", mode: 'copy', overwrite: false
    
    input:
    path(indexed_bams)
    // path vcf
    
    output:
    tuple path('*pval_out.txt*'), path('*.vcf.*.p')
    
    script:
    """
    /yerkes-cifs/runs/tools/SMaSH-master/SMaSH.py \
        -i $params.vcf_path/$params.vcf_file \
        -bam 'ALL' \
        -output_dir .
    """
}

process PLOT_HEATMAP {
    publishDir "$params.outdir/smash_out", mode: 'move', overwrite: false

    input:
    path 'pval_out.txt'

    output: 
    path 'SMaSH_pval_heatmap.png'

    script:
    """
    Rscript plot_SMaSH_heatmap.R ${params.heatmapWidth} ${params.heatmapTextSize}
    """
}
