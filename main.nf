#!/usr/bin/env nextflow
nextflow.enable.dsl=2
def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run SMaSH.nf -c nf_SMaSH.config

    Edit the nextflow.config file to add run parameters
    """.stripIndent()
}

// Show help message
params.cleanup = true  
params.help = ""
if (params.help) {
    helpMessage()
    exit 0
}

// ////////////////////////////////////////////////////
// /* --          VALIDATE INPUTS                 -- */
// ////////////////////////////////////////////////////

def validations() {
    try {
        file("${params.outdir}/", checkIfExists:true)
    } catch (e) {
        println "\nOutdir '$params.outdir' does not exist, creating path"
        file("${params.outdir}/").mkdirs()
    }
}
validations()

// ////////////////////////////////////////////////////
// /* --            PROCESSES                     -- */
// ////////////////////////////////////////////////////

process SUBSAMPLE_BAM {
    
    input:
    tuple val sample_id, path full_bam
    
    output:
    tuple val sample_id, path 'sub.bam'
    
    script:
    """
    filesize=\$(/yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools view $full_bam --no-header -c)
    frac=\$(echo "scale=3;$params.n_reads_smash / \$filesize" | bc)
    /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools view full.bam --with-header --subsample \$frac -b > sub.bam
    """
}

process SORT_BAM {
    publishDir "$params.outdir/subset_bams", mode: 'copy', overwrite: false
    cpus "$params.cpus"
    maxForks "$params.maxForks"
    
    input:
    tuple val sample_id, path subsampled_bam
    
    output:
    tuple val sample_id, path 'sorted.bam'
    
    script:
    """
    /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools sort $subsampled_bam -@ 32 -o sorted.bam
    """
}

process INDEX_BAM {
    publishDir "$params.outdir/subset_bams", mode: 'copy', overwrite: false
    cpus "$params.cpus"
    maxForks "$params.maxForks"
    
    input:
    tuple val sample_id, path sorted_bam
    
    output:
    tuple path "${sample_id}_smash.bam", path "${sample_id}_smash.bai"
    
    script:
    """
    /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools index $sorted_bam -@ 32 -o ${sample_id}_smash.bai
    """
}

process SMASH_BAM {
    container = 'ff13285f3b4f'
    publishDir "$params.outdir/smash_out", mode: 'copy', overwrite: false
    
    input:
    path sorted_bams
    path sorted_bais
    path vcf
    
    output:
    tuple path '*pval_out.txt*', path '*.vcf.*.p'
    
    script:
    """
    /yerkes-cifs/runs/tools/SMaSH-master/SMaSH.py \
        -i $vcf \
        -bam 'ALL' \
        -output_dir smash_out/
    """

}

// process PLOT_HEATMAP {
//     // container = base R container with dplyr and ggplot2
//     publishDir "$params.outdir/smash_out", mode: 'move', overwrite: false

//     input:
//     path 'pval_out.txt'

//     output: 
//     path 'pval_heatmap.png'

//     script:
//     // R script

// }

// ////////////////////////////////////////////////////
// /* --              WORKFLOW                    -- */
// ////////////////////////////////////////////////////

workflow {
    println "Workflow start: $workflow.start"
    
    // Instantiate input bams channel
    Channel
        .fromPath("${params.bam_dir}/**.bam")
        .take ( params.dev ? params.dev_n_samples : -1 )
        .map { file -> tuple(file.baseName, file) }
        set { input_bams_ch }
    
    // Prepare bams for SMaSH
    SUBSAMPLE_BAM(input_bams_ch) | SORT_BAM() | INDEX_BAM() // Publishes subsampled + sorted bams used in SMaSH
        .collect() // And passes bams and bais to SMaSH
        .set { processed_bams_ch }
    
    // Instantiate vcf channel for 
    Channel
        .fromPath("$params.vcf_path/$params.vcf_file")
        .set { vcf_ch }

    // Run SMaSH (and send 'p_val_out.txt' to heatmap process)
    SMASH(processed_bams_ch, vcf_ch) // Publishes SMaSH output
        .filter( it[1].name == 'pval_out.txt' ) // And passes 'pval_out.txt' to the heatmap process
        .set { smash_p_val_ch }

    // Plot p-value heatmap
    // PLOT_HEATMAP(smash_p_val_ch) // Publishes a heatmap representation of SMaSH output
}

workflow.onComplete {
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        publishDir  : ${params.outdir}
        """
        .stripIndent()
    if (params.emails?.trim()){
        sendMail(to: "${params.emails}", subject: 'SMaSH Run Complete', body: msg)
    } 
}
