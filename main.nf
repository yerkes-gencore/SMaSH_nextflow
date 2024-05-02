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

// ////////////////////////////////////////////////////
// /* --              WORKFLOW                    -- */
// ////////////////////////////////////////////////////

// Three alternative workflows: (1) subsample_bams_wf, (2) sort_bams_wf, (3) smash_bams_wf

// Instantiate input bams channel
Channel
    .fromPath("${params.bam_dir}/**.bam")
    .take ( params.dev ? params.dev_n_samples : -1 )
    .map { file -> [file.simpleName, file] }
    .set { input_bams_ch } 

workflow START_BY_SMASHING {
    Channel
        .fromPath("${params.bam_dir}/**.bai")
        .take ( params.dev ? params.dev_n_samples : -1 )
        .map { file -> [file.simpleName, file] }
        .set { input_bais_ch }

    input_bams_ch
        .mix(input_bais_ch)
        .groupTuple()
        .map { baseName, file -> file }
        .collect()
        .set { processed_bams_ch }

    SMASH(processed_bams_ch)
        .filter { it.name == 'pval_out.txt' } // And passes 'pval_out.txt' to the heatmap process
        .set { smash_p_val_ch }
}

workflow START_BY_SORTING {
    input_bams_ch
        .set { subsampled_bam_ch }

    SORT_BAM(subsampled_bam_ch) 
        .set { sorted_bam_ch }

    INDEX_BAM(sorted_bam_ch) // Publishes subsampled + sorted bams used in SMaSH
        .collect() // And passes bams and bais to SMaSH
        .set { processed_bams_ch }

    SMASH(processed_bams_ch)
        .filter { it.name == 'pval_out.txt' } // And passes 'pval_out.txt' to the heatmap process
        .set { smash_p_val_ch }          
}

workflow START_BY_SUBSAMPLING {
    SUBSAMPLE_BAM(input_bams_ch)
        .set { subsampled_bam_ch }

    SORT_BAM(subsampled_bam_ch) 
        .set { sorted_bam_ch }

    INDEX_BAM(sorted_bam_ch) // Publishes subsampled + sorted bams used in SMaSH
        .collect() // And passes bams and bais to SMaSH
        .set { processed_bams_ch }

    SMASH(processed_bams_ch)
        .filter { it.name == 'pval_out.txt' } // And passes 'pval_out.txt' to the heatmap process
        .set { smash_p_val_ch }          
}

workflow {
    println "Workflow start: $workflow.start"

    // If start_by_smashing, create a bai channel, mix it with the bams and skip SUBSAMPLE_BAM, SORT_BAM, and INDEX_BAM
    if ( params.start_by == "smashing" ) {
        START_BY_SMASHING()
    }
    // If start-by-sorting, skip SORT_BAM, and INDEX_BAM
    else if ( params.start_by == "sorting" ) {
        START_BY_SORTING()
    }
    // If start-by-sorting, skip SORT_BAM, and INDEX_BAM
    else if ( params.start_by == "subsampling" ) {
        START_BY_SUBSAMPLING()
    }

    // Plot p-value heatmap
    //PLOT_HEATMAP(smash_p_val_ch) // Publishes a heatmap representation of SMaSH output
}

workflow.onComplete {
    def msg = """\
        Workflow complete! The key output is ${params.outdir}/smash_out/pval_out.txt

        To visualize these data more easily, open this directory in RStudio Server and modify plot_SMaSH_heatmap.R to plot a heatmap using pval_out.txt given the sampling design for these data.

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
    if (params.emails?.trim() && !params.dev){
        sendMail(to: "${params.emails}", subject: 'SMaSH Run Complete', body: msg)
    } 
}
