#!/usr/bin/env nextflow

// ////////////////////////////////////////////////////
// /* --            PROCESSES                     -- */
// ////////////////////////////////////////////////////

include { SUBSAMPLE_BAM } from './modules.nf'
include { SORT_BAM      } from './modules.nf'
include { INDEX_BAM     } from './modules.nf'
include { SMASH         } from './modules.nf'

// ////////////////////////////////////////////////////
// /* --              WORKFLOW                    -- */
// ////////////////////////////////////////////////////

// Instantiate input bams channel
Channel
    .fromPath("${params.bam_dir}/**.bam", checkIfExists: true)
    .take ( params.dev ? params.dev_n_samples : -1 )
    .map { file -> [file.simpleName, file] }
    .set { input_bam_ch } 

// Main workflow
workflow {
    println "Workflow start: $workflow.start"

    // Subsample bams only if params.subsample == true, otherwise pass input bams to next process
    if ( params.subsample ) {
        SUBSAMPLE_BAM(input_bam_ch)
            .set { subsampled_bam_ch }
    } 
    else {
        input_bam_ch
            .set { subsampled_bam_ch }
    }

    SORTING_WORKFLOW(subsampled_bam_ch)
        .set { indexed_bam_ch }
    
    SMASH(indexed_bam_ch)
        .filter { it.name == 'pval_out.txt' } // And passes 'pval_out.txt' to the heatmap process
        .set { smash_p_val_ch }
}

// Sub-workflow to handle sorting and indexing, if necessary
workflow SORTING_WORKFLOW {
    take: 
    subsampled_bam_ch

    main: 
    // If bams are not already indexed (default), sort and index bams
    if ( !params.input_bais_exist ) {
        SORT_BAM(subsampled_bam_ch) 
            .set { sorted_bam_ch }

        INDEX_BAM(sorted_bam_ch) // Publishes subsampled + sorted bams used in SMaSH
            .collect() // And passes bams and bais to SMaSH
            .set { indexed_bam_ch }     
    }
    // If input bams are indexed (i.e. bams and bais exist in ), bypass sorting and indexing, find the bais and create a channel with both
    // But if you are also subsampling, you need to re-sort and re-index anyway
    if ( params.input_bais_exist & !params.subsample ) {
        Channel
            .fromPath("${params.bam_dir}/**.bai", checkIfExists: true)
            .take ( params.dev ? params.dev_n_samples : -1 )
            .map { file -> [file.simpleName, file] }
            .set { input_bai_ch }
        
        input_bam_ch
            .mix(input_bai_ch)
            .groupTuple()
            .map { baseName, file -> file }
            .collect()
            .set { indexed_bam_ch }
    }

    emit:
    indexed_bam_ch
}

workflow.onComplete {
    def msg = """\
        Workflow complete! The key output is ${params.outdir}/smash_out/pval_out.txt

        To visualize these data more easily, open this directory with R and modify plot_SMaSH_heatmap.R to plot a heatmap using pval_out.txt given the sampling design for these data.

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

// ////////////////////////////////////////////////////
// /* --              HELP MSG                    -- */
// ////////////////////////////////////////////////////

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