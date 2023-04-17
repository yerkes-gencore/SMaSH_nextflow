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

process process_bam_smash {
    input:
        path bam
    output:
        path "${sample}_smash.bam", emit: bam
    script:
    """
    filesize=\$(/yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools view $bam --no-header -c)
    frac=\$(echo "scale=3;${params.n_reads_smash} / \$filesize" | bc)
    /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools view $bam --with-header --subsample \$frac -b > ${sample}_smash.bam
    /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools sort $bam -o ${sample}_sorted.bam
    /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools index -o ${sample}_sorted.bam.bai
    """
}

process smash {
    input:
        path bam
        //path bai
    output:
        stdout
    script:
    """
    echo "execute smash with samples $bam"
    """
}


// ////////////////////////////////////////////////////
// /* --              WORKFLOW                    -- */
// ////////////////////////////////////////////////////

workflow {
    println "Workflow start: $workflow.start"
    Channel
        .fromPath("${params.bam_dir}/*/*.bam")
        .set { files }
        
    //process_bam_smash(files)
    //smash(process_bam_smash.out.collect()) | view()
    smash(files.collect()) | view()
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

