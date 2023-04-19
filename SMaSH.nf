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
        path("*smash.bam"), emit: bam
        path("*smash.bam.bai"), emit: bai
    script:
    // """
    // bam_basename=`basename ${bam}`
    // echo "process bam, output \${bam_basename}_sorted.bam.bai"
    // """

   """
   bam_basename=`basename ${bam}`
   filesize=\$(/yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools view $bam --no-header -c)
   frac=\$(echo "scale=3;${params.n_reads_smash} / \$filesize" | bc)
   /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools view $bam --with-header --subsample \$frac -b > \${bam_basename}_smash.bam
   /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools sort $bam -o \${bam_basename}_smash.bam
   /yerkes-cifs/runs/tools/samtools/samtools-1.17/samtools index $bam -o \${bam_basename}_smash.bam.bai
   """
}

process smash {
    input:
        path bam
        path bai

    output:
        stdout
    
    script:
    """"
    echo ${bam}
    echo ${bai}

    # zcat /yerkes-cifs/runs/tools/docker/images/smash.tar.gz | docker load

    # my_pwd=`pwd`

    # docker run \
    #     --rm \
    #     -v /yerkes-cifs/runs/tools/SMaSH-master:/yerkes-cifs/runs/tools/SMaSH-master \
    #     -v ${my_pwd}:/smash/ \
    #     -w /smash/ \
    #     ff13285f3b4f \
    #     /yerkes-cifs/runs/tools/SMaSH-master/SMaSH.py \
    #             -i /yerkes-cifs/runs/tools/SMaSH-master/snps_GRCh38.vcf \
    #             -bam 'ALL'
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
        
    process_bam_smash(files)
    smash(process_bam_smash.out.bam.collect(), process_bam_smash.out.bai.collect()) | view()
    // smash(process_bam_smash.out.collect()) | view()
    //smash(files.collect()) | view()
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

