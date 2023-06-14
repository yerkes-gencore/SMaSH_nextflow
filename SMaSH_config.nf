/* 
 */

params {      
    // Where FASTQs are stored
    bam_dir = '/yerkes-cifs/runs/Analysis/2023_Analyses/p23025_Prabhu/p23025_Prabhu_Processing/second_seq_alignments'
    
    // Where you want the results saved
    outdir = '/yerkes-cifs/runs/analyst/micah/automation/SMaSH_nextflow'  //'s3://cellranger-nextflow/test' 

    // Separate multiple email addresses by commas if desired
    emails = null
    // emails = "mpfletc@emory.edu"
    
    // Resources (cpus = 32, maxForks = 25 works ok for sblab03)
    cpus = 32
    maxForks = 25
   
    // Number of reads per sample to use for SMaSH
    n_reads_smash = 5000000

    // VCF file
    vcf_file = 'rhesus.indian.af.0.3-0.7.onTarget.noIntronStreamStop.synonymous.vcf'

    // Path to VCF dir
    vcf_path = '/yerkes-cifs/runs/Genome_references/macaca_mulatta/mmul10/SNPs'
    
    // Run in development mode (run on only a small subset of samples)
    dev = false
    dev_n_samples = 4

    // Pval heatmap params
    heatmapWidth = 30 // cm
    heatmapTextSize = 4

}

// Nextflow automatically mounts the task workdir, so it only causes problems if you try to do that here.
docker {
    enabled = true
    runOptions = "-v /yerkes-cifs/runs/tools/SMaSH-master:/yerkes-cifs/runs/tools/SMaSH-master \
                  -v $params.vcf_path:$params.vcf_path"
    // runOptions = "-v /yerkes-cifs/runs/tools/SMaSH-master:/yerkes-cifs/runs/tools/SMaSH-master \
    //               -v $PWD/subset_bams:$PWD/subset_bams/ \
    //               -v $params.vcf_path:$params.vcf_path \
    //               -w $PWD/subset_bams/"
}

profiles { 
    dryrun {

    }
}

