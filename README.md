# SMaSH_nextflow

Start by cloning this repo into the desired working dir. 

Then modify `SMaSH_config.nf` to specify at least three key parameters:
 - `bam_dir` - the directory with the bam files output by the STAR run
 - `outdir` - the directory to publish the results in (smash output will actually be published in `${outdir}/smash_out`
 - `vcf_file` - name of vcf file used by SMaSH
 - `vcf_path` - path to vcf file used by SMaSH

If you are working on a server other than sblab-03 you may want to adjust `cpus` and `maxForks` in the config too, as they are optimized for that server.

For running on a new dataset, test out on the first few (default: 4) bams using using:
`/yerkes-cifs/runs/tools/nextflow/nextflow main.nf -c SMaSH_config.nf --dev`

This test run should finish in a few minutes. If that works, run on the whole dataset by dropping the `--dev` tag:
`/yerkes-cifs/runs/tools/nextflow/nextflow main.nf -c SMaSH_config.nf`

As of 2023.06.14 the PLOT_HEATMAP process is commented out, so if you want a heatmap you will need to run `bin/plot_SMaSH_heatmap.R` on a machine with R and tidyverse installed. See issue https://github.com/yerkes-gencore/SMaSH_nextflow/issues/2#issue-1757315834
