# read_depth_genotyper
Get genotypes for regions using WSSD and SUNK

## Install
`git clone --recursive https://github.com/EichlerLab/read_depth_genotyper.git`

`alias snakesub='mkdir -p log; snakemake -p --ri --jobname "{rulename}.{jobid}" --drmaa " -P eichlerlab -q eichler-short.q -l h_rt=48:00:00  -l mfree={resources.mem}G -pe serial {threads} -R y -V -cwd -S /bin/bash " --drmaa-log-dir $PWD/log --latency-wait 60'`

## Genotype
`snakesub -j [max_jobs] get_tables`

## Plot
`snakesub -j [max_jobs]`
