# read_depth_genotyper
Get genotypes for regions using WSSD and SUNK

## Install
`git clone --recursive https://github.com/EichlerLab/read_depth_genotyper.git`

`alias snakesub_sge='mkdir -p log; snakemake --ri -k --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -j y -o ./log -e ./log {params.sge_opts} -w n -S /bin/bash" -w 60'`

## Genotype
`snakesub_sge -j [max_jobs] -T get_tables`

## Plot
`snakesub_sge -T -j [max_jobs]`
