import csv
import os
import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

configfile: "config.yaml"

shell.prefix("source {0}/env.cfg; export PYTHONPATH=$PYTHONPATH:{0}/wssd_sunk:{0}/ssf_DTS_caller; ".format(SNAKEMAKE_DIR))

COORDS = list(config.get("bedfiles").values())
REGION_NAMES = list(config.get("bedfiles").keys())

TABLE_DIR = config["table_dir"]
PLOT_DIR = config["plot_dir"]

REFERENCE = config["reference"]
DATASETS = config["datasets"]
DATATYPES = config["data_types"]
MASTER_MANIFEST = config["master_manifest"]
POP_CODES = config["pop_codes"]

ds_manifest = pd.read_table(config["datasets_file"], header=0)
ds_manifest = ds_manifest.ix[ds_manifest.reference == REFERENCE, :]
ds_manifest.index = ds_manifest.dataset

DIRS_TO_MAKE = ["log", TABLE_DIR, PLOT_DIR]

genotyper = config.get("genotyper")

if genotyper == "wssd_cc":
	include: "workflows/wssd_cc_genotyper.snake"
elif genotyper == "gglob":
	include: "workflows/gglob_genotyper.snake"

for folder in DIRS_TO_MAKE:
	if not os.path.exists(folder):
		os.makedirs(folder)

GENE_GRAM_SETTINGS = config.get("gene_gram_settings", "")
SPP = config.get("spp", 500)

def get_region_names(region_names):
	names = []
	for region in region_names:
		with open(config.get("bedfiles")[region], 'r') as reader:
			for line in reader:
				names.append("%s.%s" % (region, line.rstrip().split()[3]))
	return names

def get_coords_and_size_from_name(name, coord_files):
	 for coord_file in coord_files:
		with open(coord_file, 'r') as reader:
			for line in reader:
				dat = line.rstrip().split()
				if dat[3] == name:
					chr, start, end = dat[0], int(dat[1]), int(dat[2])
					return ("%s:%d-%d" % (chr, start, end), str(end - start))

def get_family_from_name(name, coord_files):
	for coord_file in coord_files:
		with open(coord_file, 'r') as reader:
			for line in reader:
				test_name = line.rstrip().split()[3]
				if test_name == name:
					return os.path.basename(coord_file).split(".")[0]

localrules: all, get_combined_pdfs, get_tables

rule all:   
	input:  
		"{TABLE_DIR}/num_sunks.table.tab",
		expand("%s/gene_grams/{fam}_{dataset}_{datatype}.0.{file_type}" % (PLOT_DIR), fam = REGION_NAMES, dataset = DATASETS, datatype = DATATYPES, file_type = config["plot_file_type"]),
		expand("%s/violin/{fam_name}.{dataset}_violin_{datatype}.{file_type}" % (PLOT_DIR), fam_name = get_region_names(REGION_NAMES), dataset = config["main_dataset"], datatype = DATATYPES, file_type = config["plot_file_type"]),
		expand("%s/{fam}.{plottype}_{datatype}.pdf" % PLOT_DIR, fam = REGION_NAMES, plottype=["violin", "scatter", "superpop"], datatype = DATATYPES),
		expand("{fam}/{fam}.{ds}.combined.{dt}.GMM.bed", fam = REGION_NAMES, ds = DATASETS, dt = DATATYPES)

rule get_cn_wssd_variance:
	input:  
		gts = expand("{fam}/{dataset}/{dataset}.wssd.genotypes.tab", fam = REGION_NAMES, dataset = DATASETS),
		sunks = "%s/num_sunks.table.tab" % (TABLE_DIR)
	output: 
		tab = "%s/wssd_stats_by_family.tab" % (TABLE_DIR)
	params: sge_opts = "-l mfree=2G -N get_var -l h_rt=1:00:00", families = REGION_NAMES
	run:
		wssd_stats = pd.DataFrame(columns=["name"])
		datasets = DATASETS
		size = []
		for i, dataset in enumerate(datasets):
			wssd_mean, wssd_std, cn_two = [], [], []
			sample_dat = pd.read_csv(ds_manifest.loc[dataset]["manifest"], sep='\t', header=0)
			samples = sample_dat["sample"].tolist()
			for j, fam in enumerate(REGION_NAMES):
				data = pd.read_csv("%s/%s/%s_wssd_genotypes.tab" % (fam, dataset, dataset), sep='\t', header=0)
				common_samples = [sample for sample in samples if sample in data.columns]
				cns = data[common_samples]
				wssd_mean.append(cns.mean(axis=1).mean(axis=0))
				wssd_std.append(cns.std(axis=1).mean(axis=0))
				cn_two.append(cns.applymap(lambda x: x < 2.5).sum(axis=1).sum(axis=0))
				if i == 0:
					size.append(int((data["end"] - data["start"]).mean()))
			wssd_stats[dataset + "_mean"] = wssd_mean
			wssd_stats[dataset + "_std"] = wssd_std
			wssd_stats[dataset + "_cn2"] = cn_two
		wssd_stats["name"] = list(params.families)
		wssd_stats["size"] = size
		wssd_stats.to_csv(output.tab, index=False, sep='\t')

rule get_cn_sunk_variance:
	input:  gts = expand("{fam}/{dataset}/{dataset}.sunk.genotypes.tab", fam = REGION_NAMES, dataset = DATASETS), 
			sunks = "%s/num_sunks.table.tab" % (TABLE_DIR)
	output: "%s/sunk_stats_by_region.tab" % TABLE_DIR
	params: sge_opts = "-l mfree=2G -N get_var -l h_rt=1:00:00", families = REGION_NAMES
	run:
		sunk_stats = pd.DataFrame(columns = ["name"])
		names = []
		for i, dataset in enumerate(DATASETS):
			sample_dat = pd.read_csv(ds_manifest.loc[dataset]["manifest"], header=0, sep='\t')
			samples = sample_dat["sample"].tolist()
			sunk_mean, sunk_std, cn_zero = [], [], []
			for j, fam in enumerate(params.families):
				data = pd.read_csv("%s/%s/%s_sunk_genotypes.tab" % (fam, dataset, dataset), sep='\t', header=0)
				common_samples = [sample for sample in samples if sample in data.columns]
				cns = data[common_samples]
				sunk_mean.extend(cns.mean(axis=1).tolist())
				sunk_std.extend(cns.std(axis=1).tolist())
				cn_zero.extend(cns.applymap(lambda x: x < 0.5).sum(axis=1).tolist())
				if i == 0:
					names.extend(data["name"].tolist())
			sunk_stats[dataset + "_mean"] = sunk_mean
			sunk_stats[dataset + "_std"] = sunk_std
			sunk_stats[dataset + "_cn<0.5"] = cn_zero
		sunk_stats["name"] = names
		num_sunks = pd.read_csv(input.sunks, sep='\t', header=0)
		sunk_stats = num_sunks.merge(sunk_stats, on = "name")
		sunk_stats.to_csv(output[0], index=False, sep='\t')

rule get_sunks:
	input: COORDS
	output: "%s/num_sunks.table.tab" % (TABLE_DIR)
	params: sge_opts = "-l mfree=2G -N get_SUNKs -l h_rt=0:30:00", sunks = config["ref_files"][REFERENCE]["sunk_bed"]
	run:
		for i, coords in enumerate(COORDS):
			if i == 0:
				shell("""module load bedtools/2.21.0; bedtools intersect -a {coords} -b {params.sunks} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
						 awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' > {output[0]}""")
			else:
				shell("""module load bedtools/2.21.0; bedtools intersect -a {coords} -b {params.sunks} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
						 awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' >> {output[0]}""")
		shell("""sed -i '1ichr\tstart\tend\tname\tsize\tnSUNKs' {output[0]}""")

rule plot_gene_grams:
	input: 
		bed = "{fam}/{fam}.{dataset}.combined.{datatype}.bed",
	output: 
		"%s/gene_grams/{fam}_{dataset}_{datatype}.0.{file_type}" % (PLOT_DIR)
	params: sge_opts = "-l mfree=8G -N gene_grams -l h_rt=0:30:00"
	run:
		manifest = ds_manifest.loc[wildcards.dataset]["manifest"]
		shell("""python scripts/gene_gram.py {input.bed}{manifest} {PLOT_DIR}/gene_grams/{wildcards.fam}_{wildcards.dataset}_{wildcards.datatype} --plot_type {wildcards.file_type} --spp {SPP} {GENE_GRAM_SETTINGS}""")

rule get_combined_pdfs:
	input: expand("%s/{fam}/violin_{datatype}.pdf" % PLOT_DIR, fam = REGION_NAMES, datatype = DATATYPES)

rule combine_violin_pdfs:
	input: expand("%s/{plottype}/{fam_name}.{dataset}_{plottype}_{datatype}.pdf" % (PLOT_DIR), plottype = ["violin", "scatter", "superpop"], fam_name = get_region_names(REGION_NAMES), dataset = config["main_dataset"], datatype = DATATYPES)
	output: "%s/{fam}.violin_{datatype}.pdf" % PLOT_DIR, "%s/{fam}.scatter_{datatype}.pdf" % PLOT_DIR, "%s/{fam}.superpop_{datatype}.pdf" % PLOT_DIR
	params: sge_opts = "-l mfree=8G -N pdf_combine -l h_rt=0:30:00"
	run:
		for pt in ["violin", "scatter", "superpop"]:
			outfile = [file for file in output if pt in file][0]
			shell("""gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={outfile} plots/{pt}/{wildcards.fam}*{pt}_{wildcards.datatype}.pdf""")

rule plot_violins:
	input: 
		df = expand("%s/{fam}.{{dataset}}.{{datatype}}.genotypes.df" % (TABLE_DIR), fam = REGION_NAMES)
	output: 
		violin = "%s/violin/{fam_name}.{dataset}_violin_{datatype}.{file_type}" % (PLOT_DIR), 
		scatter = "%s/scatter/{fam_name}.{dataset}_scatter_{datatype}.{file_type}" % (PLOT_DIR), 
		superpop = "%s/superpop/{fam_name}.{dataset}_superpop_{datatype}.{file_type}" % (PLOT_DIR)
	params: 
		sge_opts = "-l mfree=8G -N plot_violins -l h_rt=0:10:00",
		max_cp = "7"
	run:
		fam, name = wildcards.fam_name.split(".")[0], ".".join(wildcards.fam_name.split(".")[1:])
		input_table = [file for file in input.df if fam in file and wildcards.dataset in file][0]
		dat = pd.read_csv(input_table, sep='\t')
		(coords, size) = get_coords_and_size_from_name(name, COORDS)
		title = "_".join([name, coords, size, config["reference"], wildcards.dataset, wildcards.datatype])
		shell("""Rscript scripts/genotype_violin.R {input_table} {output.violin} {name} {wildcards.file_type} {title} 3 violin super_pop_only --max_cp {params.max_cp}""")
		shell("""Rscript scripts/genotype_violin.R {input_table} {output.scatter} {name} {wildcards.file_type} {title} 3 --max_cp {params.max_cp}""")
		shell("""Rscript scripts/genotype_violin.R {input_table} {output.superpop} {name} {wildcards.file_type} {title} 3 super_pop_only --max_cp {params.max_cp};""")

rule get_tables:
	input: expand("%s/{fam}.{dataset}.{datatype}.genotypes.df" % (TABLE_DIR), fam = REGION_NAMES, dataset = DATASETS, datatype = DATATYPES)

rule get_long_table:
	input: 
		regions = "{fam}/{fam}.{dataset}.combined.{datatype}.bed"
	output: 
		tab = "{TABLE_DIR}/{fam}.{dataset}.{datatype}.genotypes.df"
	params: 
		sge_opts = "-l mfree=8G -N make_long_table -l h_rt=0:30:00"
	resources:
		mem=8
	threads: 1
	shell:
		'''
		Rscript scripts/transform_genotypes.R {input.regions} {MASTER_MANIFEST} {POP_CODES} {wildcards.dataset} {output.tab}")
		'''

rule get_combined_GMM_genotypes:
	input: 
		bed = "{fam}/{fam}.{dataset}.combined.{datatype}.bed"
	output: 
		bed = "{fam}/{fam}.{dataset}.combined.{datatype}.GMM.bed"
	params: 
		sge_opts = "-N GMM -l h_rt=0:30:00", 
		max_cp = "12"
	resources:
		mem = 4
	threads: 1
	shell:
		'''
		python scripts/get_GMM_genotypes.py {input.bed} {output.bed} --max_cp {params.max_cp}
		'''

rule combine_genotypes:
	input: 
		tab = expand("{{fam}}/{{fam}}.{ds}.{{datatype}}.genotypes.tab", ds = DATASETS)
	output: 
		bed = "{fam}/{fam}.{dataset}.combined.{datatype,\w+}.bed"
	params: 
		sge_opts="-l mfree=2G -N combine_gt -l h_rt=0:30:00"
	resources:
		mem = 2
	threads: 1
	run:
		fn_main = [x for x in input.tab if wildcards.dataset in x][0]
		main_ds = pd.read_csv(fn_main, na_values='NA', sep='\t', header=0)
		if config["append_dataset"] is not None:
			for app_ds in config["append_dataset"]:
				if app_ds != wildcards.dataset:
					append_dataset = pd.read_csv("{fam}/{fam}.{app_ds}.{dt}.genotypes.tab".format(fam=wildcards.fam, app_ds=app_ds, dt=wildcards.datatype), header=0, sep='\t', index_col=False)
					main_ds = main_ds.merge(append_dataset, on=["chr", "start", "end", "name"])
		main_ds.to_csv(output.bed, index=False, sep='\t', na_rep='NA')
