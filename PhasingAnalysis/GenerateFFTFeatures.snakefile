# GenerateFFTFeatures.snakefile
# Robert Patton, rpatton@fredhutch.org (Ha Lab)
# v1.0, adapted from Anna-Lisa Doebley's Griffin pipeline

"""
# before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Python/3.7.4-foss-2019b-fh1
PATH="$PATH:/fh/fast/ha_g/user/rpatton/scripts/"

# command to run snakemake (remove -np at end when done validating):
snakemake -s GenerateFFTFeatures.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np

# command to run snakemake on restart (remove -np at end when done validating):
#need to add '-q restart-new' according to scicomp. This should be fixed at some point.
snakemake -s GenerateFFTFeatures.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p restart-new -q restart-new --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName} --requeue" -j 40 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"
configfile: "config/cluster_slurm.yaml"

rule all:
	input:
		expand("{results_dir}/PhasingFM.tsv",results_dir=config['results_dir']),

rule calc_phasing:
	input:
		bam = lambda wildcards: config["samples"][wildcards.samples]['bam'],
		bias = lambda wildcards: config["samples"][wildcards.samples]['GC_bias']
	output:
		fm_file = temp(config['results_dir']+"/{samples}_PhasingFM.tsv")
	params:
		sample_name = "{samples}",
		results_dir=config['results_dir'],
		sites_bed = config['sites_bed'],

		reference_genome = config['reference_genome'],
		size_range=config['size_range'],
		map_q=config['map_quality'],
		plot_list=config['plot_list'],

		cpus = config['calc_phasing']['ncpus']

	shell:
		"time GenerateFFTFeatures.py --sample_name {params.sample_name} \
		--input {input.bam} --bias {input.bias} --annotation {params.sites_bed} \
		--reference_genome {params.reference_genome} --results_dir {params.results_dir} \
		--size_range {params.size_range} --map_quality {params.map_q} \
		--plot_list {params.plot_list} --cpu {params.cpus} "

rule combine_fms:
	input:
		fm_files = expand("{results_dir}/{samples}_PhasingFM.tsv", results_dir = config['results_dir'], samples = config['samples'])
	output:
		final = expand("{results_dir}/PhasingFM.tsv", results_dir = config['results_dir'])
	params:
		results_dir=config['results_dir']
	shell:
		'CombinePhasingFM.py --inputs {input.fm_files} --results_dir {params.results_dir}'
