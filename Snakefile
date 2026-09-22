### requires modules snakemake, star, subread, trimmomatic, fastqc, samtools, py-rseqc
print(sys.version)
import pandas as pd
import subprocess
import os
from collections import defaultdict

configfile: 'project_config_file.yaml'
localrules: aggregate_fastqc, aggregate_bams, aggregate_rseqc

config = defaultdict(str, config)

sample_sheet = os.path.join(os.getcwd(), config['sample_sheet'])
samples = pd.read_table(sample_sheet, sep = ',').set_index("Sample_Name", drop = False)

QC_name = config['QC_name']
results_dir = os.path.join('results', QC_name) + '/'
### Read QC summary file (input to batch_eval checkpoint) and creates dictionary to hold the file paths held within
QC_file = os.path.join(os.getcwd(), config['QC_file'])

### This was all taken from the SC pipeline but I didn't complete the upgrade. For now, just let it be.
qcdat = pd.read_csv(QC_file)
### Add missing steps to qcdat with '' file column, including step 0 (no QC done yet)
steps = list(set(list(range(2))) - set(qcdat.step))
d = {'step':steps, 'file':['' for s in steps]}
qcdat = pd.concat([qcdat, pd.DataFrame(d)])
### Reorder rows
qcdat['step_num'] = ['step' + str(q) for q in qcdat.step]
### Add and use columns with post-step file name and whether the step should be skipped
#qcdat['post_name'] = [results_dir + 'PostQC' + str(q) + '.RDS' for q in qcdat.step]
#qcdat = qcdat.set_index('step')
#qcdat = qcdat.sort_index(ascending = True)
### Add step_name
#qcdat['step_name'] = ['Base data', 'Sample filter']
### Add original RDS file
#qcdat.iloc[0, qcdat.columns.get_loc('post_name')] = 'data/All_data.RDS'
### Add skip information
#qcdat['skip'] = [f == '' for f in qcdat.file]
#qcdat.iloc[0, qcdat.columns.get_loc('skip')] = False
#qcdat.to_csv(results_dir + 'QC_steps.csv', sep = ',')
QC_specs = defaultdict(str, zip(qcdat.step_num, qcdat.file))
print(qcdat)

### For dynamic resources requests
size_step = [8000,15000,64000, 240000]
def get_mem_mb(wildcards, attempt):
  return size_step[attempt]

### Required columns in sample_sheet are: Sample_Name, S_N, Lane, 
### List of fastq files
samples['R1'] = samples['Sample_Name'].map(str) + '_' + samples['S_N'].map(str) + '_L00' + samples['Lane'].map(str) + '_R1_001'
samples['R2'] = samples['Sample_Name'].map(str) + '_' + samples['S_N'].map(str) + '_L00' + samples['Lane'].map(str) + '_R2_001'
samples['samp'] = samples['Sample_Name'].map(str) + '_' + samples['S_N'].map(str) + '_L00' + samples['Lane'].map(str)
#samples['R1'] = samples['Sample_Name'].map(str) + '_' + samples['Sample_ID'].map(str) + '_L00' + samples['Lane'].map(str) + '_R1_001'
#samples['R2'] = samples['Sample_Name'].map(str) + '_' + samples['Sample_ID'].map(str) + '_L00' + samples['Lane'].map(str) + '_R2_001'
#samples['samp'] = samples['Sample_Name'].map(str) + '_' + samples['Sample_ID'].map(str) + '_L00' + samples['Lane'].map(str)
samples['R1_fullpath'] = list([os.path.join(config['runs_dir'], samples.FC_Name[i], 'demultiplexed', config['project'], samples['Sample_ID'][i], samples['R1'][i] + '.fastq.gz') for i in range(0,len(samples['R1']))])
samples['R2_fullpath'] = list([os.path.join(config['runs_dir'], samples.FC_Name[i], 'demultiplexed', config['project'], samples['Sample_ID'][i], samples['R2'][i] + '.fastq.gz') for i in range(0,len(samples['R2']))])
runs_fastqs = list(samples['R1_fullpath']) + list(samples['R2_fullpath'])
samples['R1_projpath'] = list([os.path.join('data', 'fastq', 'raw', samples['R1'][i] + '.fastq.gz') for i in range(0,len(samples['R1']))])
samples['R2_projpath'] = list([os.path.join('data', 'fastq', 'raw', samples['R2'][i] + '.fastq.gz') for i in range(0,len(samples['R2']))])
proj_fastqs = list(samples['R1_projpath']) + list(samples['R2_projpath'])
samples['R1_subdir'] = list([os.path.join(samples['Sample_ID'][i], samples['R1'][i] + '.fastq.gz') for i in range(0,len(samples['R1']))])
samples['R2_subdir'] = list([os.path.join(samples['Sample_ID'][i], samples['R2'][i] + '.fastq.gz') for i in range(0,len(samples['R2']))])
fastq_files = list(samples['R1']) + list(samples['R2'])

### List of test variable names from the configuration file
tests = [x.strip() for x in config['test'].split(',')]
print('Tests:')
print(tests)
### Needs testing but I think this accomadates only one stratification category but with multiple factors for that category, e.g. Timepoint:Time1,Time2,Time3
strats = [x.strip() for x in config['stratifications'].split(':')]
strat_list = ['All']
if(len(strats) > 1):
  strat_list = strat_list + [strats[0] + '-' + v.strip() for v in strats[1].split(',')]
print('Stratifications: ')
print(strat_list)

def get_test_txt(wildcards):
  return [os.path.join('test/', ID + '_test.txt') for ID in list(samples[samples['FC_Name'] == wildcards.FC].Sample_ID)]

refbed = config['refbed']
if config['refbed'] == '':
  gtfbase = os.path.splitext(os.path.basename(config['refgtf']))[0]
  refbed = 'resources/' + gtfbase + '.bed'

rule all:
  input:
    #"QC_summaries.txt",
    [results_dir + "{test}_GSEA.xls".format(test = test) for test in tests],
    [results_dir + "{test}_DESeq2.xls".format(test = test) for test in tests],
    [results_dir + "{test}_GSEA.pdf".format(test = test) for test in tests],
    [results_dir + "{test}_DE.pdf".format(test = test) for test in tests]

print(runs_fastqs)
### FCs is a list
rule aggregate_fc_and_copy_to_project:
  input:  ### ancient is required though not ideal, because 'cat' apparently modifies the timestamp of the inputs :(
    ancient([os.path.join(config['runs_dir'], FC, "raw_fastq_files.txt") for FC in list(flowcells.keys())])
  output:
    #'data/fastq/raw/copy_fastqs.sh'
    proj_fastqs
  params:
    sh_file = 'data/fastq/raw/copy_fastqs.sh'
  shell:  ### Remove it in case in exists, then create and run an sh file that copies raw fastq files from runs directory to project's data directory
    """
    mkdir -p data/fastq/raw/
    rm -f {params.sh_file}
    cat {input} > {params.sh_file}
    sed -i 's/^/cp /' {params.sh_file}
    sed -i 's/$/ data\/fastq\/raw\/ /' {params.sh_file}
    sh {params.sh_file}
    """

rule run_fastqc_once:
  input:
    "data/fastq/{sdir}/{fastq}.fastq.gz"
  output:
    "data/fastq/{sdir}/fastqc/{fastq}_fastqc.zip"
  shell:
    "fastqc --extract -o data/fastq/{wildcards.sdir}/fastqc/ {input}"

### wildcard sdir cannot go in the expand function with the non-wildcard fastq
rule aggregate_fastqc:
  input:
    ["data/fastq/{sdir}/fastqc/" + second for second in expand("{fastq}_fastqc.zip", fastq = fastq_files)]
  output:
    "data/fastq/{sdir}/multiqc_input_files.txt"
  shell:
    "ls {input} > {output}"

rule aggregate_fastq:
  input:
    ["data/fastq/{sdir}/" + second for second in expand("{fastq}.fastq.gz", fastq = fastq_files)]
  output:
    "data/fastq/{sdir}/file_list.txt"
  shell:
    "ls {input} > {output}"

rule aggregate_bams:
  input:
    expand('data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam', sample = samples['samp'])
  output:
    'data/bam/bam_files.txt'
  shell:
    'ls {input} > {output}'

rule aggregate_rseqc:
  input:
    expand('data/bam/{sample}/RSeQC/multiqc_input_files.txt', sample = samples['samp'])
  output:
    'data/bam/multiqc_input_files.txt'
  shell:
    "cat {input} > {output}"

rule multiqc:
  input:
    "data/{sdir}/multiqc_input_files.txt"
  output:
    "data/{sdir}/multiqc.html"
  #conda:
  #  "/data/vrc_his/douek_lab/snakemakes/BulkRNASeq/multiqc.yml"
  shell:
    "multiqc --file-list {input} -f --filename {output}"

rule list_multiqcs:
  input:
    "data/fastq/raw/multiqc.html",
    "data/fastq/clipped/paired/multiqc.html",
    "data/bam/multiqc.html"
  output:
    "QC_summaries.txt"
  shell:
    "ls {input} > {output}"

### Only necessary because of some trimmomatic unzipping error. May remove once that is solved or circumvented. 
rule gunzip_once:
  input: 
    R1 = "data/fastq/raw/{sample}_R1_001.fastq.gz",
    R2 = "data/fastq/raw/{sample}_R2_001.fastq.gz"
  output: 
    R1 = temp("data/fastq/raw/{sample}_R1_001.fastq"),
    R2 = temp("data/fastq/raw/{sample}_R2_001.fastq")
  shell:
    """
    gunzip -c {input.R1} > {output.R1}
    gunzip -c {input.R2} > {output.R2}
    """

rule run_trimmomatic_once:
  input:
    R1 = "data/fastq/raw/{sample}_R1_001.fastq",
    R2 = "data/fastq/raw/{sample}_R2_001.fastq"
  output:
    R1_paired = "data/fastq/clipped/paired/{sample}_R1_001.fastq.gz",
    R1_unpaired = "data/fastq/clipped/unpaired/{sample}_R1_001.fastq.gz",
    R2_paired = "data/fastq/clipped/paired/{sample}_R2_001.fastq.gz",
    R2_unpaired = "data/fastq/clipped/unpaired/{sample}_R2_001.fastq.gz"
  shell:
    #"data/bam/multiqc.html"
    #"java -Xmx7G -jar /data/vrc_his/douek_lab/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 1 {input.R1} {input.R2} {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} ILLUMINACLIP:/data/vrc_his/douek_lab/programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:4:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50"
    "java -Xmx7G -jar trimmomatic-0.39.jar PE -phred33 -threads 1 {input.R1} {input.R2} {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:4:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50"

rule create_STAR_index:
  input:
    reffa = config['reffa'],
    refgtf = config['refgtf']
  output:
    directory('resources/STAR/'),
    'resources/STAR/geneInfo.tab' 
  shell:
    """
    STAR --runMode genomeGenerate --runThreadN 1 --genomeDir resources/STAR/ --genomeFastaFiles {inputs.reffa} --sjdbGTFfile {inputs.refgtf} sjdbOverhang 150
    """

rule run_STAR:
  input:
    "data/fastq/clipped/paired/{sample}_R1_001.fastq.gz",
    "data/fastq/clipped/paired/{sample}_R2_001.fastq.gz"
  params:
    starDB = config['starDB'],
    prefix = "data/bam/{sample}/{sample}_"
  output:
    "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    "data/bam/{sample}/{sample}_Log.final.out"
  resources: mem_mb=30000
  shell:
    """
    mkdir -p "data/bam/{wildcards.sample}/"
    STAR --genomeDir {params.starDB} --readFilesIn {input} --outFileNamePrefix {params.prefix} --readFilesCommand gzip -dc --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 --outFilterMismatchNoverReadLmax 0.04 --limitSjdbInsertNsj 1200000
   """

rule run_index:
  input:
    "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
  output:
    "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
  shell:
    "samtools index {input}"

rule sort_by_name:
  input:
    bam = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    bai = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
  output:
    "data/bam/{sample}/{sample}_sortn.bam"
  shell:
    "samtools sort -n -O bam {input.bam} -o {output}"

rule run_fixmate:
  input:
    "data/bam/{sample}/{sample}_sortn.bam",
  output:
    "data/bam/{sample}/{sample}_fixmate.bam"
  shell:
    "samtools fixmate -m -O bam {input} {output}"
 
rule sort_by_coord:
  input:
    "data/bam/{sample}/{sample}_fixmate.bam",
  output:
    "data/bam/{sample}/{sample}_sortc.bam"
  shell:
    "samtools sort -O bam {input} -o {output}"

rule run_markdup:
  input:
    "data/bam/{sample}/{sample}_sortc.bam"
  output:
    "data/bam/{sample}/{sample}_dedup.bam"
  shell:
    "samtools markdup -r {input} {output}"

rule bam_to_fastq:
  input:
    "data/bam/{sample}/{sample}_dedup.bam"
  output:
    f1 = "data/fastq/deduplicated/{sample}_R1_001.fastq.gz",
    f2 = "data/fastq/deduplicated/{sample}_R2_001.fastq.gz"
  shell:
    "samtools fastq -1 {output.f1} -2 {output.f2} {input}"

###  These two RSeQC functions take so long that they will run separately
rule run_RSeQC_tin:
  input:
    bam = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    bai = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
  params:
    anno = refbed,
    inter1 = '{sample}_Aligned.sortedByCoord.out.summary.txt',
    inter2 = '{sample}_Aligned.sortedByCoord.out.tin.xls'
  output: 
    tin = 'data/bam/{sample}/RSeQC/{sample}.out.summary.txt',
    xls = 'data/bam/{sample}/RSeQC/{sample}.out.tin.xls'
  shell:
    """
    tin.py -i {input.bam} -r {params.anno} > {output.tin}
    mv {params.inter1} {output.tin}
    mv {params.inter2} {output.xls}
    """

rule run_RSeQC_gbc:
  input: 
    bam = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    bai = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
  params:
    anno = refbed,
    out_pref = 'data/bam/{sample}/RSeQC/{sample}'
  output:
    geneBody_cov = 'data/bam/{sample}/RSeQC/{sample}.geneBodyCoverage.txt'
  shell:
    "geneBody_coverage.py -i {input.bam} -o {params.out_pref} -r {params.anno}"

rule run_RSeQC_once:
  input:
    bam = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    bai = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
  params:
    anno = refbed,
    out_pref = 'data/bam/{sample}/RSeQC/{sample}'
  output:
    bam_stat = 'data/bam/{sample}/RSeQC/{sample}.bam_stat.txt',
    read_dist = 'data/bam/{sample}/RSeQC/{sample}.read_distribution.txt',
    infer_exp = 'data/bam/{sample}/RSeQC/{sample}.infer_experiment.txt',
    inner_dist = 'data/bam/{sample}/RSeQC/{sample}.inner_distance_freq.txt',
    junction_anno = 'data/bam/{sample}/RSeQC/{sample}_junction_annotation.txt',
    junction_sat = 'data/bam/{sample}/RSeQC/{sample}.junctionSaturation_plot.r',
    read_GC = 'data/bam/{sample}/RSeQC/{sample}.GC.xls',
    dup_rate = 'data/bam/{sample}/RSeQC/{sample}.pos.DupRate.xls',
    RNA_frag_size = 'data/bam/{sample}/RSeQC/{sample}.RNA_fragment_size.txt'
  shell:
    """
    mkdir -p data/bam/{wildcards.sample}/RSeQC/
    bam_stat.py -i {input.bam} > {output.bam_stat}
    read_distribution.py -i {input.bam} -r {params.anno} > {output.read_dist}
    infer_experiment.py -i {input.bam} -r {params.anno} > {output.infer_exp}
    inner_distance.py -i {input.bam} -o {params.out_pref} -r {params.anno}
    junction_annotation.py -i {input.bam} -o {params.out_pref} -r {params.anno} 2>{output.junction_anno}
    junction_saturation.py -i {input.bam} -o {params.out_pref} -r {params.anno}
    read_GC.py -i {input.bam} -o {params.out_pref}
    read_duplication.py -i {input.bam} -o {params.out_pref}
    #read_NVC.py -i {input} -o {params.out_pref}
    #stread_quality.py -i {input} -o {params.out_pref}
    RNA_fragment_size.py -i {input} -r {params.anno} > {output.RNA_frag_size}
    #RPKM_saturation.py -i {input} -o {params.out_pref} -r {params.anno}
    #mismatch_profile.py -i {input} -o {params.out_pref}
    """

rule RSeQC_aggregate_sample:
  input:
    bam_stat = 'data/bam/{sample}/RSeQC/{sample}.bam_stat.txt',
    read_dist = 'data/bam/{sample}/RSeQC/{sample}.read_distribution.txt',
    infer_exp = 'data/bam/{sample}/RSeQC/{sample}.infer_experiment.txt',
    inner_dist = 'data/bam/{sample}/RSeQC/{sample}.inner_distance_freq.txt',
    junction_anno = 'data/bam/{sample}/RSeQC/{sample}_junction_annotation.txt',
    junction_sat = 'data/bam/{sample}/RSeQC/{sample}.junctionSaturation_plot.r',
    read_GC = 'data/bam/{sample}/RSeQC/{sample}.GC.xls',
    dup_rate = 'data/bam/{sample}/RSeQC/{sample}.pos.DupRate.xls',
    RNA_frag_size = 'data/bam/{sample}/RSeQC/{sample}.RNA_fragment_size.txt',
    geneBody_cov = 'data/bam/{sample}/RSeQC/{sample}.geneBodyCoverage.txt',
    tin = 'data/bam/{sample}/RSeQC/{sample}.out.summary.txt',
    starlog = 'data/bam/{sample}/{sample}_Log.final.out'
  output:
    'data/bam/{sample}/RSeQC/multiqc_input_files.txt'
  shell:
    'ls {input} > {output}'

rule make_gtf_R_object:
  input:
    config['refgtf']
  params:
    scripts = config['scripts_dir']
  output:
    'resources/gtf.RDS'
  conda:
    'envs/make_gtf_R_object.yaml'
  shell:
    'Rscript {params.scripts}/gtf.R {input} {output}'

rule make_bed_file:
  input:
    config['refgtf']
  output:
    refbed
  conda:
    'envs/make_bed_file.yaml'
  shell:
    'gtf2bed < {input} > {output}'

rule run_featureCounts:
  input:
    bams = expand('data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam', sample = samples['samp']),
    bais = expand('data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai', sample = samples['samp'])
  output:
    'data/counts/featureCounts.txt'
  params:
    ref = config['refgtf']
  shell:
    """
    mkdir -p data/counts/
    featureCounts -p --primary --largestOverlap -a {params.ref} -o {output} {input.bams}
    """

rule compile_QC_metrics:
  input:
    counts = 'data/counts/featureCounts.txt'
  params:
    scripts = config['scripts_dir'],
    covs = sample_sheet,
    bam_dir = 'data/bam/'
  output:
    'data/Covariates_QC_metrics.csv'
  shell:
    'Rscript {params.scripts}/Compile_QC_metrics.R {params.covs} {params.bam_dir} {input} {output}'

rule filter_samples_raw:
  input:
    'data/counts/featureCounts.txt',
    'data/Covariates_QC_metrics.csv',
    'data/gtf.RDS',
    QC_specs['step1'],
  params:
    scripts = config['scripts_dir'],
    test = config['test'],
    strat = config['stratifications'],
    adjust = config['adjust'],
    batch = config['batch_candidates'],
    results_dir = results_dir,
    qc_file = config['QC_file']
  output:
    results_dir + 'QC/Covariates_QC_metrics_intermediate.csv',
    results_dir + 'counts/filteredCounts.txt',
    results_dir + 'QC/sample_filters.pdf',
    results_dir + 'QC/sample_filters.txt',
    results_dir + 'QC/feature_filters.txt',
    results_dir + 'QC/chr_dist.pdf'
  shell:
    """
    mkdir -p {params.results_dir}/QC/
    mkdir -p {params.results_dir}/counts/
    cp {params.qc_file} {params.results_dir}
    Rscript {params.scripts}/Filters.R {input} {output} '{params.test}' '{params.strat}' '{params.adjust}' '{params.batch}'
    """

rule Normalization:
  input:
    results_dir + 'counts/filteredCounts.txt',
    results_dir + 'QC/Covariates_QC_metrics_intermediate.csv'
  output:
    results_dir + 'counts/normalizedCounts.txt',
    temp(results_dir + 'counts/normalizedCounts.RDS'),
    results_dir + 'QC/normalization.pdf'
  params:
    scripts = config['scripts_dir'],
    test = config['test'],
    adjust = config['adjust'],
    batch = config['batch_candidates']
  shell:
    """
    Rscript {params.scripts}/DESeq2_normalization.R {input} '{params.test}' '{params.adjust}' '{params.batch}' {output}
    """

rule filter_samples_norm:
  input:
    results_dir + 'counts/normalizedCounts.txt',
    results_dir + 'QC/Covariates_QC_metrics_intermediate.csv',
    QC_specs['step1']
  params:
    scripts = config['scripts_dir']
  output:
    results_dir + 'QC/Covariates_QC_metrics_filter.csv',
    results_dir + 'QC/PCA_filter.pdf',
    results_dir + 'QC/PCA_filter.csv'
  shell:
    "Rscript {params.scripts}/PCA_filter.R {input} {output}"

rule Evaluate_batch:
  input:
    results_dir + 'counts/normalizedCounts.txt',
    results_dir + 'QC/Covariates_QC_metrics_filter.csv',
  params:
    scripts = config['scripts_dir'],
    candidates = config['batch_candidates'],
    method = config['batch_method']
  output:
    results_dir + 'QC/batch_evaluation.txt',
    results_dir + 'QC/batch_evaluation.pdf'
  shell:
    """
    Rscript {params.scripts}/Batch_evaluation.R {input} '{params.candidates}' '{params.method}' {output}
    """

rule ComBat:
  input:
    results_dir + 'counts/normalizedCounts.txt',
    results_dir + 'counts/normalizedCounts.RDS',
    results_dir + 'QC/Covariates_QC_metrics_filter.csv',
    results_dir + 'QC/batch_evaluation.txt'
  params:
    scripts = config['scripts_dir'],
  output:
    results_dir + 'counts/finalCounts.txt',
    results_dir + 'counts/finalCounts.RDS'
  shell:
    'Rscript {params.scripts}/ComBat.R {input} {output}'

#rule AppData:
#  input:
#    'data/counts/finalCounts.txt',
#    'data/Covariates_QC_metrics_filter.csv'
#  params:
#    scripts = config['scripts_dir'],
#  output:
#    'results/App.RData'
#  shell:
#    'Rscript {params.scripts}/Compile_AppData.R'

rule DESeq:
  input:
    results_dir + 'counts/finalCounts.txt',
    results_dir + 'counts/finalCounts.RDS',
    results_dir + 'QC/Covariates_QC_metrics_filter.csv',
    results_dir + 'QC/batch_evaluation.txt',
    'data/gtf.RDS'
  params:
    scripts = config['scripts_dir'],
    adjust = config['adjust']
  output:
    results_dir + '{test}/{strat}/DESeq2_results.txt',
    results_dir + '{test}/{strat}/DESeq2_results.pdf'
  log: 
    results_dir + '{test}/{strat}/DESeq2_results.log'
  shell:
    """
    #mkdir -p 'results/{wildcards.test}/{wildcards.strat}/'
    Rscript {params.scripts}/DESeq2.R {input} '{wildcards.test}' '{params.adjust}' '{wildcards.strat}' {output} > {log}
    #Rscript {params.scripts}/DESeq2.R {input} '{wildcards.test}' '{params.adjust}' '{wildcards.strat}' {output}
    """

strats_str = ';'.join(strat_list),
rule DE_plots:
  input:
    results_dir + 'counts/finalCounts.txt',
    results_dir + 'QC/Covariates_QC_metrics_filter.csv',
    'data/gtf.RDS',
    [results_dir + "{test}/{strat}/DESeq2_results.txt".format(test = test, strat = strat) for test in tests for strat in strat_list] 
  params:
    scripts = config['scripts_dir'],
    strats = strats_str
  output:
    results_dir + '{test}_DE.pdf'
  shell:
    "Rscript {params.scripts}/DESeq2_plots.R {output} '{wildcards.test}' '{params.strats}' {input}"

rule make_DE_excel:
  input:
    [results_dir + "{test}/{strat}/DESeq2_results.txt".format(test = test, strat = strat) for test in tests for strat in strat_list]
  params:
    scripts = config['scripts_dir']
  output:
    results_dir + "{test}_DESeq2_p{pval}.xls",
    results_dir + "{test}_DESeq2_sig_p{pval}.xls"
  shell:
    "Rscript {params.scripts}/Create_excel_output.R {output} {wildcards.pval} {input}"

rule GSEA:
  input:
    results_dir + '{test}/{strat}/DESeq2_results.txt',
    results_dir + 'counts/finalCounts.txt',
    "data/gtf.RDS"
  params:
    scripts = config['scripts_dir'],
    gmt = config['pathways'],
    species = config['species'],
    custom_sets = config['custom_sets']
  output:
    results_dir + '{test}/{strat}/fgsea_results.txt',
    results_dir + '{test}/{strat}/fgsea_results.RDS',
    results_dir + '{test}/{strat}/fgsea_results.pdf'
  shell:
    "Rscript {params.scripts}/fgsea.R {input} '{params.gmt}' '{params.species}' {output} {params.custom_sets}"

rule GSEA_plots:
  input:
    results_dir + 'counts/finalCounts.txt',
    results_dir + 'QC/Covariates_QC_metrics_filter.csv',
    'data/gtf.RDS',
    [results_dir + "{test}/{strat}/DESeq2_results.txt".format(test = test, strat = strat) for test in tests for strat in strat_list],
    [results_dir + "{test}/{strat}/fgsea_results.txt".format(test = test, strat = strat) for test in tests for strat in strat_list] 
  params:
    scripts = config['scripts_dir'],
    adjust = config['adjust'],
    gmt = config['pathways'],
    custom_sets = config['custom_sets'],
    strats = strats_str
  output:
    results_dir + '{test}_GSEA.pdf'
  shell:
    "Rscript {params.scripts}/GSEA_plots.R {output} {wildcards.test} '{params.gmt}' '{params.custom_sets}' '{params.strats}' {input}"

rule make_GSEA_excel:
  input:
    [results_dir + "{test}/{strat}/fgsea_results.txt".format(test = test, strat = strat) for test in tests for strat in strat_list]
  params:
    scripts = config['scripts_dir']
  output:
    results_dir + "{test}_GSEA.xls",
    results_dir + "{test}_GSEA_sig.xls"
  shell:
    "Rscript {params.scripts}/Create_excel_output.R {output} 0.05 {input}"

