import pandas as pd
import subprocess
import os
from collections import defaultdict

configfile: 'project_config_file.yaml'
localrules: aggregate_fastqc, aggregate_bams, aggregate_rseqc, mk_cp_sh

config = defaultdict(str, config)
config['filter'] = defaultdict(defaultdict, config['filter'])
config['filter']['sample'] = defaultdict(defaultdict, config['filter']['sample'])
config['filter']['feature'] = defaultdict(defaultdict, config['filter']['feature'])

sample_sheet = os.path.join(os.getcwd(), config['sample_sheet'])
samples = pd.read_table(sample_sheet, sep = ',').set_index("Sample_Name", drop = False)

### Filter information
sample_filter_keys = ';'.join(map(str,config['filter']['sample'].keys()))
sample_filter_values = ';'.join(map(str, config['filter']['sample'].values()))

def get_mem_h_vmem(wildcards, attempt):
  return attempt * 50
def get_fc_csv(wildcards):
  #return os.path.join(config['runs_dir'], samples_to_fc[wildcards.samp], 'SampleSheet.csv')
  return os.path.join(config['runs_dir'], wildcards.FC, 'SampleSheet.csv')
def get_tiles(wildcards):
  return lane_regex[wildcards.FC]
def get_flowcell_fastqs(wildcards):
  FC = samples_to_fc[wildcards.samp]
  return fc_fqs[FC]
def get_fqs_from_FC(wildcards):
  return ['data/fastq/raw/' + x + '.fastq.gz' for x in proj_fqs[wildcards.FC]]
def convert_fq_paths(wildcards):
  return path_dict[wildcards.fastq + '.fastq.gz']
def fc_from_sample(wildcards):
  return id_to_fc[wildcards.samp]
def fc_file_from_sample(wildcards):
  return os.path.join(config['runs_dir'], id_to_fc[wildcards.samp], 'SampleSheet.csv')
def samples_from_sample(wildcards):
  fc = id_to_fc[wildcards.samp]
  return [fc + '_test/' + s + '_test.txt' for s in list(samples[samples['flowcell_full'] == fc].samp)]
def samples_from_fc(wildcards): ### subdir is Z1/samplename_Z1_etc.fastq.gz
  sub = samples[samples['flowcell_full'] == wildcards.FC]
  return list(sub.R1_subdir) + list(sub.R2_subdir)

### Required columns in sample_sheet are: Sample_Name, S_N, Lane, 
### List of fastq files
#samples['R1'] = samples['Sample_Name'].map(str) + '_S' + samples['S_N'].map(str) + '_L00' + samples['Lane'].map(str) + '_R1_001'
#samples['R2'] = samples['Sample_Name'].map(str) + '_S' + samples['S_N'].map(str) + '_L00' + samples['Lane'].map(str) + '_R2_001'
#samples['samp'] = samples['Sample_Name'].map(str) + '_S' + samples['S_N'].map(str) + '_L00' + samples['Lane'].map(str)
samples['R1'] = samples['Sample_Name'].map(str) + '_' + samples['Sample_ID'].map(str) + '_L00' + samples['Lane'].map(str) + '_R1_001'
samples['R2'] = samples['Sample_Name'].map(str) + '_' + samples['Sample_ID'].map(str) + '_L00' + samples['Lane'].map(str) + '_R2_001'
samples['samp'] = samples['Sample_Name'].map(str) + '_' + samples['Sample_ID'].map(str) + '_L00' + samples['Lane'].map(str)
samples['R1_fullpath'] = list([os.path.join(config['runs_dir'], samples.flowcell_full[i], 'demultiplexed', config['project'], samples['Sample_ID'][i], samples['R1'][i] + '.fastq.gz') for i in range(0,len(samples['R1']))])
samples['R2_fullpath'] = list([os.path.join(config['runs_dir'], samples.flowcell_full[i], 'demultiplexed', config['project'], samples['Sample_ID'][i], samples['R2'][i] + '.fastq.gz') for i in range(0,len(samples['R2']))])
runs_fastqs = list(samples['R1_fullpath']) + list(samples['R2_fullpath'])
samples['R1_projpath'] = list([os.path.join('data', 'fastq', 'raw', samples['R1'][i] + '.fastq.gz') for i in range(0,len(samples['R1']))])
samples['R2_projpath'] = list([os.path.join('data', 'fastq', 'raw', samples['R2'][i] + '.fastq.gz') for i in range(0,len(samples['R2']))])
proj_fastqs = list(samples['R1_projpath']) + list(samples['R2_projpath'])
samples['R1_subdir'] = list([os.path.join(samples['Sample_ID'][i], samples['R1'][i] + '.fastq.gz') for i in range(0,len(samples['R1']))])
samples['R2_subdir'] = list([os.path.join(samples['Sample_ID'][i], samples['R2'][i] + '.fastq.gz') for i in range(0,len(samples['R2']))])
fastq_files = list(samples['R1']) + list(samples['R2'])

### For testing only
Lanes = list(set(samples.Lane))
print(Lanes)
### dict of first_path/second_path for fastqs
path_dict = dict(zip(list(samples.R1_projpath) + list(samples.R2_projpath), list(samples.R1_fullpath) + list(samples.R2_fullpath)))
#samples_to_fc = dict(zip(samples.samp, samples.flowcell_full))
id_to_fc = dict(zip(samples.samp, samples.flowcell_full))
### Dictionary of flowcells, FullName:AbbreviatedName
flowcells = dict(zip(samples.flowcell_full, samples.FCID))
### Dictionary whose values are lane regular expressions for each flowcell
lane_regex = dict(zip(flowcells.keys(), ['s_[' + ''.join(map(str, list(set(samples[samples['FCID'] == fc].Lane)))) + ']'  for fc in list(flowcells.values())]))
### Dict whose values are lists of raw fastq files
fc_fqs = dict(zip(flowcells.keys(), [list(samples[samples['FCID'] == fc].R1_fullpath) + list(samples[samples['FCID'] == fc].R2_fullpath) for fc in list(flowcells.values())]))
proj_fqs = dict(zip(flowcells.keys(), [list(samples[samples['FCID'] == fc].R1) + list(samples[samples['FCID'] == fc].R2) for fc in list(flowcells.values())]))
#print([os.path.isfile(f) for f in list(samples['R2_fullpath'])])

### List of test variable names from the configuration file
tests = [x.strip() for x in config['test'].split(',')]
print('Tests:')
print(tests)
strats = [x.strip() for x in config['stratifications'].split(':')]
strat_list = ['All']
if(len(strats) > 1):
  strat_list = strat_list + [strats[0] + '-' + v.strip() for v in strats[1].split(',')]
print('Stratifications: ')
print(strat_list)

def get_test_txt(wildcards):
  return [os.path.join('test/', ID + '_test.txt') for ID in list(samples[samples['flowcell_full'] == wildcards.FC].Sample_ID)]

### Input wildcard is flowcell_full, which is converted to a list of the associated samples
def check_bcl2fastq(wildcards):
  checkpoint_output = checkpoints.run_bcltofastq.get(**wildcards).output[0] 
  samps = samples_from_fc(wildcards)
  return [os.path.join(checkpoint_output, samp) for samp in samps]

rule all:
  input:
    ["results/{test}_GSEA.xls".format(test = test) for test in tests],
    ["results/{test}_DESeq2.xls".format(test = test) for test in tests]

### Split from single file -> FC
rule bcltofastq_sample_sheet:
  input:
    ancient(sample_sheet)
  params:
    scripts = config['scripts_dir'],
    proj = config['project'],
    investigator = config['investigator'],
    ref = config['species'],
    runs_dir = config['runs_dir']
  output:
    [os.path.join(config['runs_dir'], FC, 'SampleSheet.csv') for FC in list(flowcells.keys())]
  shell:
    "Rscript {params.scripts}/FC_sample_sheets.R '{params.proj}' '{params.investigator}' '{params.ref}' '{params.runs_dir}' {input}"

### flowcell -> sample
### Once per flowcell (so FC is a WC), do an split (output is a list)
### Note: it would be best to rearrange such that output is not a directory, as it will be deleted when the rule is rerun
checkpoint run_bcltofastq:
  input:
    ### Requires wildcards called 'FC' and returns the SampleSheet.csv for that FC
    get_fc_csv
  params:
    run_dir = os.path.join(config['runs_dir'], '{FC}'),
    intensities_dir = os.path.join(config['runs_dir'], '{FC}', 'Data', 'Intensities'),
    basecall_dir = os.path.join(config['runs_dir'], '{FC}', 'Data', 'Intensities', 'BaseCalls'),
    out_dir = os.path.join(config['runs_dir'], '{FC}', 'demultiplexed'),
    project = config['project'],
    tiles  = get_tiles
    #flowcell_fastqs = get_flowcell_fastqs
  output:
    directory(config['runs_dir'] + "{FC}/demultiplexed/" + config['project'] + '/')
  shell:
    """
    mkdir -p {params.out_dir}
    bcl2fastq -r 8 -w 8 --tiles {params.tiles} --sample-sheet {input} --fastq-compression-level 9 --create-fastq-for-index-reads -R {params.run_dir} -i {params.basecall_dir} --intensities-dir {params.intensities_dir} -o {params.out_dir}
    """

### For a FC (wildcard) check that bcl2fastq results are done, and return (as input) the list of fastq files
rule FC_bcltofastq:
  input: ### list of full paths to runs_dir fastq files from flowcell FC
    check_bcl2fastq
  output:
    os.path.join(config['runs_dir'], "{FC}", "demultiplexed", config['proj'], "raw_fastq_files.txt")
  shell:
    "ls {input} > {output}"

### FCs is a list
rule aggregate_fc:
  input: 
    [os.path.join(config['runs_dir'], FC, "demultiplexed", config['proj'], "raw_fastq_files.txt") for FC in list(flowcells.keys())]
  output:
    'data/fastq/raw/fastq_files.txt'
  shell:
    """
    mkdir -p data/fastq/raw/
    cat {input} > {output}
    """

rule mk_cp_sh:
  input:
    'data/fastq/raw/fastq_files.txt'
  output:
    'data/fastq/raw/copy_fastqs.sh'
  shell:
    """
    sed 's/^/cp /' {input} > {output}
    sed -i 's/$/ data\/fastq\/raw\/ /' {output}
    """

rule cp_fq:
  input:
    'data/fastq/raw/copy_fastqs.sh'
  output:
    temp(proj_fastqs)
  shell:
    """
    sh {input}
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
  shell:
    "multiqc --file-list {input} -f --filename {output}"

rule run_trimmomatic_once:
  input:
    R1 = "data/fastq/raw/{sample}_R1_001.fastq.gz",
    R2 = "data/fastq/raw/{sample}_R2_001.fastq.gz"
    #R1 = os.path.join(config['runs_dir'], 'demultiplexed', config['project'], 
  output:
    R1_paired = "data/fastq/clipped/paired/{sample}_R1_001.fastq.gz",
    R1_unpaired = "data/fastq/clipped/unpaired/{sample}_R1_001.fastq.gz",
    R2_paired = "data/fastq/clipped/paired/{sample}_R2_001.fastq.gz",
    R2_unpaired = "data/fastq/clipped/unpaired/{sample}_R2_001.fastq.gz"
  shell:
    #"data/bam/multiqc.html"
    "java -Xmx7G -jar /hpcdata/vrc/vrc1_data/douek_lab/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 1 {input.R1} {input.R2} {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} ILLUMINACLIP:/hpcdata/vrc/vrc1_data/douek_lab/programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:4:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50"

rule run_STAR_once:
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

rule run_index_once:
  input:
    "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
  output:
    "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
  shell:
    "samtools index {input}"

###  These two RSeQC functions take so long that they will run separately
rule run_RSeQC_tin:
  input:
    bam = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    bai = "data/bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
  params:
    anno = config['refbed'],
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
    anno = config['refbed'],
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
    anno = config['refbed'],
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
    'data/Covariates_QC_metrics.csv'
  params:
    scripts = config['scripts_dir'],
    test = config['test'],
    strat = config['stratifications'],
    adjust = config['adjust'],
    batch = config['batch_candidates'],
    n_reads = config['filter']['feature']['counts'],
    frac = config['filter']['feature']['fraction_samples'],
    filter_name = sample_filter_keys,
    filter_value = sample_filter_values
  output:
    temp('data/Covariates_QC_metrics_intermediate.csv'),
    'data/counts/filteredCounts.txt',
    'data/sample_filters.pdf',
    'data/sample_filters.txt'
  shell:
    """
    Rscript {params.scripts}/Filters.R {input} {output} '{params.test}' '{params.strat}' '{params.adjust}' '{params.batch}' '{params.n_reads}' '{params.frac}' '{params.filter_name}' '{params.filter_value}'
    """

rule filter_samples_norm:
  input:
    'data/counts/normalizedCounts.txt',
    'data/Covariates_QC_metrics_intermediate.csv'
  params:
    scripts = config['scripts_dir'],
    pca = config['filter']['sample']['PCA']
  output:
    'data/Covariates_QC_metrics_filter.csv',
    'data/PCA_filter.pdf'
  shell:
    "Rscript {params.scripts}/PCA_filter.R {input} '{params.pca}' {output}"

rule Normalization:
  input:
    'data/counts/filteredCounts.txt',
    'data/Covariates_QC_metrics_intermediate.csv'
  output:
    'data/counts/normalizedCounts.txt',
    temp('data/counts/normalizedCounts.RDS'),
    'data/counts/normalization.pdf'
  params:
    scripts = config['scripts_dir'],
    test = config['test'],
    adjust = config['adjust'],
    batch = config['batch_candidates']
  shell:
    """
    Rscript {params.scripts}/DESeq2_normalization.R {input} '{params.test}' '{params.adjust}' '{params.batch}' {output}
    """

rule Evaluate_batch:
  input:
    'data/counts/normalizedCounts.txt',
    'data/Covariates_QC_metrics_filter.csv',
  params:
    scripts = config['scripts_dir'],
    candidates = config['batch_candidates'],
    method = config['batch_method']
  output:
    'data/batch_evaluation.txt',
    'data/batch_evaluation.pdf'
  shell:
    """
    Rscript {params.scripts}/Batch_evaluation.R {input} '{params.candidates}' '{params.method}' {output}
    """

rule ComBat:
  input:
    'data/counts/normalizedCounts.txt',
    'data/counts/normalizedCounts.RDS',
    'data/Covariates_QC_metrics_filter.csv',
    'data/batch_evaluation.txt'
  params:
    scripts = config['scripts_dir'],
  output:
    'data/counts/finalCounts.txt',
    'data/counts/finalCounts.RDS'
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
    'data/counts/finalCounts.txt',
    'data/counts/finalCounts.RDS',
    'data/Covariates_QC_metrics_filter.csv',
    'data/batch_evaluation.txt'
  params:
    scripts = config['scripts_dir'],
    adjust = config['adjust'],
    gtf = config['refgtf']
  output:
    'results/{test}/{strat}/DESeq2_results.txt',
    'results/{test}/{strat}/DESeq2_results.pdf'
  shell:
    """
    mkdir -p 'results/{wildcards.test}/{wildcards.strat}/'
    Rscript {params.scripts}/DESeq2.R {input} '{wildcards.test}' '{params.adjust}' '{params.gtf}' '{wildcards.strat}' {output}
    """
  
rule make_DE_excel:
  input:
    ["results/{test}/{strat}/DESeq2_results.txt".format(test = test, strat = strat) for test in tests for strat in strat_list]
  params:
    scripts = config['scripts_dir']
  output:
    "results/{test}_DESeq2.xls"
  shell:
    "Rscript {params.scripts}/Create_excel_output.R {output} {input}"

rule GSEA:
  input:
    'results/{test}/{strat}/DESeq2_results.txt',
    'data/counts/finalCounts.txt'
  params:
    scripts = config['scripts_dir'],
    gmt = config['pathways'],
    species = config['species']
  output:
    'results/{test}/{strat}/fgsea_results.txt',
    'results/{test}/{strat}/fgsea_results.pdf'
  shell:
    "Rscript {params.scripts}/fgsea.R {input} '{params.gmt}' '{params.species}' {output}"

rule make_GSEA_excel:
  input:
    ["results/{test}/{strat}/fgsea_results.txt".format(test = test, strat = strat) for test in tests for strat in strat_list]
  params:
    scripts = config['scripts_dir']
  output:
    "results/{test}_GSEA.xls"
  shell:
    "Rscript {params.scripts}/Create_excel_output.R {output} {input}"

