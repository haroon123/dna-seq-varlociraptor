WORKING_DIR = "tmp/gridssworkingdir"

# rule griddss:
#     input:
#         bams=lambda wildcards: expand("results/recal/{sample}.sorted.bam", sample=get_group_samples(wildcards)),
#         bais=lambda wildcards: expand("results/recal/{sample}.sorted.bai", sample=get_group_samples(wildcards)),
#         ref="resources/genome.fasta",
#         idx=rules.bwa_index.output,
#     output:
#         vcf="results/gridss_vcf/{group}.vcf",
#         assembly="results/gridss_assembly/{group}.bam",
#     conda:
#         "../envs/gridss.yaml"
#     threads:
#         8
#     shell:
#         "./workflow/scripts/gridss.sh --reference {input.ref} --output {output.vcf} --assembly {output.assembly} --threads {threads} --jar workflow/scripts/gridss-2.9.4-gridss-jar-with-dependencies.jar --workingdir {WORKING_DIR} {input.bams}"


rule gridss_setupreference:
    input:
        reference="resources/genome.fasta",
        dictionary="resources/genome.dict",
        indices=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        multiext("resources/genome.fasta", ".gridsscache", ".img")
    log:
        "log/gridss/setupreference.log"
    wrapper:
        "file:/vol/huge/christo/snakemake-wrappers/bio/gridss/setupreference"


preprocess_endings = (".cigar_metrics", ".coverage.blacklist.bed", ".idsv_metrics", ".insert_size_histogram.pdf", ".insert_size_metrics", ".mapq_metrics", ".sv.bam", ".sv.bam.bai", ".sv_metrics", ".tag_metrics")

rule gridss_preprocess:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bai",
        reference="resources/genome.fasta",
        dictionary="resources/genome.dict",
        refindex=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa", ".dict", ".gridsscache", ".img")
    output:
        multiext("{WORKING_DIR}/{sample}.bam.gridss.working/{sample}.bam", *preprocess_endings)
    params:
        workingdir=WORKING_DIR
    log:
        "log/gridss/preprocess/{WORKING_DIR}/{sample}.preprocess.log"
    threads:
        8
    wrapper:
        "file:/vol/huge/christo/snakemake-wrappers/bio/gridss/preprocess"


rule gridss_assemble:
    input:
        bams=lambda wildcards: expand("results/recal/{sample}.sorted.bam", sample=get_group_samples(wildcards)),
        bais=lambda wildcards: expand("results/recal/{sample}.sorted.bai", sample=get_group_samples(wildcards)),
        reference="resources/genome.fasta",
        dictionary="resources/genome.dict",
        indices=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa", ".dict", ".gridsscache", ".img"),
        preprocess=lambda wildcards: expand("{working_dir}/{sample}.bam.gridss.working/{sample}.bam{ending}", working_dir=[WORKING_DIR], sample=get_group_samples(wildcards), ending=preprocess_endings)
    output:
        assembly="results/gridss_assembly/{group}.bam"
    params:
        workingdir=WORKING_DIR
    log:
        "log/gridss/assemble/{group}.log"
    threads:
        100
    wrapper:
        "file:/vol/huge/christo/snakemake-wrappers/bio/gridss/assemble"


rule gridss_call:
    input:
        bams=lambda wildcards: expand("results/recal/{sample}.sorted.bam", sample=get_group_samples(wildcards)),
        bais=lambda wildcards: expand("results/recal/{sample}.sorted.bai", sample=get_group_samples(wildcards)),
        reference="resources/genome.fasta",
        dictionary="resources/genome.dict",
        indices=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa", ".dict", ".gridsscache", ".img"),
        preprocess=lambda wildcards: expand("{working_dir}/{sample}.bam.gridss.working/{sample}.bam{ending}", working_dir=[WORKING_DIR], sample=get_group_samples(wildcards), ending=preprocess_endings),
        assembly="results/assembly/{group}.bam"
    output:
        vcf="results/gridss_vcf/{group}.vcf"
    params:
        workingdir=WORKING_DIR
    log:
        "log/gridss/call/{group}.log"
    threads:
        100
    wrapper:
        "file:/vol/huge/christo/snakemake-wrappers/bio/gridss/call"