rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi"
    log:
        "logs/bcf-index/{prefix}.log"
    conda:
        "../envs/bcftools.yaml"
    threads:
        8
    shell:
        "bcftools index --threads {threads} {input} 2> {log}"


    # wildcard_constraints:
    #     prefix="^results"

rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "logs/bam-index/{prefix}.log"
    wildcard_constraints:
        prefix=r"^results/.*"
    wrapper:
        "0.59.2/bio/samtools/index"


rule bam_index_temp:
    input:
        "{prefix}.bam"
    output:
        temp("{prefix}.bam.bai")
    log:
        "logs/bam-index/{prefix}.log"
    wildcard_constraints:
        prefix=r"^tmp/.*"
    wrapper:
        "0.59.2/bio/samtools/index"


rule tabix_known_variants:
    input:
        "resources/{prefix}.{format}.gz"
    output:
        "resources/{prefix}.{format}.gz.tbi"
    log:
        "logs/tabix/{prefix}.{format}.log"
    params:
        get_tabix_params
    cache: True
    wrapper:
        "0.59.2/bio/tabix"
