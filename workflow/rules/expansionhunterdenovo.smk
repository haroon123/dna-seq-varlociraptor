rule ehdn_merge:
    input:
        bam=expand("results/strling/extract/{sample}.bin", sample=samples["sample_name"].values),
        ref="resources/genome.fasta",
    output:
        merged="results/strling/merge/merged-bounds.txt"
    params:
        prefix="results/strling/merge/merged"
    log:
        "logs/ehdn/merge.log"
    shell:
        """
        /path/to/ExpansionHunterDenovo merge \
        --reference reference.fasta \
        --manifest manifest.tsv \
        --output-prefix example_dataset
        """


rule ehdn_profile:
    input:
        bam="results/recal/{sample}.sorted.bam",
        ref="resources/genome.fasta",
    output:
        json="results/ehdn/profile/{sample}.str_profile.json",
    params:
        prefix="results/ehdn/profile/{sample}"
    log:
        "logs/ehdn/profile/{sample}.extract.log"
    shell:
        """
        ExpansionHunterDenovo profile \
        --reads {bam} \
        --reference {ref} \
        --output-prefix {params.prefix} \
        --min-anchor-mapq 50 \
        --max-irr-mapq 40 2> {log}
        """