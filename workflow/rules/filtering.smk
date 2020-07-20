def get_filter_expression(w):
    expression = config["calling"]["filter"][w.filter].get("expression", None)
    if expression is None:
        return ""
    return f"vembrane - \"{expression}\" |"

def get_filter_region(w):
    region = config["calling"]["filter"][w.filter].get("region", None)
    if region is None:
        return ""
    return f"-T \"{region}\""


rule filter_by_annotation:
    input:
        get_annotated_bcf
    output:
        "results/calls/{group}.{filter}.filtered.bcf"
    log:
        "logs/filter-calls/{group}.{filter}.log"
    params:
        filter_expression=get_filter_expression,
        region=get_filter_region
    conda:
        "../envs/vembrane.yaml"
    shell:
        "(bcftools view {input} {params.region} | sed 's/#CHROM/##INFO=<ID=END,Number=1,Type=Integer>\\n#CHROM/' | {params.filter_expression} bcftools view -Ob > {output}) 2> {log}"


def pre_fdr_command(wc):
    if config["calling"]["fdr-control"]["events"][wc.event].get("compound_het", False):
        return "| python workflow/scripts/compound_heterozygous.py - PROB_COMPOUND_HETEROZYGOUS_CANDIDATE_MOTHER PROB_COMPOUND_HETEROZYGOUS_CANDIDATE_FATHER | bcftools view -Ob"
    else:
        return ""

rule pre_fdr:
    input:
        "results/calls/{group}.{filter}.filtered.bcf"
    output:
        "tmp/pre_fdr/{group}.{vartype}.{event}.{filter}.fdr-controlled.bcf"
    params:
        command=pre_fdr_command
    conda:
        "../envs/pre_fdr.yaml"
    shell:
        "cat {input} {params.command} > {output}"
        #"bcftools view {input} | {params.compound}"


def control_fdr_events(wc):
    if config["calling"]["fdr-control"]["events"][wc.event].get("compound_het", False):
        return "COMPOUND_HETEROZYGOUS"
    else:
        return config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]

rule control_fdr:
    input:
        "tmp/pre_fdr/{group}.{vartype}.{event}.{filter}.fdr-controlled.bcf"
    output:
        "results/calls/{group}.{vartype}.{event}.{filter}.fdr-controlled.bcf"
    log:
        "logs/control-fdr/{group}.{vartype}.{event}.{filter}.log"
    params:
        threshold=config["calling"]["fdr-control"]["threshold"],
        events=control_fdr_events
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events {params.events} --fdr {params.threshold} > {output} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_calls_input(".bcf"),
        idx=get_merge_calls_input(".bcf.csi")
    output:
        "results/merged-calls/{group}.{event}.fdr-controlled.bcf"
    log:
        "logs/merge-calls/{group}.{event}.log"
    params:
        "-a -Ob"
    wrapper:
        "0.59.2/bio/bcftools/concat"
