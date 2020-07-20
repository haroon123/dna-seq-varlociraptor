def comp_het_event(wc):
    if config["calling"]["fdr-control"]["events"][wc.event].get("compound_het", False):
        return ["COMPOUND_HETEROZYGOUS"]
    return []

rule vcf2xl:
    input:
        "results/merged-calls/{group}.{event}.fdr-controlled.bcf"
    output:
        "results/merged-calls/{group}.{event}.fdr-controlled.xlsx"
    params:
        events=lambda wc: ["PROB_" + x.upper() for x in comp_het_event(wc) + config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]],
        eventnames=lambda wc: ["p_" + x.lower() for x in comp_het_event(wc) + config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]]
    conda:
        "../envs/tsv2xl.yaml"
    shell:
        """
        bcftools view {input} | python workflow/scripts/vcf2xl.py /dev/stdin --fields ANN[*][1] ANN[*][2] ANN[*][3] ANN[*][6] {params.events} gnomad_AF --fieldnames region impact gene transcript {params.eventnames} gnomad_af depth_ est_gen_ cnt_var_ --genotype DP AF OBS[V+,V-] -out {output}
        """
        #"bcftools view {input} | python workflow/scripts/vcf2tsv.py /dev/stdin --fields ANN[*][1] --fields ANN[*][1] ANN[*][3] --fieldnames region gene dp_ af_ obs_ --genotype DP AF OBS[V+,V-] > {output}"