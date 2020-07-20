# rule download_snpeff_db:
#     output:
#         directory("resources/snpeff/{ref}")
#     log:
#         "logs/download-snpeff-db/{ref}.log"
#     params:
#         db_dir=lambda _, output: str(Path(output[0]).parent.resolve()),
#         ref="{ref}"
#     conda:
#         "../envs/snpeff.yaml"
#     cache: True
#     shell:
#         "snpEff download -dataDir {params.db_dir} {params.ref} 2> {log}"

# rule snpeff:
#     input:
#         calls="results/calls/{group}.bcf",
#         db="resources/snpeff/{build}.{snpeff_release}".format(**config["ref"])
#     output:
#         calls="results/calls/{group}.annotated.bcf",
#         stats="results/snpeff/{group}.html",
#         csvstats="results/snpeff/{group}.csv"
#     log:
#         "logs/snpeff/{group}.log"
#     params:
#         reference="{build}.{snpeff_release}".format(**config["ref"]),
#         data_dir=lambda _, input: Path(input.db).parent.resolve(),
#         extra="-Xmx4g -nodownload"
#     resources:
#         mem_mb=4000
#     wrapper:
#         "0.50.4/bio/snpeff"

rule annotate_variants:
    threads: 8
    input:
        calls="results/calls/{group}.bcf",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins"
    output:
        calls="results/calls/{group}.annotated.bcf",
        stats=report("results/calls/{group}.stats.html", caption="../report/stats.rst", category="QC")
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["annotations"]["vep"]["plugins"],
        extra="{} --vcf_info_field ANN".format(config["annotations"]["vep"]["params"])
    log:
        "logs/vep/{group}.annotate.log"
    wrapper:
        "0.60.0/bio/vep/annotate"


# TODO What about multiple ID Fields?
rule annotate_vcfs:
    threads:
        100
    input:
        bcf="results/calls/{prefix}.bcf",
        csi="results/calls/{prefix}.bcf.csi",
        annotations=get_annotation_vcfs(),
        idx=get_annotation_vcfs(idx=True)
    output:
        "results/calls/{prefix}.db-annotated.bcf"
    log:
        "logs/annotate-vcfs/{prefix}.log"
    params:
        extra="-Xmx4g",
        pipes=get_annotation_pipes
    conda:
        "../envs/snpsift.yaml"
    shell:
        "(python workflow/scripts/parallel_vcf.py --threads {threads} {input.bcf} '{params.pipes}' | bcftools view --threads {threads} -Ob > {output}) 2> {log}"


rule annotate_dgidb:
    threads:
        4
    input:
        "results/calls/{prefix}.bcf"
    output:
        "results/calls/{prefix}.dgidb.bcf"
    log:
        "logs/annotate-dgidb/{prefix}.log"
    conda:
        "../envs/annotate_dgidb.yaml"
    resources:
        dgidb_requests=1
    shell:
        "rbt vcf-annotate-dgidb {input} > {output} 2> {log}"
