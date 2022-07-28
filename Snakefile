pepfile: "sample-info/hca-downstream-pep.yaml"

rule target:
    input:
        "rtest.log"
        # expand("{sample}.txt", sample=pep.sample_table["sample_id"])

# Dummy rule used for building conda & renv environment
rule build_renv:
    input: "renv.lock"
    output: "renv/.snakemake_timestamp"
    conda: "envs/scpca-renv.yaml"
    shell:
      """
      Rscript -e "renv::restore()"
      date -u -Iseconds  > {output}
      """



rule test_R_integration:
    conda: "envs/scpca-renv.yaml"
    input:
        script = workflow.source_path("scripts/utils/test_integration-functions.R")
    output:
        "rtest.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        "Rscript {input.script} > {output}"
