# I realized that this is not used at all, but I'm leaving it for thinking
pepfile: "sample-info/hca-downstream-pep.yaml"

rule target:
    input:
        "results/human_cell_atlas/merged-sce-objects"
        # expand("{sample}.txt", sample=pep.sample_table["sample_id"])

# Rule used for building conda & renv environment
rule build_renv:
    input: workflow.source_path("renv.lock")
    output: "renv/.snakemake_timestamp"
    conda: "envs/scpca-renv.yaml"
    shell:
      """
      Rscript -e "renv::restore('{input}')"
      date -u -Iseconds  > {output}
      """



rule merge_sces:
    conda: "envs/scpca-renv.yaml"
    input:
        processed_tsv = "sample-info/hca-processed-libraries.tsv",
        sce_dir = "results/human_cell_atlas/scpca-downstream-analyses"
    output:
        directory("results/human_cell_atlas/merged-sce-objects")
    params:
        grouping_var = "project_name"
    shell:
        """
        Rscript scripts/02-prepare-merged-sce.R \
          --library_file {input.processed_tsv} \
          --grouping_var {params.grouping_var} \
          --sce_dir {input.sce_dir} \
          --merged_sce_dir {output}
        """
