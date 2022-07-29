# I realized that this is not used at all, but I'm leaving it for thinking
pepfile: "sample-info/hca-downstream-pep.yaml"

rule target:
    input:
        "results/human_cell_atlas/anndata/merged_anndata_objects",
        "results/human_cell_atlas/integrated-sce-objects/1M_Immune_Cells_integrated_fastmnn_sce.rds",
        "results/human_cell_atlas/integrated-sce-objects/1M_Immune_Cells_integrated_harmony_sce.rds"
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
        sce_dir = "{basedir}/scpca-downstream-analyses"
    output:
        directory("{basedir}/merged-sce-objects")
    params:
        grouping_var = "project_name"
    shell:
        """
        Rscript scripts/02-prepare-merged-sce.R \
          --library_file "{input.processed_tsv}" \
          --sce_dir "{input.sce_dir}" \
          --grouping_var {params.grouping_var} \
          --merged_sce_dir "{output}"
        """

rule convert_sce_anndata:
    input:
        processed_tsv = "sample-info/hca-processed-libraries.tsv",
        merged_sce_dir = "{basedir}/merged-sce-objects"
    output:
        directory("{basedir}/anndata/merged_anndata_objects")
    params:
        grouping_var = "project_name"
    shell:
        """
        Rscript scripts/02a-convert-sce-to-anndata.R \
          --library_file "{input.processed_tsv}" \
          --merged_sce_dir "{input.merged_sce_dir}" \
          --grouping_var {params.grouping_var} \
          --anndata_output_dir "{output}"
        """


rule integrate_fastmnn:
    conda: "envs/scpca-renv.yaml"
    input:
        merged_sce_dir = "{basedir}/merged-sce-objects"
    output:
        "{basedir}/integrated-sce-objects/{project}_integrated_fastmnn_sce.rds"
    params:
        merged_sce_file = "{project}_merged_sce.rds",
        seed = 2022
    shell:
       """
       Rscript scripts/03-integrate-sce.R \
         --input_sce_file "{input.merged_sce_dir}/{params.merged_sce_file}" \
         --output_sce_file "{output}" \
         --method fastMNN \
         --seed {params.seed}
       """

rule integrate_harmony:
    conda: "envs/scpca-renv.yaml"
    input:
        merged_sce_dir = "{basedir}/merged-sce-objects"
    output:
        "{basedir}/integrated-sce-objects/{project}_integrated_harmony_sce.rds"
    params:
        merged_sce_file = "{project}_merged_sce.rds",
        seed = 2022
    shell:
       """
       Rscript scripts/03-integrate-sce.R \
         --input_sce_file "{input.merged_sce_dir}/{params.merged_sce_file}" \
         --output_sce_file "{output}" \
         --method harmony \
         --seed {params.seed}
       """
