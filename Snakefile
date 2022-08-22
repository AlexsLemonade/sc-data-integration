pepfile: "sample-info/hca-project-pep.yaml"

rule target:
    input:
        expand("results/human_cell_atlas/integrated_sce/{project}_integrated_{sce_method}_sce.rds",
               project = pep.sample_table["project_name"],
               sce_method = ["fastmnn", "harmony", "seurat-cca", "seurat-rpca", "scanorama", "scvi"])

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
        sce_dir = "{basedir}/scpca-downstream-analyses",
        celltype_file = "sample-info/hca-celltype-info.tsv"
    output:
        directory("{basedir}/merged_sce")
    params:
        grouping_var = "project_name",
        num_hvg = 5000,
    shell:
        """
        Rscript scripts/02-prepare-merged-sce.R \
          --library_file "{input.processed_tsv}" \
          --sce_dir "{input.sce_dir}" \
          --celltype_info "{input.celltype_file}" \
          --grouping_var {params.grouping_var} \
          --merged_sce_dir "{output}" \
          --num_hvg {params.num_hvg} \
          --subset_hvg
        """

rule convert_sce_anndata:
    input:
        processed_tsv = "sample-info/hca-processed-libraries.tsv",
        merged_sce_dir = "{basedir}/merged_sce"
    output:
        directory("{basedir}/merged_anndata")
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
        # The input has to be the merged directory so snakemake can find it.
        # We will add the file name with params.
        merged_sce_dir = "{basedir}/merged_sce"
    output:
        "{basedir}/integrated_sce/{project}_integrated_fastmnn_sce.rds"
    params:
        merged_sce_file = "{project}_merged_sce.rds",
        seed = 2022
    shell:
        """
        Rscript scripts/03a-integrate-sce.R \
          --input_sce_file "{input.merged_sce_dir}/{params.merged_sce_file}" \
          --output_sce_file "{output}" \
          --method fastMNN \
          --seed {params.seed} \
          --corrected_only
        """

rule integrate_harmony:
    conda: "envs/scpca-renv.yaml"
    input:
        # The input has to be the merged directory so snakemake can find it.
        # We will add the file name with params.
        merged_sce_dir = "{basedir}/merged_sce"
    output:
        "{basedir}/integrated_sce/{project}_integrated_harmony_sce.rds"
    params:
        merged_sce_file = "{project}_merged_sce.rds",
        seed = 2022
    shell:
        """
        Rscript scripts/03a-integrate-sce.R \
          --input_sce_file "{input.merged_sce_dir}/{params.merged_sce_file}" \
          --output_sce_file "{output}" \
          --method harmony \
          --seed {params.seed} \
          --corrected_only
        """

rule integrate_seurat:
    conda: "envs/scpca-renv.yaml"
    input:
        # The input has to be the merged directory so snakemake can find it.
        # We will add the file name with params.
        merged_sce_dir = "{basedir}/merged_sce"
    output:
        "{basedir}/integrated_sce/{project}_integrated_seurat-{method}_sce.rds"
    wildcard_constraints:
        method = "cca|rpca"
    params:
        merged_sce_file = "{project}_merged_sce.rds",
        num_genes = 2000
    shell:
        """
        Rscript scripts/03a-integrate-sce.R \
          --input_sce_file "{input.merged_sce_dir}/{params.merged_sce_file}" \
          --output_sce_file "{output}" \
          --method seurat \
          --seurat_reduction_method {wildcards.method} \
          --num_genes {params.num_genes} \
          --corrected_only
        """

rule integrate_scanorama:
    conda: "envs/scanorama.yaml"
    input:
        merged_anndata_dir = "{basedir}/merged_anndata"
    output:
        integrated_anndata = temp("{basedir}/integrated_anndata/{project}_integrated_scanorama.h5"),
        integrated_sce = "{basedir}/integrated_sce/{project}_integrated_scanorama_sce.rds"
    params:
        merged_anndata_file = "{project}_anndata.h5",
        seed = 2022
    shell:
        """
        python scripts/03b-integrate-scanorama.py \
          --input_anndata "{input.merged_anndata_dir}/{params.merged_anndata_file}" \
          --output_anndata "{output.integrated_anndata}" \
          --seed {params.seed} \
          --use_hvg \
          --corrected_only

        Rscript scripts/04-post-process-anndata.R \
            --input_anndata_file "{output.integrated_anndata}" \
            --output_sce_file "{output.integrated_sce}" \
            --method "scanorama" \
            --seed {params.seed} \
            --corrected_only
        """

rule integrate_scvi:
    conda: "envs/scvi.yaml"
    input:
        merged_anndata_dir = "{basedir}/merged_anndata"
    output:
        integrated_anndata = temp("{basedir}/integrated_anndata/{project}_integrated_scvi.h5"),
        integrated_sce = "{basedir}/integrated_sce/{project}_integrated_scvi_sce.rds"
    params:
        merged_anndata_file = "{project}_anndata.h5",
        seed = 2022
    shell:
        """
        python scripts/03c-integrate-scvi.py \
          --input_anndata "{input.merged_anndata_dir}/{params.merged_anndata_file}" \
          --output_anndata "{output.integrated_anndata}" \
          --seed {params.seed} \
          --use_hvg \
          --corrected_only

        Rscript scripts/04-post-process-anndata.R \
            --input_anndata_file "{output.integrated_anndata}" \
            --output_sce_file "{output.integrated_sce}" \
            --method "scvi" \
            --seed {params.seed} \
            --corrected_only
        """
