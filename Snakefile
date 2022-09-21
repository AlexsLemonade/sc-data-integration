configfile: "config/config.yaml"
pepfile: config['pepfile']
INTEGRATION_METHODS = ["fastmnn", "harmony", "seurat-cca", "seurat-rpca", "scanorama", "scvi"]

rule target:
    input:
        expand(os.path.join(config['results_dir'], "integrated_sce/{project}_integrated_{integration_methods}_sce.rds"),
               project = pep.sample_table["project_name"],
               integration_methods = INTEGRATION_METHODS),
        expand(os.path.join(config['results_dir'], "analysis_reports/{project}_integration_report.html"),
               project = pep.sample_table["project_name"])


# Rule used for building conda & renv environment
rule build_renv:
    input: workflow.source_path("renv.lock")
    output: "renv/.snakemake_timestamp"
    log: "logs/build_renv.log"
    conda: "envs/scpca-renv.yaml"
    shell:
        """
        Rscript -e "renv::restore(lockfile = '{input}')" &> {log}
        date -u -Iseconds  > {output}
        """

rule merge_sces:
    conda: "envs/scpca-renv.yaml"
    input:
        processed_tsv = config["processed_tsv"],
        sce_dir = config["sce_dir"]
    output:
        directory(os.path.join(config["results_dir"], "merged_sce"))
    log:
        "logs/merge_sce.log"
    shell:
        """
        Rscript scripts/02-prepare-merged-sce.R \
          --library_file "{input.processed_tsv}" \
          --sce_dir "{input.sce_dir}" \
          --add_celltype {config[add_celltype]} \
          --celltype_info "{config[celltype_file]}" \
          --grouping_var {config[grouping_var]} \
          --merged_sce_dir "{output}" \
          --num_hvg {config[num_hvg]} \
          --subset_hvg \
          &> {log}
        """

rule convert_sce_anndata:
    conda: "envs/scpca-renv.yaml"
    input:
        processed_tsv = config["processed_tsv"],
        merged_sce_dir = "{basedir}/merged_sce"
    output:
        directory("{basedir}/merged_anndata")
    log:
        "logs/{basedir}/convert_sce_anndata.log"
    shell:
        """
        Rscript scripts/02a-convert-sce-to-anndata.R \
          --library_file "{input.processed_tsv}" \
          --merged_sce_dir "{input.merged_sce_dir}" \
          --grouping_var {config[grouping_var]} \
          --anndata_output_dir "{output}" \
          &> {log}
        """


rule integrate_fastmnn:
    conda: "envs/scpca-renv.yaml"
    input:
        # The input has to be the merged directory so snakemake can find it.
        # We will add the file name with params.
        merged_sce_dir = "{basedir}/merged_sce"
    output:
        "{basedir}/integrated_sce/{project}_integrated_fastmnn_sce.rds"
    log:
        "logs/{basedir}/integrate_fastmnn-{project}.log"
    params:
        merged_sce_file = "{project}_merged_sce.rds",
    shell:
        """
        Rscript scripts/03a-integrate-sce.R \
          --input_sce_file "{input.merged_sce_dir}/{params.merged_sce_file}" \
          --output_sce_file "{output}" \
          --method fastMNN \
          --seed {config[seed]} \
          --corrected_only \
          &> {log}
        """

rule integrate_harmony:
    conda: "envs/scpca-renv.yaml"
    input:
        # The input has to be the merged directory so snakemake can find it.
        # We will add the file name with params.
        merged_sce_dir = "{basedir}/merged_sce"
    output:
        "{basedir}/integrated_sce/{project}_integrated_harmony_sce.rds"
    log:
        "logs/{basedir}/integrate_harmony-{project}.log"
    params:
        merged_sce_file = "{project}_merged_sce.rds",
    shell:
        """
        Rscript scripts/03a-integrate-sce.R \
          --input_sce_file "{input.merged_sce_dir}/{params.merged_sce_file}" \
          --output_sce_file "{output}" \
          --method harmony \
          --seed {config[seed]} \
          --corrected_only \
          &> {log}
        """

rule integrate_seurat:
    conda: "envs/scpca-renv.yaml"
    input:
        # The input has to be the merged directory so snakemake can find it.
        # We will add the file name with params.
        merged_sce_dir = "{basedir}/merged_sce"
    output:
        "{basedir}/integrated_sce/{project}_integrated_seurat-{method}_sce.rds"
    log:
        "logs/{basedir}/integrate_seurat-{method}-{project}.log"
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
          --corrected_only \
          &> {log}
        """

rule integrate_scanorama:
    conda: "envs/scanorama.yaml"
    input:
        merged_anndata_dir = "{basedir}/merged_anndata"
    output:
        temp("{basedir}/integrated_anndata/{project}_integrated_scanorama.h5")
    log:
        "logs/{basedir}/integrate_scanorama-{project}.log"
    params:
        merged_anndata_file = "{project}_anndata.h5",
    shell:
        """
        python scripts/03b-integrate-scanorama.py \
          --input_anndata "{input.merged_anndata_dir}/{params.merged_anndata_file}" \
          --output_anndata "{output}" \
          --seed {config[seed]} \
          --use_hvg \
          --corrected_only \
          &> {log}
        """

rule integrate_scvi:
    conda: "envs/scvi.yaml"
    input:
        merged_anndata_dir = "{basedir}/merged_anndata"
    output:
        temp("{basedir}/integrated_anndata/{project}_integrated_scvi.h5")
    log:
        "logs/{basedir}/integrate_scvi-{project}.log"
    params:
        merged_anndata_file = "{project}_anndata.h5",
    shell:
        """
        python scripts/03c-integrate-scvi.py \
          --input_anndata "{input.merged_anndata_dir}/{params.merged_anndata_file}" \
          --output_anndata "{output}" \
          --continuous_covariates {config[continuous_covariates]} \
          --num_latent {config[num_latent]} \
          --seed {config[seed]} \
          --use_hvg \
          --corrected_only \
          &> {log}
        """

rule convert_anndata_sce:
    conda: "envs/scpca-renv.yaml"
    input:
        "{basedir}/integrated_anndata/{project}_integrated_{method}.h5"
    output:
        "{basedir}/integrated_sce/{project}_integrated_{method}_sce.rds"
    log:
        "logs/{basedir}/convert_anndata_sce-{method}-{project}.log"
    shell:
        """
        Rscript scripts/04-post-process-anndata.R \
          --input_anndata_file "{input}" \
          --output_sce_file "{output}" \
          --method "{wildcards.method}" \
          --seed {config[seed]} \
          --corrected_only \
          &> {log}
        """


rule generate_report:
    conda: "envs/scpca-renv.yaml"
    input:
        merged_sce_dir = rules.merge_sces.output,
        integrated_sce_files = expand("{{basedir}}/integrated_sce/{{project}}_integrated_{integration_method}_sce.rds",
                                      integration_method = INTEGRATION_METHODS)
    output:
        "{basedir}/analysis_reports/{project}_integration_report.html"
    log:
        "logs/{basedir}/generate_report-{project}.log"
    params:
        integrated_sce_dir = "{basedir}/integrated_sce"
    shell:
        """
        Rscript -e \
        "rmarkdown::render('analysis_templates/01-single-group-integration-check-template.Rmd', \
                            clean = TRUE, \
                            output_file = '{output}', \
                            output_dir = dirname('{output}'), \
                            params = list(group_name = '{wildcards.project}', \
                                          merged_sce_dir = '{workflow.basedir}/{input.merged_sce_dir}', \
                                          integrated_sce_dir = '{workflow.basedir}/{params.integrated_sce_dir}'))" \
          &> {log}
        """
