"""
Improvements:
- create_seurat can be run in parallel, and then merge?
- link dimensionality_reduction with celltype_validation 
"""
import os
from re import search
import getpass


############
## Config ##
############
host = os.uname()[1]
configfile: "/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/code/snakemake/config_bart.yaml"

# validate(config, schema="schemas/config.schema.yaml")

###########
## Rules ##
###########

rule all:
    input: 
        expand("%s/rna/qc/sample_metadata_after_qc.txt.gz" % config["directories"]["results"]),
        expand("%s/rna/doublet_detection/sample_metadata_after_doublets.txt.gz" % config["directories"]["results"]),
        config["directories"]["results"]+"/rna/mapping/sample_metadata_after_mapping.txt.gz",
        expand("%s/rna/mapping/mapping_mnn_all_samples.txt.gz" % config["directories"]["results"]),
        expand("%s/rna/mapping/mapping_mnn_{sample}.txt.gz" % config["directories"]["results"], sample=config["samples"]),
        config["directories"]["results"]+"/rna/mapping/pdf/umap_mapped_allcells.pdf",
        expand("%s/rna/dimensionality_reduction/sce/pca_features{dimred_sce_features}_pcs{dimred_sce_npcs}.txt.gz" % (config["directories"]["results"]), 
            dimred_sce_features = config["dimensionality_reduction_sce"]["features"], 
            dimred_sce_npcs = config["dimensionality_reduction_sce"]["npcs"]
            ),
        expand("%s/rna/dimensionality_reduction/sce/no_ExE/pca_features{dimred_sce_features}_pcs{dimred_sce_npcs}.txt.gz" % (config["directories"]["results"]), 
            dimred_sce_features = config["dimensionality_reduction_sce"]["features"], 
            dimred_sce_npcs = config["dimensionality_reduction_sce"]["npcs"]
            ),
        expand("%s/rna/celltype_proportions/celltype_proportions_{stage}_horizontal_barplots_per_sample.pdf"  % (config["directories"]["results"]), 
            stage = config["stages"]),
        expand("%s/rna/MiloR/{stage}_Milo_features{features}_pcs{npcs}_tdTomcorr{tdTom_corr}.rds" % (config["directories"]["results"]), 
            stage=config["run_MiloR"]["stage"],
            features=config["run_MiloR"]["features"], 
            npcs=config["run_MiloR"]["npcs"],
            tdTom_corr=config["run_MiloR"]["tdTom_corr"]),
        expand("%s/rna/differential/DEresults_stage{stage}_tdTom_corr{tdTom_corr}_{method}.txt.gz" % (config["directories"]["results"]), 
            stage=config["find_DEGs"]["stage"], 
            tdTom_corr=config["run_MiloR"]["tdTom_corr"],
            method=config["find_DEGs"]["method"])
 
##################################################
## Load count matrices and create Seurat object ##
##################################################

rule create_seurat:
    input:
        script=config["scripts"]["create_seurat"],
        input_dir=config["directories"]["original_data"]
    output:
        seurat=config["directories"]["processed_data"]+"/seurat.rds",
        metadata=config["directories"]["processed_data"]+"/metadata.txt.gz"
    params:
        outdir=config["directories"]["processed_data"],
        sample=expand("{sample}", sample=config["samples"])
    log: 
        "logs/create_seurat.log"
    threads: 
        config["slurm"]["create_seurat"]["threads"]
    resources:
        mem_mb = config["slurm"]["create_seurat"]["memory"]
    shell:
        "Rscript {input.script} --inputdir {input.input_dir} --outdir {params.outdir} --samples {params.sample} > {log}"

#####################
## Quality control ##
#####################

rule qc:
    input:
        metadata=rules.create_seurat.output.metadata,
        script=config["scripts"]["qc"]
    output:
        # qc_metrics_boxplot=config["directories"]["results"]+"/rna/qc/qc_metrics_boxplot.pdf",
        qc_metrics_histogram=config["directories"]["results"]+"/rna/qc/qc_metrics_histogram.pdf",
        metadata=config["directories"]["results"]+"/rna/qc/sample_metadata_after_qc.txt.gz"
    params:
        sample = expand("{sample}", sample=config["samples"]),
        min_nCount_RNA = config["qc"]["min_nCount_RNA"],
        max_nCount_RNA = config["qc"]["max_nCount_RNA"],
        min_nFeature_RNA = config["qc"]["min_nFeature_RNA"],
        max_nFeature_RNA = config["qc"]["max_nFeature_RNA"],
        percent_mt = config["qc"]["percent_mt"],
        percent_rib = config["qc"]["percent_rib"],
        outdir=config["directories"]["results"]+"/rna/qc"
    log: 
        "logs/qc.log"
    threads: 
        config["slurm"]["qc"]["threads"]
    resources:
        mem_mb = config["slurm"]["qc"]["memory"]
    shell:
        "Rscript {input.script} --metadata {input.metadata} --outputdir {params.outdir} --samples {params.sample} --min_nCount_RNA {params.min_nCount_RNA} \
        --max_nCount_RNA {params.max_nCount_RNA} --min_nFeature_RNA {params.min_nFeature_RNA} --max_nFeature_RNA {params.max_nFeature_RNA} \
        --ribosomal_percent_RNA {params.percent_rib} --mitochondrial_percent_RNA {params.percent_mt} > {log}"


###################################################
## Convert Seurat object to SingleCellExperiment ##
###################################################

rule seurat_to_sce:
    input:
        seurat = rules.create_seurat.output.seurat,
    	metadata = rules.qc.output.metadata,
        script = config["scripts"]["seurat_to_sce"],
    output:
        config["directories"]["processed_data"]+"/SingleCellExperiment.rds",
    params:
        sample = expand("{sample}", sample=config["samples"])
    log: 
        "logs/seurat_to_sce.log"
    threads: 
        config["slurm"]["seurat_to_sce"]["threads"]
    resources:
        mem_mb = config["slurm"]["seurat_to_sce"]["memory"]
    shell:
        "Rscript {input.script} --samples {params.sample} --seurat {input.seurat} \
        --metadata {input.metadata} --outfile {output} > {log}"


#######################
## Doublet detection ##
#######################

rule doublet_detection:
    input:
        sce = rules.seurat_to_sce.output,
        metadata = rules.qc.output.metadata,
        script = config["scripts"]["doublet_detection"]
    output:
        outfile=config["directories"]["results"]+"/rna/doublet_detection/doublets_{sample}_{doublet_score_threshold}.txt.gz"
        #metadata=config["directories"]["results"]+"/rna/doublet_detection/sample_metadata_after_doublet_detection.txt.gz"
    log: 
        "logs/doublet_detection_{sample}_{doublet_score_threshold}.log"
    threads: 
        config["slurm"]["doublet_detection"]["threads"]
    resources:
        mem_mb = config["slurm"]["doublet_detection"]["memory"]
    shell:
        "Rscript {input.script} --metadata {input.metadata} --sce {input.sce} --samples {wildcards.sample} \
        --hybrid_score_threshold {wildcards.doublet_score_threshold}  --outfile {output} > {log}"

rule parse_doublet_results:
    input:
        metadata = rules.qc.output.metadata,
        script = config["scripts"]["parse_doublets"],
        doublet_files=expand(config["directories"]["results"]+"/rna/doublet_detection/doublets_{sample}_{doublet_score_threshold}.txt.gz", 
            sample=config["samples"], doublet_score_threshold=config["doublet_detection"]["doublet_score_threshold"])
        # doublet_files = rules.doublet_detection.output
    output:
        config["directories"]["results"]+"/rna/doublet_detection/sample_metadata_after_doublets.txt.gz"
    log: 
        "logs/parse_doublet_results.log"
    threads: 
        config["slurm"]["parse_doublet_results"]["threads"]
    resources:
        mem_mb = config["slurm"]["parse_doublet_results"]["memory"]
    shell:
        "Rscript {input.script} --metadata {input.metadata} --doublet_files {input.doublet_files} --outfile {output} > {log}"

##########################
## Mapping to the atlas ##
##########################

rule mapping_mnn:
    input:
        atlas_sce = config["directories"]["atlas"]+"/SingleCellExperiment.rds",
        atlas_metadata = config["directories"]["atlas"]+"/sample_metadata.txt.gz",
        query_sce = rules.seurat_to_sce.output,
        query_metadata = rules.parse_doublet_results.output,
        script = config["scripts"]["mapping_mnn"]
    output:
        config["directories"]["results"]+"/rna/mapping/mapping_mnn_{sample}.txt.gz"
    params:
        # sample = expand("{sample}", sample=config["samples"]),
        atlas_stages=config["mapping_mnn"]["atlas_stages"],
        npcs = config["mapping_mnn"]["npcs"],
        n_neighbours = config["mapping_mnn"]["n_neighbours"]
    log: 
        "logs/mapping_mnn_{sample}.log"
    threads: 
        config["slurm"]["mapping_mnn"]["threads"]
    resources:
        mem_mb = config["slurm"]["mapping_mnn"]["memory"]
    shell:
        "Rscript {input.script} --query_samples {wildcards.sample} --atlas_stages {params.atlas_stages} --query_sce {input.query_sce} \
        --atlas_sce {input.atlas_sce} --atlas_metadata {input.atlas_metadata} --query_metadata {input.query_metadata} \
        --npcs {params.npcs} --n_neighbours {params.n_neighbours} --outfile {output} > {log}"

rule mapping_mnn_all_samples:
    input:
        atlas_sce = config["directories"]["atlas"]+"/SingleCellExperiment.rds",
        atlas_metadata = config["directories"]["atlas"]+"/sample_metadata.txt.gz",
        query_sce = rules.seurat_to_sce.output,
        # query_metadata=rules.qc.output.metadata,
        query_metadata = rules.parse_doublet_results.output,
        script = config["scripts"]["mapping_mnn"]
        # outdir=config["directories"]["results"]+"/rna/mapping"
    output:
        config["directories"]["results"]+"/rna/mapping/mapping_mnn_all_samples.txt.gz"
    params:
        samples = expand("{sample}", sample=config["samples"]),
        atlas_stages=config["mapping_mnn"]["atlas_stages"],
        npcs = config["mapping_mnn"]["npcs"],
        n_neighbours = config["mapping_mnn"]["n_neighbours"]
    log: 
        "logs/mapping_mnn_all_samples.log"
    threads: 
        config["slurm"]["mapping_mnn"]["threads"]
    resources:
        mem_mb = config["slurm"]["mapping_mnn_all_samples"]["memory"]
    shell:
        "Rscript {input.script} --query_samples {params.samples} --atlas_stages {params.atlas_stages} --query_sce {input.query_sce} \
        --atlas_sce {input.atlas_sce} --atlas_metadata {input.atlas_metadata} --query_metadata {input.query_metadata} \
        --npcs {params.npcs} --n_neighbours {params.n_neighbours} --outfile {output} > {log}"

rule parse_mapping_results:
    input:
        query_metadata = rules.parse_doublet_results.output,
        sce = rules.seurat_to_sce.output,
        mapping_mnn = expand(config["directories"]["results"]+"/rna/mapping/mapping_mnn_{sample}.txt.gz", sample=config["samples"]),
        script = config["scripts"]["parse_mapping"]
    output:
        config["directories"]["results"]+"/rna/mapping/sample_metadata_after_mapping.txt.gz"
    log: 
        "logs/parse_mapping_results.log"
    threads: 
        config["slurm"]["parse_mapping_results"]["threads"]
    resources:
        mem_mb = config["slurm"]["parse_mapping_results"]["memory"]
    shell:
        "Rscript {input.script} --metadata {input.query_metadata} --sce {input.sce} --mapping_mnn {input.mapping_mnn} --outfile {output} > {log}"


##############################
## Dimensionality reduction ##
##############################

rule dimensionality_reduction_sce: 
    input:
        script=config["scripts"]["dimensionality_reduction_sce"],
        sce=rules.seurat_to_sce.output,
        metadata=rules.parse_mapping_results.output
    output:
        # config["directories"]["results"]+"/rna/dimensionality_reduction/umap_features{dimred_sce_features}_pcs{dimred_sce_npcs}_neigh{n_neighbors}_dist{min_dist}.txt.gz",
        config["directories"]["results"]+"/rna/dimensionality_reduction/sce/pca_features{dimred_sce_features}_pcs{dimred_sce_npcs}.txt.gz"
    params:
        outdir = config["directories"]["results"]+"/rna/dimensionality_reduction/sce",
        n_neighbors = config["dimensionality_reduction_sce"]["n_neighbors"],
        min_dist = config["dimensionality_reduction_sce"]["min_dist"],
        vars_to_regress = config["dimensionality_reduction_sce"]["vars_to_regress"],
        batch_correction = config["dimensionality_reduction_sce"]["batch_correction"],
        colour_by = config["dimensionality_reduction_sce"]["colour_by"],
        sample = expand("{sample}", sample=config["samples"]),
    log: 
        "logs/dimensionality_reduction_features{dimred_sce_features}_pcs{dimred_sce_npcs}.log"
        # "logs/dimensionality_reduction_features{dimred_sce_features}_pcs{dimred_sce_npcs}_neigh{n_neighbors}_dist{min_dist}.log"
    threads: 
        config["slurm"]["dimensionality_reduction_sce"]["threads"]
    resources:
        mem_mb = config["slurm"]["dimensionality_reduction_sce"]["memory"]
    shell:
        "Rscript {input.script} --sce {input.sce} --metadata {input.metadata} --npcs {wildcards.dimred_sce_npcs} --features {wildcards.dimred_sce_features} \
        --vars_to_regress {params.vars_to_regress} --samples {params.sample} --batch_correction  {params.batch_correction} \
        --n_neighbors {params.n_neighbors} --min_dist {params.min_dist} --colour_by {params.colour_by} --outdir {params.outdir} > {log}"

rule dimensionality_reduction_sce_remove_ExE_cells: 
    input:
        script=config["scripts"]["dimensionality_reduction_sce"],
        sce=rules.seurat_to_sce.output,
        metadata=rules.parse_mapping_results.output
    output:
        config["directories"]["results"]+"/rna/dimensionality_reduction/sce/no_ExE/pca_features{dimred_sce_features}_pcs{dimred_sce_npcs}.txt.gz"
    params:
        outdir = config["directories"]["results"]+"/rna/dimensionality_reduction/sce/no_ExE",
        n_neighbors = config["dimensionality_reduction_sce"]["n_neighbors"],
        min_dist = config["dimensionality_reduction_sce"]["min_dist"],
        vars_to_regress = config["dimensionality_reduction_sce"]["vars_to_regress"],
        batch_correction = config["dimensionality_reduction_sce"]["batch_correction"],
        colour_by = config["dimensionality_reduction_sce"]["colour_by"],
        sample = expand("{sample}", sample=config["samples"]),
    log: 
        "logs/dimensionality_reduction_features{dimred_sce_features}_pcs{dimred_sce_npcs}_noExE.log"
        # "logs/dimensionality_reduction_features{dimred_sce_features}_pcs{dimred_sce_npcs}_neigh{n_neighbors}_dist{min_dist}.log"
    threads: 
        config["slurm"]["dimensionality_reduction_sce"]["threads"]
    resources:
        mem_mb = config["slurm"]["dimensionality_reduction_sce"]["memory"]
    shell:
        "Rscript {input.script} --sce {input.sce} --metadata {input.metadata} --npcs {wildcards.dimred_sce_npcs} --features {wildcards.dimred_sce_features} \
        --vars_to_regress {params.vars_to_regress} --remove_ExE_cells --samples {params.sample} --batch_correction  {params.batch_correction}  \
        --n_neighbors {params.n_neighbors} --min_dist {params.min_dist} --colour_by {params.colour_by} --outdir {params.outdir} > {log}"
        #--vars_to_regress {params.vars_to_regress} 
        
##########################
## Plot mapping results ##
##########################

rule plot_mapping_results: 
    input:
        script = config["scripts"]["plot_mapping_results"],
        query_metadata=rules.parse_mapping_results.output,
        atlas_metadata = config["directories"]["atlas"]+"/sample_metadata.txt.gz"
    output:
        config["directories"]["results"]+"/rna/mapping/pdf/umap_mapped_allcells.pdf"
    params:
        samples = expand("{sample}", sample=config["samples"]),
        outdir = config["directories"]["results"]+"/rna/mapping/pdf"
    log: 
        "logs/plot_mapping_results.log"
    threads: 
        config["slurm"]["plot_mapping_results"]["threads"]
    resources:
        mem_mb = config["slurm"]["plot_mapping_results"]["memory"]        
    shell:
        "Rscript {input.script} --query_metadata {input.query_metadata} --atlas_metadata {input.atlas_metadata} \
        --samples {params.samples} --outdir {params.outdir} > {log}"
        
################################
## Plot cell type proportions ##
################################

rule plot_celltype_proportions: 
    input:
        script = config["scripts"]["plot_celltype_proportions"],
        metadata=rules.parse_mapping_results.output
    output:
        expand(config["directories"]["results"]+"/rna/celltype_proportions/celltype_proportions_{stage}_horizontal_barplots_per_sample.pdf", stage=config["stages"])
    params:
        samples = expand("{sample}", sample=config["samples"]),
        celltype_label = config["plot_celltype_proportions"]["celltype_label"],
        outdir = config["directories"]["results"]+"/rna/celltype_proportions"
#    conda:
#        "environment.yaml"
    log: 
        "logs/plot_celltype_proportions.log"
    threads: 
        config["slurm"]["plot_celltype_proportions"]["threads"]
    resources:
        mem_mb = config["slurm"]["plot_celltype_proportions"]["memory"]        
    shell:
        "Rscript {input.script} --metadata {input.metadata} --celltype_label {params.celltype_label} \
        --samples {params.samples} --outdir {params.outdir} > {log}"

###################
## Perform MiloR ##
###################

rule run_MiloR:
    input:
        script = config["scripts"]["run_MILO"],
        metadata = rules.parse_mapping_results.output,
        sce = rules.seurat_to_sce.output,
        atlas_metadata = config["directories"]["atlas"]+"/sample_metadata.txt.gz",
        dimred = expand("%s/rna/dimensionality_reduction/sce/pca_features{dimred_sce_features}_pcs{dimred_sce_npcs}.txt.gz" % (config["directories"]["results"]), 
            dimred_sce_features = config["dimensionality_reduction_sce"]["features"], 
            dimred_sce_npcs = config["dimensionality_reduction_sce"]["npcs"]
            )
    output:
        #expand(config["directories"]["results"]+"/rna/MiloR/{stage}_Milo_features{features}_pcs{npcs}.rds", 
        #    stage=config["run_MiloR"]["stage"],
        #    features = config["run_MiloR"]["features"],
        #    npcs = config["run_MiloR"]["npcs"])
        config["directories"]["results"]+"/rna/MiloR/{stage}_Milo_features{features}_pcs{npcs}_tdTomcorr{tdTom_corr}.rds"
    params:
        outdir = config["directories"]["results"]+"/rna/MiloR",
        n_neighbors = config["run_MiloR"]["n_neighbors"],
        prop = config["run_MiloR"]["prop"],
        tdTom_corr = config["run_MiloR"]["tdTom_corr"]
    log: 
        "logs/run_MiloR_{stage}_{features}_{npcs}_{tdTom_corr}.log"
    threads: 
        config["slurm"]["run_MiloR"]["threads"]
    resources:
        mem_mb = config["slurm"]["run_MiloR"]["memory"]        
    conda:
        "env_rna"
    shell:
        "Rscript {input.script} --metadata {input.metadata} --sce {input.sce} --atlas_metadata {input.atlas_metadata} --batch_correction --regression \
        --n_neighbors {params.n_neighbors} --npcs {wildcards.npcs} --features {wildcards.features} --stage {wildcards.stage}  --prop {params.prop} --tdTom_corr {wildcards.tdTom_corr} \
         --outdir {params.outdir} > {log}"  

###############
## find DEGs ##
###############
rule find_DEGs:
    input:
        script = config["scripts"]["find_DEGs"],
        metadata = rules.parse_mapping_results.output,
        sce = rules.seurat_to_sce.output
    output:
        #expand(config["directories"]["results"]+"/rna/MiloR/{stage}_Milo_features{features}_pcs{npcs}.rds", 
        #    stage=config["run_MiloR"]["stage"],
        #    features = config["run_MiloR"]["features"],
        #    npcs = config["run_MiloR"]["npcs"])
        config["directories"]["results"]+"/rna/differential/DEresults_stage{stage}_tdTom_corr{tdTom_corr}_{method}.txt.gz"
    params:
        outdir = config["directories"]["results"]+"/rna/differential"
    log: 
        "logs/findDEGs_{stage}_{tdTom_corr}_{method}.log"
    threads: 
        config["slurm"]["find_DEGs"]["threads"]
    resources:
        mem_mb = config["slurm"]["find_DEGs"]["memory"]        
    conda:
        "env_rna"
    shell:
        "Rscript {input.script} --metadata {input.metadata} --sce {input.sce} \
        --stage {wildcards.stage} --tdTom_corr {wildcards.tdTom_corr} --method {wildcards.method} \
         --outdir {params.outdir} > {log}"  