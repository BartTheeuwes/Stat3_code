suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))

#########
## I/O ##
#########

io <- list()
io$basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA"
io$atlas.basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/"
io$gene_metadata <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/Mmusculus_genes_BioMart.87.txt"


io$metadata <- file.path(io$basedir,"processed/metadata.txt.gz")

# TFs
io$TFs <- file.path(io$basedir,"results_new/TFs.txt")

# RNA
io$rna.anndata <- file.path(io$basedir,"processed/anndata.h5ad")
io$rna.seurat <- file.path(io$basedir,"processed/seurat.rds")
io$rna.sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$rna.differential <- file.path(io$basedir,"results/rna/differential")
io$rna.pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment.rds")

# RNA atlas (Pijuan-Sala2019)
io$rna.atlas.metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
io$rna.atlas.marker_genes.up <- file.path(io$atlas.basedir,"results/marker_genes/all_stages/marker_genes.txt.gz")
# io$rna.atlas.marker_genes.all <- file.path(io$atlas.basedir,"results/marker_genes/all_stages/marker_genes_all.txt.gz")
io$rna.atlas.marker_TFs.up <- file.path(io$atlas.basedir,"results/differential/celltypes/TFs/TF_markers/marker_TFs_up.txt.gz")
# io$rna.atlas.marker_TFs.all <- file.path(io$atlas.basedir,"results/differential/celltypes/TFs/TF_markers/marker_TFs_all.txt.gz")
io$rna.atlas.differential <- file.path(io$atlas.basedir,"results/differential")
# io$rna.atlas.average_expression_per_celltype <- file.path(io$atlas.basedir,"results/marker_genes/all_stages/avg_expr_per_celltype_and_gene.txt.gz")
io$rna.atlas.sce.pseudobulk <- file.path(io$atlas.basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk.rds")
io$rna.atlas.sce <- file.path(io$atlas.basedir,"final/SingleCellExperiment.rds")
io$rna.atlas.celltype_proportions <- file.path(io$atlas.basedir,"results/celltype_proportions/celltype_proportions.txt.gz")

#############
## Options ##
#############

opts <- list()

opts$stages <- c(
  "E7.5",
  "E8.5",
  "E9.5"
)

opts$samples <- c(
'SLX-21143_SITTA2_HTJH3DSX2',
'SLX-21143_SITTA4_HTJH3DSX2',
'SLX-21143_SITTB2_HTJH3DSX2',
'SLX-21143_SITTB4_HTJH3DSX2',
'SLX-21143_SITTC4_HTJH3DSX2',
'SLX-21143_SITTD4_HTJH3DSX2',
'SLX-21143_SITTD5_HTJH3DSX2',
'SLX-21143_SITTE3_HTJH3DSX2',
'SLX-21143_SITTE5_HTJH3DSX2',
'SLX-21143_SITTF3_HTJH3DSX2',
'SLX-21143_SITTF5_HTJH3DSX2',
'SLX-21143_SITTG1_HTJH3DSX2',
'SLX-21143_SITTG3_HTJH3DSX2',
'SLX-21143_SITTG5_HTJH3DSX2',
'SLX-21143_SITTH1_HTJH3DSX2',
'SLX-21143_SITTH3_HTJH3DSX2'
)

opts$rename.samples <- c(
'SLX-21143_SITTA2_HTJH3DSX2' = 'Sample_5',
'SLX-21143_SITTA4_HTJH3DSX2' = 'Sample_9',
'SLX-21143_SITTB2_HTJH3DSX2' = 'Sample_6',
'SLX-21143_SITTB4_HTJH3DSX2' = 'Sample_10',
'SLX-21143_SITTC4_HTJH3DSX2' = 'Sample_11',
'SLX-21143_SITTD4_HTJH3DSX2' = 'Sample_12',
'SLX-21143_SITTD5_HTJH3DSX2' = 'Sample_1',
'SLX-21143_SITTE3_HTJH3DSX2' = 'Sample_13',
'SLX-21143_SITTE5_HTJH3DSX2' = 'Sample_2',
'SLX-21143_SITTF3_HTJH3DSX2' = 'Sample_14',
'SLX-21143_SITTF5_HTJH3DSX2' = 'Sample_3',
'SLX-21143_SITTG1_HTJH3DSX2' = 'Sample_7',
'SLX-21143_SITTG3_HTJH3DSX2' = 'Sample_15',
'SLX-21143_SITTG5_HTJH3DSX2' = 'Sample_4',
'SLX-21143_SITTH1_HTJH3DSX2' = 'Sample_8',
'SLX-21143_SITTH3_HTJH3DSX2' = 'Sample_16'
)

opts$sample2stage <- c(
'SLX-21143_SITTA2_HTJH3DSX2' = 'E8.5',
'SLX-21143_SITTA4_HTJH3DSX2' = 'E9.5',
'SLX-21143_SITTB2_HTJH3DSX2' = 'E8.5',
'SLX-21143_SITTB4_HTJH3DSX2' = 'E9.5',
'SLX-21143_SITTC4_HTJH3DSX2' = 'E9.5',
'SLX-21143_SITTD4_HTJH3DSX2' = 'E9.5',
'SLX-21143_SITTD5_HTJH3DSX2' = 'E7.5',
'SLX-21143_SITTE3_HTJH3DSX2' = 'E9.5',
'SLX-21143_SITTE5_HTJH3DSX2' = 'E7.5',
'SLX-21143_SITTF3_HTJH3DSX2' = 'E9.5',
'SLX-21143_SITTF5_HTJH3DSX2' = 'E7.5',
'SLX-21143_SITTG1_HTJH3DSX2' = 'E8.5',
'SLX-21143_SITTG3_HTJH3DSX2' = 'E9.5',
'SLX-21143_SITTG5_HTJH3DSX2' = 'E7.5',
'SLX-21143_SITTH1_HTJH3DSX2' = 'E8.5',
'SLX-21143_SITTH3_HTJH3DSX2' = 'E9.5'
)


opts$sample2tomato = c(
'SLX-21143_SITTA2_HTJH3DSX2' = TRUE,
'SLX-21143_SITTA4_HTJH3DSX2' = TRUE,
'SLX-21143_SITTB2_HTJH3DSX2' = FALSE,
'SLX-21143_SITTB4_HTJH3DSX2' = FALSE,
'SLX-21143_SITTC4_HTJH3DSX2' = TRUE,
'SLX-21143_SITTD4_HTJH3DSX2' = FALSE,
'SLX-21143_SITTD5_HTJH3DSX2' = TRUE,
'SLX-21143_SITTE3_HTJH3DSX2' = TRUE,
'SLX-21143_SITTE5_HTJH3DSX2' = FALSE,
'SLX-21143_SITTF3_HTJH3DSX2' = FALSE,
'SLX-21143_SITTF5_HTJH3DSX2' = TRUE,
'SLX-21143_SITTG1_HTJH3DSX2' = TRUE,
'SLX-21143_SITTG3_HTJH3DSX2' = TRUE,
'SLX-21143_SITTG5_HTJH3DSX2' = FALSE,
'SLX-21143_SITTH1_HTJH3DSX2' = FALSE,
'SLX-21143_SITTH3_HTJH3DSX2' = FALSE
)


opts$sample2pool <- c(
'SLX-21143_SITTA2_HTJH3DSX2' = 'E8.5_rep1',
'SLX-21143_SITTA4_HTJH3DSX2' = 'E9.5_rep3',
'SLX-21143_SITTB2_HTJH3DSX2' = 'E8.5_rep1',
'SLX-21143_SITTB4_HTJH3DSX2' = 'E9.5_rep3',
'SLX-21143_SITTC4_HTJH3DSX2' = 'E9.5_rep4',
'SLX-21143_SITTD4_HTJH3DSX2' = 'E9.5_rep4',
'SLX-21143_SITTD5_HTJH3DSX2' = 'E7.5_rep1',
'SLX-21143_SITTE3_HTJH3DSX2' = 'E9.5_rep1',
'SLX-21143_SITTE5_HTJH3DSX2' = 'E7.5_rep1',
'SLX-21143_SITTF3_HTJH3DSX2' = 'E9.5_rep1',
'SLX-21143_SITTF5_HTJH3DSX2' = 'E7.5_rep2',
'SLX-21143_SITTG1_HTJH3DSX2' = 'E8.5_rep2',
'SLX-21143_SITTG3_HTJH3DSX2' = 'E9.5_rep2',
'SLX-21143_SITTG5_HTJH3DSX2' = 'E7.5_rep2',
'SLX-21143_SITTH1_HTJH3DSX2' = 'E8.5_rep2',
'SLX-21143_SITTH3_HTJH3DSX2' = 'E9.5_rep2'
)


opts$celltypes <- c(
'Epiblast', 
'Primitive Streak', 
'ExE ectoderm', 
'Visceral endoderm', 
'ExE endoderm', 
'Non-neural ectoderm', 
'Nascent mesoderm', 
'Parietal endoderm', 
'Ectoderm', 
'Anterior Primitive Streak', 
'Haematoendothelial progenitors', 
'Caudal epiblast', 
'Blood progenitors', 
'Intermediate mesoderm', 
'Paraxial mesoderm', 
'Lateral plate mesoderm', 
'Mesenchyme', 
'PGC', 
'Node', 
'Gut tube', 
'Embryo proper endothelium', 
'Cardiopharyngeal progenitors SHF', 
'Notochord', 
'Amniotic ectoderm', 
'Venous endothelium', 
'Presomitic mesoderm', 
'Cardiomyocytes FHF 1', 
'Allantois', 
'Cranial mesoderm', 
'EMP', 
'Limb mesoderm', 
'Anterior somitic tissues', 
'Pharyngeal mesoderm', 
'Allantois endothelium', 
'Thyroid primordium', 
'Erythroid', 
'Hindbrain neural progenitors', 
'Cardiomyocytes SHF 1', 
'NMPs', 
'Pharyngeal endoderm', 
'Dorsal spinal cord progenitors', 
'Anterior cardiopharyngeal progenitors', 
'Placodal ectoderm', 
'Optic vesicle', 
'Ventral forebrain progenitors', 
'Spinal cord progenitors', 
'Hindgut', 
'Caudal mesoderm', 
'Embryo proper mesothelium', 
'Neural tube', 
'Midbrain/Hindbrain boundary', 
'Posterior somitic tissues', 
'Midgut', 
'Migratory neural crest', 
'Ventral hindbrain progenitors', 
'Surface ectoderm', 
'YS mesothelium', 
'Limb ectoderm', 
'Somitic mesoderm', 
'NMPs/Mesoderm-biased', 
'Cardiomyocytes FHF 2', 
'Foregut', 
'Dermomyotome', 
'Kidney primordium', 
'Otic placode', 
'Cardiomyocytes SHF 2', 
'Midbrain progenitors', 
'Cardiopharyngeal progenitors FHF', 
'Epicardium', 
'Hindbrain floor plate', 
'Late dorsal forebrain progenitors', 
'Dorsal hindbrain progenitors', 
'Sclerotome', 
'YS endothelium', 
'Endocardium', 
'MEP', 
'Epidermis', 
'Megakaryocyte progenitors', 
'Early dorsal forebrain progenitors', 
'Cardiopharyngeal progenitors', 
'Chorioallantoic-derived erythroid progenitors', 
'YS mesothelium-derived endothelial progenitors', 
'Dorsal midbrain neurons', 
'Branchial arch neural crest', 
'Forelimb', 
'Frontonasal mesenchyme', 
'Otic neural progenitors'
)

opts$stage.colors = c(
                  "mixed_gastrulation" = "#A9A9A9",
                  "E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD",
                    "E8.75" = '#3C1ACE',
                    "E9.0" = '#9A28F7' ,
                    "E9.25" =  '#E228F7',
                    "E9.5" = '#F728B8' 
 )
opts$stage.colors <- viridis::viridis(n=length(opts$stages))
names(opts$stage.colors) <- rev(opts$stages)

opts$celltype.colors = c(
 "Epiblast" = "#635547",
"Primitive Streak" = "#DABE99",
"Caudal epiblast" = "#9e6762",

"PGC" = "#FACB12",

"Anterior Primitive Streak" = "#c19f70",
"Node"="#153b3d",
"Notochord" = "#0F4A9C",



"Gut tube" = "#EF5A9D",
"Hindgut" = "#F397C0",
"Midgut" = "#ff00b2",
"Foregut" = "#ffb7ff",
"Pharyngeal endoderm"="#95e1ff",
"Thyroid primordium"="#97bad3",

"Nascent mesoderm" = "#C594BF",
"Intermediate mesoderm" = "#139992",
"Caudal mesoderm" = "#3F84AA",
"Lateral plate mesoderm" = "#F9DFE6",
"Limb mesoderm" = "#e35f82",
"Forelimb" = "#d02d75",
"Kidney primordium" = "#e85639",
"Presomitic mesoderm"="#5581ca",#"#0000ff",#blue
"Somitic mesoderm" = "#005579",
"Posterior somitic tissues" = "#5adbe4",#"#40e0d0",#turquoise



"Paraxial mesoderm" = "#8DB5CE",
"Cranial mesoderm" = "#456722",#"#006400",#darkgreen
"Anterior somitic tissues"= "#d5e839",
"Sclerotome" = "#e3cb3a",#"#ffff00",#yellow
"Dermomyotome" = "#00BFC4",#"#a52a2a",#brown



"Pharyngeal mesoderm" = "#C9EBFB",
"Cardiopharyngeal progenitors" = "#556789",
"Anterior cardiopharyngeal progenitors"="#683ed8",



"Allantois" = "#532C8A",
"Mesenchyme" = "#cc7818",
"YS mesothelium" = "#ff7f9c",
"Epicardium"="#f79083",
"Embryo proper mesothelium" = "#ff487d",



"Cardiopharyngeal progenitors FHF"="#d780b0",
"Cardiomyocytes FHF 1"="#a64d7e",
"Cardiomyocytes FHF 2"="#B51D8D",



"Cardiopharyngeal progenitors SHF"="#4b7193",
"Cardiomyocytes SHF 1"="#5d70dc",
"Cardiomyocytes SHF 2"="#332c6c",



"Haematoendothelial progenitors" = "#FBBE92",
"Blood progenitors" = "#6c4b4c",
"Erythroid" = "#C72228",
"Chorioallantoic-derived erythroid progenitors"="#E50000",
"Megakaryocyte progenitors"="#e3cb3a",
"MEP"="#EF4E22",
"EMP"="#7c2a47",



"YS endothelium"="#ff891c",
"YS mesothelium-derived endothelial progenitors"="#AE3F3F",
"Allantois endothelium"="#2f4a60",
"Embryo proper endothelium"="#90e3bf",
"Venous endothelium"="#bd3400",
"Endocardium"="#9d0049",

"NMPs/Mesoderm-biased" = "#89c1f5",
"NMPs" = "#8EC792",

"Ectoderm" = "#ff675c",



"Optic vesicle" = "#bd7300",

"Ventral forebrain progenitors"="#a0b689",
"Early dorsal forebrain progenitors"="#0f8073",
"Late dorsal forebrain progenitors"="#7a9941",
"Midbrain/Hindbrain boundary"="#8ab3b5",
"Midbrain progenitors"="#9bf981",
"Dorsal midbrain neurons"="#12ed4c",
"Ventral hindbrain progenitors"="#7e907a",
"Dorsal hindbrain progenitors"="#2c6521",
"Hindbrain floor plate"="#bf9da8",
"Hindbrain neural progenitors"="#59b545",



"Neural tube"="#233629",



"Migratory neural crest"="#4a6798",
"Branchial arch neural crest"="#bd84b0",
"Frontonasal mesenchyme"="#d3b1b1",



"Spinal cord progenitors"="#6b2035",
"Dorsal spinal cord progenitors"="#e273d6",

"Non-neural ectoderm" = "#f7f79e",
"Surface ectoderm" = "#fcff00",
"Epidermis" = "#fff335",
"Limb ectoderm" = "#ffd731",
"Amniotic ectoderm" = "#dbb400",



"Placodal ectoderm" = "#ff5c00",



"Otic placode"="#f1a262",
"Otic neural progenitors"="#00b000",

"Visceral endoderm" = "#F6BFCB",
"ExE endoderm" = "#7F6874",
"ExE ectoderm" = "#989898",
"Parietal endoderm" = "#1A1A1A"
)

opts$tdTom.color = c(
    'TRUE' = "red", 
    'FALSE' = "black"
)

opts$chr <- paste0("chr",c(1:19,"X","Y"))
