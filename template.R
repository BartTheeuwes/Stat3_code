
####################
## Load libraries ##
####################

#####################
## Define settings ##
#####################

source("/rds/project/rds-SDzz0CATGms/users/bt392/multiome_practice/settings.R")
source("/rds/project/rds-SDzz0CATGms/users/bt392/multiome_practice/utils.R")


# Define I/O
io$outdir <- paste0(io$basedir,"/results/results/differential")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

#####################
## Load metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples]


###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, normalise = TRUE, cells = sample_metadata$cell)
