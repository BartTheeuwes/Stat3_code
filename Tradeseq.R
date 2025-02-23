## Load packages
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(BiocParallel))

here::i_am("blood/Tradeseq.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers = 32
set.seed(6)

## I/O
args = list()
args$outdir <- paste0(io$basedir,"/results/rna/blood")
sce_filt = readRDS(file.path(args$outdir, 'sce_slingshot_filt.rds'))

# Set all cells to belong to one trajectory
cellWeights <- rep(1,ncol(sce_filt))

# create a model matrix -> Comparing within Embryo pools
U <- model.matrix(~sce_filt$pool)
print('Run Trade-seq')
sceGAM <- fitGAM(counts = counts(sce_filt),
                 conditions = as.factor(sce_filt$tdTom_corr),
                 U = U,
                 pseudotime= sce_filt$pseudotime,
                 cellWeights=cellWeights,
                 nknots=6, # put at 6!
                 parallel = TRUE,
                 BPPARAM=BPPARAM,
                 verbose=T)

print('Saving Trade-seq results')
saveRDS(sceGAM, file.path(args$outdir,"sceGAM.rds"))