suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(argparse))

here::i_am("mapping/run/mnn/mapping_mnn.R")

# Load mapping functions
source(here::here("mapping/run/mnn/mapping_functions.R"))

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--query_sce',       type="character",               help='SingleCellExperiment file for the query')
p$add_argument('--atlas_sce',       type="character",               help='SingleCellExperiment file for the atlas')
p$add_argument('--query_metadata',  type="character",               help='metadata file for the query')
p$add_argument('--atlas_metadata',  type="character",               help='metadata file for the atlas')
p$add_argument('--npcs',            type="integer",                 help='Number of principal components')
p$add_argument('--n_neighbours',    type="integer",                 help='Number of neighbours')
p$add_argument('--use_marker_genes',action = "store_true",          help='Use marker genes?')
p$add_argument('--cosine_normalisation',      action = "store_true",          help='Use cosine normalisation?')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# I/O
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir

## START TEST ##
# args$atlas_stages <- c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")
# args$query_samples <- opts$samples
# args$query_sce <- io$rna.sce
# args$query_sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
# args$atlas_sce <- io$rna.atlas.sce
# args$query_metadata <- paste0(io$basedir,"/results/rna/qc/sample_metadata_after_qc.txt.gz")
# args$atlas_metadata <- io$rna.atlas.metadata
# args$test <- FALSE
# args$npcs <- 50
# args$n_neighbours <- 25
# args$use_marker_genes <- FALSE
# args$cosine_normalisation <- FALSE
# args$outdir <- paste0(io$basedir,"/results/rna/mapping/test")
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(args$query_metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$query_samples]
if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)

# Load SingleCellExperiment
sce_query <- load_SingleCellExperiment(args$query_sce, cells = meta_query$cell, remove_non_expressed_genes = TRUE)

# Update colData
tmp <- meta_query %>% .[cell%in%colnames(sce_query)] %>% setkey(cell) %>% .[colnames(sce_query)]
stopifnot(tmp$cell == colnames(sce_query))
colData(sce_query) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_query),] %>% DataFrame()

################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  .[stage%in%args$atlas_stages] %>%
  .[,sample:=factor(sample)] %>%
  .[,celltype:=celltype_extended_atlas]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)

# Load SingleCellExperiment
sce_atlas <- load_SingleCellExperiment(args$atlas_sce, normalise = TRUE, cells = meta_atlas$cell, remove_non_expressed_genes = TRUE)

# Update colData
tmp <- meta_atlas %>% .[cell%in%colnames(sce_atlas)] %>% setkey(cell) %>% .[colnames(sce_atlas)]
stopifnot(tmp$cell == colnames(sce_atlas))
colData(sce_atlas) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_atlas),] %>% DataFrame()

# Sanity cehcks
stopifnot(sum(is.na(rownames(sce_atlas)))==0)
stopifnot(sum(duplicated(rownames(sce_atlas)))==0)

#####################
## Define gene set ##
#####################

# Get gene metadata
gene_metadata <- fread(io$gene_metadata) %>% .[,c("chr","ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce_atlas)] %>%
  .[!duplicated(symbol)]

# Imprinted genes
imprint = gene_metadata[c(grep('maternally', gene_metadata$description),
                       grep('paternally', gene_metadata$description)), symbol]
#Other imprinted genes: 
#- Nnat (https://www.genecards.org/cgi-bin/carddisp.pl?gene=NNAT)
#- Grb10 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=GRB10)

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))

# Filter some genes manually
genes.intersect <- genes.intersect[grep("^Rik|Rik$|^mt-|^Rps|^Rpl|^Gm",genes.intersect,invert=T)] # filter out non-informative genes
genes.intersect <- genes.intersect[grep("^Hbb|^Hba",genes.intersect,invert=T)] # test removing Haem genes 
genes.intersect <- genes.intersect[!genes.intersect %in% c(imprint, 'Grb10', 'Nnat')] # remove imprinted genes
genes.intersect <- genes.intersect[!genes.intersect %in% c("Xist", "Tsix")] # remove Xist & Tsix
genes.intersect <- genes.intersect[!genes.intersect=="tomato-td"] # remove tomato itself
genes.intersect <- genes.intersect[!genes.intersect %in% gene_metadata[chr=="chrY",symbol]] # no genes on y-chr 

# Subset SingleCellExperiment objects
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

#######################
## Feature selection ##
#######################

if (args$use_marker_genes) {
  # Load marker genes
  marker_genes.dt <- fread(io$rna.atlas.marker_genes)
  genes_to_use <- genes.intersect[genes.intersect %in% unique(marker_genes.dt$gene)]
} else {
  # # Load gene statistics from the atlas
  # gene_stats.dt <- fread(paste0(io$atlas.basedir,"/results/gene_statistics/gene_statistics.txt.gz")) %>%
  #   .[gene%in%genes.intersect]
  # genes_to_use <- gene_stats.dt %>% setorder(-var_pseudobulk, na.last = T) %>% head(n=2500) %>% .$gene  
  
  # Calculate mean-variance relationship and extract HVGs
  decomp <- modelGeneVar(sce_atlas, block=sce_atlas$sample)
 # genes_to_use <- rownames(decomp)[decomp$p.value<=0.01 & decomp$mean>0.1]
  genes_to_use <- decomp[order(decomp$FDR),] %>% head(n=2500) %>% rownames
}

stopifnot(genes_to_use%in%rownames(sce_atlas))
stopifnot(genes_to_use%in%rownames(sce_query))

sce_query = sce_query[genes_to_use,]
sce_atlas = sce_atlas[genes_to_use,]

#########
## Map ##
#########

# TO-DO: TRY COSINE NORMALISATION
mapping  <- mapWrap(
  sce_atlas = sce_atlas,
  meta_atlas = meta_atlas,
  sce_query = sce_query,
  meta_query = meta_query,
  genes = genes_to_use,
  npcs = args$npcs,
  k = args$n_neighbours,
  cosineNorm = args$cosine_normalisation,
  order = NULL
)


##########
## Save ##
##########

mapping.dt <- mapping$mapping %>% 
  .[,c("cell","celltype.mapped","celltype.score","stage.mapped", "cellstage.score", "closest.cell")] %>% 
  as.data.table

fwrite(mapping.dt, args$outfile, sep="\t")
