suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(miloR))
suppressPackageStartupMessages(library(patchwork))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--atlas_metadata',  type="character",                               help='Atlas metadata file')
p$add_argument('--stage',           type="character",                               help='Stages to include')
p$add_argument('--features',        type="integer",    default=2500,                help='Number of features')
p$add_argument('--regression',      action="store_true",                            help='Regress out variable')
p$add_argument('--batch_correction',      action="store_true",                            help='Batch correction by sample')
p$add_argument('--npcs',            type="integer",    default=50,                  help='Number of PCs')
p$add_argument('--n_neighbors',     type="integer",    default=30,                  help='Number of neighbours')
p$add_argument('--tdTom_corr',      type="character",                               help='Keep tdTom+ cells in WT samples (otherwise removed)')
p$add_argument('--prop',            type="double",     default=0.3,                 help='Proportion of cells to sample for MILO')
#p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################
here::i_am("processing/1_create_seurat_rna.R")
source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("mapping/run/mnn/mapping_functions.R"))
test = FALSE

if(test){
## START TEST ##
    args = list()
args$sce <- io$rna.sce
args$metadata <- io$metadata
args$stage <- c('E7.5')#, 'E8.5', 'E9.5')
args$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$atlas_metadata <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/sample_metadata.txt.gz"
args$features <- 2500
args$npcs <- 50 
args$n_neighbors = 45 
args$prop = 0.15 # put at 0.15 - 0.2
args$remove_ExE_cells <- FALSE
args$tdTom_corr = NULL
args$outdir <- paste0(io$basedir,"/results/rna/MiloR/test")
## END TEST ##
}

# If passing multiple timepoints split in vector
args$stage = strsplit(args$stage, "_")[[1]] 

dir.create(args$outdir, recursive=TRUE, showWarnings = FALSE)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
   .[pass_rnaQC==TRUE & doublet_call==FALSE & stage %in% args$stage] %>%
    .[,pool:=stringr::str_replace_all(sample,opts$sample2pool)]


if(args$tdTom_corr!='True'){
    sample_metadata = sample_metadata[tdTom==tdTom_corr]
}

if(test){
    sample_metadata = sample_metadata[sample(1:nrow(sample_metadata), nrow(sample_metadata)/10)]
}

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame


###############
## Get PCA   ##
###############

# Get Reduced Dims
# PCA
if(length(args$stage)>1){
    message('loading precomputed PCA')
    pca_file = sprintf("%s/results/rna/dimensionality_reduction/sce/pca_features%d_pcs%d.txt.gz", io$basedir, args$features, args$npcs)
    pca = fread(pca_file)[match(colnames(sce), cell)] %>% tibble::column_to_rownames("cell") %>% as.matrix
} else {
    ## Find HVGs - detection on WT samples only
    # Get gene metadata
    gene_metadata <- fread(io$gene_metadata) %>% .[,c("chr","ens_id","symbol")] %>%
      .[symbol!="" & ens_id%in%rownames(sce)] %>%
      .[!duplicated(symbol)]

    # Imprinted genes
    imprint = gene_metadata[c(grep('maternally', gene_metadata$description),
                           grep('paternally', gene_metadata$description)), symbol]
    #Other imprinted genes: 
    #- Nnat (https://www.genecards.org/cgi-bin/carddisp.pl?gene=NNAT)
    #- Grb10 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=GRB10)

    decomp <- modelGeneVar(sce[,sample_metadata[tdTom==FALSE, cell]], block=colData(sce[,sample_metadata[tdTom==FALSE, cell]])$sample) # Only detect HVGs from WT samples
    decomp <- decomp[decomp$mean > 0.01,]
    hvgs <- decomp[order(decomp$FDR),] %>% 
        as.data.table(., keep.rownames=T) %>% 
        .[grep("^Rik|Rik$|^mt-|^Rps|^Rpl|^Gm",rn,invert=T)] %>% # filter out non-informative genes
        .[grep("^Hbb|^Hba",rn,invert=T)] %>% # test removing Haem genes 
        .[!rn %in% c(imprint, 'Grb10', 'Nnat')] %>%  # remove imprinted genes
        .[!rn %in% c("Xist", "Tsix")] %>%  # remove Xist & Tsix
        .[!rn=="tomato-td"] %>% # remove tomato itself
        .[!rn%in%gene_metadata[chr=="chrY",symbol]] %>%
         head(n=args$features) %>% .$rn

    sce_filt <- sce[hvgs,]
    
    # Regress Variables
    if(args$regression){
       args$vars_to_regress = c('nFeature_RNA', 'nCount_RNA') # , 'mitochondrial_percent_RNA', 'ribosomal_percent_RNA'
       print(sprintf("Regressing out variables: %s", paste(args$vars_to_regress,collapse=" ")))
       logcounts(sce_filt) <- RegressOutMatrix(
         mtx = logcounts(sce_filt), 
         covariates = colData(sce_filt)[,args$vars_to_regress,drop=F]
        )
    }
    
    # Run PCA
    if (args$batch_correction==TRUE) {
        library(batchelor)
        print('batch correcting by sample')
        pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)$sample, d = args$npcs)
        pca.corrected <- reducedMNN(pca)$corrected
        colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
        pca.corrected = pca.corrected[match(colnames(sce_filt), rownames(pca.corrected)),]
        reducedDim(sce_filt, "PCA") <- pca.corrected
    } else {
        sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)  
    }
    
    # Save PCA coordinates
    outfile <- sprintf("%s/%s_pca_features%d_pcs%d.txt.gz",args$outdir, paste(args$stage, collapse ='_'), args$features, args$npcs)
    pca.dt <- reducedDim(sce_filt,"PCA") %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
    fwrite(pca.dt, outfile)
}

# Umap -> doesn't currently give option to chose NN & dist of final umap
umap_file = sprintf("%s/results/rna/dimensionality_reduction/sce/umap_features%d_pcs%d_neigh30_dist0.3.txt.gz", io$basedir, args$features, args$npcs)
umap = fread(umap_file)[match(colnames(sce), cell)] %>% tibble::column_to_rownames("cell") %>% as.matrix

reducedDim(sce_filt,"UMAP") = umap

########################
## Create MILO object ##
########################

# create object
Milo_sce <- Milo(sce_filt)

colData(Milo_sce)$sample_ko = paste0(colData(Milo_sce)$sample, '_',  colData(Milo_sce)$tdTom_corr)

# Build graph
Milo_sce <- buildGraph(Milo_sce, 
                       k = args$n_neighbors, 
                       d = args$npcs, 
                       reduced.dim = "PCA")

# Identify neighbourhoods from NN-cells
# lower prop for larger datasets
Milo_sce <- makeNhoods(Milo_sce, 
                       prop = args$prop, 
                       k = args$n_neighbors, 
                       d = args$npcs, 
                       refined = TRUE, 
                       reduced_dims = "PCA")


# count cells of different samples per neighbourhood
Milo_sce <- countCells(Milo_sce, meta.data = as.data.frame(colData(Milo_sce)), sample="sample_ko")

# Calculate Neighbourhood Distances (most time consuming step)
Milo_sce <- calcNhoodDistance(Milo_sce, d=args$npcs, reduced.dim = "PCA")

# Build neighbourhood graph
Milo_sce <- buildNhoodGraph(Milo_sce)

# Save MILO object
outfile = sprintf("%s/%s_Milo_features%d_pcs%d_tdTomcorr%s.rds",args$outdir, paste(args$stage, collapse ='_'), args$features, args$npcs, args$tdTom_corr)
saveRDS(Milo_sce, outfile)

# Embryo design table
embryo_design <- data.frame(colData(Milo_sce))[,c("sample_ko", "sample", 'tdTom_corr', 'tdTom', 'stage', 'pool')]

embryo_design <- distinct(embryo_design)
rownames(embryo_design) <- embryo_design$sample_ko

# Test differential abundance per hood
da_results <- testNhoods(Milo_sce, design = ~ pool + tdTom_corr, design.df = embryo_design, reduced.dim="PCA")

# Save MILO results
outfile = sprintf("%s/%s_Milo_features%d_pcs%d_tdTomcorr%s_DAresults.csv",args$outdir, paste(args$stage, collapse ='_'), args$features, args$npcs, args$tdTom_corr)
write.csv(da_results, outfile, row.names=FALSE)


##################
## Plot results ##
##################

p1 = ggplot(da_results, aes(PValue)) + 
    geom_histogram(bins=50) + 
    theme_bw()

p2 = ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1) + ## Mark significance threshold (10% FDR)
    theme_bw()


## Plot single-cell UMAP
p3 = plotReducedDim(Milo_sce, dimred = "UMAP", colour_by="tdTom_corr", text_by = "celltype.mapped_mnn", 
                          text_size = 3, point_size=0.3, text_colour = "grey30") +
    scale_color_manual(values=opts$tdTom.color, name = 'tdTom') + 
    theme_void() + 
    theme(legend.position='right') +
    guides(fill="none")

## Plot neighbourhood graph
p4 = plotNhoodGraphDA(Milo_sce, da_results, layout="UMAP",alpha=0.05) 

# plot stats + graphs
outfile = sprintf("%s/%s_Milo_features%d_pcs%d_tdTomcorr%s_plots.pdf",args$outdir, paste(args$stage, collapse ='_'), args$features, args$npcs, args$tdTom_corr)
pdf(outfile, width=14, height = 7)
    p1 + p2
    p3 + p4 +
      plot_layout(guides="collect")
dev.off()

# If significant results create beeswarm plot
if(TRUE %in% unique(da_results$SpatialFDR < 0.05)){
    # Find Neighbourhood groups
    da_results <- groupNhoods(Milo_sce, da_results, max.lfc.delta = 0.5, overlap = 20)

    # Annotate hoods by celltype
    da_results <- annotateNhoods(Milo_sce, da_results, coldata_col = "celltype.mapped_mnn")

    # Annotate 'mixed' hoods when there is not one main celltype
    da_results$celltype.mapped_mnn <- ifelse(da_results$celltype.mapped_mnn_fraction < 0.4, "Mixed", da_results$celltype.mapped_mnn)

    # Save MILO object
    outfile = sprintf("%s/%s_Milo_features%d_pcs%d_tdTomcorr%s_DAresults.csv",args$outdir, paste(args$stage, collapse ='_'), args$features, args$npcs, args$tdTom_corr)
    write.csv(da_results, outfile, row.names=FALSE)
    
    # plot DA beeswarm
    p6 = plotDAbeeswarm(da_results, group.by = "celltype.mapped_mnn")

    # # save beeswarm
    outfile = sprintf("%s/%s_Milo_features%d_pcs%d_tdTomcorr%s_beeswarm.pdf",args$outdir, paste(args$stage, collapse ='_'), args$features, args$npcs, args$tdTom_corr)
    pdf(outfile, width=14, height = length(unique(da_results$celltype.mapped_mnn))/2.5)
        print(p6)
    dev.off()
}







