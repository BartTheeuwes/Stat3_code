suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################
here::i_am("processing/1_create_seurat_rna.R")
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results_new/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$samples <- opts$samples
# args$celltype_label <- "celltype.mapped_mnn"
# args$outdir <- file.path(io$basedir,"results_new/rna/celltype_proportions")
## END TEST ##

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples & !is.na(eval(as.name(args$celltype_label)))] %>%
  setnames(args$celltype_label,"celltype")

# if (!"stage"%in%colnames(sample_metadata)) {
#   sample_metadata[,stage:=stringr::str_replace_all(sample,opts$sample2stage)]
# }
sample_metadata[,stage:=stringr::str_replace_all(sample,opts$sample2stage)]

table(sample_metadata$stage)

################################################
## Calculate cell type proportions per sample ##
################################################

to.plot <- sample_metadata %>%
  .[,N:=.N,by="sample"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","stage","celltype")] %>%
  setorder(sample)  %>% .[,sample:=factor(sample,levels=args$samples)]


#########################
## Horizontal barplots ##
#########################

# Rename "_" to " " in cell types
# to.plot[,celltype:=stringr::str_replace_all(celltype,"_"," ")]
# names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=names(opts$celltype.colors))]

for (i in unique(to.plot$stage)) {
  p <- ggplot(to.plot[stage==i], aes(x=celltype, y=celltype_proportion)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~sample, nrow=1, scales="free_x") +
    coord_flip() +
    labs(y="Proportion of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots_per_sample.pdf",args$outdir,i), width=10, height=8)
  print(p)
  dev.off()
}


################################################
## Calculate cell type proportions by TdTom ##
################################################

to.plot <- sample_metadata %>%
  .[,N:=.N,by=c("tdTom_corr","stage")] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("tdTom_corr","stage","celltype")] 

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=names(opts$celltype.colors))]

for (i in unique(to.plot$stage)) {
  p <- ggplot(to.plot[stage==i], aes(x=celltype, y=celltype_proportion)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~tdTom_corr, nrow=1, scales="free_x") +
    coord_flip() +
    labs(y="Proportion of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    ) 

  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots_per_tdTom.pdf",args$outdir,i), width=10, height=8)
  print(p)
  dev.off()
  
  p <- ggplot(to.plot[stage==i], aes(x=celltype, y=celltype_proportion)) +
    geom_bar(aes(fill=tdTom_corr), stat="identity", color="black",  position = "dodge") +
    scale_fill_manual(values=c('black', 'red')) +
    coord_flip() +
    labs(y="Proportion of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    ) 
  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots_per_tdTom_dodged.pdf",args$outdir,i), width=10, height=8)
  print(p)
  dev.off()
}


#####################################################
## Calculate cell type proportions by TdTom No ExE ##
#####################################################

to.plot <- sample_metadata %>%
    .[!celltype%in%c("Visceral endoderm","ExE endoderm","ExE ectoderm","Parietal endoderm")] %>%
  .[,N:=.N,by=c("tdTom_corr","stage")] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("tdTom_corr","stage","celltype")] 

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=names(opts$celltype.colors))]

for (i in unique(to.plot$stage)) {
  p <- ggplot(to.plot[stage==i], aes(x=celltype, y=celltype_proportion)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~tdTom_corr, nrow=1, scales="free_x") +
    coord_flip() +
    labs(y="Proportion of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    ) 

  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots_per_tdTom_no_ExE.pdf",args$outdir,i), width=10, height=8)
  print(p)
  dev.off()
  
  p <- ggplot(to.plot[stage==i], aes(x=celltype, y=celltype_proportion)) +
    geom_bar(aes(fill=tdTom_corr), stat="identity", color="black",  position = "dodge") +
    scale_fill_manual(values=c('black', 'red')) +
    coord_flip() +
    labs(y="Proportion of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    ) 
  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots_per_tdTom_no_ExE_dodged.pdf",args$outdir,i), width=10, height=8)
  print(p)
  dev.off()
}



