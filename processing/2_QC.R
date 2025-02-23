suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(argparse))

#####################
## Define arguments ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--metadata',       type="character",                    help='Metadata')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--min_nCount_RNA',       type="double",                    help='Minimum number of log10(reads)')
p$add_argument('--max_nCount_RNA',       type="double",                    help='Maximum number of log10(reads)')
p$add_argument('--min_nFeature_RNA',       type="integer",                    help='Minimum number of expressed genes')
p$add_argument('--max_nFeature_RNA',       type="integer",                    help='Maximum number of expressed genes')
p$add_argument('--mitochondrial_percent_RNA',       type="integer",                    help='Maximum percentage of mitochondrial reads')
p$add_argument('--ribosomal_percent_RNA',       type="integer",                    help='Maximum percentage of ribosomal reads')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE)) 

message(args$min_nCount_RNA) 
message(args$max_nCount_RNA) 


message(args$min_nFeature_RNA) 
message(args$max_nFeature_RNA) 

message(args$mitochondrial_percent_RNA) 
message(args$ribosomal_percent_RNA) 

#####################
## Define settings ##
#####################
here::i_am("processing/1_create_seurat_rna.R")
source(here::here("settings.R"))

# Sanity checks
stopifnot(args$samples%in%opts$samples)

###############
## Load data ##
###############

metadata <- fread(args$metadata) %>% 
    .[sample%in%args$samples] %>%
    # .[,pass_rnaQC:=nFeature_RNA>args$min_nFeature_RNA & nCount_RNA>2**args$log_nCount_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA]
    .[,pass_rnaQC := nFeature_RNA<=args$max_nFeature_RNA & nFeature_RNA>=args$min_nFeature_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA & ribosomal_percent_RNA<args$ribosomal_percent_RNA & log10(nCount_RNA)<=args$max_nCount_RNA & log10(nCount_RNA)>=args$min_nCount_RNA]

table(metadata$pass_rnaQC)

#####################
## Plot QC metrics ##
#####################
# Pre QC plot
to.plot <- metadata %>% 
    .[nFeature_RNA<=15000 & mitochondrial_percent_RNA<=60 & ribosomal_percent_RNA<=60] %>% # remove massive outliers for plotting
    .[,nCount_RNA := log10(nCount_RNA)] %>% 
    melt(id.vars=c("sample","cell","stage"), measure.vars=c("nCount_RNA", "nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))

facet.labels <- c("nCount_RNA" = "log10 num. of reads", "nFeature_RNA" = "Num. of genes", "mitochondrial_percent_RNA" = "Mitochondrial %", "ribosomal_percent_RNA" = "Ribosomal %")

hline = data.frame(variable = c('nCount_RNA', 'nCount_RNA', 'nFeature_RNA', 'nFeature_RNA', 'mitochondrial_percent_RNA', 'ribosomal_percent_RNA'),
                  value = c(args$min_nCount_RNA, args$max_nCount_RNA, args$min_nFeature_RNA, args$max_nFeature_RNA, args$mitochondrial_percent_RNA, args$ribosomal_percent_RNA))

## Box plot 
p <- ggplot(to.plot, aes_string(x="sample", y="value", fill="stage")) +
    geom_violin() +
    geom_boxplot(outlier.shape=NA, coef=1.5, fill='white', width=0.2) +
    geom_hline(data=hline, aes(yintercept=value), color='red', linetype='longdash', size=1) + 
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    scale_fill_manual(values=opts$stage.colors) +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
        axis.title.x = element_blank()
    )

pdf(sprintf("%s/qc_metrics_boxplot_preQC.pdf",args$outputdir), width=13, height=5)
    print(p)
dev.off()

# Post QC plot
to.plot <- metadata %>% .[pass_rnaQC==TRUE] %>%
    .[nFeature_RNA<=15000 & mitochondrial_percent_RNA<=60 & ribosomal_percent_RNA<=60] %>% # remove massive outliers for plotting
    .[,nCount_RNA := log10(nCount_RNA)] %>% 
    melt(id.vars=c("sample","cell","stage"), measure.vars=c("nCount_RNA", "nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))

facet.labels <- c("nCount_RNA" = "Num. of reads", "nFeature_RNA" = "Num. of genes", "mitochondrial_percent_RNA" = "Mitochondrial %", "ribosomal_percent_RNA" = "Ribosomal %")
    
## Box plot 
p <- ggplot(to.plot, aes_string(x="sample", y="value", fill="stage")) +
    geom_violin() +
    geom_boxplot(outlier.shape=NA, coef=1.5, fill='white', width=0.2) +
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    scale_fill_manual(values=opts$stage.colors) +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
        axis.title.x = element_blank()
    )

pdf(sprintf("%s/qc_metrics_boxplot_postQC.pdf",args$outputdir), width=13, height=5)
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir))
    print(p)
dev.off()

## histogram 

tmp <- data.table(
    variable = c("nCount_RNA", "nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA"),
    value = c(0, args$min_nFeature_RNA, args$mitochondrial_percent_RNA, args$ribosomal_percent_RNA)
)
# tmp <- data.table(
#     variable = c("nFeature_RNA"),
#     value = c(args$min_nFeature_RNA)
# )

p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
# p <- ggdensity(to.plot, x="value", fill="sample") +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free", nrow=1) +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.75))
    )
    
pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outputdir), width=17, height=6)
print(p)
dev.off()


#########################################################
## Plot fraction of cells that pass QC for each sample ##
#########################################################

to.plot <- metadata %>%
    .[,mean(pass_rnaQC,na.rm=T),by=c("sample","stage")]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="stage") +
    scale_fill_manual(values=opts$stage.colors) +
    labs(x="", y="Fraction of cells that pass QC (RNA)") +
    # facet_wrap(~stage)
    theme(
        legend.position = "none",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
    )

pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outputdir), width=8, height=5)
    print(p)
dev.off()

##############################
## Plot nFeatures vs nGenes ##
##############################

to.plot <- metadata 

p <- ggplot(to.plot, aes(nCount_RNA, nFeature_RNA, col=pass_rnaQC)) +
    geom_point() +
    geom_hline(yintercept = c(args$min_nFeature_RNA, args$max_nFeature_RNA), color='red', linetype='longdash', size=1) +
    geom_vline(xintercept = c(10^args$min_nCount_RNA, 10^args$max_nCount_RNA), color='red', linetype='longdash', size=1) +
    # facet_wrap(~stage)
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
    )

pdf(sprintf("%s/qc_nFeatures_nGenes.pdf",args$outputdir), width=8, height=5)
    print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, paste0(args$outputdir,"/sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")

