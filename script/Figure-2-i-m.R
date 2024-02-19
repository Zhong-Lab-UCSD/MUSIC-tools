# Figure-2-i-m.R

require(GenomicRanges)
require(tidyverse)
require(tibble)
require(data.table)
require(magrittr)
require(InteractionSet)
library(parallel)
require(karyoploteR)
require(regioneR)
require(ggpubr)

`%ni%` = Negate(`%in%`)
options(stringsAsFactors=FALSE)

gtfgr36 <- readRDS("/database/Human/genome/Robj/gtf36_gene_granges.rds")
source("MUSIC_utils2.R") # import Range_from_Rle function from utils
MASTER_OUTDIR = "/outputs/"

## ---- Figure 2-i,j  RNA ends expression ----------


H1_cells <- readRDS("human_cells_reads_95purity_filtered.rds")

RNA_gr <- readRDS("/data/1_reads_gr_h5/human_purity95_readsfilter_RNA_gr.rds")

hg38_chr_size <- fread("/database/Human/genome/hg38_carno_chrom_size.txt")
hg38_chr_gr <- hg38_chr_size %>% dplyr::mutate(V3=1) %>% filter(V1 %ni% c("chrM")) %>%
  GenomicRanges::makeGRangesFromDataFrame(.,keep.extra.columns = T,  seqnames.field = "V1", start.field = "V3", end.field = "V2")

RNA_gr_plus <- RNA_gr %>% filter(strand == "+")
RNA_gr_neg <- RNA_gr %>% filter(strand == "-")

# MUSIC.K562.pos_cov <- coverage(c(RNA_gr_plus, hg38_chr_gr))-1 # in order to keep the full length of the genome.
# export.bw(MUSIC.K562.pos_cov, "/data/8_RNA_H1_compare/MUSIC_H1_bulk.pos.bw")

# MUSIC.K562.neg_cov <- coverage(c(RNA_gr_neg, hg38_chr_gr))-1 # in order to keep the full length of the genome.
# export.bw(MUSIC.K562.neg_cov, "/data/8_RNA_H1_compare/MUSIC_H1_bulk.neg.bw")


## 2. compare gene expression using bw 

setClass("SampleBWPath", slots=list(sample_name="character", pos_bw_path="character", neg_bw_path="character"))

iMARGI <- new("SampleBWPath", sample_name="iMARGI",
              pos_bw_path="/nascentTX/data/iMARGI_cis_trans_bw/iMARGI_H1_ctrl_gi.hg38.pos.bw",
              neg_bw_path="/nascentTX/data/iMARGI_cis_trans_bw/iMARGI_H1_ctrl_gi.hg38.neg.bw")

MUSIC <- new("SampleBWPath", sample_name="MUSIC",
             pos_bw_path="/data/8_RNA_H1_compare/MUSIC_H1_bulk.pos.bw",
             neg_bw_path="/data/8_RNA_H1_compare/MUSIC_H1_bulk.neg.bw")

RNAseq <- new("SampleBWPath", sample_name="RNA-seq",
              pos_bw_path="/database/RNA_SEQ/H1_RNAseq/H1-hESC_RNAseq_hg38/GSM2400153_ENCFF959NLR_plus_strand_signal_of_unique_reads_GRCh38.bigWig",
              neg_bw_path="/database/RNA_SEQ/H1_RNAseq/H1-hESC_RNAseq_hg38/GSM2400153_ENCFF829TDN_minus_strand_signal_of_unique_reads_GRCh38.bigWig")


gene_expr <- function(SampleBWPath_obj, TU){
  
  # path to coverage Rle
  
  pos_cov <- rtracklayer::import.bw(SampleBWPath_obj@pos_bw_path) %>%
    keepSeqlevels(., hg38_chr_size$V1, pruning.mode="coarse") %>% coverage(., weight = "score")
  
  neg_cov <- rtracklayer::import.bw(SampleBWPath_obj@neg_bw_path) %>%
    keepSeqlevels(., hg38_chr_size$V1, pruning.mode="coarse") %>% coverage(., weight = "score")
  
  if(sum(runValue(neg_cov)[[1]])<0){neg_cov <- neg_cov*(-1)}
  
  TU_pos <- TU[which(strand(TU)=="+")]
  TU_neg <- TU[which(strand(TU)=="-")]
  
  intron_pos_sum <- Range_from_Rle(RANGES = TU_pos, cov_rle = pos_cov, type = "sum")
  intron_neg_sum <- Range_from_Rle(RANGES = TU_neg, cov_rle = neg_cov, type = "sum")
  
  TU_pos$cov_sum <- intron_pos_sum
  TU_neg$cov_sum <- intron_neg_sum
  
  res <- c(TU_pos, TU_neg) %>% as_tibble %>% dplyr::select(width, cov_sum, gene_type, TU_ID) %>%
    mutate(mean_cov_perK = cov_sum*1000/width) %>%
    mutate(sample=SampleBWPath_obj@sample_name)
  
  return(res)
}

iMARGI_expr <- gene_expr(iMARGI, TU = gtfgr36 %>% mutate(TU_ID = seq(1, length(.))))
MUSIC_expr <- gene_expr(MUSIC, TU = gtfgr36 %>% mutate(TU_ID = seq(1, length(.))))
RNAseq_expr <- gene_expr(RNAseq, TU = gtfgr36 %>% mutate(TU_ID = seq(1, length(.))))



plt_compare_TU_expr <- function(x_tech, y_tech){
  
  # x_tech and y_tech is the name of compared technology that will be plotted and measure on x or y axis.
  # test is the name of the other compared technology
  x_expr <- eval(parse(text = paste0(x_tech,"_expr")))
  y_expr <- eval(parse(text = paste0(y_tech,"_expr")))
  
  plt_df <- data.frame(x_tech = x_expr$mean_cov_perK, 
                       y_tech = y_expr$mean_cov_perK, 
                       type = factor(x_expr$gene_type)) %>% 
    set_colnames(c(x_tech, y_tech, "type"))
  
  print(paste0("number of compared genes:", nrow(plt_df)))
  
  g <- ggscatter(plt_df, x=x_tech, y=y_tech, color="#435c61",
                 size=0.1, pch='.',alpha=0.1,
                 xlab=paste0(x_tech, ' RNA end (CPK)'),
                 ylab=paste0(y_tech, ' RNA end (CPK)'),
                 add = "reg.line",
                 conf.int = TRUE,  
                 cor.method = "spearman",
                 legend.title = "Density",
                 title = "Gene expression") + 
    xscale("log10", .format = T) +  yscale("log10", .format = T) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon")+
    gradient_fill("jama") + 
    stat_cor(aes(label = ..r.label..), size=2) +
    theme(
      axis.text = element_text(size=6),
      axis.title = element_text(size=8),
      legend.title=element_text(size=6),
      legend.text =element_text(size=6),
      legend.position = "None",
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0),
      legend.key.width = unit(0.1,"in"),
      legend.key.height = unit(0.1,"in"),
      legend.background = element_blank(),
      plot.title = element_text(hjust=0.5, size=8),
      axis.line = element_line(size=0.2),
      plot.margin = unit(c(0,0,0,0),'in')
    ) 
  
  return(g)
}
require(patchwork)
g_MUSIC_iMARGI <- 
  plt_compare_TU_expr(x_tech = "MUSIC", y_tech = "iMARGI") + 
  ggtitle("H1 gene expression")+xlab("ensemble MUSIC RNA reads")+ylab("iMARGI RNA reads")+
  theme(plot.title = element_blank())

g_MUSIC_RNAseq <- 
  plt_compare_TU_expr(x_tech = "MUSIC", y_tech = "RNAseq") + 
  ggtitle("H1 gene expression")+xlab("ensemble MUSIC RNA reads")+ylab("RNA-seq")+
  theme(axis.title.x = element_blank())



png("/Figures/Figure-2ij.png",
    width = 1.7, height = 3, units = 'in', res = 300)

g_MUSIC_RNAseq %>% aplot::insert_bottom(g_MUSIC_iMARGI) 

dev.off()


## ---- Figure K gene track ----------

setClass("SampleBWPath", slots=list(sample_name="character", pos_bw_path="character", neg_bw_path="character"))

iMARGI <- new("SampleBWPath", sample_name="iMARGI",
              pos_bw_path="/nascentTX/data/iMARGI_cis_trans_bw/iMARGI_H1_ctrl_gi.hg38.pos.bw",
              neg_bw_path="/nascentTX/data/iMARGI_cis_trans_bw/iMARGI_H1_ctrl_gi.hg38.neg.bw")

MUSIC <- new("SampleBWPath", sample_name="MUSIC",
             pos_bw_path="/data/8_RNA_H1_compare/MUSIC_H1_bulk.pos.bw",
             neg_bw_path="/data/8_RNA_H1_compare/MUSIC_H1_bulk.neg.bw")

RNAseq <- new("SampleBWPath", sample_name="RNA-seq",
              pos_bw_path="/database/RNA_SEQ/H1_RNAseq/H1-hESC_RNAseq_hg38/GSM2400153_ENCFF959NLR_plus_strand_signal_of_unique_reads_GRCh38.bigWig",
              neg_bw_path="/database/RNA_SEQ/H1_RNAseq/H1-hESC_RNAseq_hg38/GSM2400153_ENCFF829TDN_minus_strand_signal_of_unique_reads_GRCh38.bigWig")



top3_sc <- c("BC3_91-BC2_93-BC1_85", "BC3_42-BC2_50-BC1_17", "BC3_36-BC2_36-BC1_35")
SC_1 <- new("SampleBWPath", sample_name=top3_sc[1],
            pos_bw_path = paste0("/data/8_RNA_H1_compare/",top3_sc[1],"_MUSIC.pos.bw"),
            neg_bw_path = paste0("/data/8_RNA_H1_compare/",top3_sc[1],"_MUSIC.neg.bw"))

SC_2 <- new("SampleBWPath", sample_name=top3_sc[2], 
            pos_bw_path= paste0("/data/8_RNA_H1_compare/",top3_sc[2],"_MUSIC.pos.bw"),
            neg_bw_path= paste0("/data/8_RNA_H1_compare/",top3_sc[2],"_MUSIC.neg.bw"))

SC_3 <- new("SampleBWPath", sample_name=top3_sc[3],
            pos_bw_path= paste0("/data/8_RNA_H1_compare/",top3_sc[3],"_MUSIC.pos.bw"),
            neg_bw_path= paste0("/data/8_RNA_H1_compare/",top3_sc[3],"_MUSIC.neg.bw"))


require(karyoploteR)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)
require(zoo)

genometrackplot <- function(chr, start, end, log_=TRUE){ 
  
  chr.region <- toGRanges(paste0(chr,":",start,"-",end))
  # highlight_TU <- TU.gr[queryHits(findOverlaps(TU.gr, chr.region))]
  # plot params
  plot_params <- getDefaultPlotParams(plot.type=1)
  plot_params$data1inmargin <- 5
  plot_params$ideogramheight <- 10
  plot_params$bottommargin <- 0.05
  plot_params$topmargin <- 0.05
  plot_params$leftmargin <- 0.35
  plot_params$rightmargin <- 0.05
  
  genenamecex <- 0.85
  namecex <- 0.8 # track label fontsize 
  axiscex <- 0.8 # fontsize of axis labels 
  tickcex <- 0.5
  labelmargin <- 0.08 # distance from label to axis
  MARGIN=0.5 # distance between tracks 
  
  kp <- plotKaryotype(genome="hg38", zoom = chr.region, use.cache = F, plot.params = plot_params, 
                      ideogram.plotter=NULL, labels.plotter = NULL)
  
  # plot gene structure
  genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                      karyoplot=kp,
                                      plot.transcripts = TRUE, 
                                      plot.transcripts.structure = TRUE)
  genes.data <- addGeneNames(genes.data)
  genes.data.merge <- mergeTranscripts(genes.data)
  
  kp <- kpPlotGenes(kp, data=genes.data.merge, r0=0, r1=0.16, gene.name.cex = genenamecex, 
                    gene.name.position = "left",transcript.name.position = "left")
  
  
  histone.marks <- list(H3K27ac = "/database/ENCODE/H1_chipseq/ENCFF423TVA_H3K27ac.bw",
                        H3K4me3 = "/database/ENCODE/H1_chipseq/ENCFF629XVB_H3K4me3.bw"
  )
  
  TOTAL_track <- length(histone.marks)
  starts <- sapply(seq(1, TOTAL_track), function(x){autotrack(x,TOTAL_track,r0=0.18, r1=0.28, margin = MARGIN)$r0})
  ends <- sapply(seq(1, TOTAL_track), function(x){autotrack(x,TOTAL_track,r0=0.18, r1=0.28, margin = MARGIN)$r1})
  TRACK_colors <- c("grey50","grey50","grey50")
  
  for(i in seq(1,length(histone.marks))){
    
    kp <- kpPlotBigWig(kp, data=histone.marks[[i]], data.panel=1, 
                       col=TRACK_colors[i], border=TRACK_colors[i], 
                       r0=starts[i], r1=ends[i], ymax="visible.region")
    #col=TRACK_colors[i], r0=starts[i], r1=ends[i], ymax=7)
    computed.ymax <- kp$latest.plot$computed.values$ymax
    kpAxis(kp, ymin=0, ymax=round(computed.ymax,digits = 0), r0=starts[i], r1=ends[i], numticks = 2, cex=axiscex)
    kp <- kpAddLabels(kp, labels=names(histone.marks)[i], r0=starts[i], r1=ends[i], data.panel = 1, label.margin = labelmargin,
                      cex = namecex)
  }
  
  
  # plot negative strand 
  neg.tracks <- list(`MUSIC (Single cell 1)` =   SC_1@neg_bw_path,
                     `MUSIC (Single cell 2)` =   SC_2@neg_bw_path,
                     `MUSIC (Single cell 3)` =   SC_3@neg_bw_path,
                     `MUSIC (Ensemble)`      = MUSIC@neg_bw_path,
                     `iMARGI (bulk)`         = iMARGI@neg_bw_path,
                     `RNA-seq (bulk)`        = RNAseq@neg_bw_path
  )
  
  TOTAL_track <- length(neg.tracks)
  starts <- sapply(seq(1,TOTAL_track), function(x){autotrack(x,TOTAL_track,r0=0.3,r1=0.65,margin = MARGIN)$r0})
  ends <- sapply(seq(1,TOTAL_track), function(x){autotrack(x,TOTAL_track,r0=0.3,r1=0.65,margin = MARGIN)$r1})
  TRACK_colors <- rep("#aed7e6",length(neg.tracks))
  
  for(i in seq(1,length(neg.tracks))){
    
    bw.gr <- import.bw(neg.tracks[[i]])
    bw.cov <- coverage(bw.gr, weight = "score")
    y <- bw.cov[[chr]][start:end] %>% as.vector() 
    x <- seq(start, end)
    
    if(min(y)<0){y <- y*(-1)}
    
    if(log_ == T){
      y_log10 <- log(y+1)
    } else {y_log10 <- y}
    
    y_zoo <- zoo::zoo(y_log10) # use sliding window
    y_window <- zoo::rollapply(y_zoo, width = as.integer((end-start)/1000), by = 1, FUN = mean, align = "center")
    y_smooth <- y_log10
    y_smooth[index(y_window)] <- y_window
    
    
    kp <- kpArea(kp, chr=chr, x=x, y=y_smooth,  ymin=0, ymax=max(y_smooth)+0.001,
                 col=TRACK_colors[i], border=TRACK_colors[i], 
                 data.panel=1,
                 r0=ends[i], r1=starts[i])
    
    kpAxis(kp, ymin=0, ymax=round(max(y_log10),digits = 1)+0.1, r0=ends[i], r1=starts[i], numticks = 2, cex=axiscex)
    kp <- kpAddLabels(kp, labels=names(neg.tracks)[i], r0=starts[i], r1=ends[i], data.panel = 1, 
                      label.margin = labelmargin, cex = namecex)
    
  }
  
  
  
  # plot positive strand 
  pos.tracks <- list(`MUSIC (Single cell 1)` =   SC_1@pos_bw_path,
                     `MUSIC (Single cell 2)` =   SC_2@pos_bw_path,
                     `MUSIC (Single cell 3)` =   SC_3@pos_bw_path,
                     `MUSIC (Ensemble)`      = MUSIC@pos_bw_path,
                     `iMARGI (bulk)`         =   iMARGI@pos_bw_path,
                     `RNA-seq (bulk)`        = RNAseq@pos_bw_path
  )
  
  TOTAL_track <- length(pos.tracks)
  starts <- sapply(seq(1,TOTAL_track), function(x){autotrack(x,TOTAL_track,r0=0.66,r1=1,margin = MARGIN)$r0})
  ends <- sapply(seq(1,TOTAL_track), function(x){autotrack(x,TOTAL_track,r0=0.66,r1=1,margin = MARGIN)$r1})
  TRACK_colors <- rep("#c99575", length(pos.tracks))
  
  for(i in seq(1,length(pos.tracks))){
    
    bw.gr <- import.bw(pos.tracks[[i]])
    bw.cov <- coverage(bw.gr, weight = "score")
    y <- bw.cov[[chr]][start:end] %>% as.vector() 
    x <- seq(start, end)
    
    if(log_ == T){
      y_log10 <- log(y+1)
    } else {y_log10 <- y}
    
    
    y_zoo <- zoo::zoo(y_log10) # use sliding window
    y_window <- zoo::rollapply(y_zoo, width = as.integer((end-start)/1000), by = 1, FUN = mean, align = "center")
    y_smooth <- y_log10
    y_smooth[zoo::index(y_window)] <- y_window
    
    
    kp <- kpArea(kp, chr=chr, x=x, y=y_smooth, ymin=0, ymax=max(y_smooth)+0.001,
                 col=TRACK_colors[i], border =TRACK_colors[i], 
                 data.panel=1, r0=starts[i], r1=ends[i])
    
    kpAxis(kp, ymin=0, ymax=round(max(y_log10), digits = 1)+0.1, r0=starts[i], r1=ends[i], numticks = 2, cex=axiscex)
    kp <- kpAddLabels(kp, labels=names(pos.tracks)[i], r0=starts[i], r1=ends[i], data.panel = 1, 
                      label.margin = labelmargin, cex = namecex)
  }
  # add highlighter
  #kpPlotRegions(kp, data=highlight_TU, col="#AAAAAA10", avoid.overlapping=F)
  
  # dev.off()
}



png("/Figures/Figure-2k-XIST-sc-track.png",
    width = 5, height=8, units = 'in', res = 300)
genometrackplot(chr = "chrX", start = 73733735, end = 74360684, log_ = T) # XIST and JPX 
dev.off()




## ---- Figure 2i nsaRNA RAL ----------


DNA_human <- readRDS("/data/1_reads_gr_h5/human_purity95_readsfilter_DNA_gr.rds")
RNA_human <- readRDS("/data/1_reads_gr_h5/human_purity95_readsfilter_RNA_gr.rds")

nsaRNA_gr <- gtfgr36 %>%
  filter(gene_name %in% c("7SK", "MALAT1", "U1","U2","RNU4ATAC","U4","U5","U6","RNU6ATAC","U11","U12"))

nsa_RNA_gr <- RNA_human[subjectHits(findOverlaps(nsaRNA_gr, RNA_human, ignore.strand=F)) %>% unique]
nsa_DNA_gr <- DNA_human %>% filter(CBMB %in% nsa_RNA_gr$CBMB)

# K562 compartment
H1_compartment <- 
  fread("/HiC/compartments/H1_control/H1_control_1000000.txt") %>% 
  #mutate(eigen = (-1)*eigen) %>%  # score is reversed 
  mutate(compt = ifelse(eigen>0, "A", "B")) %>% mutate(end = coord + 1000000) %>% rename(start = coord) %>%  # always be careful if AB is not correspond to positive and negative
  makeGRangesFromDataFrame(., start.field = "start", end.field = "end", seqnames.field = "chr", keep.extra.columns = T) %>% 
  plyranges::mutate(start = start+1)

H1_compartment_reduce <- unlist(GenomicRanges::reduce(split(H1_compartment, H1_compartment$compt)))
H1_compartment_reduce$cmpt <- names(H1_compartment_reduce)
H1_compartment_reduce %<>% mutate(color = ifelse(cmpt=="A", "red", "steelblue"))

# export to bigwig
hg38_chr_size <- fread("/database/Human/genome/hg38_carno_chrom_size.txt")
hg38_chr_gr <- hg38_chr_size %>% dplyr::mutate(V3=1) %>% dplyr::filter(!grepl("chrM",V1)) %>%
  GenomicRanges::makeGRangesFromDataFrame(.,keep.extra.columns = T,seqnames.field = "V1",start.field = "V3",end.field = "V2")



RES = 1e6
hg38_1mb_tiles <- gen_windows(window.size = RES, species = 'hg38')

human_clusters <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_human.sort.cluster_rna_dna.csv"))


## MUSIC nsaRNA RAL 

nsa_clusizedf <- subset_RD_gen_clusizedf(nsa_DNA_gr, nsa_RNA_gr,
                                         DNA_subregion = NULL, RNA_subregion = NULL,
                                         CB_selected = NULL, CBMB_selected = NULL,
                                         clu_size_min = 2, clu_size_max = 1000)

nsa_RAL <- gen_RAL_1d(nsa_clusizedf$dna_gr, nsa_clusizedf$clu_size, hg38_1mb_tiles, single_cell=F)


## MUSIC pre-mRNA RAL, intron or exon-intron junction
require(GenomicFeatures)
txdb <- GenomicFeatures::makeTxDbFromGFF("/database/Human/genome/gencode.v36.primary_assembly.annotation.gff3",
                                         format = "gff")
introns <- intronsByTranscript(txdb) %>% unlist


mRNA <-  RNA_human[queryHits(findOverlaps(RNA_human, gtfgr36 %>% filter(gene_type=="protein_coding"), ignore.strand=F))]

intron_RNA_human <- mRNA[queryHits(findOverlaps(mRNA, introns, minoverlap = 15))]
intron_DNA_human <- DNA_human %>% filter(CBMB %in% intron_RNA_human$CBMB)

premRNA_clusizedf <- subset_RD_gen_clusizedf(intron_DNA_human, intron_RNA_human,
                                             DNA_subregion = NULL, RNA_subregion = NULL,
                                             CB_selected = NULL, CBMB_selected = NULL,
                                             clu_size_min = 2, clu_size_max = 1000)


premRNA_RAL <- gen_RAL_1d(premRNA_clusizedf$dna_gr, premRNA_clusizedf$clu_size, hg38_1mb_tiles, single_cell=F)




library(karyoploteR)

H1_SPIN <- fread("/database/4DN/SPIN/H1_new.SPIN.JAWG.25kb.9_state.bed")
H1_SPIN <- H1_SPIN %>% dplyr::filter(!grepl("NAN",V4)) %>%
  makeGRangesFromDataFrame(.,seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T) %>%
  filter(V4 %in% c("Speckle"))
#filter(V4 %in% c("Speckle"))

H1_compartment <- 
  fread("/HiC/compartments/H1_control/H1_control_1000000.txt") %>% 
  #mutate(eigen = (-1)*eigen) %>%  # score is reversed 
  mutate(compt = ifelse(eigen>0, "A", "B")) %>% mutate(end = coord + 1000000) %>% rename(start = coord) %>%  # always be careful if AB is not correspond to positive and negative
  makeGRangesFromDataFrame(., start.field = "start", end.field = "end", seqnames.field = "chr", keep.extra.columns = T) %>% 
  plyranges::mutate(start = start+1)

H1_compartment_reduce <- unlist(GenomicRanges::reduce(split(H1_compartment, H1_compartment$compt)))
H1_compartment_reduce$cmpt <- names(H1_compartment_reduce)
H1_compartment_reduce %<>% mutate(color = ifelse(cmpt=="A", "red", "steelblue"))


preRNA_gr <- hg38_1mb_tiles
preRNA_gr$y <- premRNA_RAL$window_weight
preRNA_gr %<>%
  filter(seqnames %in% c(paste0("chr",c(1:22)))) %>% group_by(seqnames) %>%
  mutate(y = y-min(y)) %>%
  mutate(chr_max = max(y)) %>% mutate(y = y/chr_max) %>% mutate(y = ifelse(is.nan(y), 0, y)) %>%
  ungroup()

nsaRNA_gr <- hg38_1mb_tiles
nsaRNA_gr$y <- nsa_RAL$window_weight
nsaRNA_gr %<>%
  filter(seqnames %in% c(paste0("chr",c(1:22)))) %>% group_by(seqnames) %>%
  mutate(y = y-min(y)) %>%
  mutate(chr_max = max(y)) %>% mutate(y = y/chr_max) %>% mutate(y = ifelse(is.nan(y), 0, y)) %>%
  ungroup()

saveRDS(preRNA_gr, "/data/Figure-2l-preRNA_gr.rds")
saveRDS(nsaRNA_gr, "/data/Figure-2l-nsaRNA_gr.rds")

preRNA_gr <- readRDS("/data/Figure-2l-preRNA_gr.rds")
nsaRNA_gr <- readRDS("/data/Figure-2l-nsaRNA_gr.rds")


# plot without TSA signals and iMARGI signals
plot_params <- getDefaultPlotParams(plot.type=1)
plot_params$data1inmargin <- 5
plot_params$data1height <- 200
plot_params$ideogramheight <- 20
plot_params$bottommargin <- 0.05
plot_params$topmargin <- 50
plot_params$leftmargin <- 0.3
plot_params$rightmargin <- 0.05

genenamecex <- 0.85
namecex <- 1 # track label fontsize
axiscex <- 0.8 # fontsize of axis labels
tickcex <- 0.5
labelmargin <- 0.01 # distance from label to axis


pbmcapply::pbmclapply(c("chr1"), function(CHR){
  
  
  # png(paste0("/Figures/4_RNA_DNA_val/nsaRNA_RAL_chr/",CHR,"_nsaRNA_preRNA_TSA_compartment_lines.png"),
  #     width = 8, height = 5.5, units = 'in', res = 400)
  pdf(paste0("/Figures/Figure-2l-nsaRNA_preRNA_compartment_lines.pdf"),
      width = 5, height = 5)
  
  kp <- plotKaryotype(genome = 'hg38', chromosomes = CHR, plot.params = plot_params)
  
  
  # MUSIC nsaRNA
  kpLines(kp, r0 = 0.3, r1 = 0.45, lwd =2,
          data = nsaRNA_gr, col="#75cac4")
  kpAddLabels(kp, labels="nsaRNA",  r0 = 0.32, r1 = 0.42, data.panel = 1,
              label.margin = labelmargin, cex = namecex)
  
  # MUSIC pre-mRNA
  kpLines(kp, r0 = 0.5, r1 = 0.65, lwd =2,
          data = preRNA_gr, col="#ec756f")
  kpAddLabels(kp, labels="pre-mRNA", r0 = 0.52, r1 = 0.62, data.panel = 1,
              label.margin = labelmargin, cex = namecex)
  
  
  
  # SPIN states
  kpPlotRegions(kp, data=H1_SPIN,
                r0=0.04, r1=0.10, col="#842b2d")
  kpAddLabels(kp, labels="SPIN-Speckle",  r0=0.0, r1=0.12, data.panel = 1,
              label.margin = labelmargin, cex = namecex)
  
  # H1 compartment score
  kpPlotRegions(kp, data=H1_compartment_reduce,
                r0=0.16, r1=0.22, col=H1_compartment_reduce$color)
  kpAddLabels(kp, labels="A/B",  r0=0.14, r1=0.24, data.panel = 1,
              label.margin = labelmargin, cex = namecex)
  
  kpAddMainTitle(kp, main = CHR, cex = namecex*1.5)
  
  dev.off()
  
})


## ---- Figure-2m-scatter-nsaRNA-premRNA ----------

g_nsa_pRNA <- ggpubr::ggscatter(data.frame(nsa_RAL=nsa_RAL$window_weight, premRNA_RAL=premRNA_RAL$window_weight),
                                x='nsa_RAL', y='premRNA_RAL', conf.int = TRUE,  cor.method = 'spearman', add = 'reg.line',
                                xlab="nsaRNA normalized counts",
                                ylab="pre-mRNA\nnormalized counts",
                                colour="black", pch=21, alpha=0.2, size=0.2)+
  xscale('log2', .format = T) + yscale('log2', .format = T)+
  stat_cor(aes(label = ..r.label..), size=2) +
  theme(
        axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        legend.title=element_text(size=6),
        legend.text =element_text(size=6),
        legend.position = "None",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(size=0.2),
        plot.margin = unit(c(0,0,0,0),'in'))


require(patchwork)
pdf("/Figures/Figure-2m-nsaRAL_preRNA_TSA_cor.pdf",
    width = 1.8, height = 1.45)
g_nsa_pRNA
dev.off()

