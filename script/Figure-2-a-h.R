# Figure-2.R


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
options(scipen=999999999)


MASTER_OUTDIR = "/outputs/"
source("MUSIC_utils2.R")

gtfgr36 <- readRDS("/mnt/extraids/SDSC_NFS/wenxingzhao/database/Human/genome/Robj/gtf36_gene_granges.rds")

human_cells <- readRDS("/data/human_cells_reads_95purity_filtered.rds")

DNA_human <- readRDS("/data/1_reads_gr_h5/human_purity95_readsfilter_DNA_gr.rds")

RES = 1e6
hg38_1mb_tiles <- gen_windows(window.size = RES, species = 'hg38')

human_clusters <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_human.sort.cluster_rna_dna.csv"))


## ---- Figure 2a ----------

# whole chromosome contact matrix of MUSIC and Micro-C

MUSIC_HiC_comb_heatmap <- function(CHR){
  
  hg38_chr_gr <- get_hg38_chrgr()
  hg38_CHR_gr <- hg38_1mb_tiles %>% plyranges::filter(seqnames == CHR)
  
  CHR_bins_ID <- hg38_CHR_gr %>% as.data.frame() %>%
    mutate(binID = paste0(seqnames,"_",start-1)) %>% mutate(binID_num = as.character(seq(1, nrow(.)))) %>%
    dplyr::select(binID, binID_num)
  
  DNA_human_chr <- DNA_human %>% plyranges::filter(seqnames == CHR)
  DNA_chr_cluster_size <- DNA_human_chr %>% plyranges::group_by(CB, CBMB) %>%
    plyranges::summarize(dna_reads = n()) %>% as.data.frame() %>%
    dplyr::filter(dna_reads > 1)
  
  onechr <- gen_DD_2d_diag_contact(dna_gr = DNA_human_chr,
                                   clu_size_df = DNA_chr_cluster_size %>% filter(dna_reads<1000),
                                   tiles = hg38_CHR_gr,
                                   CB_selected = NULL, clu_selected=NULL,
                                   return_sparse=F, diag_rm = T)
  saveRDS(onechr, paste0(
    "/data/Figure-2a-",CHR,"-music-raw-mtx.rds"))
  
  
  
  E14_HiC_chr <- strawr::straw("NONE", "/H1/4DNFI2TK7L2F.hic",
                               gsub("chr","",CHR), gsub("chr","",CHR), "BP", 1000000) %>%
    mutate(x = paste0(CHR, "_", x),
           y = paste0(CHR, "_", y)) %>% filter(x!=y) %>%
    left_join(CHR_bins_ID, by=c("x"="binID")) %>% rename(x_binID_num = binID_num) %>%
    left_join(CHR_bins_ID, by=c("y"="binID")) %>% rename(y_binID_num = binID_num)
  
  contact_matrix_full = matrix(0, nrow = nrow(CHR_bins_ID), ncol = nrow(CHR_bins_ID))
  rownames(contact_matrix_full) <- CHR_bins_ID$binID_num
  colnames(contact_matrix_full) <- CHR_bins_ID$binID_num
  contact_matrix_full[as.matrix(E14_HiC_chr[,c("x_binID_num","y_binID_num")])] = E14_HiC_chr$counts
  contact_matrix_full[as.matrix(E14_HiC_chr[,c("y_binID_num","x_binID_num")])] = E14_HiC_chr$counts
  diag(contact_matrix_full) <- 0
  
  saveRDS(contact_matrix_full, 
          paste0("/data/Figure-2a-",CHR,"-MicroC-raw-mtx.rds"))
  
  
}

MUSIC_HiC_comb_heatmap("chr1")


# plot
CHR = "chr1"
CHR_bins_ID <-gen_windows(window.size = RES, species = 'hg38') %>% filter(seqnames==CHR)  %>% as.data.frame()


onechr <- readRDS(paste0("/data/Figure-2a-",CHR,"-music-raw-mtx.rds"))
onechr <- log(onechr+1)
onechr <-onechr/quantile(onechr, 1)

contact_matrix_full <- readRDS(paste0("/data/Figure-2a-",CHR,"-MicroC-raw-mtx.rds"))
contact_matrix_full <- log(contact_matrix_full+1)
contact_matrix_full <- contact_matrix_full/quantile(contact_matrix_full,1)

new <- matrix(0, nrow = nrow(CHR_bins_ID), ncol = nrow(CHR_bins_ID))
new[upper.tri(new)] <- contact_matrix_full[upper.tri(contact_matrix_full)]
new[lower.tri(new)] <- onechr[lower.tri(onechr)]


pdf("/Figures/Figure-2a-chr-MUSIC-MICROC-contact.pdf",
    width = 2.4, height = 2.1)
g_all_chrs_map <- single_hm_from_full_mtx(contact_matrix_full = new,
                                          log_values = F, legend_position = 'right',
                                          scale_factor = 1)+
  ggtitle("chromatin interactions\nMicro-C vs. ensemble MUSIC") + 
  theme(plot.title = element_text(hjust=0.5, size=8),
        legend.text = element_text(size=5),
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.key.width = unit(0.1,'in')) +
  annotate("text", x=nrow(CHR_bins_ID)-45, y=nrow(CHR_bins_ID)-45, label="Micro-C", size=2.5) +
  annotate("text", x=45, y=45, label="MUSIC", size=2.5)
print(g_all_chrs_map)
dev.off()


## ---- Figure 2e ----------

# genomic distance contact frequency separated by cluster size in H1 compared with MicroC


## DD pair distance 

human_cells <- readRDS("human_cells_reads_95purity_filtered.rds")

# # function cal dist bin freq for each cell
# df2dist <- function(df){
#   
#   # df is from h5 DNA group
#   clu_2_100 <- df %>% add_count(CBMB) %>% filter(n<=100 & n>=1) %>% mutate(ID = seq(1, nrow(.)))
#   res <- clu_2_100 %>% left_join(clu_2_100, by=c("CBMB"="CBMB")) %>% 
#     filter(ID.x < ID.y, seqnames.x == seqnames.y) %>% 
#     dplyr::rename(clu_size = `n.x`) %>% 
#     mutate(dist = abs(start.x - start.y)) %>%
#     dplyr::filter(dist>10) %>% 
#     mutate(dist_cut = cut(dist, dist_bins)) %>%
#     select(clu_size, dist_cut) %>% 
#     dplyr::group_by(clu_size, dist_cut) %>% dplyr::summarize(ct = length(clu_size))
#   return(res)
# }
# 
# human_DD_df <- lapply(human_cells, function(CB){
#   print(paste0("DD-",CB))
#   DNA_h5 <- H5Fopen("/data/1_reads_gr_h5/human_purity95_readsfilter_reads.h5")
#   group <- paste0(CB,"/DNA")
#   df <- h5read(file = DNA_h5, name = group)
#   h5closeAll()
#   df_res <- df2dist(df)
#   return(df_res)
# }) %>% reduce(full_join, by=c("clu_size"="clu_size", "dist_cut"="dist_cut")) %>% 
#   mutate_if(is.numeric , replace_na, replace = 0) %>% 
#   filter(!is.na(dist_cut)) %>% 
#   rowwise() %>% mutate(count = sum(c_across(starts_with("ct"))))
#   
# 
# saveRDS(human_DD_df, "./data/6_DD_RD_RR_genomic_dist_contact_freq/human_DD_intra_clu_dist.rds")
# rm(human_cell_split)
# rm(human_DD_df)


human_DD_df <- readRDS("/data/6_DD_RD_RR_genomic_dist_contact_freq/human_DD_intra_clu_dist.rds") %>%
  dplyr::select(clu_size, dist_cut, count)

DIST_BINS <- seq(10, 2.5*10e8, length.out = 50000)
#DIST_BINS <- c(0, DIST_BINS)
# only plot upto chromosome long arm
hg38_arm_size <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
                       col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
hg38_arm_size %>%
  mutate(width = chromEnd - chromStart) %>%
  mutate(arm = substring(name, 1,1)) %>%
  group_by(chrom, arm) %>% summarize(size = sum(width)) %>%
  arrange(size) %>% filter(chrom %in% c(paste0("chr",1:22),"chrX"))
# longest chrom arm is 125mbp.
DIST_BINS <- DIST_BINS[DIST_BINS<15000000]


range_bin_summary <- function(clu_size_min,
                              clu_size_max,
                              dist_bins,
                              size_label){
  
  human_DD_clu2_intra_aggr <- human_DD_df %>%
    dplyr::filter(clu_size >= clu_size_min & clu_size <= clu_size_max) %>%
    mutate(dist_cut = as.character(dist_cut)) %>%
    dplyr::filter(dist_cut != "(0,10]") %>%
    mutate(dist_end = gsub("^.*,", "", dist_cut)) %>%
    mutate(dist_end = gsub("]","", dist_end)) %>%
    mutate(dist_end = as.numeric(dist_end)) %>%
    dplyr::filter(dist_end < 150000000) %>%
    dplyr::group_by(dist_cut) %>% summarize(ct = sum(count)) %>%
    mutate(freq_ratio = ct / sum(ct)) %>%
    mutate(size = size_label)
  
  return(human_DD_clu2_intra_aggr)
}

human_DD_clu2 <- range_bin_summary(clu_size_min = 2, clu_size_max = 2, dist_bins = DIST_BINS, size_label = "2")
human_DD_clu3_5 <- range_bin_summary(clu_size_min = 3, clu_size_max = 5, dist_bins = DIST_BINS, size_label = "[3,5]")
human_DD_clu6_10 <- range_bin_summary(clu_size_min = 6, clu_size_max = 10, dist_bins = DIST_BINS, size_label = "[6,10]")
human_DD_clu10_20 <- range_bin_summary(clu_size_min = 11, clu_size_max = 20, dist_bins = DIST_BINS, size_label = "[11,20]")
human_DD_clu20_50 <- range_bin_summary(clu_size_min = 21, clu_size_max = 50, dist_bins = DIST_BINS, size_label = "[21,50]")
human_DD_clu50_100 <- range_bin_summary(clu_size_min = 51, clu_size_max = 100, dist_bins = DIST_BINS, size_label = "[51,100]")


# the Pc(S) function of micro-c data is computed using cooltools 

dist_bins <- seq(10, 1.5*10e8, length.out = 50000)

microC_GDCF <- fread("/H1/4DNFI9GMP2J8_mcool_1k_genome_distance_contact_freq.txt")

resolution = 1000

res <- 
  microC_GDCF %>% mutate(s_bp = dist*resolution) %>% 
  group_by(s_bp) %>% 
  summarize(ncount = sum(count.sum)) %>% 
  mutate(dist_cut = cut(s_bp, dist_bins)) %>% 
  dplyr::filter(dist_cut != "(10,30010.6]") %>%
  mutate(dist_end = gsub("^.*,", "", dist_cut)) %>%
  mutate(dist_end = gsub("]","", dist_end)) %>%
  mutate(dist_end = as.numeric(dist_end)) %>%
  dplyr::filter(dist_end < 150000000) %>%
  dplyr::group_by(dist_cut) %>% summarize(ct = sum(ncount)) %>%
  mutate(freq_ratio = ct / sum(ct)) %>% mutate(size="Micro-C")



g_sciMARGI_microc_GDCF <- 
  rbind(
    human_DD_clu2,
    human_DD_clu3_5,
    human_DD_clu6_10,
    human_DD_clu10_20,
    human_DD_clu20_50,
    human_DD_clu50_100,
    res) %>%
  mutate(size = factor(size, levels = c("2", "[3,5]","[6,10]","[11,20]","[21,50]","[51,100]", "Micro-C"),
                       labels=c("2","3-5","6-10","11-20","21-50","51-100","Micro-C"))) %>%
  mutate(dist_end = gsub("^.*,", "", dist_cut)) %>%
  mutate(dist_end = gsub("]","", dist_end)) %>%
  mutate(dist_end = as.numeric(dist_end)) %>%
  ggplot() +
  #  annotate("rect", xmin=10^(6), xmax = 2*10^7, ymin=10^(-4), ymax=0.1, alpha=0.1)+
  geom_line(aes(x=dist_end, y=freq_ratio,
                color = size),
            size=0.3, alpha = 0.7) +
  scale_x_log10(labels=scales::trans_format('log10',math_format(10^.x)),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(4.5e4, 2*1e7)) +
  scale_y_log10(labels=scales::trans_format('log10',math_format(10^.x)),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(10^(-4.2), 0.2),
                expand = c(0, 0)) +
  annotation_logticks(size=0.2) +
  scale_color_manual(name = "DNA cluster size", values = c(RColorBrewer::brewer.pal(7,"Blues")[2:7], "black"))+
  labs(x = expression(paste("Genomic distance, ",italic("s")," (bp)")),
       y = expression(paste("Contact probability, ",italic("Pc"))))+
  annotate("text", x=10^(5.4), y=10^(-1), label="Micro-C",size=2.5)+
  annotate("text", x=10^(5), y=10^(-1.9), label="MUSIC",size=2.5, color ="#6e8cb9")+
  ggtitle("Chromatin interactions\n ") +
  
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill="white"),
    axis.text = element_text(size=6),
    axis.title = element_text(size=8),
    legend.position = c(0.8,0.8),
    plot.title = element_text(size=8, hjust=0.5),
    legend.text = element_text(size=6),
    legend.title = element_blank(),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0),
    legend.key.width = unit(0.1,"in"),
    legend.key.height = unit(0.1,"in"),
    legend.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "in"),
    axis.ticks = element_line(size=0.2),
    axis.line= element_line(size = 0.2, linetype = "solid", colour = "black")
  )


pdf("/Figures/Figure-2e-GDCF.pdf",
    width = 2.2, height = 2.2)
g_sciMARGI_microc_GDCF
dev.off()


## ---- Figure-2 b-d ----------

# TAD level cluster sizes

# normalize to the max value and when plot log scale and scale to 20
# require(BSgenome.Mmusculus.UCSC.hg38)

RES2 = 50000 
cluster_size <- 1000


compute_TAD_region_with_different_clu_size <- function(ROI_CHR, ROI_start, ROI_end,
                                                    cluster_size_min, cluster_size_max, pre_log_scale,
                                                    RES2 = 50000){
  
  ROI <- toGRanges(paste0(ROI_CHR,":",ROI_start,"-",ROI_end), genome = 'hg38')
  
  ROI_gr_bins <- slidingWindows(ROI, width = RES2, step = RES2) %>% as_granges() %>% plyranges::filter(width==RES2)
  
  ROI_DNA_gr <- DNA_human[queryHits(findOverlaps(DNA_human, ROI)) %>% unique]
  
  ROI_bins_ID <- ROI_gr_bins %>% as.data.frame() %>% 
    mutate(binID = paste0("bin",start)) %>% mutate(binID_num = as.character(seq(1, nrow(.)))) %>% 
    dplyr::select(binID, binID_num)
  
  DNA_chr_cluster_size_raw <- ROI_DNA_gr %>% plyranges::group_by(CB, CBMB) %>%
    plyranges::summarize(dna_reads = n()) %>% as.data.frame() 
  
  DNA_chr_cluster_size <- DNA_chr_cluster_size_raw %>%
    dplyr::filter(dna_reads >= cluster_size_min & dna_reads <= cluster_size_max)
  
  ROI_mtx <- gen_DD_2d_diag_contact(dna_gr = ROI_DNA_gr,
                                    clu_size_df = DNA_chr_cluster_size, # change cluster size within the ROI is more useful because only serveral mb off diag should be largely influence by large clusters.
                                    tiles = ROI_gr_bins,
                                    CB_selected = NULL, clu_selected=NULL,
                                    return_sparse=F, diag_rm = T)
  
  saveRDS(ROI_mtx, paste0(
    "/data/Figure2-bcd-MUSIC-clumin", cluster_size_min,
    "-max",cluster_size_max,".rds"))
  
  
  E14_HiC_chr <- strawr::straw("NONE", "/H1/4DNFI2TK7L2F.hic",
                               paste0(gsub("chr","",ROI_CHR),":",ROI_start,":",ROI_end-1),
                               paste0(gsub("chr","",ROI_CHR),":",ROI_start,":",ROI_end-1), "BP", RES2) %>%
    mutate(x = paste0("bin", x),
           y = paste0("bin", y)) %>% filter(x!=y) %>% 
    left_join(ROI_bins_ID, by=c("x"="binID")) %>% rename(x_binID_num = binID_num) %>% 
    left_join(ROI_bins_ID, by=c("y"="binID")) %>% rename(y_binID_num = binID_num)
  
  contact_matrix_full = matrix(0, nrow = nrow(ROI_bins_ID), ncol = nrow(ROI_bins_ID))
  rownames(contact_matrix_full) <- ROI_bins_ID$binID_num
  colnames(contact_matrix_full) <- ROI_bins_ID$binID_num
  contact_matrix_full[as.matrix(E14_HiC_chr[,c("x_binID_num","y_binID_num")])] = E14_HiC_chr$counts
  contact_matrix_full[as.matrix(E14_HiC_chr[,c("y_binID_num","x_binID_num")])] = E14_HiC_chr$counts
  
  saveRDS(contact_matrix_full, paste0(
    "/data/Figure2-bcd-MICROC.rds"))

}

purrr::pmap(list(c(1, 10, 50), c(10, 50, 100), c(10,1,1)), function(x, y, z){
  compute_TAD_region_with_different_clu_size(ROI_CHR = "chr12",
                                          ROI_start = 112000000, ROI_end = 122000000,
                                          cluster_size_min = x, cluster_size_max = y, pre_log_scale = z)
  
  
  
plot_TAD_region_with_different_clu_size <- function(ROI_CHR, ROI_start, ROI_end,
                                                    cluster_size_min, cluster_size_max, pre_log_scale,
                                                    RES2 = 50000, title){
  
  ROI <- toGRanges(paste0(ROI_CHR,":",ROI_start,"-",ROI_end), genome = 'hg38')
  
  ROI_gr_bins <- slidingWindows(ROI, width = RES2, step = RES2) %>% as_granges() %>% plyranges::filter(width==RES2)

  ROI_bins_ID <- ROI_gr_bins %>% as.data.frame() %>% 
    mutate(binID = paste0("bin",start)) %>% mutate(binID_num = as.character(seq(1, nrow(.)))) %>% 
    dplyr::select(binID, binID_num)
  
  # plot 
  ROI_mtx <- readRDS(paste0(
    "/data/Figure2-bcd-MUSIC-clumin", cluster_size_min,
    "-max",cluster_size_max,".rds"))
  
  contact_matrix_full <- readRDS(paste0(
    "/data/Figure2-bcd-MICROC.rds"))
  
  
  ROI_mtx_log <- log(ROI_mtx*pre_log_scale+1)
  
  ROI_mtx_norm <- ROI_mtx_log/quantile(ROI_mtx_log,1)
  contact_matrix_full_log <- log(contact_matrix_full+1)
  contact_matrix_full_norm <- contact_matrix_full_log/quantile(contact_matrix_full_log,1)
  
  new <- matrix(0, nrow = nrow(ROI_bins_ID), ncol = nrow(ROI_bins_ID))
  new[upper.tri(new)] <- contact_matrix_full_norm[upper.tri(contact_matrix_full_norm)]
  new[lower.tri(new)] <- ROI_mtx_norm[lower.tri(ROI_mtx_norm)]
  
  
  g_all_chrs_map <- single_hm_from_full_mtx(contact_matrix_full = new, 
                                            log_values = F, legend_position = 'None',
                                            scale_factor = 1)+
    ggtitle(paste0(title,"\ncluster size: ",cluster_size_min+1,"-",cluster_size_max)) +
    theme(plot.title = element_text(hjust=0.5, size=8),
          legend.text = element_text(size=6),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    annotate("text", x=nrow(ROI_bins_ID)-35, y=nrow(ROI_bins_ID)-35, label="Micro-C", size=2) +
    annotate("text", x=35, y=30, label=paste0("MUSIC"), size=2)
  
  g_all_chrs_map
  
  
}



pdf("/Figures/Figure2-b-TAD_region_diff_clusize_chr12_112000000_122000000.pdf",
    width = 4.2, height = 1.8) # directly use this 

region3_lists <- purrr::pmap(list(c(1, 10, 50), c(10, 50, 100), c(10,1,1), c("Small clusters","Median clusters","Large clusters")), function(x, y, z, t){
  plot_TAD_region_with_different_clu_size(ROI_CHR = "chr12",
                                          ROI_start = 112000000, ROI_end = 122000000,
                                          cluster_size_min = x, cluster_size_max = y, pre_log_scale = z, title=t)
})

cowplot::plot_grid(plotlist = region3_lists, nrow = 1)
dev.off()

## ---- Figure 2f plot split clusters ----------
# define plotting region

CHR = "chr12"
START = 114000000
END = 117000000
RES = 10000
hg38_tiles <- gen_windows(window.size = RES, species = 'hg38') %>% filter(seqnames %in% paste0("chr",c(1:22,"X")))
ROI_gr <- toGRanges(paste0(CHR,":",START,"-",END))



DNA_gr_ROI <- DNA_human[queryHits(findOverlaps(DNA_human, ROI_gr, ignore.strand=T)) %>% unique]

local_clu_size <- DNA_gr_ROI$CBMB %>% table %>% as.data.frame() %>% set_colnames(c("CBMB", "dna_reads"))

saveRDS(DNA_gr_ROI, "/data/Figure-2-fgh-DNA_gr_ROI.rds")
saveRDS(local_clu_size, "/data/Figure-2-fgh-local_clu_size.rds")

DNA_gr_ROI <- readRDS("/data/Figure-2-fgh-DNA_gr_ROI.rds")
local_clu_size <- readRDS("/data/Figure-2-fgh-local_clu_size.rds")

plot_cluster_split <- function(local_clu_size, DNA_gr_ROI, clu_min, clu_max){
  
  filtered_clusters <- local_clu_size %>% filter(dna_reads>=clu_min, dna_reads<=clu_max)
  DNA_ROI_clu_gr <- DNA_gr_ROI %>% filter(CBMB %in% filtered_clusters$CBMB)
  reads_df <- DNA_ROI_clu_gr %>% mutate(frag_mid = (start+end)/2) %>% mcols %>% as.data.frame()
  CBMB_integer <- as.integer(factor(reads_df$CBMB, levels = unique(reads_df$CBMB)))
  reads_df %<>% mutate(CBMB_int = CBMB_integer)
  
  g_clu_split_CBMB <- ggplot(reads_df, aes(x=frag_mid, y=CBMB_int))+
    geom_point(size=0.01, color = "#f0cc45", pch='.')+
    geom_line(aes(group=CBMB), size=0.01, color="#acada2")+
    scale_x_continuous(expand = c(0,0), limits = c(START, END))+
    scale_y_continuous(expand = c(0,0))+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill="white"),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = unit(c(0,0,0,0),"in"),
      axis.ticks = element_blank()
    )
  return(g_clu_split_CBMB)
}

g_clu_2 <- plot_cluster_split(local_clu_size = local_clu_size, DNA_gr_ROI = DNA_gr_ROI, clu_min = 2, clu_max = 2)
g_clu_3_10 <- plot_cluster_split(local_clu_size = local_clu_size, DNA_gr_ROI = DNA_gr_ROI, clu_min = 3, clu_max = 10)


ROI_tiles <- hg38_tiles[queryHits(findOverlaps(hg38_tiles, ROI_gr)) %>% unique()]

# plot triangle using plotgardener engine
AD_chrcontact <- gen_DD_2d_diag_contact(dna_gr = DNA_gr_ROI,
                                        clu_size_df = local_clu_size %>% filter(dna_reads>=2, dna_reads<=10),
                                        tiles = ROI_tiles,
                                        clu_selected=NULL,
                                        return_sparse=T, diag_rm = T) %>%
  set_colnames(c("bin1", "bin2", 'counts')) %>% filter(bin1<=bin2) %>%
  mutate(bin1 = start(ROI_tiles)[bin1], bin2 = start(ROI_tiles)[bin2]) %>%
  set_colnames(c(CHR, CHR, "counts")) %>% as.data.frame()


AD_chrcontact$counts <- log2(AD_chrcontact$counts*1000+1)


params_genetrack <- pgParams(
  chrom = CHR,
  chromstart = START, chromend = END,
  just = c("left", "top"), x=0.2, y0 = -0.1,
  width = 1.8, default.units = "inches",
  assembly = "hg38"
)


png("/Figures/Figure-2-fgh.png",
    width = 2, height = 4, res = 300, units = 'in')

pageCreate(width = 2, height = 4, default.units = "inches")

hicPlot <- plotHicTriangle(
  data = AD_chrcontact, params = params_genetrack,
  palette = colorRampPalette(colors = c(c("#ffffff","#ffc69d","#d4745d","#a42121"))),
  just = c("left", "top"),
  zrange = quantile(AD_chrcontact$counts, c(0.02, 0.98)),
  y = params_genetrack$y0 , colorTrans="linear",
  height = 1.2)

annoGenomeLabel(
  plot = hicPlot, params = params_genetrack,
  scale = "Mb", y = params_genetrack$y0 + 1.25,
  just = c("left", "top"), fontsize = 6
)

plotGG(g_clu_2,
       params = params_genetrack,
       y=params_genetrack$y0+1.5, height = 1.1,
       x=0.18, width = 1.82
       )

annoGenomeLabel(
  plot = hicPlot, params = params_genetrack,
  scale = "Mb", y = params_genetrack$y0 + 2.6,
  just = c("left", "top"),fontsize = 6
)

plotGG(g_clu_3_10,
       params = params_genetrack,
       x=0.18, width = 1.82,
       y=params_genetrack$y0+2.85, height = 1.1,
       )

annoGenomeLabel(
  plot = hicPlot, params = params_genetrack,
  scale = "Mb", y = params_genetrack$y0 + 3.95,
  just = c("left", "top"),fontsize = 6
)

plotText(label = "2D contact map(ensemble MUSIC)",
         fontsize = 8, x=0.25, y=0,
         just = c("left", "top")
)

plotText(label = "Pairwise interactions",
         fontsize = 8, x=0.25, y=params_genetrack$y0 + 1.5,
         just = c("left", "top")
)

plotText(label = "Multiple interactions",
         fontsize = 8, x=0.25, y=params_genetrack$y0 + 2.85,
         just = c("left", "top")
)

pageGuideHide()

dev.off()


