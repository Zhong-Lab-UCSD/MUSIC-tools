# Figure-4.R

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
require(patchwork)
require(wesanderson)
require(ComplexHeatmap)
require(WGCNA)
require(hdWGCNA)
require(Seurat)
require(plotgardener)
require(RColorBrewer)
require(TxDb.Hsapiens.UCSC.hg38.refGene)
require(png)

`%ni%` = Negate(`%in%`)
options(stringsAsFactors=FALSE)

source("MUSIC_utils2.R")


# sample meta info, copy pasted from google sheet: https://docs.google.com/spreadsheets/d/1q7RR6euJ0MTmcdclTPZnpLb8yqBsSNa9/edit#gid=1163861362
sample_info <- fread("brain_sample_info.tab") %>% 
  mutate(lib = paste0("lib", Order)) %>% 
  mutate(expired_age = gsub("≥", "", expired_age) %>% as.numeric) 

brain_list <- paste0("lib",c(1:13, 15))


cell_type_color <- list(#Neu = "#3d76b0", 
  Ex = "#7eb5c6", #L23 = "#8da94b", L5 = "#70b672", L6 = "#72ba96",
  In = "#bebada", #VIP = "#818e6c", PVALB = "#9ca58e", LAMP5 = "#90b061", SST = "#869077",
  Oli = "#f4b461", Ast = "#ee7e72", Opc = "#f7cce5",
  Mic = "#448369", Vas = "#bcde91") 



gtfgr36 <- readRDS("/database/Human/genome/Robj/gtf36_gene_granges.rds")
GOI_gr <- gtfgr36 %>% plyranges::filter(gene_name == "XIST")
# RD contact map split by disease, cell type 

windows.1mb <- gen_windows(1e6, 'hg38') %>% plyranges::filter(seqnames %in% paste0("chr", c(1:22, "X")))


# per cell total RNA
CB_total_rna <- lapply(brain_list, function(LIB){
  
  dna_clu_size <- fread(
    paste0("/Brain_tissue_",LIB,"_deepseq/outputs/6_stats/merge_RNA_human.sort.cell_reads.csv")) %>%
    set_colnames(c("CB","rna_reads")) %>%
    mutate(CB_new = paste0(LIB,"_",CB)) %>% dplyr::select(-CB)
  
}) %>% do.call(rbind, .)

cell_anno <- readRDS("/data/brain_merged_subcluster_seurat_obj.rds")@meta.data %>% 
  rownames_to_column('CB') %>% 
  left_join(CB_total_rna, by=c("CB"="CB_new"))


## ---- Figure-4a SEX ----------

merged_brain <- readRDS("/data/brain_merged_subcluster_seurat_obj.rds")

DefaultAssay(merged_brain) <- "RNA"
merged_brain %<>% NormalizeData() %>% ScaleData
# XIST raw counts in cell
DAT <- #merged_brain@assays$SCT@scale.data["XIST",] %>% as.data.frame() %>%
  merged_brain@assays$RNA@data["XIST",] %>% as.data.frame() %>%
  rownames_to_column('CB') %>% set_colnames(c("CB", "Expr")) %>%
  cbind(., cell_anno %>% dplyr::select(-CB)) %>%
  dplyr::select(CB, Expr, cell_type, type, sex) %>% 
  mutate(cell_type = gsub("Ex","ExN",cell_type)) %>% 
  mutate(cell_type = gsub("In", "InN", cell_type)) %>% 
  mutate(cell_type = factor(cell_type, levels=c("ExN", "InN", "Oli","Ast","Mic","Opc","Vas")))

g_expr_F_M_XIST <- DAT %>%
  ggboxplot(., x = "cell_type", y = "Expr",
            color = "sex", palette = c("#d49299","#27647c"),
            add = c("jitter"),
            ylab="Normalized expression",
            xlab = "Cell type", size=0.2,
            add.params = list(size=0.05, shape='.', alpha=0.2),
            legend = "right",
            title = "XIST RNA levels in single cells") +
  stat_compare_means(aes(group = sex), label = "p.signif", method = "t.test", hide.ns = T) +
  geom_text(data = DAT %>% group_by(cell_type, sex) %>% dplyr::summarize(len = length(Expr)),
            aes(x = factor(cell_type),
                y = min(DAT$Expr) - 0.6,
                label = paste(len,"\n"),
                group = sex),
            position = position_dodge(.8), size=1.8, angle=60)+
  scale_x_discrete(expand = c(0.2,0))+
  theme( plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_rect(fill="white"),
         axis.text = element_text(size=6),
         axis.title = element_text(size=8),
         axis.title.x = element_blank(),
         legend.position = "right",
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
         axis.line= element_line(size = 0.2, linetype = "solid", colour = "black"))

g_expr_F_M_XIST

pdf("/Figures/Figure-4a-XIST-sex-expr.pdf",
    width = 2.2, height = 1.7)
g_expr_F_M_XIST
dev.off()

## ---- Figure 4b RD contact 2D map chrX ----------

# RD heatmap for all cells combined

# pbmcapply::pbmclapply(names(cell_type_color), function(CTYPE){
#   
#   CLU_SIZE = 5000
#   
#   remaining_cells <- cell_anno %>% 
#     filter(cell_type == CTYPE, sex=='F')
#   
#   RNA_gr <- import_grs(CTYPE = CTYPE, mole = "RNA", cond = "all", sex = "F", import_clu = F) %>% filter(CB %in% remaining_cells$CB)
#   DNA_gr <- import_grs(CTYPE = CTYPE, mole = "DNA", cond = "all", sex = "F", import_clu = F) %>% filter(CB %in% remaining_cells$CB)
#   
#   clusize <- import_clusize(CTYPE = CTYPE, mole = "RNA", cond = "all", sex = "F") %>% 
#     filter( (rna_reads+dna_reads) < CLU_SIZE, rna_reads>=1, dna_reads>=1)
#   
#   chrX_RD_mtx <- gen_RD_2d_contact(dna_gr = DNA_gr,
#                                    rna_gr = RNA_gr,
#                                    clu_size_df = clusize,
#                                    DNA_tiles = chrX_win1mb,
#                                    RNA_tiles = chrX_win1mb,
#                                    return_sparse=F, diag_rm = F) 
#   saveRDS(chrX_RD_mtx,
#           paste0("/data/new_5_chrX_RD_DD_mtx/nocutoff_chrX_",CTYPE,"_RD_mtx.rds"))
#   
# }, mc.cores = 10)

chrX_RD_mtx_all <- 
  lapply(names(cell_type_color), function(CTYPE){
    
    RD_mtx <- readRDS(paste0("/data/new_5_chrX_RD_DD_mtx/nocutoff_chrX_",CTYPE,"_RD_mtx.rds"))
    tol_cell <- cell_anno %>% filter(sex=="F", cell_type == CTYPE) %>% nrow()
    avg_mtx <- RD_mtx / tol_cell
    return(avg_mtx)
  }) %>% Reduce('+', .)

g_female <- 
  single_hm_from_full_mtx(chrX_RD_mtx_all, log_values = T, scale_factor = 1000)+
  ggtitle("RNA-chromatin association map\nChrX, female cells")+
  theme(legend.position = 'right',
        plot.title = element_text(hjust=0.5, size=8),
        legend.text = element_text(size=6),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.width =unit(0.1, 'in'),
        legend.key.height =unit(0.1, 'in'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        plot.margin = unit(c(0,0,0,0),'in'))


pdf("/Figures/Figure-4b-RD2Dmap-chrX.pdf",
    width = 1.8,height = 1.8)
g_female
dev.off()


## ---- Figure 4c XIST genome coverage 1D ----------

all_cell_info <- readRDS("/data/new_4_XIST/all_cell_XIST_chrXDNA_tol_reads_info.rds")

CTYPES = c("Ex", "In", "Vas", "Ast", "Opc", "Oli", "Mic")
CLU_SIZE = 5000

CB_total_rna <- lapply(brain_list, function(LIB){
  
  dna_clu_size <- fread(
    paste0("/Brain_tissue_",LIB,"_deepseq/outputs/6_stats/merge_RNA_human.sort.cell_reads.csv")) %>%
    set_colnames(c("CB","rna_reads")) %>%
    mutate(CB_new = paste0(LIB,"_",CB)) %>% dplyr::select(-CB)
  
}) %>% do.call(rbind, .)

cell_anno <- readRDS("/data/brain_merged_subcluster_seurat_obj.rds")@meta.data %>% 
  rownames_to_column('CB') %>% 
  left_join(CB_total_rna, by=c("CB"="CB_new"))


info_df <- lapply(CTYPES, function(CTYPE){
  readRDS(paste0("/Figures/12_1_XIST/cell_type_sc_XIST_RD/",CTYPE,"_",CLU_SIZE,"_info_df.rds"))
}) %>% do.call(rbind,.) %>% 
  left_join(CB_total_rna, by=c("CB"="CB_new")) 


GOI_RD_1d <- lapply(CTYPES, function(CTYPE){
  readRDS(paste0("/Figures/12_1_XIST/cell_type_sc_XIST_RD/",CTYPE,"_",CLU_SIZE,"_RD_1D_sc_mtx.rds"))
}) %>% do.call(rbind,.)

## for XIST zero cells, use zero to make it full
info_df_all <- info_df %>% 
  dplyr::select(CB, chrX_sign_prop, chrX_nonzero, class, CTYPE, chrX_DNA_ct) %>% # notice that the initially chrX_DNA_ct is XIST associated cluster DNA reads count
  rename(XIST_clu_DNA_ct = chrX_DNA_ct) %>% 
  full_join(all_cell_info, by=c("CB"="CB", "CTYPE"="cell_type")) %>% 
  mutate(chrX_sign_prop = replace_na(chrX_sign_prop, 0),
         chrX_nonzero = replace_na(chrX_nonzero, 0),
         class = replace_na(class, "X_depleted"),
         XIST_clu_DNA_ct = replace_na(XIST_clu_DNA_ct, 0)) %>% 
  left_join(cell_anno %>% dplyr::select(CB, type), by=c("CB"="CB"))


empty_RD_matrix <- matrix(0, nrow = nrow(info_df_all)-nrow(GOI_RD_1d), ncol = ncol(GOI_RD_1d))
colnames(empty_RD_matrix) <- colnames(GOI_RD_1d)
rownames(empty_RD_matrix) <- setdiff(info_df_all$CB, rownames(GOI_RD_1d))

GOI_RD_1d_all <- rbind(GOI_RD_1d, empty_RD_matrix)


submtx_cellnum_norm <- colSums(GOI_RD_1d_all)

all_ctype_avg <- data.frame(x = 1:length(submtx_cellnum_norm),
                            y = submtx_cellnum_norm,
                            nz_bin = apply(GOI_RD_1d_all,2,function(x){sum(x>0)}))

# sum RAL (ensemble RAL)
g_bar <- ggbarplot(all_ctype_avg, x="x", y="y",
                   xlab="genomic bins", ylab = "RAL",
                   title = "Genome-wide distribution of chromatin-associated XIST RNA in female",
                   fill = "#b7726b", color = "#b7726b")+
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color=NA),
    panel.background = element_rect(fill=NA),
    axis.text = element_text(size=6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=8),
    axis.title.x = element_blank(),
    legend.position = "right",
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


chr_colors <- data.frame(chr = paste0("chr",c(1:22,"X")),
                         color = Seurat::DiscretePalette(23, palette = "stepped"))
bin_colors_vec <- windows.1mb %>% as.data.frame() %>% left_join(chr_colors, by=c("seqnames"="chr")) %>% .[["color"]]

# each chr middle point loc 
chr_mid_point <- as.integer((seqnames(windows.1mb) %>% table)/2)[1:23]
chr_pre_size <- (c(0, seqnames(windows.1mb) %>% table) %>% cumsum())[1:23]
chr_name_marker <- chr_mid_point+chr_pre_size

g_chr <- 
  ggplot(data.frame(x = seq(1, length(submtx_cellnum_norm)), y = 1))+
  geom_col(aes(x = x, y=y), fill = bin_colors_vec, color=bin_colors_vec)+
  xlab("Chromosomes")+
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits = c(-0.8, 1))+
  annotate("text", x=chr_name_marker, label=c(1:22,"X"),y=(-0.8),size=1.8)+
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=8),
        panel.border = element_rect(fill = NA, color=NA),
        panel.background = element_rect(fill=NA),
        plot.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
g_chr

pdf("/Figures/Figure-4c-XIST-RD-1D-genome.pdf",
    width = 4.4, height = 1)
# pdf("/Figures/12_1_XIST/Figure_allcelltype_XIST_genome_RD_barplot.pdf",
#     width = 10, height = 1.6)
g_bar %>% aplot::insert_bottom(g_chr, height = 0.5) 
dev.off()


## ---- Figure 4d XIST+ XIST- contrast DD contact ----------

# remain_DNA_gr <- readRDS("/data/new_4_XIST/sc_GDCF/RNAreads_filtered_remain_sc_chrX_DNA_gr.rds")
# 
# XIST_linked_DNA_gr <- remain_DNA_gr %>% filter(CBMB %in% XIST_RNA_gr$CBMB)
# 
# XIST_no_linked_DNA_gr <- remain_DNA_gr %>% filter(CBMB %ni% XIST_RNA_gr$CBMB)
# 
# windows.1mb <- gen_windows(1e6, 'hg38') %>% plyranges::filter(seqnames %in% paste0("chr", c(1:22, "X")))
# 
# chrX_win1mb <- windows.1mb %>% plyranges::filter(seqnames=="chrX")
# 
# # remaining_cells <- unique(c(remain_DNA_gr$CB, remain_RNA_gr$CB))
# # 
# # clusize_allctype <- pbmcapply::pbmclapply(names(cell_type_color), function(CTYPE){
# #   clusize <- import_clusize(CTYPE = CTYPE, mole = "RNA", cond = "all", sex = "F") %>% 
# #     mutate(CB = gsub("^(lib.*BC3.*BC2.*BC1.*)-[ATCGN]{16}-lib.*$", "\\1", CBMB)) %>% 
# #     filter(CB %in% remaining_cells)
# #   return(clusize)
# # }, mc.cores=10) %>% do.call(rbind,.)
# # saveRDS(clusize_allctype, "/data/new_4_XIST/sc_GDCF/remaining_cells_global_clusize.rds")
# clusize_allctype <- readRDS("/data/new_4_XIST/sc_GDCF/remaining_cells_global_clusize.rds")
# 
# 
# DD_Xlinked_mtx <- gen_DD_2d_diag_contact(dna_gr = XIST_linked_DNA_gr,
#                                          clu_size_df = clusize_allctype %>% filter(dna_reads<1000),
#                                          tiles = chrX_win1mb,
#                                          return_sparse = F,diag_rm = T)
# 
# DD_no_Xlinked_mtx <- gen_DD_2d_diag_contact(dna_gr = XIST_no_linked_DNA_gr,
#                                             clu_size_df = clusize_allctype %>% filter(dna_reads<1000),
#                                             tiles = chrX_win1mb,
#                                             return_sparse = F,diag_rm = T)
# saveRDS(DD_Xlinked_mtx, "/data/new_4_XIST/All_XIST_linked_DD_mtx.rds")
# saveRDS(DD_no_Xlinked_mtx, "/data/new_4_XIST/All_non_XIST_linked_DD_mtx.rds")

DD_Xlinked_mtx <- readRDS("/data/new_4_XIST/All_XIST_linked_DD_mtx.rds")
DD_no_Xlinked_mtx <- readRDS("/data/new_4_XIST/All_non_XIST_linked_DD_mtx.rds")

DD_Xlinked_mtx_nrom <- DD_Xlinked_mtx/max(DD_Xlinked_mtx)
DD_no_Xlinked_mtx_nrom <- DD_no_Xlinked_mtx/max(DD_no_Xlinked_mtx)


linked_contrast <- log2(DD_Xlinked_mtx_nrom/DD_no_Xlinked_mtx_nrom)
linked_contrast[is.infinite(linked_contrast)] <- NA

g_DD_Xlinked_mtx <- single_hm_from_full_mtx(DD_Xlinked_mtx_nrom,
                                            legend_position = "bottom", scale_factor = 10, log_values = T) +
  ggtitle("XIST+  chromatin")+
  theme(legend.position = 'None',
        plot.title = element_text(hjust=0.5, size=7),
        legend.text = element_text(size=6),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.width =unit(0.1, 'in'),
        legend.key.height =unit(0.1, 'in'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        plot.margin = unit(c(0,0,0,0),'in'))


g_DD_nolinked_chrX <- single_hm_from_full_mtx(DD_no_Xlinked_mtx_nrom,
                                              legend_position = "bottom", scale_factor = 10, log_values = T) +
  ggtitle("XIST-  chromatin")+
  theme(legend.position = 'None',
        plot.title = element_text(hjust=0.5, size=7),
        legend.text = element_text(size=6),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.width =unit(0.1, 'in'),
        legend.key.height =unit(0.1, 'in'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        plot.margin = unit(c(0,0,0,0),'in'))


g_contrast <- single_hm_from_full_mtx(linked_contrast-mean(linked_contrast, na.rm=T),
                                      palette = c("#476ba8", "#f1f7e1", "#c9462e"),
                                      mybreaks = c(-3,0,3),
                                      legend_position = "bottom", scale_factor = 1, log_values = F) +
  ggtitle("XIST+ vs. XIST-")+
  theme(legend.position = 'None',
        plot.title = element_text(hjust=0.5, size=7),
        legend.text = element_text(size=6),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.width =unit(0.1, 'in'),
        legend.key.height =unit(0.1, 'in'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        plot.margin = unit(c(0,0,0,0),'in'))



pdf("/Figures/Figure-4d-XIST-linked-nolinked-contrast-DD.pdf",
    width = 2.8, height = 1.2)
(g_DD_Xlinked_mtx|g_DD_nolinked_chrX|g_contrast)
dev.off()

# the legend need to add manually (plot with legend first)


## ---- Figure 4g XIST linked contrast curve ----------

# import genomic bins for GDCF
hg38_arm_size <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
                       col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
hg38_arm_size %>%
  mutate(width = chromEnd - chromStart) %>%
  mutate(arm = substring(name, 1,1)) %>%
  group_by(chrom, arm) %>% summarize(size = sum(width)) %>%
  arrange(size) %>% filter(chrom %in% c(paste0("chr",1:22),"chrX"))
# longest chrom arm is 125mbp.
dist_bins <- seq(10, 1.5*10e8, length.out = 50000)

# function 
one_chr_DNAgr_dist_df <- function(DNA_gr){
  
  onechr <- DNA_gr %>% as.data.frame() %>% 
    filter(seqnames == "chrX") %>% 
    mutate(ID = seq(1, nrow(.))) %>% dplyr::select(CB, start, CBMB, ID)
  
  clusize_dist <- onechr %>%
    left_join(onechr %>% dplyr::select(-CB), by=c("CBMB" = "CBMB"), multiple = "all") %>%
    filter(ID.x < ID.y) %>%
    mutate(dist = abs(start.x - start.y)) %>%
    dplyr::filter(dist>10) %>%
    mutate(dist_cut = cut(dist, dist_bins)) %>%
    dplyr::group_by(dist_cut) %>%
    dplyr::summarize(ct = length(ID.x)) %>% ungroup
  
  clusize_dist %<>% mutate(freq_ratio = ct/sum(ct))
  
  return(clusize_dist)
  
}

remaining_cells <- info_df_all %>% filter(tol_RNA>=5000) %>% .[['CB']]



remain_DNA_gr <- readRDS("/data/new_4_XIST/sc_GDCF/RNAreads_filtered_remain_sc_chrX_DNA_gr.rds")

remain_RNA_gr <- pbmcapply::pbmclapply( c("Ex", "In", "Oli", "Mic", "Vas", "Opc", "Ast"), function(CTYPE){
  
  RNA_gr <- import_grs(CTYPE = CTYPE, mole = "RNA", cond = "all", sex="F", import_clu = F) %>% 
    filter(CB %in% remaining_cells) 
  
  return(RNA_gr)
  
}, mc.cores = 10) %>% do.call(c,.)


# remaining cells, no XIST DNA reads
DNA_gr_nonXIST <- remain_DNA_gr %>% filter(CB %in% (info_df_all %>% filter(CB %in% remaining_cells, chrX_nonzero == 0) %>% .[["CB"]]))
DNA_gr_nonXIST_df <- one_chr_DNAgr_dist_df(DNA_gr_nonXIST) %>% mutate(type="no_coat_DNA")


# remaining cells, no XIST, but involved in RD clusters' DNA
gtfgr36 <- readRDS("/database/Human/genome/Robj/gtf36_gene_granges.rds")
GOI_gr <- gtfgr36 %>% plyranges::filter(gene_name == "XIST")
XIST_RNA_gr <- remain_RNA_gr[queryHits(findOverlaps(remain_RNA_gr, GOI_gr, ignore.strand=F))] # XIST RNA grs and clusters

DNA_gr_nonXIST_in_RD_clu <- DNA_gr_nonXIST %>% filter(CBMB %in% remain_RNA_gr$CBMB) %>% filter(CBMB %ni% XIST_RNA_gr$CBMB)
DNA_gr_nonXIST_in_RD_clu_df <- one_chr_DNAgr_dist_df(DNA_gr_nonXIST_in_RD_clu) %>% mutate(type="no_coat_RNA_asso_DNA")


# remaining cells, split into chunck, calculate XIST linked and non-linked separately

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

N_chunk_cells <- chunk2((info_df_all %>% 
                           arrange(desc(chrX_nonzero)) %>% 
                           filter(CB %in% remaining_cells, chrX_nonzero > 0))$CB,
                        3) # split into 3 chuncks 


chuncks_df_all <- lapply(1:length(N_chunk_cells), function(idx){
  
  chunck = N_chunk_cells[[idx]]
  
  XIST_linked_DNA_gr <- remain_DNA_gr %>% filter(CB %in% chunck, CBMB %in% XIST_RNA_gr$CBMB)
  XIST_no_linked_DNA_gr <- remain_DNA_gr %>% filter(CB %in% chunck, CBMB %ni% XIST_RNA_gr$CBMB)
  
  XIST_linked_DNA_gr_df <- one_chr_DNAgr_dist_df(XIST_linked_DNA_gr) %>% mutate(type = paste0("Coated_group",idx,"_XIST_linked"))
  XIST_no_linked_DNA_gr_df <- one_chr_DNAgr_dist_df(XIST_no_linked_DNA_gr) %>% mutate(type = paste0("Coated_group",idx,"_XIST_notlinked"))
  
  return(rbind(XIST_linked_DNA_gr_df, XIST_no_linked_DNA_gr_df))
}) %>% do.call(rbind,.)


# combine all and plot
chuncks_group <- chuncks_df_all %>% mutate(group = gsub("Coated_(group.*)_XIST.*$", '\\1', type))
saveRDS(chuncks_group, "/data/Figure-4-XIST-GDCF-curve-chuncks_group.rds")

non_coated_group <- rbind(DNA_gr_nonXIST_df, DNA_gr_nonXIST_in_RD_clu_df ) %>% mutate(group="non_coated")
saveRDS(non_coated_group, "/data/Figure-4-XIST-GDCF-curve-non_coated_group.rds")

all_no_XIST_linked <- rbind(chuncks_group, non_coated_group) %>% 
  filter(grepl("no_coat_DNA|notlinked",type)) %>% 
  group_by(dist_cut) %>% summarize(ct = sum(ct)) %>% 
  mutate(freq_ratio = ct / sum(ct)) %>% mutate(type = 'no_XIST_linked')

XIST_linked_onecurve <- rbind(chuncks_group, non_coated_group) %>% 
  filter(grepl("_XIST_linked",type)) %>% dplyr::select(-group) %>% 
  group_by(dist_cut) %>% summarize(ct = sum(ct)) %>% 
  mutate(freq_ratio = ct/sum(ct)) %>% mutate(type = "XIST_linked")

otherRNA_linked <- non_coated_group %>% 
  filter(type == "no_coat_RNA_asso_DNA") %>% dplyr::select(-group) %>% 
  mutate(type = "otherRNA_linked")

# rbind(all_no_XIST_linked, XIST_linked_onecurve, otherRNA_linked)  %>% 
g_xlinked_nolinked <- 
  rbind(all_no_XIST_linked, XIST_linked_onecurve)  %>% 
  # mutate(type = factor(type, levels = c("enriched_coat", "not_enriched", "male"))) %>%
  mutate(dist_end = gsub("^.*,", "", dist_cut)) %>%
  mutate(dist_end = gsub("]","", dist_end)) %>%
  mutate(dist_end = as.numeric(dist_end)) %>%
  ggplot() +
  geom_line(aes(x=dist_end, y=freq_ratio,color = type),
            size=0.3, alpha = 0.7) +
  scale_x_log10(labels=scales::trans_format('log10',scales::math_format(10^.x)),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(4.5e4, 1*1e8)) +
  scale_y_log10(labels=scales::trans_format('log10',scales::math_format(10^.x)),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(10^(-4), 10^(-3)),
                expand = c(0, 0)) +
  annotation_logticks(size=0.2,
                      short = unit(0.01, "in"),
                      mid = unit(0.02, "in"),
                      long = unit(0.03, "in")) +
  ggsci::scale_color_jama()+
  labs(x = expression(paste("Genomic distance, ",italic("s")," (bp)")),
       y = expression(paste("Contact probability, ",italic("Pc"))))+
  ggtitle("chrX chromatin interactions")+
  annotate("text", x=10^(5.1), y=10^(-3.4), label="XIST+", color="#b1895e", size=2.5)+
  annotate("text", x=10^(6), y=10^(-3.2), label="XIST-", color="#6a8593", size=2.5)+
  theme(
    legend.position = "None",
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill="white"),
    axis.text = element_text(size=5),
    axis.title = element_text(size=7),
    plot.title = element_text(size=7, hjust = 0.5),
    legend.text = element_text(size=6),
    legend.title = element_text(size=7, hjust=0.5),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0),
    legend.key.width = unit(0.1,"in"),
    legend.key.height = unit(0.1,"in"),
    legend.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "in"),
    axis.ticks = element_line(size=0.2),
    axis.line= element_line(size = 0.2, linetype = "solid", colour = "black")
  )

pdf("/Figures/Figure-4g-XIST-linked-not-2curves.pdf",
    width = 1.3, height = 1.3)
g_xlinked_nolinked
dev.off()


# chisquare test test if otherRNA-linked and noXIST-linked are statistically different
chisq2_merged_df <- XIST_linked_onecurve %>% left_join(all_no_XIST_linked, by=c("dist_cut"="dist_cut"))
chisq.test(data.frame(one = chisq2_merged_df$ct.x, two = chisq2_merged_df$ct.y) %>% as.matrix())

## ---- Figure 4h GDCF XIST 4 groups ----------

chuncks_group <- readRDS("/data/Figure-4-XIST-GDCF-curve-chuncks_group.rds")
non_coated_group <- readRDS("/data/Figure-4-XIST-GDCF-curve-non_coated_group.rds")


g_GDCF_groups <- rbind(chuncks_group, non_coated_group)  %>% 
  mutate(group = factor(group, levels = c("non_coated","group3","group2","group1"), 
                        labels=c("XAL=0, Group1","Low XAL, Group2","Median XAL, Group3","High XAL, Group4"))) %>%
  mutate(dist_end = gsub("^.*,", "", dist_cut)) %>%
  mutate(dist_end = gsub("]","", dist_end)) %>%
  mutate(dist_end = as.numeric(dist_end)) %>%
  ggplot() +
  geom_line(aes(x=dist_end, y=freq_ratio,color = type),
            size=0.3, alpha = 0.7) +
  scale_x_log10(labels=scales::trans_format('log10',scales::math_format(10^.x)),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(4.5e4, 1*1e8)) +
  scale_y_log10(labels=scales::trans_format('log10',scales::math_format(10^.x)),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(10^(-4.5), 10^(-2.5)),
                expand = c(0, 0)) +
  annotation_logticks(size=0.2,
                      short = unit(0.01, "in"),
                      mid = unit(0.02, "in"),
                      long = unit(0.03, "in")) +
  facet_wrap(~group, nrow=1)+
  scale_color_manual(name = "chrX-XIST enrich type",
                     values = c("#dd5a50","grey70","#e5a246","grey70","#c57eb9","grey70","#3a8076","grey70"),
                     labels = c("XIST+","XIST-","XIST+","XIST-","XIST+","XIST-","Any_RNA+","Any_RNA-"))+
  labs(x = expression(paste("Genomic distance, ",italic("s")," (bp)")),
       y = expression(paste("Contact probability, ",italic("Pc"))))+
  ggtitle("chrX, chromatin interactions")+
  theme(
    legend.position = "None",
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill="white"),
    axis.text = element_text(size=5),
    axis.title = element_text(size=7),
    plot.title = element_text(size=7, hjust = 0.5),
    legend.text = element_text(size=6),
    legend.title = element_text(size=7, hjust=0.5),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0),
    legend.key.width = unit(0.1,"in"),
    legend.key.height = unit(0.1,"in"),
    legend.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "in"),
    axis.ticks = element_line(size=0.2),
    axis.line= element_line(size = 0.2, linetype = "solid", colour = "black"),
    strip.text.x = element_text(size = 6,margin = unit(c(0.02,0,0.02,0), "in"))
  )


pdf("/Figures/Figure-4h-4groups-GDCF-curves.pdf",
    width = 4.2, height = 1.4)
g_GDCF_groups
dev.off()


## ---- Figure 4i ExN sc XIST RD heatmap ----------

CTYPE = "Ex"
CLU_SIZE = 5000
cutoff = 1000

Ex_info_df <- info_df_all %>% filter(CTYPE == "Ex") %>% filter(tol_RNA>=cutoff, chrX_DNA>1000)
Ex_AD_GOI_RD_1d <- GOI_RD_1d_all[Ex_info_df$CB,]


Ex_row_orders <- Ex_info_df %>% arrange(desc(num_XIST_RD_clu)) %>% .[["CB"]]

Ex_col_anno <- HeatmapAnnotation(chr=as.character(seqnames(windows.1mb)),
                                 `cRAL` = anno_barplot(apply(Ex_AD_GOI_RD_1d, 2, function(x){sum(x)}),
                                                       gp=gpar(col="#b7726b",fill="#b7726b"),
                                                       height = unit(0.15, "in"),
                                                       axis_param=list(gp=gpar(fontsize = 5))),
                                 col = list(chr = setNames(chr_colors$color, chr_colors$chr)),
                                 annotation_name_side = c("left",'right'),
                                 annotation_name_gp = gpar(fontsize=6),
                                 simple_anno_size = unit(0.1, "in"),
                                 gap = unit(0.01, "in"),
                                 annotation_name_rot = 0,
                                 show_legend = F) # add chr marker

Ex_row_anno <- rowAnnotation(`XIST RNA` = anno_barplot(setNames(Ex_info_df$XIST_ct, Ex_info_df$CB)[Ex_row_orders],
                                                       gp=gpar(fill="#dec790",col="#dec790"),
                                                       width = unit(0.15, "in"),
                                                       axis_param=list(gp=gpar(fontsize = 5))),
                             `log2(chrX DNA)` = anno_barplot(setNames(log2(Ex_info_df$chrX_DNA+1), Ex_info_df$CB)[Ex_row_orders],
                                                             gp=gpar(fill="#cedfc1",col="#cedfc1"),
                                                             width = unit(0.15, "in"),
                                                             axis_param=list(gp=gpar(fontsize = 5))),
                             annotation_name_rot = 30,
                             annotation_name_side = "top",
                             annotation_name_gp = gpar(fontsize=6),
                             col = list(XIST_ct = "grey70",
                                        X_DNA_log = "#cfb39d"),
                             gap = unit(0.01, "in"))


pdf("/Figures/Figure-4h-Ex-sc-hm-XIST-RAL.pdf",
    width = 3.6, height = 2.6)
draw(
  ComplexHeatmap::Heatmap(Ex_AD_GOI_RD_1d[Ex_row_orders,],
                          name="RAL",
                          cluster_rows = F, cluster_columns = F,
                          show_row_names = F, show_column_names = F,
                          bottom_annotation = Ex_col_anno,
                          right_annotation = Ex_row_anno,
                          row_title = "Individual ExNs",
                          column_split = factor(seqnames(windows.1mb) %>% as.character() %>% gsub("chr","",.), levels=c(1:22,"X")),
                          column_gap = unit(0.2, "mm"),
                          use_raster = T,
                          raster_quality=10, raster_by_magick=TRUE,
                          column_title_gp = gpar(fontsize=5),
                          row_title_gp = gpar(fontsize=7),
                          col = circlize::colorRamp2(c(0, 0.0005, 0.001, 0.003),
                                                     c("#4192bf", "#d1e5ef", "white", "#dc5537")),
                          heatmap_legend_param = list(title = "RAL",
                                                      title_gp = gpar(fontsize=6),
                                                      legend_width = unit(0.02, "in"),
                                                      legend_height = unit(0.8, 'in'),
                                                      labels_gp=gpar(fontsize=5))
  ),
  padding = unit(c(0.1, 0.1, 0.1, 0.1), "in")
)
dev.off()

## ---- Figure 4j ExN, InN, Ast GDCF curve ----------

celltype_GDCF_df <- readRDS("/data/new_4_XIST/Ex_In_Ast_XIST_linked_GDCF.rds")

g_cell_type_xlinked_nolinked <- 
  celltype_GDCF_df %>% mutate(cell_type = gsub("()_XIST.*$","\\1", type)) %>% 
  filter(cell_type %in% c("Ex", "In", "Ast")) %>%
  mutate(cell_type = factor(cell_type, levels = c("Ex", "In", "Ast"), labels=c("ExN","InN","Ast"))) %>%
  mutate(dist_end = gsub("^.*,", "", dist_cut)) %>%
  mutate(dist_end = gsub("]","", dist_end)) %>%
  mutate(dist_end = as.numeric(dist_end)) %>%
  ggplot() +
  geom_line(aes(x=dist_end, y=freq_ratio,color = type),
            size=0.3, alpha = 0.7) +
  scale_x_log10(labels=scales::trans_format('log10',scales::math_format(10^.x)),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(4.5e4, 1.5*1e8)) +
  scale_y_log10(labels=scales::trans_format('log10',scales::math_format(10^.x)),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(10^(-5), 10^(-3)),
                expand = c(0, 0)) +
  facet_wrap(~cell_type, nrow=1)+
  annotation_logticks(size=0.2,
                      short = unit(0.01, "in"),
                      mid = unit(0.02, "in"),
                      long = unit(0.03, "in")) +
  scale_color_manual(values=c("#e27262","#ebb6af", "#7eb5c6","#abbdc4", "#928dcc","#d0ceef"))+
  # scale_color_manual(values = c("#7eb5c6", "#CCCCCC", "#bebada", "#CCCCCC", "#f4b461", "#CCCCCC", "#ee7e72", "#CCCCCC", "#f7cce5", "#CCCCCC", "#448369", "#CCCCCC", "#bcde91", "#CCCCCC"))+
  labs(x = expression(paste("Genomic distance, ",italic("s")," (bp)")),
       y = expression(paste("Contact probability, ",italic("Pc"))))+
  ggtitle("chrX chromatin interactions")+
  theme(
    legend.position = "None",
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill="white"),
    axis.text = element_text(size=5),
    axis.title = element_text(size=7),
    plot.title = element_text(size=7, hjust = 0.5),
    legend.text = element_text(size=6),
    legend.title = element_text(size=7, hjust=0.5),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0),
    legend.key.width = unit(0.1,"in"),
    legend.key.height = unit(0.1,"in"),
    legend.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "in"),
    axis.ticks = element_line(size=0.2),
    axis.line= element_line(size = 0.2, linetype = "solid", colour = "black"),
    strip.text.x = element_text(size = 6,margin = unit(c(0.02,0,0.02,0), "in"))
  )

pdf("/Figures/Figure-4j-ExN-InN-Ast-XIST-curve.pdf",
    width = 2.9, height = 1.3)
g_cell_type_xlinked_nolinked
dev.off()



## ---- Figure 4k ExN InN Ast contact contrast matrix ----------


pbmcapply::pbmclapply(c("Ex", "In", "Ast"), function(CTYPE){
  
  CLU_SIZE = 5000
  cutoff = 5000
  
  remaining_cells <- cell_anno %>%
    # left_join(CB_total_rna, by=c("CB"="CB_new")) %>%
    dplyr::filter(rna_reads >= cutoff) %>%
    dplyr::filter(cell_type == CTYPE)
  
  RNA_gr <- import_grs(CTYPE = CTYPE, mole = "RNA", cond = "all", sex = "F", import_clu = F) %>% filter(CB %in% remaining_cells$CB)
  
  DNA_gr <- import_grs(CTYPE = CTYPE, mole = "DNA", cond = "all", sex = "F", import_clu = F) %>% filter(CB %in% remaining_cells$CB)
  
  GOI_gr <- gtfgr36 %>% plyranges::filter(gene_name == "XIST")
  
  XIST_clu <- RNA_gr[queryHits(findOverlaps(RNA_gr, GOI_gr)) %>% unique]$CBMB %>% unique
  
  clusize <- import_clusize(CTYPE = CTYPE, mole = 'DNA', cond = "all", sex = "F") %>% 
    filter( (rna_reads+dna_reads) < CLU_SIZE, dna_reads>=2)
  
  XIST_clusize <- clusize %>% filter(CBMB %in% XIST_clu)
  not_XIST_clusize <- clusize %>% filter(CBMB %ni% XIST_clu)
  
  XIST_linked_DNA_gr <- DNA_gr %>% filter(CBMB %in% XIST_clu) %>% filter(seqnames == "chrX")
  XIST_notlinked_DNA_gr <- DNA_gr %>% filter(CBMB %ni% XIST_clu) %>% filter(seqnames == "chrX")
  
  
  Xlinked_chrX_DD_mtx <- gen_DD_2d_diag_contact(dna_gr = XIST_linked_DNA_gr,
                                                clu_size_df = XIST_clusize,
                                                tiles = chrX_win1mb,
                                                return_sparse=F, diag_rm = T) 
  
  
  
  saveRDS(Xlinked_chrX_DD_mtx,
          paste0("/data/new_5_chrX_RD_DD_mtx/cutoff5000_",
                 CTYPE, "_Xlinked_chrX_DD_mtx.rds"))
  # Xlinked_chrX_DD_mtx_norm <- Xlinked_chrX_DD_mtx/(length(XIST_linked_DNA_gr)/1e6)
  
  
  non_Xlinked_chrX_DD_mtx <- gen_DD_2d_diag_contact(dna_gr = XIST_notlinked_DNA_gr,
                                                    clu_size_df = not_XIST_clusize,
                                                    tiles = chrX_win1mb,
                                                    return_sparse=F, diag_rm = T) 
  saveRDS(non_Xlinked_chrX_DD_mtx,
          paste0("/data/new_5_chrX_RD_DD_mtx/cutoff5000_",
                 CTYPE, "_noXlinked_chrX_DD_mtx.rds"))
  
}, mc.cores = 3)


three_ctypes_contrast_matrix <- 
  lapply(c("Ex", "In", "Ast"), function(CTYPE){
    
    # CTYPE = "In"
    CTYPE_Xlinked_DD <- readRDS(paste0("/data/new_5_chrX_RD_DD_mtx/cutoff5000_",CTYPE,"_Xlinked_chrX_DD_mtx.rds"))
    CTYPE_noXlinked_DD <- readRDS(paste0("/data/new_5_chrX_RD_DD_mtx/cutoff5000_",CTYPE,"_noXlinked_chrX_DD_mtx.rds"))
    
    DD_Xlinked_mtx_nrom <- CTYPE_Xlinked_DD/quantile(CTYPE_Xlinked_DD, 0.95)
    DD_no_Xlinked_mtx_nrom <- CTYPE_noXlinked_DD/quantile(CTYPE_noXlinked_DD, 0.95)
    
    linked_contrast <- log2(DD_Xlinked_mtx_nrom/DD_no_Xlinked_mtx_nrom)
    linked_contrast[is.infinite(linked_contrast)] <- NA
    
    # 
    # g_DD_Xlinked_mtx <- single_hm_from_full_mtx(DD_Xlinked_mtx_nrom,
    #                                             legend_position = "bottom", scale_factor = 1, log_values = T) +
    #   ggtitle("XIST linked DD contact")
    # 
    # 
    # g_DD_nolinked_chrX <- single_hm_from_full_mtx(DD_no_Xlinked_mtx_nrom,
    #                                               legend_position = "bottom", scale_factor = 1, log_values = T) +
    #   ggtitle("No XIST linked DD contact")
    # 
    # 
    g_contrast <- single_hm_from_full_mtx(linked_contrast,
                                          palette = c("#476ba8", "#f1f7e1", "#c9462e"),
                                          mybreaks = c(-4,0,4),
                                          legend_position = "bottom", scale_factor = 1, log_values = F) +
      ggtitle(CTYPE) +
      theme(
        legend.position = 'right',
        plot.title = element_text(hjust=0.5, size=7),
        legend.text = element_text(size=6),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.key.width =unit(0.1, 'in'),
        legend.key.height =unit(0.1, 'in'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        plot.margin = unit(c(0,0,0,0),'in')
      )
    
    return(g_contrast)
    
  })

pdf("/Figures/Figure-4k-ExN-InN-Ast-XIST-contrast-DD.pdf",
    width = 2.9, height = 1.3)
three_ctypes_contrast_matrix[[1]] %>% aplot::insert_right(three_ctypes_contrast_matrix[[2]]) %>% 
  aplot::insert_right(three_ctypes_contrast_matrix[[3]])
dev.off()