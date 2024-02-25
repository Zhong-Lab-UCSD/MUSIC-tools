# Figure-3.R
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
require(Seurat)
require(harmony)
`%ni%` = Negate(`%in%`)
options(stringsAsFactors=FALSE)
options(scipen=3)

theme_no_axis <- theme(axis.line = element_blank(),
                       axis.ticks = element_blank(),
                       axis.text = element_blank(),
                       axis.title = element_blank())

source("MUSIC_utils2.R")
gtfgr36 <- readRDS("/database/Human/genome/Robj/gtf36_gene_granges.rds")

# sample meta info, copy pasted from google sheet: https://docs.google.com/spreadsheets/d/1q7RR6euJ0MTmcdclTPZnpLb8yqBsSNa9/edit#gid=1163861362
sample_info <- fread("brain_sample_info.tab") %>% 
  mutate(lib = paste0("lib", Order))

brain_list <- paste0("lib",c(1:13, 15))


cell_type_color <- list(#Neu = "#3d76b0", 
  Ex = "#7eb5c6", #L23 = "#8da94b", L5 = "#70b672", L6 = "#72ba96",
  In = "#bebada", #VIP = "#818e6c", PVALB = "#9ca58e", LAMP5 = "#90b061", SST = "#869077",
  Oli = "#f4b461", Ast = "#ee7e72", Opc = "#f7cce5",
  Mic = "#448369", Vas = "#bcde91") 

## ---- Figure 3a, brain UMAP plot ----------

brains <- readRDS("/data/brain_merged_subcluster_seurat_obj.rds")

# load brain markers
Ex <- c("SLC17A7","CAMK2A","NRGN")
In <- c("GAD1","GAD2")
Ast <- c("AQP4", "GFAP")
Mic <- c("CD74", "CSF1R", "C3")
OPC <- c("PDGFRA","VCAN", "CSPG4")
Oli <- c("MBP", "PLP1","MOBP")
Vas <- c("FLT1",  "EMCN", "VWF")


marker_list <- list(
  Ex = Ex, In = In, 
  Oli = Oli, Ast = Ast, Opc = OPC, Mic = Mic, Vas = Vas
)

# prepare all sub celltype color

sub_celltype_color <- list(
  Ast = "#ee7e72", 
  Opc = "#f7cce5",
  Vas = "#bcde91",
  Mic_MS4A = "#7ac4ac",
  Mic = "#448369",
  OLI1 = "#dfbb3f",
  OLI2 = "#cb8a31",
  OLI3 = "#f5e24d",
  OLI4 = "#673d11",
  Pvalb = "#933861",
  Sst = "#c67eb8",
  Lamp5 = "#9598bf",
  Vip = "#3d4383",
  `L2/3_IT` = "#193d62",
  L5_IT = "#5f869e",
  `L5/6_NP` = "#c2f3fc",
  L6_CT = "#7eb5c6",
  L6_IT = "#82a9fa"    
) 

clabel <- brains@meta.data$sub_celltype %>% table %>% as.data.frame() %>% 
  set_names(c("ctype", "n")) %>% mutate(label = paste0(ctype," (n=",n,")")) 

cl <- setNames(clabel$label, clabel$ctype)

# 1. Cluster plot
as <- -15 # arrow start
al <- 3 # arrow length

g_clu <- DimPlot(brains, reduction = "umap", pt.size = 1,
                 group.by = "sub_celltype",
                 cols = sub_celltype_color,
                 label = T, label.size = 2, repel = T, raster = T) + 
  theme_no_axis + theme(plot.title = element_blank(),
                        aspect.ratio = 1,
                        plot.margin = unit(c(0,0,0,0), "in"),
                        legend.position = "None")+
  annotate("segment", x = as, xend = as+al, y = as, yend = as,size=0.2,
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  annotate("segment", x = as, xend = as, y = as, yend = as+al,size=0.2,
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  annotate("text", x = as+al/2, y=as-0.8, label="UMAP1",size=2)+
  annotate("text", x = as-0.8, y=as+al/2, label="UMAP2", angle=90,size=2) + 
  scale_color_manual(breaks = names(sub_celltype_color),
                     values = sub_celltype_color,
                     labels = cl)
pdf("/Figures/Figure-3a-UMAP-brain-single-cell-merged.pdf",
    width = 3, height = 3)
g_clu
dev.off()


## ---- Figure-3b cell type GDCF curves ----------
# export_clusize_dist <- function(cell_type, sex, cond){
#   
#   opc_clusize <- readRDS(paste0("/data/2_celltype_sex_disease_DNA_RNA_gr_clusize/",
#                             cell_type, "_", sex, "_", cond, "_clusize.rds")) %>% 
#     filter(dna_reads<=100 & dna_reads>=2) %>% mutate(CBMB = paste0(lib,"_",CBMB))
#   
#   opc_dna <- readRDS(paste0("/data/2_celltype_sex_disease_DNA_RNA_gr_clusize/",
#                             cell_type, "_", sex, "_", cond, "_DNA_gr.rds")) %>% filter(seqnames %ni% c("chrY", "chrM")) %>% 
#     mutate(CBMB = paste0(lib, "_", CBMB)) %>% filter(CBMB %in% opc_clusize$CBMB)
#   seqlevels(opc_dna) <- paste0("chr", c(1:22,"X"))
#   
#   all_chr <- split(opc_dna, seqnames(opc_dna)) %>% pbmcapply::pbmclapply(., function(x){
#     
#     print(head(x))
#     onechr <- x %>% as.data.frame() %>% mutate(CB = paste0(lib, "_", CB)) %>% mutate(ID = seq(1, nrow(.)))
#     
#     clusize_dist <- onechr %>%  
#       left_join(onechr %>% dplyr::select(start, CBMB, ID), by=c("CBMB" = "CBMB"),multiple = "all") %>% 
#       filter(ID.x < ID.y) %>% 
#       mutate(dist = abs(start.x - start.y)) %>%
#       dplyr::filter(dist>10) %>%
#       mutate(dist_cut = cut(dist, dist_bins)) %>%
#       dplyr::group_by(CB, CBMB, dist_cut) %>% 
#       dplyr::summarize(ct = length(ID.x)) %>% ungroup %>% 
#       left_join(opc_clusize %>% dplyr::select(CBMB, dna_reads), by=c("CBMB"='CBMB')) %>% 
#       group_by(CB, CBMB, dist_cut, dna_reads) %>% 
#       summarize(Freq = sum(ct)) %>% ungroup %>% dplyr::select(-CBMB)
# 
#     return(clusize_dist)
#     
#   }, mc.cores = 6) %>% do.call(rbind,.)
#   
#   
#   saveRDS(all_chr,
#           paste0("/data/new_9_celltype_clusize_genomic_dist_contact_freq/",
#                  cell_type, "_", sex, "_", cond, "_contact_freq_genomic_dist.rds")
#           )
#   
# }
# 
# expand.grid(names(cell_type_color),
#             c("F", "M"),
#             c("AD", "CTRL")) %>% 
#   apply(., 1, function(x){
#     print(x[[1]])
#     print(x[[2]])
#     print(x[[3]])
#     export_clusize_dist(cell_type = x[[1]], sex=x[[2]], cond=x[[3]])
#   })


import_all <- function(CTYPE, sex = "all", cond = 'all'){
  
  if(sex == 'all'){SEX = c("F", "M")} else {SEX = sex}
  if(cond == 'all'){COND = c("AD", "CTRL")} else {COND = cond}
  
  res <- lapply(SEX, function(se){
    lapply(COND, function(cd){
      
      df <- readRDS(paste0("/data/new_9_celltype_clusize_genomic_dist_contact_freq/",
                           CTYPE,"_",se,"_",cd,"_contact_freq_genomic_dist.rds"))
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>% mutate(CTYPE = CTYPE) %>% 
    mutate(size = cut(dna_reads, c(0, 2, 5, 10, 20, 50, 100))) %>% 
    mutate(size = factor(size, levels = c("(0,2]", "(2,5]","(5,10]","(10,20]","(20,50]","(50,100]"))) 
  
  return(res)
  
}



combined_curve_df <- pbmcapply::pbmclapply(names(cell_type_color), function(ctype){
  
  print(ctype)
  import_all(ctype)
  
}, mc.cores = 10) %>% do.call(rbind,.) 


plot_curve <- function(dist_df, group_by_vars = c("size"), color_by_vars, 
                       line_legend_title, palette, TITLE = "DNA-DNA interactions"){
  
  # group_by_vars should be a list of variable to group_by 
  # color_by_vars should be a list of variable to alpha by 
  # main_color should be one color to vary based upon
  
  g <- dist_df %>% 
    mutate(dist_end = gsub("^.*,", "", dist_cut)) %>%
    mutate(dist_end = gsub("]","", dist_end)) %>%
    mutate(dist_end = as.numeric(dist_end)) %>%
    group_by_(.dots = c(group_by_vars, "dist_end")) %>% 
    summarize(Freq = sum(Freq)) %>% 
    mutate(freq_ratio = Freq/sum(Freq)) %>% 
    mutate(CTYPE=factor(CTYPE, levels=c("Ex","In","Opc","Ast","Mic","Oli","Vas"))) %>% 
    ggplot() +
    #annotate("rect", xmin=10^(6), xmax = 2*10^7, ymin=10^(-4), ymax=0.1, alpha=0.1)+
    geom_line(aes_string(x="dist_end", y="freq_ratio",
                         color = color_by_vars),
              size=0.3, alpha = 0.7) +
    scale_x_log10(labels=scales::trans_format('log10', scales::math_format(10^.x)),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  limits = c(4.5e4, 1.5e8)) + # 4.5e4, 2*1e7
    scale_y_log10(labels=scales::trans_format('log10', scales::math_format(10^.x)),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  limits = c(10^(-5.5), 0.05),
                  expand = c(0, 0)) +
    annotation_logticks(size=0.2,
                        short = unit(0.01, "in"),
                        mid = unit(0.02, "in"),
                        long = unit(0.03, "in")) +
    # scale_color_manual(name = "DNA cluster size", values = RColorBrewer::brewer.pal(7,"Blues")[2:7])+
    scale_color_manual(name = line_legend_title,
                       values = palette)+ # RColorBrewer::brewer.pal(7,"Blues")[2:7]
    labs(x = expression(paste("Genomic distance, ",italic("s")," (bp)")),
         y = expression(paste("Contact probability, ",italic("Pc"))))+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill="white"),
      axis.text = element_text(size=6),
      axis.title = element_text(size=8),
      legend.position = c(0.25,0.4),
      plot.title = element_blank(),
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
  
  
  return(g)
}
pdf("/Figures/Figure-3b-GDCF-ctype.pdf",
    width = 1.8, height = 1.8)

combined_curve_df %>% 
  plot_curve(., group_by_vars = c("CTYPE"), color_by_vars = c("CTYPE"),
             line_legend_title = "Cell type", palette = unlist(cell_type_color))

dev.off()


## ---- Figure 3c Ex scGDCF heatmap barplot ----------

# combine the heatmap with barplot showing the proportion of each category
# use the proportion to indicating the contact frequency
# proportion also needs to be log scaled


average_top10_avg_index <-  
  
  apply(sc_dist_freq_mtx, 1, function(x){
    
    top_indices <- order(x, decreasing = TRUE)[1:10]
    return(mean(top_indices))
    
  })


result <- as.integer(sub("^\\((\\d+\\.\\d+e\\+\\d+).*", "\\1", colnames(sc_dist_freq_mtx)[average_top10_avg_index] ))

# merged two bins into one from dist_bins
dist_bins <- lseq(5000, 1.5 * 1e8, length.out = 150)
dist_bins_odd <- dist_bins[seq(1,length(dist_bins), 5)]
length(dist_bins_odd)

g_bar <- data.frame(CB = names(average_top10_avg_index),
                    start = result) %>% 
  mutate(bin_cut = cut(start, c(0,1e5, 1e6, 1e7, 1.5e8))) %>% 
  left_join(cell_anno %>% dplyr::select(CB, cell_type) %>% 
              mutate(cell_type=gsub("Ex","ExN",cell_type)) %>% 
              mutate(cell_type=gsub("In","InN",cell_type))) %>%
  dplyr::filter(cell_type == "ExN") %>% 
  group_by(cell_type, bin_cut) %>% summarize(ct = dplyr::n()) %>% 
  ungroup %>% group_by(cell_type) %>% mutate(ctype_tol = sum(ct)) %>% 
  mutate(prop = ct / ctype_tol) %>% mutate(prop_int = prop*1e5) %>% 
  complete(bin_cut) %>% mutate(prop_int = ifelse(is.na(prop_int), 1, prop_int)) %>% 
  ggbarplot(., x='cell_type', y='prop_int', fill="bin_cut",
            position = position_dodge(width=0.7), ylab="Proportion",size=0)+
  scale_fill_manual(name="Genomic Distance",
                    labels=c("0-100kb","100kb-1Mb",'1Mb-10Mb','10Mb-150Mb'),
                    values=c("grey60","grey70","grey80","grey90"))+
  scale_y_log10(expand=c(0,0),
                breaks = c(1e1,1e2,1e3,1e4),
                labels = c('1e-4', "1e-3","1e-2", '1e-1'))+
  scale_x_discrete(expand = c(0,0))+
  theme(legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust=0.5, size=8),
        plot.margin = unit(c(0,0,0,0), "in"),
        axis.ticks = element_line(size=0.2),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.x = element_blank())


# heatmap 
Ex_cells <- intersect(intersect_cells, cell_anno2 %>% filter(cell_type=="ExN") %>% .[["CB"]])

horizontal_mtx <- sc_dist_freq_mtx[Ex_cells,2:ncol(sc_dist_freq_mtx)]


# add row distance annotation

idx2 <- sapply(distance_near, function(x) {
  which.min(abs(dist_bins - x))
})
labels2 <- distance_near %>% sapply(., format_bp)

new_col_anno = HeatmapAnnotation(foo = anno_mark(at = idx2,
                                                 side = "bottom",
                                                 labels = labels,
                                                 labels_gp = gpar(fontsize=6),
                                                 link_gp = gpar(size=0.2),
                                                 link_height = unit(0.1,'in')))


grob_hm = grid.grabExpr(
  draw(
    Heatmap(
      horizontal_mtx,
      name = "Frequency",
      cluster_rows = F,
      cluster_columns = F,
      col = circlize::colorRamp2(c(0, 0.01, 0.05), c("#e0f0f3", "#81acce", "#343d96")),
      row_title_gp = gpar(fontsize=8),
      show_row_names = F,
      show_column_names = F,
      bottom_annotation = new_col_anno,
      row_title = "Individual cells\nExN",
      show_heatmap_legend = F,
      use_raster = T
    ),
    padding = unit(c(0, 0, 0, 0), "in")
  )
)

pdf("/Figures/Figure-3c-ExN-sc-bar-hm.pdf",
    width = 1.5, height = 2.2)

(g_bar/grob_hm)+plot_layout(heights = c(0.4,1))

dev.off()


## ---- Figure 3d scGDCF heatmap ----------

CLUSIZE <- 100

# lapply(names(cell_type_color), function(ctype) {
#   
#   opc_dna <- import_grs(CTYPE = ctype,
#                         mole = "DNA",
#                         import_clu = F)
#   
#   opc_clusize <- import_clusize(CTYPE = ctype, mole = "DNA")
#   
#   clu_cutoff <- opc_clusize %>% filter(dna_reads>1, dna_reads < CLUSIZE)
#   
#   all_chr <-
#     split(opc_dna, seqnames(opc_dna)) %>% pbmcapply::pbmclapply(., function(x) {
#       print(head(x))
#       onechr <-
#         x %>% as.data.frame() %>% mutate(ID = seq(1, nrow(.))) %>%
#         filter(CBMB %in% clu_cutoff$CBMB)
#       
#       clusize_dist <- onechr %>%
#         left_join(
#           onechr %>% dplyr::select(start, CBMB, ID),
#           by = c("CBMB" = "CBMB"),
#           multiple = "all"
#         ) %>%
#         filter(ID.x < ID.y) %>%
#         mutate(dist = abs(start.x - start.y)) %>%
#         dplyr::filter(dist > 10) %>%
#         mutate(dist_cut = cut(dist, dist_bins)) %>%
#         dplyr::group_by(CB, dist_cut) %>%
#         dplyr::summarise(ct = dplyr::n(), .groups = "drop")
#       
#       return(clusize_dist)
#       
#     }, mc.cores = 12) %>% do.call(rbind, .) %>% 
#     group_by(CB, dist_cut) %>% dplyr::summarize(CT = sum(ct), .groups = "drop") %>% # for all chromosomes
#     group_by(CB) %>% mutate(sum_ct = sum(CT)) %>% 
#     mutate(freq = CT / sum_ct) %>%
#     dplyr::select(CB, dist_cut, freq)
#   
#   
#   saveRDS(
#     all_chr,
#     paste0(
#       "/data/new_9_1_SC_DD_GDCF_heatmap/",
#       ctype,
#       "_SC_contact_freq_genomic_dist.rds"
#     )
#   )
#   
# })


## ---- 2. export single cell contact frequency cluster size weighted version ----------

# normalize the weight of contacts by cluster size (1/dna_reads)

# CLUSIZE <- 100
# 
# lapply(names(cell_type_color), function(ctype) {
#   
#   opc_dna <- import_grs(CTYPE = ctype,
#                         mole = "DNA",
#                         import_clu = F)
#   
#   opc_clusize <- import_clusize(CTYPE = ctype, mole = "DNA")
#   
#   clu_cutoff <- opc_clusize %>% filter(dna_reads < CLUSIZE, dna_reads>1)
#   
#   all_chr <-
#     split(opc_dna, seqnames(opc_dna)) %>% pbmcapply::pbmclapply(., function(x) {
#       
#       onechr <-
#         x %>% as.data.frame() %>% mutate(ID = seq(1, nrow(.))) %>%
#         filter(CBMB %in% clu_cutoff$CBMB)
#       
#       clusize_dist <- onechr %>%
#         left_join(
#           onechr %>% dplyr::select(start, CBMB, ID),
#           by = c("CBMB" = "CBMB"),
#           multiple = "all"
#         ) %>%
#         filter(ID.x < ID.y) %>%
#         mutate(dist = abs(start.x - start.y)) %>%
#         dplyr::filter(dist > 10) %>%
#         mutate(dist_cut = cut(dist, dist_bins)) %>%
#         left_join(opc_clusize %>% mutate(weight = 1/dna_reads) %>% select(CBMB, weight), by=c("CBMB"="CBMB")) %>% 
#         dplyr::group_by(CB, dist_cut) %>%
#         dplyr::summarise(ct = sum(weight), .groups = "drop")
#       
#       return(clusize_dist)
#       
#     }, mc.cores = 12) %>% do.call(rbind, .) %>% 
#     group_by(CB, dist_cut) %>% dplyr::summarize(CT = sum(ct), .groups = "drop") %>% # for all chromosomes
#     group_by(CB) %>% mutate(sum_ct = sum(CT)) %>% 
#     mutate(freq = CT / sum_ct) %>%
#     dplyr::select(CB, dist_cut, freq)
#   
#   
#   saveRDS(
#     all_chr,
#     paste0(
#       "/data/new_9_1_SC_DD_GDCF_heatmap/",
#       ctype,
#       "_SC_contact_freq_genomic_dist_clusize_weighted.rds"
#     )
#   )
#   
# })




# ---- 3. Comupte sc contact frequency mtx ----------

all_cell_all_lib_dist_freq <-
  
  pbmcapply::pbmclapply(names(cell_type_color), function(ctype) {
    dna_clu_size <- readRDS(
      paste0(
        "/data/new_9_1_SC_DD_GDCF_heatmap/",
        ctype,
        "_SC_contact_freq_genomic_dist.rds"
      )
    )
    
  }) %>% do.call(rbind, .) %>% 
  filter(!is.na(dist_cut))


all_cell_all_lib_dist_freq <-
  all_cell_all_lib_dist_freq %>%
  filter(dist_cut != "(1e+03,5e+03]") %>% # newly added
  mutate(dist_cut_char = as.character(dist_cut)) %>%
  mutate(cut_start = gsub("\\((.*),(.*)\\]$", "\\1", dist_cut_char) %>% as.numeric()) %>%
  mutate(cut_end = gsub("\\((.*),(.*)\\]$", "\\2", dist_cut_char) %>% as.numeric()) %>%
  mutate(cut_interval = (cut_end - cut_start) / 1e6) %>%
  mutate(weighted_freq = freq / cut_interval) %>%  # weighted sum normalized by interval total length
  group_by(CB) %>% dplyr::mutate(sum_ = sum(weighted_freq)) %>% 
  mutate(freq = weighted_freq / sum_)


sc_dist_freq_mtx <- all_cell_all_lib_dist_freq %>% filter(!is.na(dist_cut)) %>% ungroup %>%
  dplyr::select(CB, dist_cut, freq) %>%
  pivot_wider(., names_from = "dist_cut", values_from = "freq", values_fill = 0) %>%
  column_to_rownames('CB')

sc_dist_freq_mtx_rot <- sc_dist_freq_mtx[,rev(levels(all_cell_all_lib_dist_freq$dist_cut)[2:150])] %>% t()


# ---- 4. compute most frequent genomic bin add mean_gd ----------

# for each cell add a feature indicating the average distance (bin) for the top5 most frequent bin

calculate_average_from_range <- function(range_string) {
  # Step 1: Extract the numeric values
  range_values <- unlist(strsplit(gsub("[\\[\\]()]", "", range_string), ","))
  
  range_values <- range_values %>% sapply(., function(x){gsub("\\(|\\]", "", x) %>% as.numeric}) %>% unname()
  
  # Step 3: Calculate the average
  average_value <- mean(range_values)
  
  return(average_value)
}

mean_frequent_bin <- 
  apply(sc_dist_freq_mtx_rot, 2, function(x){
    avg_range <- 
      rownames(sc_dist_freq_mtx_rot)[order(unlist(x), decreasing = TRUE)[1:10] %>% mean %>% round] %>% 
      calculate_average_from_range
    return(avg_range)
  })
most_freq_bin_df <- mean_frequent_bin %>% as.data.frame() %>% rownames_to_column() %>% set_colnames(c("CB","mean_gd"))

if("mean_dg" %ni% colnames(cell_anno)){
  cell_anno %<>% left_join(most_freq_bin_df, by=c('CB'="CB"))
} else{print("mean_dg already exists")}


# long short threshold
THRE = 3e4
cell_anno %<>% mutate(DD_type = ifelse(mean_gd>THRE, 'long', 'short'))



## ---- 5. compute cell aging score ----------

DefaultAssay(merged_brain) <- "RNA"
merged_brain %<>% Seurat::NormalizeData() %>% Seurat::ScaleData()

expr <- merged_brain@assays$RNA@data %>% t()
sum(rownames(expr) == cell_anno$CB)


human_aging_marker <- fread("/database/aging/SCALE/Download_human_aging.csv")


aging_marker_expr <- merged_brain@assays$RNA@data[intersect(human_aging_marker$Symbol, rownames(merged_brain)),] # 714 marker genes

cell_anno$age =  as.numeric(as.character(ifelse(cell_anno$expired_age=="≥90", 90, cell_anno$expired_age)))

# correlation sign
signs <- apply(aging_marker_expr, 1, function(x){
  return(sign(cor(x, cell_anno$age))) # 
})


# calculate SCALE scores, according to https://github.com/ChengLiLab/SCALE/blob/main/demo/Calculate_SCALE_score.R
aging_gene <- rownames(aging_marker_expr)

weights = signs*
  Matrix::rowSums(as.matrix(merged_brain@assays$RNA@counts[aging_gene,])>0)/dim(merged_brain@assays$RNA@data)[2]

z_matrix = t(
  scale(
    t(as.matrix(merged_brain@assays$RNA@data[aging_gene,])),
    center = T, scale = T))

SCALE_score = t(z_matrix[aging_gene,]) %*% weights
cell_anno$aging_scores = SCALE_score


cell_anno$age =  as.numeric(as.character(ifelse(cell_anno$expired_age=="≥90", 90, cell_anno$expired_age)))

# order columns (cells) by the maxium value row index
intersect_cells <- colnames(sc_dist_freq_mtx_rot)[
  apply(sc_dist_freq_mtx_rot, 2, function(x){
    
    top_indices <- order(x, decreasing = TRUE)[1:10]
    return(mean(top_indices))
    
  }) %>% order
]



# add row distance annotation
distance_near <- c(5e3, 1e4, 1e5, 1e6, 5e6, 1e7, 5e7)

# actual cut range that closes to near
dist_bins_rev <- dist_bins %>% rev()

idx <- sapply(distance_near, function(x) {
  which.min(abs(dist_bins_rev - x))
})
# labels <- dist_bins_rev[idx] %>% sapply(., format_bp)
labels <- distance_near %>% sapply(., format_bp)



row_anno = rowAnnotation(foo = anno_mark(at = idx,
                                         side = "left",
                                         labels = labels))

col_anno = HeatmapAnnotation(

  `Pc(s) group` = setNames(cell_anno$rename_DD_type, cell_anno$CB)[intersect_cells],
  `Chronological age` = setNames(cell_anno$age, cell_anno$CB)[intersect_cells],
  `Transcriptomic age` = setNames(cell_anno$aging_scores, cell_anno$CB)[intersect_cells],
  `AD pathology` = setNames(cell_anno$type, cell_anno$CB)[intersect_cells],
  Sex = setNames(cell_anno$sex, cell_anno$CB)[intersect_cells],
  
  # spliceosome 
  RNU4ATAC = setNames(log2(expr_df$RNU4ATAC+1), cell_anno$CB)[intersect_cells],
  # `RNU2-2P` = setNames(log2(expr_df$`RNU2-2P`+1), cell_anno$CB)[intersect_cells],
  RNU12 = setNames(log2(expr_df$RNU12+1), cell_anno$CB)[intersect_cells],
  `RNU5A-1` = setNames(log2(expr_df$`RNU5A-1`+1), cell_anno$CB)[intersect_cells],
  `RNU4-2` = setNames(log2(expr_df$`RNU4-2`+1), cell_anno$CB)[intersect_cells],
  
  col = list(

    `AD pathology` = setNames(c("darkred", "grey"), nm = c("AD", "CTRL")),
    Sex = setNames(c("#e79099", "#3f6a80"), nm = c("F", "M")),
    `Pc(s) group` = setNames(c("orange", "grey80"), nm = c("LCS-erosion", "LCS-perserved")),
    `Chronological age` = circlize::colorRamp2(c(50, 70, 90), c("#f4faf9", "#beddd8", "#1b9587")),
    `Transcriptomic age` = circlize::colorRamp2(c(-15, 0, 10), c("blue", "white", "red")),
    RNU4ATAC = circlize::colorRamp2(c(0, 3), c("#f7f4f9", "#17107b")),
    RNU12 = circlize::colorRamp2(c(0, 8), c("#f7f4f9", "#17107b")),
    `RNU5A-1` = circlize::colorRamp2(c(0, 3), c("#f7f4f9", "#17107b")),
    `RNU4-2` = circlize::colorRamp2(c(0, 3), c("#f7f4f9", "#17107b")) 
    
  ),
  annotation_name_side = "left"
)


cell_anno2 <- cell_anno %>% mutate(cell_type = gsub("Ex","ExN",cell_type)) %>% mutate(cell_type = gsub("In","InN",cell_type))
pdf("/Figures/3_Dist_Freq_heatmap/DD_dist_freq_weight_norm_binsize_norm_celltype_biomarker2.pdf",
    width = 15, height = 6)

set.seed(1117)
draw(
  Heatmap(
    sc_dist_freq_mtx_rot[,intersect_cells],
    name = "Freq",
    cluster_rows = F,
    cluster_columns = F,
    col = circlize::colorRamp2(c(0, 0.01, 0.05), c("#e0f0f3", "#81acce", "#343d96")),
    # col = circlize::colorRamp2(c(0, 0.005, 0.008, 0.01, 0.03, 0.05), c("#e0f0f3", "#deedf2","#b6cde0","#7ba4ce", "#343d96", "#252981")),
    column_split = factor(setNames(cell_anno2$cell_type, cell_anno2$CB),
                          levels=c("ExN","InN", "Opc", "Ast", "Mic", "Oli", "Vas"))[intersect_cells],
    # column_split = setNames(paste0(cell_anno2$cell_type, cell_anno2$type), cell_anno2$CB)[intersect_cells],
    show_row_names = F,
    show_column_names = F,
    bottom_annotation = col_anno,
    left_annotation = row_anno,
    use_raster = T
  ),
  padding = unit(c(2, 10, 10, 10), "mm")
)
dev.off()

## ---- Figure 3e ----------
############# add sample size n number to the boxplot
DAT = cell_anno %>% mutate(cell_type=factor(cell_type, 
                                            levels=rev(c("Ex", "In", "Ast", "Opc", "Vas", "Mic", "Oli")))) %>% 
  mutate(DD_type = case_when(DD_type == "long" ~ "LCS-eroded",
                             DD_type == "short" ~ "LCS-Perserved"))

compare_means(
  aging_scores ~ rename_DD_type, group.by = "cell_type", data = DAT,
  method = "wilcox"
)

compare_means(
  age ~ rename_DD_type, group.by = "cell_type", data = DAT,
  method = "wilcox"
)

g_SCALE <- ggboxplot(DAT,
                     x='cell_type', y='aging_scores', fill = 'DD_type',
                     palette = c("#f1aa34", "grey60"),
                     size=0.2, color='grey70', outlier.size = 0.3,
                     ylab = "Transcriptomic aging score") + 

  stat_compare_means(aes(group=`DD_type`, label = sprintf("p = %2.2e", as.numeric(..p.format..))), method='wilcox', hide.ns=T, size=2,
                     method.args = list(alternative = "less")
  )+
  geom_text(data = DAT %>% group_by(cell_type, DD_type) %>% summarize(len=dplyr::n()),
            aes(x = factor(cell_type),
                y = min(DAT$aging_scores) - 1,
                label = paste(len,"\n"),
                group = DD_type),
            position = position_dodge(.8), size=2, angle=0)+
  scale_y_continuous(position = "right", limits = c(-22, 11)) +
  coord_flip() +  
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size=6),
        axis.title = element_text(size = 7),
        legend.title = element_blank(),
        # legend.position = c(0.8,0.1),
        legend.text=element_text(size=5),
        legend.key.size = unit(0.1, "in"),
        legend.margin=margin(t = 0, unit='cm')
        )


pdf("/data/Erosion/Figure-4-e-LCS-erosion-tx-age.pdf",
    width=1.7, height = 2.1)
g_SCALE

dev.off()


## ---- Figure 3f eQTL----------

OUT_DIR = "/MUSIC_eQTL_contacts/eQTL_contacts_V4/"
snp_pos = data.frame(fread("/annotation/snp_pos.txt"))

##### eQTL data
eQTL_all = data.frame(fread("/database/brain_celltype_eQTL/pval005_filtered_eQTL_brain_celltype_combined.txt")) %>% 
  filter(Nominal_pval < 10^-4 & ctype %ni% c("Pericytes","Endothelial")) %>%
  mutate(ctype = recode(ctype, 'Excitatory' = 'Ex', "Oligodendrocytes" = "Oli", "Microglia" = "Mic", "Inhibitory" = "In",
                        'Astrocytes' = 'Ast', 'OPCs' = 'Opc'))
eQTL_all$gene = sapply(eQTL_all$Gene_id, function(x){strsplit(x, "_")[[1]][1]})
eQTL_all = eQTL_all %>% mutate(eQTL = paste0(SNP_id, "_", gene),
                               eQTL_ctype = paste0(SNP_id, "_", gene, "_", ctype))

eQTL_all = merge(eQTL_all, snp_pos[,c(1,2)], by.x = "SNP_id", by.y = "SNP")
eQTL_all$chr = sapply(eQTL_all$SNP_id_hg38, function(x){strsplit(x, "\\:")[[1]][1]})
eQTL_all$pos = sapply(eQTL_all$SNP_id_hg38, function(x){strsplit(x, "\\:")[[1]][2]})

combined <- purrr::reduce(lapply(sort(unique(eQTL_all$ctype)),
                                 function(x){
                                   return(eQTL_all %>% filter(eQTL_all$ctype == x) %>% select(eQTL) %>% mutate(V2 = 1) %>% rename(!!x := "V2"))
                                 }), full_join)
combined[is.na(combined)] <- 0
colnames(combined)[3:4] = c("ExN","InN")

cell_type_color <- list(
  Ex = "#7EB5C6",
  In = "#BEBADA",
  Oli = "#F4B461",
  Ast = "#EE7E72",
  Opc = "#F7CCE5",
  Mic = "#448369",
  Vas = "#BCDE91"
)

pdf("/MUSIC_eQTL_contacts/eQTL_contacts_V4/eQTL_upset.pdf", width = 6.5, height = 5)
UpSetR::upset(combined, 
              sets = c("Mic","Opc","Oli","Ast","InN","ExN"),
              sets.x.label = "Number of DD contacts",
              order.by = "freq",
              nintersects = 18,
              sets.bar.color = cell_type_color[1:6],
              keep.order = T,
              number.angles = 45,
              text.scale = c(1.8, 1.8, 1.2, 1.2, 1.6, 0),
              mainbar.y.label = "",
              point.size = 2.5,
              line.size = 0.3)
dev.off()


eQTL_all_ct = eQTL_all %>%
  group_by(eQTL) %>%
  summarize(ct = length(unique(ctype)))

eQTL_unique = eQTL_all_ct %>% 
  filter(ct==1) %>%
  select(eQTL) %>% as.data.frame()

eQTL_all_unique = eQTL_all %>% 
  filter(eQTL %in% eQTL_unique$eQTL) %>%
  mutate(eQTL_ctype = paste0(eQTL,"_",ctype)) # eQTL-gene pairs that appear only in one cell type


### Total eQTL-gene pairs and MUSIC DD
eQTL_gene_pairs_MUSIC_ovlp_num_total <-
  mclapply(c( "Ast", "Ex", "In", "Mic", "Oli", "Opc"), function(cell_type) # cell_type refers to MUSIC
  {
    RD_pairs <- readRDS(paste0("/project/14____scMARGI/eQTL/eQTL_rds/",cell_type,"_eQTL_ovlp.rds"))
    
    ##### Remove eQTL-gene pairs that are shared across cell types
    RD_pairs <- RD_pairs %>%
      mutate(eQTL = paste0(RD_pairs$SNP, "_", RD_pairs$gene))
    RD_pairs = RD_pairs %>%
      filter(eQTL %in% eQTL_all_unique$eQTL)
    
    res <- RD_pairs %>% group_by(SNP, gene, ctype, mole.x, mole.y) %>%
      summarize(clu_num = uniqueN(CBMB.x), .groups = 'drop') %>%
      group_by(mole.x, mole.y, ctype) %>%
      summarize(ct = dplyr::n(), .groups = 'drop') %>%
      filter(ctype %ni% c("Pericytes", "Endothelial")) %>%
      mutate(ctype = recode(ctype, 'Excitatory' = 'Ex', "Oligodendrocytes" = "Oli", "Microglia" = "Mic", "Inhibitory" = "In",
                            'Astrocytes' = 'Ast', 'OPCs' = 'Opc')) %>%
      mutate(ctype = factor(ctype, levels=c("Ast", "Ex", "In", "Mic", "Oli", "Opc"))) %>%
      pivot_wider(., names_from = 'ctype', values_from = 'ct', values_fill = 0) %>%
      mutate(music_ctype = cell_type) %>%
      filter(mole.x == "DNA", mole.y == "DNA")
    return(res)
  }, mc.cores = 6) %>% do.call(rbind,.)


# Chi-squared 6x6 table
temp = data.frame(eQTL_gene_pairs_MUSIC_ovlp_num_total %>% select(c(3:9))) %>% `rownames<-`(.[,"music_ctype"]) %>% select(-7)
temp = temp[c("Ex","In","Ast","Oli","Opc","Mic"),c("Ex","In","Ast","Oli","Opc","Mic")]
colnames(temp)[1:2] = c("ExN","InN")
rownames(temp)[1:2] = c("ExN","InN")
print(chisq.test(as.table(as.matrix(temp))))

# Chi-squared 2x2 tables
df <- lapply(c("Ex","In","Ast","Oli","Opc","Mic"), 
             function(cell_type) {
               
               temp = reshape2::melt(eQTL_gene_pairs_MUSIC_ovlp_num_total[,3:9], id = "music_ctype")
               tbl = as.table(matrix(c(temp %>% filter(music_ctype == cell_type & variable == cell_type) %>% select(value) %>% as.numeric,
                                       temp %>% filter(music_ctype == cell_type & variable != cell_type) %>% select(value) %>% sum(),
                                       temp %>% filter(music_ctype != cell_type & variable == cell_type) %>% select(value) %>% sum(),
                                       temp %>% filter(music_ctype != cell_type & variable != cell_type) %>% select(value) %>% sum()), nrow = 2, ncol = 2))
               
               p_val = chisq.test(tbl)$p.value
               OR = data.frame(oddsratio.wald(tbl)$measure)$estimate[2]
               CI_lower = data.frame(oddsratio.wald(tbl)$measure)$lower[2]
               CI_upper = data.frame(oddsratio.wald(tbl)$measure)$upper[2]
               n = sum(tbl)
               
               temp = tbl
               colnames(temp) = c("in_cell_type_MUSIC", "not_in_cell_type_MUSIC")
               rownames(temp) = c("in_cell_type_eQTL", "not_in_cell_type_eQTL")
               sn = ifelse(cell_type %in% c("Ex","In"), paste0(cell_type,"N"), cell_type)
               return(data.frame(cell_type = ifelse(cell_type %in% c("Ex","In"), paste0(cell_type, "N"), cell_type),  p_val, OR, CI_lower, CI_upper, n))
               
             }) %>% do.call(rbind,.)


df$cell_type = factor(df$cell_type, levels = c("ExN","InN","Ast","Oli","Opc","Mic"))
df$y = c(6:1)

pdf("/MUSIC_eQTL_contacts/eQTL_contacts_V4/chisq_result_CTS_eQTL_vs_all.pdf", width = 2.5, height = 1.2)
ggplot(df, aes(x=OR, y=y)) +
  geom_errorbar(aes(xmin=CI_lower, xmax=CI_upper), width=.3, position=position_dodge(.9), linewidth = 0.3) +
  geom_point(aes(x=OR, y=y), color = "blue", size = 1) +
  scale_x_continuous(breaks = c(0.9,1,1.7), labels = as.character(c(0.9,1,1.7)), limits = c(0.9,1.7)) +
  scale_y_continuous(breaks = c(1:6), label = as.character(sort(df$cell_type, decreasing = T))) +
  ylab("") +
  labs(x = "Odds ratio", y = NULL) +
  theme_bw() +
  theme(text = element_text(size=6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.3)
dev.off()