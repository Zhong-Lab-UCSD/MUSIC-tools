# Figure1.R
require(tidyverse)
require(data.table)
require(magrittr)
require(ggpubr)

options(stringsAsFactors=FALSE)
options(scipen=3)

MASTER_OUTDIR = "/outputs/"

DD_col = "#5cc4ef"
RR_col = "#f3af42"
RD_col = "#bebdbe"
      
gtfgr36 <- readRDS("/mnt/extraids/SDSC_NFS/wenxingzhao/database/Human/genome/Robj/gtf36_gene_granges.rds")


### Figure 1a-c ----------------
cell_dna_human_reads <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_DNA_human.sort.cell_reads.csv")) %>% set_names(c("CB","human_dna_reads"))
cell_rna_human_reads <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_RNA_human.sort.cell_reads.csv")) %>% set_names(c("CB","human_rna_reads"))

human_cell_reads <- cell_dna_human_reads %>% left_join(cell_rna_human_reads, by=c("CB"="CB")) %>% filter(CB %in% human_cells) %>% mutate(human_rna_reads = replace_na(human_rna_reads, 0)) %>% mutate(total_reads = human_dna_reads + human_rna_reads)

human_tol_contacts <- human_clusters %>%
  mutate(DD_con = dna_reads*(dna_reads-1)/2, 
         RR_con = rna_reads*(rna_reads-1)/2, 
         RD_con = as.numeric(dna_reads)*as.numeric(rna_reads)) %>% 
  group_by(CB) %>% summarize(DD_contacts = sum(DD_con),
                             RR_contacts = sum(RR_con),
                             RD_contacts = sum(RD_con))


human_clusters_summary <- human_clusters %>% 
  group_by(CB) %>% summarize(DD_num = sum(type=="DD"),
                             RD_num = sum(type=="RD"),
                             RR_num = sum(type=="RR")) %>% 
  mutate(total_clusters = DD_num+RD_num+RR_num) %>% 
  mutate(DD_prop = DD_num/total_clusters,
         RD_prop = RD_num/total_clusters,
         RR_prop = RR_num/total_clusters)



h_cell_order <- human_cell_reads %>% arrange(desc(human_dna_reads)) %>% .[['CB']]

human_comb_df <- human_cell_reads %>% left_join(human_clusters_summary, by=("CB"="CB")) %>% 
  left_join(human_tol_contacts, by=c("CB"="CB")) %>% 
  dplyr::select(CB, human_dna_reads, human_rna_reads, DD_num, RD_num, RR_num,
                DD_contacts, RR_contacts, RD_contacts) 

total_h_cell <- h_cell_order %>% length

h_x_labels_at <- quantile(1:total_h_cell, c(0, 0.33, 0.66, 1)) %>% as.integer()
h_x_labs = rep("", total_h_cell)
for(i in h_x_labels_at){
  h_x_labs[i] <- paste0("cell #",i)
}


g_all_h <- human_comb_df %>% 
  melt() %>% mutate(group = ifelse(grepl("reads",variable), "reads", 
                                   ifelse(grepl("contacts", variable), "contact", "clusters"))) %>% 
  mutate(group=factor(group, levels=c("reads","clusters","contact"))) %>% 
  ggplot()+geom_point(aes(x=factor(CB, levels = h_cell_order),
                          y=value, color=variable), alpha=0.5, size=0.1) + 
  facet_wrap(~group, ncol = 1, strip.position = "right", scales = "free_y")+
  scale_y_log10(labels = scales::scientific) + 
  scale_color_manual(labels = c("DNA reads", "RNA reads",
                                "DD clusters", "RD clusters", "RR clusters",
                                "DD contacts", "RD contacts", "RR contacts"),
                     values=c("#354d90", "#e8552b",DD_col, RD_col, RR_col, DD_col, RD_col, RR_col),
                     breaks = c("human_dna_reads", "human_rna_reads",
                                "DD_num", "RD_num", "RR_num",
                                "DD_contacts", "RD_contacts", "RR_contacts"))+
  scale_x_discrete(labels=h_x_labs)+
  xlab("Single cells") + ylab("Number") + ggtitle("Human H1 cells") +
  guides(colour = guide_legend(override.aes = list(size=1)))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"),,
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size=0.2),
        strip.background = element_rect(fill=NA),
        axis.title.x = element_text(size=8),
        axis.title.y = element_blank(),
        axis.text = element_text(size=6),
        legend.text = element_text(size=6),
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust=0.5, size=8),
        plot.margin = unit(c(0,0,0,0), "in"),
        axis.ticks = element_line(size=0.2),
        axis.line= element_line(size = 0.2, linetype = "solid", colour = "black")) 


pdf("/Figures/Figure-1a-c.pdf",
    width=4, height = 4)
g_all_h
dev.off()