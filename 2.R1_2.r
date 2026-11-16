
########P value region select

mark_top_pvalue_regions <- function(df, psotion="Chromosome", pvalue="log10FDR",FDR=40, n_top = 100, window = 5000) {

names(df)[names(df) == psotion] <- "position2"
names(df)[names(df) == pvalue] <- "pvalue"
  # 创建新列并初始化为"ns"
  df$top <- "ns"
  # 创建未标记行的索引
  unmarked_idx <- 1:nrow(df)

  P_FDR <- which(df$pvalue >= FDR)

  for (i in 1:n_top) {
    # 检查是否还有未标记的行
    if (length(unmarked_idx) == 0 ||  length( P_FDR) == 0) {
      warning(paste("共找到TOP", i, "的-log10FDR>=",FDR, sep = ""))
      break }
    # 在未标记行中找到pvalue最小的行
    min_pval_idx <- unmarked_idx[which.max(df$pvalue[unmarked_idx])]
    center_pos <- df$position2[min_pval_idx]
    # 确定± 连锁范围内的所有行
    position_idx <- which(df$position2 >= (center_pos - window) & df$position2 <= (center_pos + window))
    # 只保留未标记的行
    position_idx <- position_idx[position_idx %in% unmarked_idx]
    
    # 如果范围内没有未标记行，则跳过
    if (length(position_idx) == 0) {next}
    # 标记这些行为TOPi
    df$top[position_idx] <- paste0("Locus", i)
    # 从未标记索引中移除已标记的行
    unmarked_idx <- setdiff(unmarked_idx, position_idx)
    P_FDR <- setdiff(P_FDR, position_idx )
  }  

names(df)[names(df) =="position2"] <- psotion
names(df)[names(df) =="pvalue"] <- pvalue
  return(df)
}

###########################pseudomonas&&&Enterobacterales
###########################pseudomonas&&&Enterobacterales
###########################pseudomonas&&&Enterobacterales

library(MetBrewer)
library(ggplot2)
library(ggrepel)
library(gg.gap)
library(tidyverse)

gwas <- readr::read_tsv("R1_2_result_pseudo_GWAS_less_three.xls")
names(gwas) <- c("chr", "position","qvalue")
gwas$tax <- "Pseudomonadales" 

gwas2 <- readr::read_tsv("R1_2_result_Enterobacterales_GWAS_less_three.xls")
names(gwas2) <- c("chr", "position","qvalue")
gwas2$tax <- "Enterobacterales" 
gwas2$qvalue[gwas2$qvalue == 0] <- 1

gwas <- bind_rows(gwas , gwas2)

gwag_result <- gwas %>% group_by(chr) %>% 
  summarise(chr_len=max(position)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(gwas,.,by=c("chr"="chr")) %>%
  arrange(chr,position) %>%
  mutate(Chromosome=position+tot)%>%
  mutate(SNP_id = paste0("chr.", chr,".", position))

gwag_result$log10FDR <- -log10(gwag_result$qvalue)

######################筛选SNP
significance_gwas <- mark_top_pvalue_regions (gwag_result, psotion="Chromosome", pvalue="log10FDR",FDR =51, n_top = 30, window = 100000)

######################保存数据
significance_gwas1_top_loci <- significance_gwas %>% group_by(top) %>%
summarize(end=max(position),start = min(position),chr = mean(chr)) 

write.csv(significance_gwas1_top_loci,"Fig1_significance_100K_top_loci.csv")

significance_gwas_top1 <- significance_gwas  %>% filter( log10FDR >= 20 ) %>% mutate(top = as.character(top))%>% filter(top != "ns" )
write.csv(significance_gwas_top1,"Fig1_significance_gwas_100k_top.csv")

significance_gwas_top1$chrs <- "chr"
significance_gwas_top1_bed <- significance_gwas_top1  %>% mutate(start = position -1 ) %>% 
unite("chr_ifo", chrs , chr, sep = "", remove = FALSE) %>% select(chr_ifo,start,position,log10FDR,top,SNP_id)
write.csv(significance_gwas_top1_bed,"Fig1_significance_100k_top10_bed.csv")


####################################Ext Fig2a

colorsF = c( chr1 = "#A8CFF1", chr2 = "#9EC5C3", Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7",  others ="#ECEAEA")

library(tidyverse)
library(ggplot2)

all_locus <- read.csv("R1_2_top_locus.csv",head = T)

p8 = ggplot(all_locus, aes(x = reorder(top, level), y = top_log10FDR, fill = top_hit)) + 
   geom_bar(stat = "identity", position = "stack", width = 0.8,colour ="white",linewidth=0.1) + 
    xlab("") + ylab("Top association (-log10(P value)") +
    scale_fill_manual(values = colorsF)+ coord_cartesian(ylim = c(51,54))+
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
          strip.text = element_text(size=8),text = element_text(size = 8),
          axis.line = element_line(color = "black",linewidth = 0.25),
          axis.ticks =  element_line(color = "black",linewidth = 0.25),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(color="black",size=8),
          axis.text.x = element_text(color="black",size=8,angle = 45),
          panel.spacing.x = unit(0,"cm"),
          panel.border = element_blank(),
          panel.spacing = unit(0,"lines"))

ggsave(paste0("Ext Fig2a.pdf"), p8, width=8, height=4, units="in")
ggsave(paste0("Ext Fig2a.png"), p8, width=8, height=4, units="in") 


######################Ext Fig2b

colorsF = c( chr1 = "#A8CFF1", chr2 = "#9EC5C3", Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7",  others ="#ECEAEA")

library(tidyverse)
library(ggplot2)

all_locus <- read.csv("R1_2_top_locus.csv",head = T)

all_locus_snp <-  all_locus  %>%  select(top,Pseudomonadales,Enterobacterales,level)  %>%   pivot_longer(
    cols = c(Pseudomonadales, Enterobacterales),
    # cols = -top,
    names_to = "Order",
    values_to = "SNP_num") 

write.csv(all_locus_snp,"Ext Fig2b_source.csv") 

p7 = ggplot(all_locus_snp , aes(x = reorder(top, level), y = SNP_num, fill = Order)) + 
   geom_bar(stat = "identity", position = "fill", width = 0.8,colour ="white",linewidth=0.1) + 
    xlab("") + ylab("Associated-SNP ratio (%)") +
    scale_fill_manual(values = colorsF)+
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
          strip.text = element_text(size=8),text = element_text(size = 8),
          axis.line = element_line(color = "black",linewidth = 0.25),
          axis.ticks =  element_line(color = "black",linewidth = 0.25),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(color="black",size=8),
          axis.text.x = element_text(color="black",size=8,angle = 45),
          panel.spacing.x = unit(0,"cm"),
          panel.border = element_blank(),
          panel.spacing = unit(0,"lines"))

ggsave(paste0("Ext Fig2b.pdf"), p7, width=8, height=4, units="in")
ggsave(paste0("Ext Fig2b.jpg"), p7, width=8, height=4, units="in") 

######################作图数据
# significance_gwas1 <- significance_gwas %>% distinct()  %>% filter( log10FDR >= 20 )

significance_gwas1 <- significance_gwas %>% distinct()  %>% filter( log10FDR >= 10 )

significance_gwas1$mark <- "others" 
significance_gwas2 <- significance_gwas1 %>%
  mutate(mark = as.character(mark),
    top = as.character(top),
    chr_num = as.numeric(str_extract(as.character(chr), "\\d+")),
    mark = ifelse(!is.na(chr_num) & (chr_num %% 2 == 0) & (chr_num > 0) & (log10FDR <= 20),"chr1", mark),
    mark = ifelse(!is.na(chr_num) & (chr_num %% 2 != 0) & (chr_num > 0) & (log10FDR <= 20),"chr2", mark),
    chr_num = NULL) %>%  mutate(top= case_when(top != "ns"~ "Locus_n",TRUE ~ top)) %>% 
    mutate( mark = case_when(top == "Locus_n" ~ "Locustop", TRUE ~ mark))

significance_gwas3  <- significance_gwas2 %>%  filter( mark == "Locustop"  ) %>%  filter( log10FDR >= 20 )
significance_gwas_chr  <- significance_gwas2 %>%  filter( mark != "Locustop" & mark != "others"  )
significance_gwas4  <- significance_gwas2 %>%  filter( mark == "others"  )

write.csv(significance_gwas3,"Fig1f plot_Enterbacter_Pseudomonadales_plot_Locustop.csv")

axis <- gwag_result %>% group_by(chr) %>%
summarize(center=(max(Chromosome) + min(Chromosome))/2)

#############################注释
loci_ann <- significance_gwas1 %>% filter (top != "ns" ) %>% group_by(top) %>%
summarize(chrpos= median(Chromosome),high = max(log10FDR)) 

#########################################################
colors = c( chr1 = "#A8CFF1", chr2 = "#9EC5C3", Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7",  others ="#ECEAEA")

theme<-theme(plot.margin=unit(c(1,1,1,1),units=,"mm"),
               strip.text = element_text(size=8),text = element_text(size = 8),
               axis.line = element_line(color = "black",linewidth = 0.25),
               axis.ticks =  element_line(color = "black",linewidth = 0.25),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(color="black",size=8),
               axis.text.x = element_text(color="black",size=8,angle = 0),
               panel.spacing.x = unit(0,"cm"),
               panel.border = element_blank(),
               panel.spacing = unit(0,"lines"),
               legend.position = "none")

##############################
##############################Ext Fig2c

p1<-ggplot() +
  geom_point(data=significance_gwas4,aes(x=Chromosome,y= log10FDR,color=as.factor(mark)),shape= 21,fill = "transparent",size=0.1,alpha=0.6,stroke = 0.01) + 
  geom_point(data=significance_gwas_chr,aes(x=Chromosome,y= log10FDR,color=as.factor(mark)),shape= 16,size=0.25,alpha=1,stroke = 0.2) +
  geom_point(data=significance_gwas3,aes(x=Chromosome,y= log10FDR,color=as.factor(tax)),shape= 16,size=0.25,alpha=0.6,stroke =0.2) +
  scale_color_manual(values= colors )  +
  scale_x_continuous(labels= axis$chr, breaks= axis$center, expand = c(0.01,0.01))+
  scale_y_continuous(expand=c(0,0))+
  geom_vline(xintercept = 70985273 ,color="black",linetype=2 ,linewidth = 0.01 )+
  labs(y="-log10(P value)") + theme +
  ggrepel::geom_text_repel(data=loci_ann,  aes(label=top,x=chrpos,y=high),  seed = 42,
  nudge_y = 5,angle = 90 ,
  direction = "x",
  segment.size = 0.1,
  box.padding = 0.01,size=2,max.overlaps = 40)

label ="Ext Fig2c"
width=10
height=3.5

ggsave(paste0(label,"labled.","mahhton_clear.png"), p1, width=width, height=height, units="in",dpi = 1200)
ggsave(paste0(label,"labled.","mahhton_clear.pdf"), p1, width=width, height=height, units="in",limitsize=FALSE)


###############################Ext Fig2d

library(ggplot2)
library(gggenes)
library(tidyverse)

colors = c("#4e3910", "#845d29", "#ae8548", "#e3c28b", "#4fb6ca", "#178f92", "#175f5d", "#054544")

gene_struc <- read.csv("R1_2_top1_genes_struc.csv",head = T)

p23 <- ggplot(data = gene_struc %>% filter(sturcture== "mRNA") , aes(y = chrm)) +
  geom_gene_arrow(aes(xmin = start, xmax = end,forward = dirict1,  fill = SYMBOL), , fill = "white",
  arrowhead_width = unit(4, "mm"),arrowhead_height=unit(6, "mm") , arrow_body_height=unit(6, "mm") ) +
  geom_subgene_arrow(data = gene_struc %>% filter(sturcture== "CDS"),
  aes(xmin = start, xmax = end, xsubmin = start, xsubmax = end,  fill = SYMBOL),
  arrowhead_width = unit(0, "mm"),arrowhead_height=unit(0, "mm") , arrow_body_height=unit(6, "mm")) +
   scale_fill_manual(values= colors ) +
  theme_genes() + 
  facet_wrap(~ locus_id, scales = "free", nrow = 1, shrink = FALSE) 


ggsave(paste0("Ext Fig2d.pdf"), p23, width=20, height=8, units="in")
ggsave(paste0("Ext Fig2d.png"), p23, width=20, height=8, units="in") 


#######################################Fig1f plot

p4<-ggplot() +
  geom_point(data=significance_gwas4,aes(x=Chromosome,y= log10FDR,color=as.factor(mark)),shape= 21,fill = "transparent",size=0.01,alpha=0.6,stroke = 0.01) + 
  geom_point(data=significance_gwas_chr,aes(x=Chromosome,y= log10FDR,color=as.factor(mark)),shape= 16,size=0.2,alpha=0.8,stroke = 0.2) +
  geom_point(data=significance_gwas3,aes(x=Chromosome,y= log10FDR,color=as.factor(tax)),shape= 16,size=0.2,alpha=0.6,stroke =0.2) +
  scale_color_manual(values= colors )  +  coord_cartesian(ylim = c(10,55))+
  scale_x_continuous(labels= axis$chr, breaks= axis$center, expand = c(0.01,0.01))+
  scale_y_continuous(expand=c(0,0))+
  labs(y="-log10(P value)") + theme 

label ="Fig1f"
width =10
height =1.5

ggsave(paste0(label,"no_lable.0.2","mahhton_clear.png"), p4, width=width, height=height, units="in",dpi = 1200)
ggsave(paste0(label,"no_lable.0.2","mahhton_clear.pdf"), p4,width=width, height=height, units="in",limitsize=FALSE)

#######################################Fig1h

library(ggplot2)
library(gggenes)
library(tidyverse)

colors = c("#4e3910", "#845d29", "#ae8548", "#e3c28b", "#4fb6ca", "#178f92", "#175f5d", "#054544")

gene_struc <- read.csv("R1_2_PAL1_STRC1.csv",head = T)

p23 <- ggplot(data = gene_struc %>% filter(sturcture== "mRNA") , aes(y = chrm)) +
  geom_gene_arrow(aes(xmin = start, xmax = end,forward = dirict1), , fill = "white",
  arrowhead_width = unit(4, "mm"),arrowhead_height=unit(6, "mm") , arrow_body_height=unit(6, "mm") ) +
  geom_subgene_arrow(data = gene_struc %>% filter(sturcture== "CDS"),
  aes(xmin = start, xmax = end, xsubmin = start, xsubmax = end),  fill =  "#4fb6ca",
  arrowhead_width = unit(0, "mm"),arrowhead_height=unit(0, "mm") , arrow_body_height=unit(6, "mm")) +
  geom_feature( data =gene_struc %>% filter(sturcture== "SNP"),aes(x = start, y = chrm), size =0.5, feature_height=unit(15, "mm"))+
  geom_feature_label(data =gene_struc %>% filter(sturcture== "SNP"),aes(x = start, y = chrm, label = SYMBOL), 
  label_height=unit(10, "mm"),feature_height=unit(15, "mm"),size=8, angle = 45)+
   scale_fill_manual(values= colors ) +
  theme_genes() 

ggsave(paste0("Fig1h.pdf"), p23, width=4, height=8, units="in")
ggsave(paste0("Fig1h.jpg"), p23, width=4, height=8, units="in") 

