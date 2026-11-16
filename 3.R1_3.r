
#####################################
#####################################
#####################################functions
###############################################
box_ANA <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,
                     colors = c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                                "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96"),level,
                                shapes=c(21,22,23,24),label="ana_"
){
  library(rstatix)
  library(ggplot2)
  library( ggpubr )
  library( dplyr )
  library( tidyr )

  dfn <- metadata[c(col_name, analysisgroup)]
  names(dfn) <- c("Content", "Group")
  dfn <- dfn%>%mutate(Group=factor(Group,levels =level))
  stat.test<- dfn %>% t_test(Content ~ Group)%>% add_xy_position(fun = "max", x = "Group")
  readr::write_delim(as.data.frame(stat.test),paste0(label,"_",analysisgroup,col_name,"-df-t-test.txt"))

  sumorder<-dfn %>%group_by(Group) %>%summarise(
    allmean = mean(Content, na.rm = TRUE),       
    sd1     = sd(Content, na.rm = TRUE),         
    SNP     = analysisgroup,                      
    n1      = sum(!is.na(Content)),               
    se1     = sd1 / sqrt(n1),                   
    .groups = "drop")

  sumorder_wide <- sumorder %>%  select( Group,SNP,allmean)  %>%
  pivot_wider(names_from = Group, values_from = c(allmean))
  write.csv(as.data.frame(sumorder_wide),paste0(label,"_",analysisgroup,".","-wider-means.csv"))

  sumorder_wide_SE<- sumorder %>%  select( Group,SNP, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
  write.csv(as.data.frame(sumorder_wide_SE),paste0(label,"_",analysisgroup,".","-wider-se.csv"))

  theme<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
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
               legend.position = "none"
  )
    
p1 <-ggplot(dfn,aes(Group,Content)) +
     geom_boxplot(aes(fill=Group,alpha = 0.8),alpha=1, linewidth=0.25, width = 0.8,outlier.shape = NA, outlier.size = 0)+
    geom_point(aes(fill=Group,alpha=0.7),colour ="black",alpha=0.7,size =2.5, shape = 21, stroke = 0.1,position = position_jitter(0.2))+    
     scale_fill_manual(values = colors)+
     stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F,tip.length = 0)+
     scale_size_continuous(range=c(1,3)) +  
     ylab(ylabs)+xlab("")+
     theme
p1
  
ggsave(paste0(label,"_",analysisgroup,col_name,".pdf"), p1, width=widths, height=heights, units="in")
ggsave(paste0(label,"_",analysisgroup,col_name,".png"), p1, width=widths, height=heights, units="in")

return(p1)
} 


#####################################

paopao_plot2 <- function (data1, cut_off_x = 0.05, cut_off_y = 1, COL = c("x", "y" ,"z","name","group"),xlab= "", ylab="", width=5, height=4, 
                          color,label="ana_",topn=5) {
  library(ggplot2)
  library(ggrepel)
  
  data = as.data.frame(data1[, COL], row.names = row.names(data1))
  colnames(data)<- c("x", "y" ,"z","name","group")
  theme1<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
               strip.text = element_text(size=8),text = element_text(size = 8),
               axis.line = element_line(color = "black",linewidth = 0.25),
               axis.ticks =  element_line(color = "black",linewidth = 0.25),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(color="black",size=8),
               axis.text.x = element_text(color="black",size=8,angle = 0),
               panel.spacing.x = unit(0,"cm"),
               panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       ) 
  # 绘制火山图
  p1 <- ggplot(data, aes(x = x, y = y, colour = group)) +
    geom_point(aes(size = z),shape= 21,size=1.8,fill = "transparent",alpha=0.8,stroke =0.5 ) +
    scale_color_manual(values = color) +
    geom_vline(xintercept = cut_off_x, lty = 4, col = "black", lwd = 0.25) +
    geom_hline(yintercept = cut_off_y, lty = 4, col = "black", lwd = 0.25) +
    labs(x = xlab, y = ylab) +
    theme1 +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank())
  p1


  top_data <- order(data$y, decreasing = TRUE)[1:topn]
  top_rows_up <- data[top_data, ]
  
  p2 <- p1 + geom_text_repel(
    data = top_rows_up,
    aes(label = name), size = 1,colour = "black",
    box.padding = unit(0.1, "lines"),max.overlaps = 100,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
 
  p2


  top_data <- order(data$x, decreasing = TRUE)[1:topn]
  top_rows_up <- data[top_data, ]
  
  p3 <- p2 + geom_text_repel(
    data = top_rows_up,
    aes(label = name), size = 1,colour = "red",
    box.padding = unit(0.1, "lines"),max.overlaps = 100,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
 
  p3  
  
  ggsave(paste0(label,".","volcano.pdf"), p1, width=width, height=height, units="in")
  ggsave(paste0(label,".","volcano.png"), p1, width=width, height=height,units="in")
  
  ggsave(paste0(label,"lable.","volcano.pdf"), p2, width=width, height=height, units="in")
  ggsave(paste0(label,"lable.","volcano.png"), p2, width=width, height=height,units="in")

  ggsave(paste0(label,"lable2.","volcano.pdf"), p3, width=width, height=height, units="in")
  ggsave(paste0(label,"lable2.","volcano.png"), p3, width=width, height=height,units="in")

  return(p2)
}


#####################################

paopao_plot4 <- function (data1, cut_off_x = 0.05, cut_off_y = 1, COL = c("v1", "v2" ,"pvalue","y","name","group"),xlab= "", ylab="", width=5, height=4, 
                          color,label="ana_",topn=5) {
  library(ggplot2)
  library(ggrepel)
  
  data = as.data.frame(data1[, COL], row.names = row.names(data1))
  colnames(data)<- c("v1", "v2" ,"FDR","y","name","group")

  data$FC <- data$v1/data$v2
  data$abslog2FC <- abs(log2(data$FC))
  data$log10FDR <- -log10(data$FDR)
  write.csv(data,paste0(label,".","foldcg.csv"))

  theme1<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
               strip.text = element_text(size=8),text = element_text(size = 8),
               axis.line = element_line(color = "black",linewidth = 0.25),
               axis.ticks =  element_line(color = "black",linewidth = 0.25),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(color="black",size=8),
               axis.text.x = element_text(color="black",size=8,angle = 0),
               panel.spacing.x = unit(0,"cm"),
               panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       ) 
  # 绘制火山图
  p1 <- ggplot(data, aes(x = y  , y = abslog2FC, colour = group)) +
    geom_point(aes(size = log10FDR),shape= 21,fill = "transparent",alpha=0.8,stroke =0.5 ) +
    scale_color_manual(values = color) +
    labs(x = xlab, y = ylab) +
    theme1 +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank())
  p1


  top_data <- order(data$abslog2FC, decreasing = TRUE)[1:topn]
  top_rows_up <- data[top_data, ]
  
  p2 <- p1 + geom_text_repel(
    data = top_rows_up,
    aes(label = name), size = 1,colour = "blue",
    box.padding = unit(0.1, "lines"),max.overlaps = 100,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
 
  p2
  
  top_data <- order(data$log10FDR, decreasing = TRUE)[1:topn]
  top_rows_up <- data[top_data, ]
  
  p3 <- p2 + geom_text_repel(
    data = top_rows_up,
    aes(label = name), size = 1,colour = "red",
    box.padding = unit(0.1, "lines"),max.overlaps = 100,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
 
  p3
  

  top_data <- order(data$y, decreasing = TRUE)[1:topn]
  top_rows_up <- data[top_data, ]
  
  p4 <- p3 + geom_text_repel(
    data = top_rows_up,
    aes(label = name), size = 1,colour = "black",
    box.padding = unit(0.1, "lines"),max.overlaps = 100,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
 
  p4

  ggsave(paste0(label,".","volcano.pdf"), p1, width=width, height=height, units="in")
  ggsave(paste0(label,".","volcano.png"), p1, width=width, height=height,units="in")
  
  ggsave(paste0(label,"lable.","volcano.pdf"), p2, width=width, height=height, units="in")
  ggsave(paste0(label,"lable.","volcano.png"), p2, width=width, height=height,units="in")
  return(p2)
}

#####################################

volcano_plot <- function (data1, cut_off_FDR = 0.05, cut_off_log2FC = 1, COL = c("v1", "v2" ,"Pvalue","name"),width=5, height=4,
                          color=c(Up="#FF690B",Down="#4F95FF",no="#FDFDFD"),label="ana_",topn=5) {
  library(ggplot2)
  library(ggrepel)
  
  data = as.data.frame(data1[, COL], row.names = row.names(data1))
  colnames(data)<- c("v1", "v2" ,"FDR","gene_name")
  data$FC <- data$v1/data$v2
  data$log2FC <- log2(data$FC)
  data$log10FDR <- log10(data$FDR)
  # 设定差异基因阈值
  
  data$Sig = ifelse(data$FDR < cut_off_FDR & abs(data$log2FC) >= cut_off_log2FC,
                    ifelse(data$log2FC > cut_off_log2FC, 'Up', 'Down'), 'no')
  
  producedata = as.data.frame(data[, !(names(data) %in% c("FDR","gene_name"))], row.names = row.names(data))
  newdata <- cbind(data1,producedata)
  
  write.csv(as.data.frame(newdata), file =paste0(label,"plot_OUT.csv") , row.names = FALSE)
  # 绘制火山图
theme1<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
               strip.text = element_text(size=8),text = element_text(size = 8),
               axis.line = element_line(color = "black",linewidth = 0.25),
               axis.ticks =  element_line(color = "black",linewidth = 0.25),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(color="black",size=8),
               axis.text.x = element_text(color="black",size=8,angle = 0),
               panel.spacing.x = unit(0,"cm"),
               panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       )


  p1 <- ggplot(data, aes(x = log2FC, y = -log10FDR)) +
    geom_point( aes( size = abs(-log10FDR), colour = Sig), alpha = 0.8, stroke = 0.3) +
    scale_color_manual(values = c(Up="#FF690B",Down="#4F95FF",no="#DFDBDB")) +
    geom_vline(xintercept = c(-cut_off_log2FC, cut_off_log2FC), lty = 4, col = "black", lwd = 0.25) +
    geom_hline(yintercept = -log10(cut_off_FDR), lty = 4, col = "black", lwd = 0.25) +
    labs(x = "log2(Flod change)", y = "-log10FDR") +
    theme1
 
  p1

  up <- data[data$Sig  == 'Up', ]
  top_5_up <- order(up$log10FDR, decreasing = FALSE)[1:topn]
  top_5_rows_up <- up[top_5_up, ]
  
  down <- data[data$Sig  == 'Down', ]
  top_5_down <- order(down$log10FDR, decreasing = FALSE)[1:topn]
  top_5_rows_down <- down[top_5_down, ]
  
  
  p2 <- p1 + geom_text_repel(
    data = top_5_rows_up,
    aes(label = gene_name), size = 1,colour = "black",
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
  p2 <-p2 + geom_text_repel(
    data = top_5_rows_down,
    aes(label = gene_name), size = 1,colour = "black",
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
  
  p2
  
  
  ggsave(paste0(label,".","volcano.pdf"), p1, width=width, height=height, units="in")
  ggsave(paste0(label,".","volcano.png"), p1, width=width, height=height, units="in")
  
  ggsave(paste0(label,"lable.","volcano.pdf"), p2, width=width, height=height, units="in")
  ggsave(paste0(label,"lable.","volcano.png"), p2, width=width, height=height, units="in")
  
  return(p2)
}
##########################################

Regression_sample2 <-function(data,ana_col,level,xlabs="x",ylabs="y", width=6.0, height=8.2,label="lab",shapes = c(21,22,23,24)
,colors=c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E")){

library(ggplot2)
sampFiletax2<-as.data.frame(data[, ana_col], row.names= row.names(data))
colnames(sampFiletax2)<-c("x","y","Group","shape")
sampFiletax2$Group = factor(sampFiletax2$Group , levels = level)
theme1<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
               strip.text = element_text(size=8),text = element_text(size = 8),
               axis.line = element_line(color = "black",linewidth = 0.25),
               axis.ticks =  element_line(color = "black",linewidth = 0.25),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(color="black",size=8),
               axis.text.x = element_text(color="black",size=8,angle = 0),
               panel.spacing.x = unit(0,"cm"),
               panel.border = element_blank(),
       panel.spacing = unit(0,"lines"),,legend.position = "top",plot.title = element_text(hjust = 0 , vjust = - 25 ,size=8),
      legend.background = element_blank(),legend.key = element_blank(),legend.title=element_blank())

library(psych)
ct <- corr.test(sampFiletax2$x,sampFiletax2$y, method = "pearson")
corr_r =ct$r
corr_p = ct$p.adj

p21 <-ggplot(sampFiletax2,aes(x,y)) +
geom_point(aes(x,y,colour = Group,fill=Group,shape = shape), alpha=0.65, size =2,  stroke = 0.2)+
scale_colour_manual(values = colors)+
scale_fill_manual(values = colors)+
scale_shape_manual(values = shapes)+
stat_ellipse(aes(color = Group),level = 0.68,linewidth=0.1) +

geom_smooth(method = "lm", formula = NULL,linewidth=0.2,color="black",linetype="solid",aes(group=1))+
ggtitle(paste( " P = ",sprintf("%.4e",corr_p),   " R = ", sprintf("%.4f",corr_r), sep = ""))+
ylab(ylabs)+xlab(xlabs)+theme1

p21
ggsave(paste(label,"_",ana_col[1],"_",ana_col[2],"_","Line_pearson_R",".pdf"), p21, width=width, height=height, units="in")
ggsave(paste(label,"_",ana_col[1],"_",ana_col[2],"_","Line_pearson_R",".jpg"),p21, width=width, height=height, units="in")
}


###############################################################
##################################################Fig plot
##################################################Fig plot
#########################################Fif1g

data1=read.csv("R1_3_top1_beta2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

COL = c("R2" ,"log10FDR", "Pvalue","snp_name","tax")
color =  c(Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7")
xlab= "Efect of microbiome"
ylab= "Gwas association"
cut_off_x = 0.05
cut_off_y = 50
label="Fif1g"
width=4
height=4
topn=1

paopao_plot2 (data1, cut_off_x, cut_off_y, COL ,xlab=xlab, ylab=ylab, 
width=width, height=height, color = color ,label,topn=topn)

#########################################Sup Fig2a

data1=read.csv("R1_3_all_cds_beta.csv", header=T,  comment.char="", stringsAsFactors=F)
COL = c("R2" ,"log10FDR", "Pvalue","snp_name","tax")
color =  c(Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7")
xlab= "Efect size of phyllosphere microbiome (R2)"
ylab= "Gwas association"
cut_off_x = 0.1
cut_off_y = 50
label="Sup Fig2a"
width=5
height=5
topn=3

paopao_plot2 (data1, cut_off_x, cut_off_y, COL ,xlab=xlab, ylab=ylab, 
width=width, height=height, color = color ,label,topn=topn)


#########################################Sup Fig2b

data1=read.csv("R1_3_all_cds_beta.csv", header=T,  comment.char="", stringsAsFactors=F)

analysisgroup = "tax"
col_name = "R2"
ylabs = "Efect size of phyllosphere microbiome (R2)"
widths = 4
heights = 5
level = c("Pseudomonadales","Enterobacterales")
colors =c(Pseudomonadales = "#D7D540",Enterobacterales="#699BD7" )
labels="Sup Fig2b"

box_ANA (data1,analysisgroup , col_name ,ylabs ,widths = widths,heights =heights,colors =colors,level=level ,label=labels)

###################Sup Fig3a

data1=read.csv("R1_3_all_chao1.csv", header=T,  comment.char="", stringsAsFactors=F)
COL = c("G1_mean",	"G2_mean","Pvalue","log10FDR", "snp_name","tax")
color =  c(Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7")
ylab= "Efect of alpha diversity (chao1)"
xlab= "Gwas association"
cut_off_x = 0.01
cut_off_y = 50

width=6
height=5
topn=3
label="Sup Fig3a"

paopao_plot4 (data1, cut_off_x, cut_off_y, COL ,xlab=xlab, ylab=ylab, 
width=width, height=height, color = color ,label,topn=topn)


#############################Sup Fig3b

data1=read.csv("Sup Fig3a.foldcg.csv", header=T,  comment.char="", stringsAsFactors=F)
analysisgroup = "group"
col_name = "abslog2FC"
label="Sup Fig3b"
ylabs = "Efect of alpha diversity (chao1)"
widths = 4
heights = 5
level = c("Pseudomonadales","Enterobacterales")
colors =c(Pseudomonadales = "#D7D540",Enterobacterales="#699BD7" )

box_ANA (data1,analysisgroup , col_name ,ylabs ,widths = widths,heights =heights,colors =colors,level=level ,label=label)

###############################Sup Fig3c

data1=read.csv("R1_3_all_shannon.csv", header=T,  comment.char="", stringsAsFactors=F)
COL = c("G1_mean",	"G2_mean","Pvalue","log10FDR", "snp_name","tax")
color =  c(Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7")
ylab= "Efect of alpha diversity (shannon)"
xlab= "Gwas association"
cut_off_x = 0.01
cut_off_y = 50
label="Sup Fig3c"
width=6
height=5
topn=3

paopao_plot4 (data1, cut_off_x, cut_off_y, COL ,xlab=xlab, ylab=ylab, 
width=width, height=height, color = color ,label,topn=topn)

######################################Sup Fig3d

data1=read.csv("Sup Fig3c.foldcg.csv", header=T,  comment.char="", stringsAsFactors=F)
analysisgroup = "group"
col_name = "abslog2FC"
label="Sup Fig3d"
ylabs = "Efect of alpha diversity (shannon)"
widths = 4
heights = 5
level = c("Pseudomonadales","Enterobacterales")
colors =c(Pseudomonadales = "#D7D540",Enterobacterales="#699BD7" )

box_ANA (data1,analysisgroup , col_name ,ylabs ,widths = widths,heights =heights,colors =colors,level=level ,label=label)

#########################################Sup Fig4a

data1=read.csv("R1_3_all_psuedo.csv", header=T,  comment.char="", stringsAsFactors=F)
COL = c("G1_mean",	"G2_mean","Pvalue","log10FDR", "snp_name","tax")
color =  c(Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7")
ylab= "Efect of Pseudomonadales"
xlab= "Gwas association"
cut_off_x = 0.01
cut_off_y = 50
label="Sup Fig4a"
width=6
height=5
topn=3

paopao_plot4 (data1, cut_off_x, cut_off_y, COL ,xlab=xlab, ylab=ylab, 
width=width, height=height, color = color ,label,topn=topn)


###########################################Sup Fig4b

data1=read.csv("Sup Fig4a.foldcg.csv", header=T,  comment.char="", stringsAsFactors=F)
analysisgroup = "group"
col_name = "abslog2FC"
label="Sup Fig4b"
ylabs = "Efect of Pseudomonadales"
widths = 4
heights = 5
level = c("Pseudomonadales","Enterobacterales")
colors =c(Pseudomonadales = "#D7D540",Enterobacterales="#699BD7" )

box_ANA (data1,analysisgroup , col_name ,ylabs ,widths = widths,heights =heights,colors =colors,level=level ,label=label)


#########################################Sup Fig4c

data1=read.csv("R1_3_all_enter.csv", header=T,  comment.char="", stringsAsFactors=F)
COL = c("G1_mean",	"G2_mean","Pvalue","log10FDR", "snp_name","tax")
color =  c(Pseudomonadales = "#D7D540",Pantoe="#c3d6ce",Enterobacterales="#699BD7")
ylab= "Efect of Enterobacterales"
xlab= "Gwas association"
cut_off_x = 0.01
cut_off_y = 50
label="Sup Fig4c"
width=6
height=5
topn=3

paopao_plot4 (data1, cut_off_x, cut_off_y, COL ,xlab=xlab, ylab=ylab, 
width=width, height=height, color = color ,label,topn=topn)


##########################################Sup Fig4d
data1=read.csv("Sup Fig4c.foldcg.csv", header=T,  comment.char="", stringsAsFactors=F)
analysisgroup = "group"
col_name = "abslog2FC"
label="Sup Fig4d"
ylabs = "Efect of Enterobacterales"
widths = 4
heights = 5
level = c("Pseudomonadales","Enterobacterales")
colors =c(Pseudomonadales = "#D7D540",Enterobacterales="#699BD7" )

box_ANA (data1,analysisgroup , col_name ,ylabs ,widths = widths,heights =heights,colors =colors,level=level ,label=label)


###################Fig1j
library(tidyverse)

data1=read.csv("R1_3_order_effect.csv", header=T,  comment.char="", stringsAsFactors=F)
data1 <- data1  %>%  mutate(meanra=(hapA+hapD)/2)
data1 <- data1 %>% filter(meanra >= 0.01  )

COL = c("hapA" ,"hapD", "Pvalue","tax")
label="Fig1j"
width=4
height=5
topn=2

volcano_plot(data1, cut_off_FDR = 0.01, cut_off_log2FC = 0, COL = COL, 
width=width, height=height,label=label,topn=topn)

####################Sup Fig5a-c
Metabolic=read.csv("R1_3_metabolite.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
Metabo <- Metabolic[, colnames(Metabolic) %in% c("tree_id","chr.2.25894232","HCA","FLA","SPA")]

Metabo <- Metabo %>% filter( chr.2.25894232 != "N"  )

analysisgroup = "chr.2.25894232"
col_name = "HCA"
ylabs = "Concentration (mg/g)"
widths = 4
heights = 4
level = c("HapA","HapD")
colors =c(HapA= "#6ED470" ,  HapD = "#9EA4E5" )

labels="Sup Fig5a"
box_ANA (Metabo,analysisgroup , col_name ,ylabs ,widths = widths,heights =heights,colors =colors,level=level ,label=labels)

labels="Sup Fig5b"
col_name = "FLA"
box_ANA (Metabo,analysisgroup , col_name ,ylabs ,widths = widths,heights =heights,colors =colors,level=level ,label=labels)

labels="Sup Fig5c"
col_name = "SPA"
box_ANA (Metabo,analysisgroup , col_name ,ylabs ,widths = widths,heights =heights,colors =colors,level=level ,label=labels)


#########Sup Fig5d

data1=read.csv("R1_3_metabolic_effect.csv", header=T,  comment.char="", stringsAsFactors=F)
COL = c("hapA" ,"hapD", "Pvalue","metabolic")
label="Sup Fig5d"
width=4
height=5
topn=3

volcano_plot(data1, cut_off_FDR = 0.01, cut_off_log2FC = 0, COL = COL, 
width=width, height=height,label=label,topn=topn)

###########################Fig1k

Corre <-read.csv("R1_3_correlation.csv", header=T, row.names=1, comment.char="")

width=3
height=5
ana_col<-c("HCA","Pseudomonadales","hap1","group3")
xlabs= "Concentration (mg/g)"
ylabs="Relative abundance (%)"

colors <-c(A= "#6ED470" ,  D = "#9EA4E5")
level = c("D","A")
shapes = c(Wild=21,JAP=21,IND=21)
label="Fig1k"

Regression_sample2 (Corre ,ana_col,level,xlabs=xlabs,ylabs=ylabs, width=width, height=height,label=label,shapes = shapes,colors=colors)










