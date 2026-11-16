
################################################################################################
##############################################
#####################################################################

box_ANA2 <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,
                     colors = c("#6ED470","#9EA4E5","#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                                "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96"),
                                shapes=c(21,22,23,24),label="ana_"
){
  library(rstatix)
  library(ggplot2)
  library( ggpubr )
  library( dplyr )
  library( tidyr )

  dfn <- metadata[c(col_name, analysisgroup)]
  names(dfn) <- c("variable", "Group")
  dfn$variable <- as.numeric(dfn$variable)
  stat.test<- dfn %>% t_test(variable ~ Group) %>% add_significance(p.col = "p") 

  test.p = paste(analysisgroup ,"G1:", stat.test$group1,"G2:",stat.test$group2,  "Pvalue:", stat.test$p,"N1:",stat.test$n1,"N2:",stat.test$n2,sep = ",")
  write.table(test.p,paste0(label,".",col_name,"ALL","-t_test.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  

  sumorder<-dfn %>%group_by(Group) %>%summarise(
    allmean = mean(variable, na.rm = TRUE),       
    sd1     = sd(variable, na.rm = TRUE),         
    SNP     = analysisgroup,                      
    n1      = sum(!is.na(variable)),               
    se1     = sd1 / sqrt(n1),                   
    .groups = "drop")


  sumorder_wide <- sumorder %>%  select( Group,SNP,allmean)  %>%
  pivot_wider(names_from = Group, values_from = c(allmean))
  names(sumorder_wide) <- c("SNP", "G1", "G2")
  write.csv(as.data.frame(sumorder_wide),paste0(label,".",analysisgroup,".","-wider-means.csv"))

  means = paste("SNP:",sumorder_wide$SNP ,"G1:", sumorder_wide$G1,"G2:",sumorder_wide$G2,sep = ",")
  write.table(means, paste0(label,".",col_name,"ALL","-mean.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

  sumorder_wide_SE<- sumorder %>%  select( Group,SNP, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
  names(sumorder_wide_SE) <- c("SNP", "G1", "G2")
  write.csv(as.data.frame(sumorder_wide_SE),paste0(label,".",analysisgroup,".","-wider-se.csv"))

  SE = paste("SNP:",sumorder_wide_SE$SNP ,"G1:", sumorder_wide_SE$G1,"G2:",sumorder_wide_SE$G2,sep = ",")
  write.table(SE, paste0(label,".",col_name,"ALL","-se.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  
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
    
p1 <-ggplot(dfn,aes(Group,variable)) +
     geom_boxplot(aes(fill=Group,alpha = 0.8),alpha=1, linewidth=0.25, width = 0.8,outlier.shape = NA, outlier.size = 0)+
    geom_point(aes(fill=Group,alpha=0.7),colour ="black",alpha=0.7,size =2.5, shape = 21, stroke = 0.1,position = position_jitter(0.2))+    
     scale_fill_manual(values = colors)+
    ggtitle(paste(analysisgroup , "Pvalue:", stat.test$p,ssep = ","))+
#     stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F,tip.length = 0)+
     scale_size_continuous(range=c(1,3)) +  
     ylab(ylabs)+xlab("")+
#     ggtitle(paste(label,"_",col_name ))+
     theme
p1
  
ggsave(paste0(label,"_",analysisgroup,col_name,".pdf"), p1, width=widths, height=heights, units="in")
ggsave(paste0(label,"_",analysisgroup,col_name,".png"), p1, width=widths, height=heights, units="in")

return(p1)
} 


#####################################################################

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
  
  write.csv(as.data.frame(newdata), file =paste0(label,"FDR_OUT.csv") , row.names = FALSE)
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

################################################################################################

paopao_plot3 <- function (data1,  COL = c("v1", "v2" ,"p","r2","name","group"),xlab= "", ylab="", width=5, height=4, 
                          color, shape, label="ana_",topn=5) {
  library(ggplot2)
  library(ggrepel)
  library( dplyr )
  library( tidyr )
  
  data = as.data.frame(data1[, COL], row.names = row.names(data1))
  colnames(data)<- c("v1", "v2" ,"p","r2","name","group")
  data <- data %>% mutate( FC = abs(log2(v1/v2))) %>% mutate( log10P =-log10(p))
  write.csv(as.data.frame( data),paste0(label, ".",analysisgroup,".","plot.csv"),append =T)

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
  p1 <- ggplot(data, aes(x = r2, y = log10P, colour = group, shape= group)) +
    geom_point(aes(size = FC),fill = "transparent",alpha=0.8,stroke =0.5 ) +
    scale_color_manual(values = color) +
    scale_shape_manual(values = shape) +
    #  xlim(c(-30, 30)) +
    labs(x = xlab, y = ylab) +
    theme1 +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank())
  p1


  top_data <- order(data$r2, decreasing = TRUE)[1:topn]
  top_rows_up <- data[top_data, ]
  
  p2 <- p1 + geom_text_repel(
    data = top_rows_up,
    aes(label = name,colour = group), size = 1,
    nudge_x = 0.2,
    direction =  "y" ,
    box.padding = unit(0.1, "lines"),max.overlaps = 100,
    point.padding = unit(0, "lines"), segment.size=0.25 ,segment.color = "black", show.legend = FALSE)
 
  p2


  top_data <- order(data$log10P, decreasing = TRUE)[1:topn]
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

#####################################################data
######################################################fig visulization
######################################################fig visulization

library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggspatial)

map_data <- st_read("R4_RiceProduction_v1.shp")
map_data <- st_set_crs(map_data, 4326)
map_data <- st_transform(map_data,4326)

worlds_rice = read.csv("R4_rice_point2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
anno_sf <- st_as_sf(worlds_rice,coords = c("lon", "lat"),crs = 4326)
ids <- row.names(worlds_rice)
world <- ne_countries(scale = "small", returnclass = "sf")

theme2<-
theme(      
  axis.text.y = element_text(color="black",size=8),
  axis.text.x = element_text(color="black",size=8,angle = 0),
  panel.background = element_rect(fill = "#81D6F0", color = NA),
  panel.grid.major = element_line(color = "#F9F9F9",alpha=0.8),
  panel.grid.minor = element_line(color = "#FCFEFE",alpha=0.8))

p2<-
ggplot()+
geom_sf(data = world,fill="black",size=0,color="black")+
geom_sf(data = map_data ,aes(fill=P_total,color=P_total),stroke=0)+
scale_fill_gradient(low ="#000000",high = "#FFF23B")+
scale_color_gradient(low ="#000000",high = "#FFF23B")+
theme2

########################################################
########################################################Ext Fig3
shapes = c(HapA=21, HapD=24)
colors <- c(HapA= "#6ED470" ,  HapD = "#9EA4E5")

p3<- p2+
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice, aes(x=lon, y = lat, fill = Hap,shape= Hap ) ,
size= 0.8,colour='black',position = position_jitter(width = 1, height = 1) , stroke=.01)+
scale_fill_manual(values = colors) +
scale_shape_manual(values = shapes ) +
labs(x="Longitude",y="Latitude") 

#########################################################Ext Fig3a
label="Ext Fig3a"

p6 <- p3 +
coord_sf( ylim = c( -55,75),xlim = c( -150,150), crs = 4326)

ggsave(paste0(label,".","rice_product_area3-1.","-maps.pdf"), p6 , width=10, height=7, units="in",limitsize = FALSE)
ggsave(paste0(label,".","rice_product_area3-1.","-maps.png"), p6 , width=10, height=7, units="in",dpi = 1200, limitsize = FALSE)

#############################
#############################procution_rufipogon_nivala
#############################Ext Fig3d

shapes = c(rufipogon=24, nivala=21)

colors <- c(rufipogon="#969bc7",JAP ="#e9a00e", IND="#3c7c3d",nivala ="#2080E8")
label="Ext Fig3d"

p7<- p2 +
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice  %>% filter(subgroup == "rufipogon" | subgroup == "nivala" ) , 
aes(x=lon, y = lat,shape= subgroup ,fill = subgroup) ,
size= 1.5, colour='black',position = position_jitter(width = 0.5, height = 0.5) , stroke=.01)+
scale_fill_manual(values = colors) +
scale_shape_manual(values = shapes ) +
labs(x="Longitude",y="Latitude") +
coord_sf( ylim = c( 0,30),xlim = c( 70,120), crs = 4326)

ggsave(paste0(label,".","rice_product_area7.","-maps.pdf"), p7 , width=6, height=5, units="in",limitsize = FALSE)
ggsave(paste0(label,".","rice_product_area7.","-maps.png"), p7 , width=6, height=5, units="in",dpi = 1200, limitsize = FALSE)


############################Ext Fig3b

library(dplyr)
trail=read.csv("R4_all_geo_climeate2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

trail[is.na(trail)]<-"N" 
new_tra <- trail %>% filter( HAP != "N"  )
analysisgroup = "HAP"
widths = 3
heights = 4

colors <-c(HapA= "#6ED470" ,  HapD = "#9EA4E5")

col_names = "Longitude"
metadata2 <- new_tra %>% filter(Longitude != "N"  )
ylabs = col_names
col_name = col_names
label<-paste("Ext Fig3b","_",col_name )
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,colors = colors,
heights =heights)

############################Ext Fig3c

library(dplyr)
trail=read.csv("R4_all_geo_climeate2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

trail[is.na(trail)]<-"N" 
new_tra <- trail %>% filter( HAP != "N"  )
analysisgroup = "HAP"
widths = 3
heights = 4

colors <-c(HapA= "#6ED470" ,  HapD = "#9EA4E5")

col_names = "Latitude"
metadata2 <- new_tra %>% filter(Latitude != "N"  )
ylabs = col_names
col_name = col_names
label<-paste("Ext Fig3c","_",col_name )
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,colors = colors,
heights =heights)


##############################################################Fig4a

colors <-c(HapA= "#6ED470" ,  HapD = "#9EA4E5")

library(tidyverse)
library(ggplot2)

all_rice_ad <- read.csv("R4_ad_structure.csv",head = T)
width=4
height=4

p7 = ggplot(all_rice_ad , aes(x = reorder(subgroup, level), y = rice, fill = HAP)) + 
   geom_bar(stat = "identity", position = "fill", width = 0.8,colour ="white",linewidth=0.1) + 
    xlab("") + ylab("Haplotype ratio of OsPAL1 (%)") +
    scale_fill_manual(values = colors)+
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

ggsave(paste0("Fig4a.","rice_ad_group.pdf"), p7, width=width, height=height,units="in")
ggsave(paste0("Fig4a.","rice_ad_group.jpg"), p7, width=width, height=height,units="in") 

#############################################################Fig4b

data1=read.csv("R4_climate_effect.csv", header=T,  comment.char="", stringsAsFactors=F)

COL = c("HapA" ,"HapD", "Pvalue","Discribe")
label="Fig4b"
width=4
height=4
topn=3

volcano_plot(data1, cut_off_FDR = 0.05, cut_off_log2FC = 0, COL = COL, 
width=width, height=height,label=label,topn=topn)

#############################################################Fig4c

data1=read.csv("R4_climate_effect2.csv", header=T,  comment.char="", stringsAsFactors=F)

COL = c("nivala" ,"rufipogon", "Pvalue","Discribe")
label="Fig4c"
width=4
height=4
topn=3

volcano_plot(data1, cut_off_FDR = 0.05, cut_off_log2FC = 0, COL = COL, 
width=width, height=height,label=label,topn=topn)

#####################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(terra)
library(raster)
library(sf)
library(tibble)

worlds_rice = read.csv("R4_rice_point2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
anno_sf <- st_as_sf(worlds_rice,coords = c("lon", "lat"),crs = 4326)

####################################

mean_df = read.csv("R4_world_average_solar_5m.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

theme3<-
theme(axis.text.y = element_text(color="black",size=10),
  axis.text.x = element_text(color="black",size=10,angle = 0),
  panel.background = element_rect(fill = "transparent", color = NA),
  panel.grid.major = element_line(color = "#000000",alpha=0.8,linewidth =0.1),
  panel.grid.minor = element_line(color = "#000000",alpha=0.8,linewidth =0.1))

p0<-
ggplot() +
  geom_raster(data = mean_df, aes(x = x, y = y, fill = solar_rad)) +
  scale_fill_gradientn(name = "kJ m⁻² day⁻¹",
  colours = c("#000000",  "#6baed6", "#fcbf49", "#e7298a") )+
  coord_equal() +
  theme3+
  labs( x = "Longitude", y = "Latitude")

shapes = c(HapA=21, HapD=24)
colors <- c(HapA= "#6ED470" ,  HapD = "#9EA4E5")

p3<- p0+
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice, aes(x=lon, y = lat, fill = Hap,shape= Hap ) ,
size= 1.2,colour='black',position = position_jitter(width = 1, height = 1) , stroke=.01)+
scale_fill_manual(values = colors) +
scale_shape_manual(values = shapes ) +
labs( x = "Longitude", y = "Latitude") 

##############################################Sup Fig9a
p6 <- p3 +
coord_sf( ylim = c( -55,75),xlim = c( -150,150), crs = 4326)

ggsave(paste0("Sup Fig9a.","rice_Solar_world_wide.","-maps.pdf"), p6 , width=10, height=7, units="in",limitsize = FALSE)
ggsave(paste0("Sup Fig9a.","rice_Solar_world_wide.","-maps.png"), p6 , width=10, height=7, units="in",dpi = 1200, limitsize = FALSE)

###############################################
###############################################Fig4d
p00<-
ggplot() +
  geom_raster(data = mean_df, aes(x = x, y = y, fill = solar_rad)) +
  scale_fill_gradientn(name = "kJ m⁻² day⁻¹",
  colours = c("#000000",  "#050505") )+
  coord_equal() +
  theme3+
  labs( x = "Longitude", y = "Latitude")

p4 <- p00 +
  ggnewscale::new_scale_fill() +
  stat_density_2d(
  data = worlds_rice %>% filter(Hap == "HapA"),
  aes(x = lon, y = lat, fill = after_stat(level), alpha = (after_stat(level)/2) ),  
#  bins = 10,
  geom = "polygon",
  contour_var = "ndensity",
  show.legend =   TRUE) +
  scale_fill_gradient(low = "#6ED470", high = "#2D722E", name = "Density (HapA)")

p4 <- p4 +
  ggnewscale::new_scale_fill() +
  stat_density_2d(
    data = worlds_rice %>% filter(Hap == "HapD"),
    aes(x = lon, y = lat, fill = after_stat(level), alpha =(after_stat(level)/2) ),
#    bins = 10,
    geom = "polygon",
    contour_var = "ndensity",
    show.legend =   TRUE) + coord_equal() + 
  scale_fill_gradient(low = "#C8CBF2", high = "#3F4BF2", name = "Density (HapD)")

p4 <- p4 + coord_sf( ylim = c( -55,75),xlim = c( -150,150), crs = 4326)
p7 <- p4 + coord_sf( ylim = c( 0,35),xlim = c( 60,130), crs = 4326)

ggsave(paste0("Fig4d.","rice_Solar_ad_intensity.","-maps.pdf"), p7 , width=10, height=7, units="in",limitsize = FALSE)
ggsave(paste0("Fig4d.","rice_Solar_ad_intensity.","-maps.png"), p7 , width=10, height=7, units="in",dpi = 1200, limitsize = FALSE)

########################################Fig4e

shapes = c(rufipogon=24, nivala=21)
colors <- c(rufipogon="#969bc7",JAP ="#e9a00e", IND="#3c7c3d",nivala ="#2080E8")

p5<- p0 +
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice  %>% filter(subgroup == "rufipogon" | subgroup == "nivala" ) , 
aes(x=lon, y = lat,shape= subgroup ,fill = subgroup) ,
size= 1.5, colour='black',position = position_jitter(width = 0.5, height = 0.5) , stroke=.01)+
scale_fill_manual(values = colors) +
scale_shape_manual(values = shapes ) +
labs(x = "Longitude", y = "Latitude") 

p6 <- p5 +
coord_sf( ylim = c( 0,35),xlim = c( 60,130), crs = 4326)

ggsave(paste0("Fig4e.","rice_Solar_nivala_rugipogon.","-maps.pdf"), p6 , width=10, height=7, units="in",limitsize = FALSE)
ggsave(paste0("Fig4e.","rice_Solar_nivala_rugipogon.","-maps.png"), p6 , width=10, height=7, units="in",dpi = 1200, limitsize = FALSE)


######################################################temprature
####################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(terra)
library(raster)
library(sf)
library(tibble)

worlds_rice = read.csv("R4_rice_point2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
anno_sf <- st_as_sf(worlds_rice,coords = c("lon", "lat"),crs = 4326)

mean_df = read.csv("R4_world_average_temp_5m.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

theme3<-
theme(axis.text.y = element_text(color="black",size=10),
  axis.text.x = element_text(color="black",size=10,angle = 0),
  panel.background = element_rect(fill = "transparent", color = NA),
  panel.grid.major = element_line(color = "#000000",alpha=0.8,linewidth =0.1),
  panel.grid.minor = element_line(color = "#000000",alpha=0.8,linewidth =0.1))

p0<-
ggplot() +
  geom_raster(data = mean_df, aes(x = x, y = y, fill = wc2.1_5m_bio_1)) +
  scale_fill_gradientn(name = "°C",
  colours = c("#000000",  "#F4F8FB", "#e7298a") )+
  coord_equal() +
  theme3+
  labs( x = "Longitude", y = "Latitude")

#########################################

shapes = c(HapA=21, HapD=24)
colors <- c(HapA= "#6ED470" ,  HapD = "#9EA4E5")

p3<- p0+
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice, aes(x=lon, y = lat, fill = Hap,shape= Hap ) ,
size= 1.2,colour='black',position = position_jitter(width = 1, height = 1) , stroke=.01)+
scale_fill_manual(values = colors) +
scale_shape_manual(values = shapes ) +
labs( x = "Longitude", y = "Latitude") 

#################################################Sup Fig9b
p6 <- p3 +
coord_sf( ylim = c( -55,75),xlim = c( -150,150), crs = 4326)
label = "Sup Fig9b."

ggsave(paste0(label,"-maps.pdf"), p6 , width=10, height=7, units="in",limitsize = FALSE)
ggsave(paste0(label,"-maps.png"), p6 , width=10, height=7, units="in",dpi = 1200, limitsize = FALSE)

#################################################Sup Fig9d
shapes = c(rufipogon=24, nivala=21)
colors <- c(rufipogon="#969bc7",JAP ="#e9a00e", IND="#3c7c3d",nivala ="#2080E8")

p5.1<- p0 +
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice  %>% filter(subgroup == "rufipogon" | subgroup == "nivala" ) , 
aes(x=lon, y = lat,shape= subgroup ,fill = subgroup) ,
size= 1.5, colour='black',position = position_jitter(width = 0.5, height = 0.5) , stroke=.01)+
scale_fill_manual(values = colors) +
scale_shape_manual(values = shapes ) +
labs(x = "Longitude", y = "Latitude") 

p6.1 <- p5.1 +
coord_sf( ylim = c( 0,35),xlim = c( 70,120), crs = 4326)
label = "Sup Fig9d."

ggsave(paste0(label,"-maps.pdf"), p6.1 , width=5, height=7, units="in",limitsize = FALSE)
ggsave(paste0(label,"-maps.png"), p6.1 , width=5, height=7, units="in",dpi = 1200, limitsize = FALSE)


######################################################precipitation precipitation
####################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(terra)
library(raster)
library(sf)
library(tibble)

worlds_rice = read.csv("R4_rice_point2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
anno_sf <- st_as_sf(worlds_rice,coords = c("lon", "lat"),crs = 4326)

mean_df = read.csv("R4_world_average_precip_5m.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

theme3<-
theme(axis.text.y = element_text(color="black",size=10),
  axis.text.x = element_text(color="black",size=10,angle = 0),
  panel.background = element_rect(fill = "transparent", color = NA),
  panel.grid.major = element_line(color = "#000000",alpha=0.8,linewidth =0.1),
  panel.grid.minor = element_line(color = "#000000",alpha=0.8,linewidth =0.1))

p0<-
ggplot() +
  geom_raster(data = mean_df, aes(x = x, y = y, fill = wc2.1_5m_bio_12)) +
#  scale_fill_viridis_c(name = "kJ m-2 day-1",option = "D" )+
  scale_fill_gradientn(name = "mm",
  colours = c("#000000", "#FBFBFB", "#48BBFE") )+
  coord_equal() +
  theme3+
  labs( x = "Longitude", y = "Latitude")

##################################################

shapes = c(HapA=21, HapD=24)
colors <- c(HapA= "#6ED470" ,  HapD = "#9EA4E5")

p3<- p0+
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice, aes(x=lon, y = lat, fill = Hap,shape= Hap ) ,
size= 1.2,colour='black',position = position_jitter(width = 1, height = 1) , stroke=.01)+
scale_fill_manual(values = colors) +
scale_shape_manual(values = shapes ) +
labs( x = "Longitude", y = "Latitude") 

#################################################Sup Fig9c
p6 <- p3 +
coord_sf( ylim = c( -55,75),xlim = c( -150,150), crs = 4326)
label = "Sup Fig9c."

ggsave(paste0(label,"world-maps.pdf"), p6 , width=10, height=7, units="in",limitsize = FALSE)
ggsave(paste0(label,"world-maps.png"), p6 , width=10, height=7, units="in",dpi = 1200, limitsize = FALSE)

################################################Sup Fig9e

shapes = c(rufipogon=24, nivala=21)
colors <- c(rufipogon="#969bc7",JAP ="#e9a00e", IND="#3c7c3d",nivala ="#2080E8")

p5.2<- p0 +
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice  %>% filter(subgroup == "rufipogon" | subgroup == "nivala" ) , 
aes(x=lon, y = lat,shape= subgroup ,fill = subgroup) ,
size= 1.5, colour='black',position = position_jitter(width = 0.5, height = 0.5) , stroke=.01)+
scale_fill_manual(values = colors) +
scale_shape_manual(values = shapes ) +
labs(x = "Longitude", y = "Latitude") 

p6.2 <- p5.2 +
coord_sf( ylim = c( 0,35),xlim = c( 70,120), crs = 4326)
label = "Sup Fig9e."

ggsave(paste0(label,"-maps.pdf"), p6.2 , width=5, height=7, units="in",limitsize = FALSE)
ggsave(paste0(label,"-maps.png"), p6.2 , width=5, height=7, units="in",dpi = 1200, limitsize = FALSE)


############################traits analysis
############################traits analysis
#############################Ext Fig4a

data1=read.csv("R4_rdp2_trial_effect.csv", header=T,  comment.char="", stringsAsFactors=F)
COL = c("HapA" ,"HapD", "Pvalue","Discribe")
label="Ext Fig4a."
width=5
height=4
topn=6

volcano_plot(data1, cut_off_FDR = 0.05, cut_off_log2FC = 0, COL = COL, 
width=width, height=height,label=label,topn=topn)

#############################Ext Fig4b

data=read.csv("R4_trial_solar_corre.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
ana_col = c("describe" ,"average_solar_r", "average_solar_p","significance")
label="Ext Fig4b"

topn=6
colors=c(Neg="#FF690B",Pos="#4F95FF",no="#FDFDFD")

library(ggplot2)
library(dplyr)

bind_result<-as.data.frame(data[, ana_col], row.names= row.names(data))
colnames(bind_result)<-c("tax","Estimate","p_adj","significance")
bind_result$tax <- factor(bind_result$tax, levels = unique(bind_result$tax))
bind_result$y <- reorder(bind_result$tax, -bind_result$Estimate)
bind_result$log10FDR <- -log10(bind_result$p_adj)

theme1<-theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),units=,"cm"),
              strip.text = element_text(size=8),text = element_text(size = 8),
              axis.line = element_line(color = "black",linewidth = 0.25),
              axis.ticks =  element_line(color = "black",linewidth = 0.25),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.text.y = element_text(color="black",size=8),
              axis.text.x = element_text(color="black",size=8,angle = 0),
              panel.spacing.x = unit(0,"cm"),
              panel.border = element_rect(color = "black", linewidth = 0.25, fill = NA),
              legend.position = "right",
      panel.spacing = unit(0,"lines"),plot.title = element_text(hjust = 0 , vjust = - 25 ,size=8),
      legend.background = element_blank(),legend.key = element_blank(),legend.title=element_blank() )

p1<- ggplot(bind_result, aes(x=Estimate, y=y, colour =significance ,fill= significance,size = log10FDR)) + 
geom_point(alpha=0.8, shape = 21, stroke = 0.1)+
scale_colour_manual(values = colors)+
scale_fill_manual(values = colors)+
ylab("")+xlab("Corrlation of solar (R)")+
theme1


width=6
height=4
ggsave(paste(label,"_","forest_plot",".pdf"), p1, width=width, height=height, units="in")
ggsave(paste(label,"_","forest_plot",".jpg"),p1, width=width, height=height, units="in")
p1

#############################trait world-wide
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggspatial)


worlds_rice = read.csv("R4_all_geo_climeate2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
anno_sf <- st_as_sf(worlds_rice,coords = c("Longitude", "Latitude"),crs = 4326)

mean_df = read.csv("R4_world_average_solar_5m.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

theme3<-
theme(axis.text.y = element_text(color="black",size=10),
  axis.text.x = element_text(color="black",size=10,angle = 0),
  panel.background = element_rect(fill = "transparent", color = NA),
  panel.grid.major = element_line(color = "#000000",linewidth =0.1),
  panel.grid.minor = element_line(color = "#000000",linewidth =0.1))

p0<-
ggplot() +
  geom_raster(data = mean_df, aes(x = x, y = y, fill = solar_rad)) +
#  scale_fill_viridis_c(name = "kJ m-2 day-1",option = "D" )+
  scale_fill_gradientn(name = "kJ m⁻² day⁻¹",
  colours = c("#000000",  "#6baed6", "#fcbf49", "#e7298a") )+
  coord_equal() +
  theme3+
  labs( x = "Longitude", y = "Latitude")

##################################################

shapes = c(HapA=21, HapD=24)
worlds_rice2 <- worlds_rice  %>% filter(  hdg_80head != "N"  )
worlds_rice2$hdg_80head <- as.numeric(worlds_rice2$hdg_80head)

p3<- p0+
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice2, aes(x=Longitude, y = Latitude, fill = hdg_80head ,shape= HAP ) ,
size= 1.0,colour='black',position = position_jitter(width = 3, height = 3) , stroke=.01)+
  scale_fill_gradientn(name = "days",
  colours = c("#FBFBFB","#0022FF") ) +
scale_shape_manual(values = shapes ) +
labs( x = "Longitude", y = "Latitude") 

########################################Ext Fig4c
p6 <- p3 +
coord_sf( ylim = c( -55,55),xlim = c( -150,150), crs = 4326)

ggsave(paste0("Ext Fig4c.","rice_Solar_heading_data.","-maps.pdf"), p6 , width=10, height=7, units="in",limitsize = FALSE)
ggsave(paste0("Ext Fig4c.","rice_Solar_heading_data.","-maps.png"), p6 , width=10, height=7, units="in",dpi = 1200, limitsize = FALSE)

#########################################
worlds_rice2 <- worlds_rice  %>% filter( liglt != "N"  ) 
worlds_rice2$liglt <- as.numeric(worlds_rice2$liglt)

p32<- p0+
ggnewscale::new_scale_fill()+
geom_point(data = worlds_rice2, aes(x=Longitude, y = Latitude, fill = liglt ,shape= HAP ) ,
size= 1.0,colour='black',position = position_jitter(width = 3, height = 3) , stroke=.01)+
  scale_fill_gradientn(name = "Ligule length (mm)",
  colours = c( "#FBFBFB", "#2A0386") ) +
scale_shape_manual(values = shapes ) +
labs( x = "Longitude", y = "Latitude") 

########################################Ext Fig4d
p62 <- p32 +
coord_sf( ylim = c( -55,55),xlim = c( -150,150), crs = 4326)

ggsave(paste0("Ext Fig4d.","rice_Solar_Ligule_length.","-maps.pdf"), p62 , width=10, height=7, units="in",limitsize = FALSE)
ggsave(paste0("Ext Fig4d.","rice_Solar_Ligule_length.","-maps.png"), p62 , width=10, height=7, units="in",dpi = 1200, limitsize = FALSE)

###############################################Fig4f

library(ggplot2)
library(gggenes)
library(tidyverse)

colors = c("#DB9E31","#4e3910", "#845d29", "#ae8548", "#e3c28b", "#4fb6ca", "#178f92", "#175f5d", "#054544")
gene_struc <- read.csv("R4_chromsome2.csv",head = T)

p23 <- ggplot(data = gene_struc %>% filter(sturcture== "Chrome") , aes(y = chrm)) +
  geom_gene_arrow(aes(xmin = start, xmax = end,forward = dirict1,  fill = SYMBOL) , fill = "white",
  arrowhead_width = unit(0, "mm"),arrowhead_height=unit(0, "mm") , arrow_body_height=unit(4, "mm") ) +
  geom_subgene_arrow(data = gene_struc %>% filter(sturcture== "500k"),
  aes(xmin = start, xmax = end, xsubmin = start, xsubmax = end,  fill = SYMBOL),
  arrowhead_width = unit(0, "mm"),arrowhead_height=unit(0, "mm") , arrow_body_height=unit(4, "mm")) +
  geom_feature( data = gene_struc %>% filter(sturcture== "SNP"),
  aes(x =  start, y = chrm), size =0.5, feature_height=unit(15, "mm"))+
  geom_feature_label(data = gene_struc %>% filter(sturcture== "SNP"),
  aes(x =  start, y = chrm, label = SYMBOL), 
  label_height=unit(10, "mm"),feature_height=unit(15, "mm"),size=8, angle = 45)+
  scale_fill_manual(values= colors ) +
  theme_genes()

ggsave(paste0("Fig4f.","Chromosome2_locus_id.pdf"), p23, width=20, height=8, units="in")
ggsave(paste0("Fig4f.","Chromosome2_locus_id.jpg"), p23, width=20, height=8, units="in") 

###############################################Fig4g

library(ggplot2)
library(gggenes)
library(tidyverse)

label="Fig4g"
colors = c(ld0.8= "#4fb6ca", other ="#DDDDDD", "#175f5d", "#054544")
snp_ld <- read.csv("R4_snp_all_ld.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

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

  p1 <- ggplot(snp_ld, aes(x = position, y = ld_r2, colour = lable,fill = lable)) +
    geom_point( shape= 21,alpha=0.8,stroke =0.1 ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "", y = "LD (R2)") +
    geom_vline(xintercept = c(25782447, 26018974), lty = 4, col = "black", lwd = 0.1) +
    geom_hline(yintercept = 0.8, lty = 4, col = "black", lwd = 0.1) +
    theme1 +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank())

  ggsave(paste0(label,".ld.","snp.pdf"), p1, width=10, height=1, units="in")
  ggsave(paste0(label,".ld.","snp.png"), p1, width=10, height=1,units="in")

#########################################Fig4h

library(ggplot2)
library(gggenes)
library(tidyverse)

colors = c("#C77CFF","#e3c28b","#175f5d","#4e3910", "#845d29", "#ae8548", "#178f92", "#4fb6ca",   "#054544")
gene_struc <- read.csv("R4_chromsome2_ld0.8.csv",head = T)

p24 <- ggplot(data = gene_struc , aes(y = chrm)) +
  geom_gene_arrow(aes(xmin = start, xmax = end,forward = dirict1,  fill = SYMBOL2) ,
  arrowhead_width = unit(4, "mm"),arrowhead_height=unit(6, "mm") , arrow_body_height=unit(6, "mm") )+
   scale_fill_manual(values= colors ) +
  theme_genes()

ggsave(paste0("Fig4h.","chromsome2_ld0.8_locus_id.pdf"), p24, width=20, height=8, units="in")
ggsave(paste0("Fig4h.","Cchromsome2_ld0.8_locus_id.jpg"), p24, width=20, height=8, units="in") 

########################################Fig4i

library(ggplot2)
library(gggenes)
library(tidyverse)

colors = c(OsCRY2="#e3c28b", OsPAL1="#4fb6ca","#4e3910", "#845d29", "#ae8548", "#178f92", "#4fb6ca")
gene_struc <- read.csv("R4_cry2_pal1_struc.csv",head = T)

p23 <- ggplot(data = gene_struc %>% filter(sturcture== "mRNA") , aes(y = chrm)) +
  geom_gene_arrow(aes(xmin = start, xmax = end,forward = dirict1,  fill = SYMBOL), , fill = "white",
  arrowhead_width = unit(4, "mm"),arrowhead_height=unit(6, "mm") , arrow_body_height=unit(6, "mm") ) +
  geom_subgene_arrow(data = gene_struc %>% filter(sturcture== "CDS"),
  aes(xmin = start, xmax = end, xsubmin = start, xsubmax = end,  fill = SYMBOL),
  arrowhead_width = unit(0, "mm"),arrowhead_height=unit(0, "mm") , arrow_body_height=unit(6, "mm")) +
  geom_feature( data = gene_struc %>% filter(sturcture== "SNP"),aes(x = start, y = chrm), size =0.5, feature_height=unit(15, "mm"))+
  geom_feature_label(data = gene_struc %>% filter(sturcture== "SNP"),aes(x = start, y = chrm, label = SYMBOL), 
  label_height=unit(10, "mm"),feature_height=unit(15, "mm"),size=8, angle = 45)+
   scale_fill_manual(values= colors ) +
  theme_genes() + 
  facet_wrap(~ locus_id, scales = "free", nrow = 1, shrink = FALSE) 

ggsave(paste0("Fig4i.","cyr2_pla1.pdf"), p23, width=20, height=8, units="in")
ggsave(paste0("Fig4i.","cyr2_pla1.jpg"), p23, width=20, height=8, units="in") 

#########################################################Ext Fig5a

data1=read.csv("R4_cry2_cds_snp_effect_trails.csv", header=T,  comment.char="", stringsAsFactors=F)

COL = c("v1", "v2" ,"Pvalue","ld","snp_name","group")

color =  c(
cult_repro = "#8487E2",
liglt = "#8487E2",
hdg_80head =  "#1FE948",
mat = "#1FE948",
grwd = "#C4EA1A",
grwt100 = "#C4EA1A")

shape = c(cult_repro = 21,
liglt = 22,
hdg_80head =  21,
mat = 22,
grwd = 21,
grwt100 = 22)
xlab=  "LD (r2)"
ylab= "Efect of trails (-log10(p value))"
label="Ext Fig5a"
width=4
height=4
topn=18

paopao_plot3 (data1,  COL =COL, xlab= xlab, ylab=ylab, width=width, height=height, color, shape, label=label,topn=topn)

#############################################################Ext Fig5b

library(dplyr)
library(ggplot2)
label = "Ext Fig5b"

FC=read.csv("R4_FC_hap.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
FC$log2fc <- log2(FC$v1/FC$v2)
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
               panel.spacing = unit(0,"lines") )

colorsF=c(HapA_cryA.vs.HapD_cryB= "#3CE753" , HapA_cryB.vs.HapD_cryA= "#8EB9ED")
FC2 <- FC %>% filter( group == "trail" ) 

write.csv (FC2, "Ext Fig5b_plot.csv")
widths = 6
heights =4

p1 <-ggplot(FC2) +
 geom_col(aes(x=log2fc ,y=factor,fill=FC),
 width = 0.8,colour ="black",linewidth=0.1	) + 
    xlab("") + ylab("") +
    geom_vline(xintercept = 0, lty = 1, col = "black", lwd = 0.1) +
    scale_fill_manual(values = colorsF)+
     ylab(ylabs)+xlab("")+
     theme
p1
  
ggsave(paste0(label,".pdf"), p1, width=widths, height=heights, units="in")
ggsave(paste0(label,".png"), p1, width=widths, height=heights, units="in")

######################################################Ext Fig5c

FC2 <- FC %>% filter( group == "micro" ) 
label = "Ext Fig5c"
widths = 3
heights =4

write.csv (FC2, "Ext Ext Fig5c_plot.csv")
p1 <-ggplot(FC2) +
 geom_col(aes(y=log2fc ,x=factor,fill=FC),
 width = 0.8,colour ="black",linewidth=0.1	) + 
    xlab("") + ylab("") +
    geom_hline(yintercept = 0, lty = 1, col = "black", lwd = 0.1) +
    scale_fill_manual(values = colorsF)+
     ylab(ylabs)+xlab("")+
     theme
p1
  
ggsave(paste0(label,".pdf"), p1, width=widths, height=heights, units="in")
ggsave(paste0(label,".png"), p1, width=widths, height=heights, units="in")



