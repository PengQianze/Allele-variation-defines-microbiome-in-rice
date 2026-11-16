
##############################################16s
##############################################

DF_ANA <- function (taxonomy,metadata,analysisgroup = "TimeGrowth",repea=6,label="ana_"
){
  library(dplyr )
  library(rstatix )
  library(tidyr)
  
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  colnames(sampFile)[1] = "Group"
  taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(sampFile)]
  
  uniqlegth <- function (x){
    y = length(unique(x))
    return (y)}
  
  uniq.num.df  <- data.frame(num = apply(taxonomy, 1, uniqlegth))
  uniq.num.df.2 <- filter(uniq.num.df, num > repea*3)
  taxonomy2 <- taxonomy[rownames(uniq.num.df.2), ]
  dat1 <- as.data.frame(merge(t(taxonomy2),sampFile,by='row.names'))
  
  df1 <-dat1[,-1]
  dfnorder = as.data.frame(reshape2::melt(df1, id.vars = c("Group")))

  sumorder<-plyr::ddply(dfnorder,c("Group","variable"),summarise,
                        allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
                        n1=sum(!is.na(value)),se1=sd1/sqrt(n1))
  
  max_mean_group<- sumorder %>% group_by(variable)  %>% slice_max(allmean, n =1) %>%
   select(variable, enrichment = Group, max_mean = allmean)
  
  sumorder2<-
    full_join(sumorder, max_mean_group, by = 'variable',relationship = "many-to-many")
  write.csv(as.data.frame(sumorder2),paste0(label,".",analysisgroup,"-SUM.csv"))
  
  sumorder_wide<- sumorder %>%  select(variable, Group, allmean)  %>%
   pivot_wider(names_from = Group, values_from = c(allmean))
  write.csv(as.data.frame(sumorder_wide),paste0(label,".",analysisgroup,"-wider-means.csv"))
  
  sumorder_wide<- sumorder %>%  select(variable, Group, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
 write.csv(as.data.frame(sumorder_wide),paste0(label,".",analysisgroup,"-wider-se.csv"))
  
  
  orderstat.test<- dfnorder %>% group_by(variable)%>% t_test(value ~ Group)%>%adjust_pvalue(method = "fdr") %>%
    add_significance(p.col = "p") # Benjamini & Yekutieli 
  write.csv(as.data.frame(orderstat.test),paste0(label,".",analysisgroup,"-df-t-test.csv"))
  
  return(sumorder)
}

###########################
###########################genus_stack

STACK_ANA2 <- function (taxnonmy,metadata,analysisgroup, label = "ana_", level = level,width=4, height=3.5,
                       colorsF = c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                                   "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                                   "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39") ,colors = colors,
                       topN = 10 
){
  library(reshape2)
  library(ggplot2)
  library(tibble)
  library( rstatix)
  idx = rownames(metadata) %in% colnames(taxonomy)
  metadata = metadata[idx, , drop = F]
  taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  colnames(sampFile)[1] = "Group"  
  mean_sort = as.data.frame(taxonomy[(order(-rowSums(taxonomy))), ])
  idx = grepl("unassigned|unclassified|unknown", rownames(mean_sort), ignore.case = T)
  mean_sort = rbind(mean_sort[!idx, ], mean_sort[idx, ])
  other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
  mean_sort = mean_sort[1:(topN - 1), ]
  mean_sort = rbind(mean_sort, other)
  rownames(mean_sort)[topN] = c("Other")
  colnames(sampFile)[1] = "Group"
  amplic2 = t(mean_sort)
  
  mat_t2 = merge(sampFile, amplic2, by = "row.names")
  mat_t2 = mat_t2[, c(-1)]
  mat_mean = aggregate(mat_t2[, -1], by = mat_t2[1], FUN = mean)
  mat_mean_final = do.call(rbind, mat_mean)[-1, ]
  colnames(mat_mean_final) = mat_mean$Group
  mean_sort2 = as.data.frame(mat_mean_final)
  mean_sort2$Taxonomy = rownames(mean_sort2)
  data_all = as.data.frame(melt(mean_sort2, id.vars = c("Taxonomy")))
  data_all$Taxonomy = factor(data_all$Taxonomy, levels = rownames(mean_sort2))
  data_all$value = as.numeric(data_all$value)
  
  write.csv(mat_mean_final,paste0(label,".",analysisgroup,"stack.csv"))

  data_all <- data_all %>%mutate(variable=factor(variable,levels =level))
  
  p7 = ggplot(data_all, aes(x = variable, y = value, fill = Taxonomy)) + 
   geom_bar(stat = "identity", position = "stack", width = 0.8,colour ="white",linewidth=0.1	) + 
    xlab("") + ylab("Relative abundance (%)") +
    scale_fill_manual(values = colorsF)+
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
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
          legend.text = element_text(face = "italic")
    )
  

  ggsave(paste0(label,".",analysisgroup,"stack.pdf"), p7, width=width, height=height, units="in")
  ggsave(paste0(label,".",analysisgroup,"stack.png"), p7, width= width, height=height, units="in")


library(dplyr )
library(rstatix )
library(tidyr)
dat1 <- as.data.frame(merge(amplic2,sampFile,by='row.names'))
row_name<-dat1[,1]
df1 <-dat1[,-1]
rownames(df1)<-row_name
df <- df1 %>%mutate(Group=factor(Group,levels =level))
xlong<-length(level)*1.5
ylong<-length(level)*2


for (col_name in colnames(df)) {
  if (col_name == "Group"){
    next
  }
filename <- paste0(label,".",analysisgroup,col_name,".pdf")
dfn <- df[c({{col_name}}, "Group")]
names(dfn) <- c("Relative_aboundance", "Group")
stat.test<- dfn %>% t_test(Relative_aboundance ~ Group)%>%
add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
readr::write_csv(as.data.frame(dfn), paste0(label,".",analysisgroup,".",col_name,".plot.csv"))

p <-ggplot(dfn,aes(Group,Relative_aboundance)) +
geom_boxplot(aes(fill=Group,alpha = 0.8),alpha=1, linewidth=0.25, width = 0.8,outlier.shape = NA, outlier.size = 0)+
geom_point(aes(fill=Group,alpha=0.7),colour ="black",alpha=0.7,size =2.5, shape = 21, stroke = 0.1,position = position_jitter(0.2))+    
scale_fill_manual(values = colors)+
ggpubr::stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F,tip.length = 0, label.size =3)+
 xlab("") + ylab("Relative abundance (%)") +
scale_size_continuous(range=c(1,3)) +  
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=8),
        axis.line = element_line(color = "black",linewidth = 0.25),
        axis.ticks =  element_line(color = "black",linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        axis.title.y  = element_text(size = 8),
        axis.text.y = element_text(color="black",size=8),
        axis.text.x = element_text(color="black",size=8,angle=0),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
        legend.position = "none",
       panel.spacing = unit(0,"lines")
       )+
      labs(title =parse(text = paste0("italic('", col_name, "')")))+
  coord_cartesian(clip = "off")

ggsave(paste0(filename), p, width=xlong, height=ylong, units="in")
ggsave(paste0(label,".",analysisgroup,".",col_name,".png"), p, width=xlong, height=ylong, units="in")
}
}

##########################################

PCOA_ANA <- function (taxonomy,metadata,analysisgroup = "TimeGrowth",level,width=8, height=8,
                      colors= c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                   "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                              "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39"),
                      shapes = c(21,22,23,24,25),label="ana_"
){
  
  library( vegan  )
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  idx = rownames(metadata) %in% colnames(taxonomy)
  metadata = metadata[idx, , drop = F]
  taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  
  data <- vegdist(t(taxonomy),method = "bray")
  write.table(as.matrix(data), paste0(label,"_",analysisgroup,"bray.txt"), sep = '\t', col.names = NA, quote = FALSE)#输出距离矩阵
  otu_matrix<-read.table(paste0(label,"_",analysisgroup,"bray.txt"), header=T, row.names=1, sep="\t", comment.char="")
  colnames(sampFile)[1] = "group"
  
  #非限制性排序及组间显著性分析
  #sampFile = as.data.frame(metadata[, "group"], row.names = row.names(metadata))
  pcoa = cmdscale(otu_matrix, k = 3, eig = T)
  points = as.data.frame(pcoa$points)
  eig = pcoa$eig
  points = cbind(points, sampFile)
  colnames(points) = c("x", "y", "z", "group")
  points$group = factor(points$group, levels = level)

  #显著性
  adonis_table = adonis2(otu_matrix ~ group, data = sampFile, permutations = 10000)
  adonis_pvalue = adonis_table$`Pr(>F)`[1]
  adonis_R2=adonis_table$`R2`[1]
  adonis_R2pvalue = paste("group","R2:", adonis_R2, "pvalue:", adonis_pvalue, sep = "\t")
  write.table(adonis_R2pvalue, file = paste0(label,".",analysisgroup,"-adonis2test.txt"), append = TRUE, 
              sep = "\t", quote = F, row.names = F, col.names = F)
  
  pcoa.fit <- envfit(pcoa ~ group, data=sampFile, perm=999)
  envfitR2 = paste("group R2:", pcoa.fit$factors$r, "Pvalue:", pcoa.fit$factors$pvals,sep = "\t")
  write.table(envfitR2,paste0(label,".",analysisgroup,"-envfit.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  
  points <- points %>% tibble::rownames_to_column(var = "lable")

  write.csv(points,paste0(label,".",analysisgroup,"pcoa_plot.csv"))
  #可视化
  p2 = ggplot(points, aes(x = x, y = y, color = group,shape = group,fill = group)) + 
    labs(x = paste("PCo 1 (", format(100 * eig[1]/sum(eig), digits = 4), "%)", sep = ""), 
         y = paste("PCo 2 (", format(100 * eig[2]/sum(eig), digits = 4), "%)", sep = ""), color = "group")+
    scale_fill_manual(values = colors)+
    scale_colour_manual(values = colors)+
    scale_shape_manual(values = shapes)+
    geom_point(alpha = 0.8, size = 2.5,color = "black",stroke = 0.1) + 
    stat_ellipse(aes(color = group),level = 0.68,linewidth=0.25) +
    ggtitle(paste( " P = ", pcoa.fit$factors$pvals,  " R2 = ", round(pcoa.fit$factors$r,digits = 4), sep = ""))+

    theme(
      strip.text = element_text(size=8),text = element_text(size = 8),
      axis.line = element_line(color = "black",linewidth = 0.25),
      axis.ticks= element_line(colour = "black",linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.y = element_text(color="black",size=8),
      axis.text.x = element_text(color="black",size=8,angle = 0),
      panel.spacing.x = unit(0,"cm"),
      panel.border = element_blank(),
      panel.spacing = unit(0,"lines"),legend.position = "top",plot.title = element_text(hjust = 0 , vjust = - 25 ,size=8),
      legend.background = element_blank(),legend.key = element_blank(),legend.title=element_blank()
    )
  
  ggsave(paste0(label,".",analysisgroup,"-PCOA12.pdf"), p2, width=width, height=height, units="in",limitsize = FALSE)
  ggsave(paste0(label,".",analysisgroup,"-PCOA12.png"), p2, width=width, height=height, units="in",limitsize = FALSE)
  
  return(p2)
}

##########################################

#############################

box_ANA <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,
                     colors = c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                                "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96"),level,
                                shapes=c(21,22,23,24),label="ana_"
){
  library(rstatix)
  library(ggplot2)
  library( ggpubr )

  dfn <- metadata[c(col_name, analysisgroup)]
  names(dfn) <- c("Content", "Group")
  dfn$Content<-as.numeric(dfn$Content)
  dfn <- dfn%>%mutate(Group=factor(Group,levels =level))
  stat.test<- dfn %>% t_test(Content ~ Group)%>%
  add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
  readr::write_csv(as.data.frame(dfn),paste0(label,"_",analysisgroup,col_name,".plot.csv"))

  sumorder<-dfn %>%group_by(Group) %>%summarise(
    allmean = mean(Content, na.rm = TRUE),       
    sd1     = sd(Content, na.rm = TRUE),         
    factor     = analysisgroup,                      
    n1      = sum(!is.na(Content)),               
    se1     = sd1 / sqrt(n1),                   
    .groups = "drop")
  write.csv(as.data.frame(sumorder),paste0(label,"_",analysisgroup,col_name,"-sumorder.csv"))

  theme<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
               strip.text = element_text(size=6),text = element_text(size = 6),
               axis.line = element_line(color = "black",linewidth = 0.25),
               axis.ticks =  element_line(color = "black",linewidth = 0.25),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(color="black",size=6),
               axis.text.x = element_text(color="black",size=6,angle = 0),
               panel.spacing.x = unit(0,"cm"),
               panel.border = element_blank(),
               panel.spacing = unit(0,"lines"),
               legend.position = "none"
  )
    
p1 <-ggplot(dfn,aes(Group,Content)) +
     geom_boxplot(aes(fill=Group,alpha = 0.8),alpha=1, linewidth=0.25, width = 0.8,outlier.shape = NA, outlier.size = 0)+
    geom_point(aes(fill=Group),colour ="black",alpha=0.7,size =2.5, shape = 21, stroke = 0.1,position = position_jitter(0.2))+    
     scale_fill_manual(values = colors)+
     stat_pvalue_manual(stat.test,label = "p = {p.adj}",hide.ns = F,tip.length = 0, label.size =2.25)+
     scale_size_continuous(range=c(1,3)) +  
     ylab(ylabs)+xlab("")+
     theme
p1
  
ggsave(paste0(label,"_",analysisgroup,col_name,".pdf"), p1, width=widths, height=heights, units="in")
ggsave(paste0(label,"_",analysisgroup,col_name,".png"), p1, width=widths, height=heights, units="in")

return(p1)
} 


#############################################

box_ANA3 <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,level=c(1,2),
                     colors = c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
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
  dfn <- dfn%>%mutate(Group=factor(Group,levels =level))
  stat.test<- dfn %>% t_test(variable ~ Group) %>% add_significance(p.col = "p") 
  readr::write_csv(as.data.frame(dfn),paste0(label,"_",analysisgroup,col_name,".plot.csv"))

  test.p = paste(analysisgroup ,"G1:", stat.test$group1,"G2:",stat.test$group2,  "Pvalue:", stat.test$p,"N1:",stat.test$n1,"N2:",stat.test$n2,ssep = ",")
  write.table(test.p,paste0(col_name,"ALL","-t_test.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  

  sumorder<-dfn %>%group_by(Group) %>%summarise(
    allmean = mean(variable, na.rm = TRUE),       
    sd1     = sd(variable, na.rm = TRUE),         
    SNP     = analysisgroup,                      
    n1      = sum(!is.na(variable)),               
    se1     = sd1 / sqrt(n1),                   
    .groups = "drop")
  write.csv(as.data.frame(sumorder),paste0(analysisgroup,".","-sumorder.csv"))

  sumorder_wide <- sumorder %>%  select( Group,SNP,allmean)  %>%
  pivot_wider(names_from = Group, values_from = c(allmean))
  names(sumorder_wide) <- c("SNP", "G1", "G2")
  write.csv(as.data.frame(sumorder_wide),paste0(analysisgroup,".","-wider-means.csv"))

  means = paste("SNP:",sumorder_wide$SNP ,"G1:", sumorder_wide$G1,"G2:",sumorder_wide$G2,ssep = ",")
  write.table(means, paste0(col_name,"ALL","-mean.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

  sumorder_wide_SE<- sumorder %>%  select( Group,SNP, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
  names(sumorder_wide_SE) <- c("SNP", "G1", "G2")
  write.csv(as.data.frame(sumorder_wide_SE),paste0(analysisgroup,".","-wider-se.csv"))

  SE = paste("SNP:",sumorder_wide_SE$SNP ,"G1:", sumorder_wide_SE$G1,"G2:",sumorder_wide_SE$G2,ssep = ",")
  write.table(SE, paste0(col_name,"ALL","-se.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  
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
    geom_point(aes(fill=Group),colour ="black",alpha=0.7,size =2.5, shape = 21, stroke = 0.1,position = position_jitter(0.2))+    
     scale_fill_manual(values = colors)+
    ggtitle(paste("p = ", stat.test$p))+
     scale_size_continuous(range=c(1,3)) +  
     ylab(ylabs)+xlab("")+
     theme
p1
  
ggsave(paste0(label,"_",analysisgroup,col_name,".pdf"), p1, width=widths, height=heights, units="in")
ggsave(paste0(label,"_",analysisgroup,col_name,".png"), p1, width=widths, height=heights, units="in")

return(p1)
} 

###################################################data
################################################

library(dplyr)
taxonomy=read.csv("R2R3_Genus RA_16s.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("R2R3_metadata_16s.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0
colorscat <-c(Pantoea="#5C8DC7",Pseudomonas="#EF7F29",Curtobacterium="#9A6731",Sphingomonas="#3EB1DB",Acinetobacter="#9B990A",
Escherichia="#6ECCDD",Methylobacterium="#298C8C",Deinococcus="#6FD8AB",Enterobacter="#AD70AE",Parahymenobacter="#C8A580",
Other="#B4B4B4","#4B86C6","#AFD8E6",Paenibacillus="#D6B28C",Exiguobacterium="#D0B937",Erwinia="#BD52BE",
"#AED8E5","#4B86C6","#FED700","#EC7722")

colors <- c(A="#3CE753", D="#8EB9ED")
level = c("A","D")
shapes<- c(A=21, D=21)
analysisgroup = "genotype"
##################################
metadata2 <- metadata %>% filter( varity == "NPB"  ) 
label = paste("Fig 2a")
width=3
height=3.5

STACK_ANA2(taxnonmy,metadata2,analysisgroup,label = label, level = level,width=width, height=height, colorsF = colorscat,colors = colors, topN = 10 )
DF_ANA (taxonomy =taxonomy,metadata2,analysisgroup,repea=0,label=label) 
#################################

taxonomy[is.na(taxonomy)]<-0
label = paste("Fig S6a")
width = 2.5
height = 4
PCOA_ANA (taxonomy,metadata2,analysisgroup ,level,width=width, height=height,colors= colors,shapes = shapes,label=label)

################################
metadata2 <- metadata %>% filter( varity == "ZH11"  ) 
label = paste("Fig 2b")
width=3
height=3.5
STACK_ANA2(taxnonmy,metadata2,analysisgroup,label = label, level = level,width=width, height=height, colorsF = colorscat,colors = colors, topN = 10 )
DF_ANA (taxonomy =taxonomy,metadata2,analysisgroup,repea=0,label=label)  

#########################

label = paste("Sup Fig6b")
width = 2.5
height = 4
PCOA_ANA (taxonomy,metadata2,analysisgroup ,level,width=width, height=height,colors= colors,shapes = shapes,label=label)

#################
metadata2 <- metadata %>% filter( varity == "9311"  ) 
label = paste("Fig 2c")
width=3
height=3.5
STACK_ANA2(taxnonmy,metadata2,analysisgroup,label = label, level = level,width=width, height=height, colorsF = colorscat,colors = colors, topN = 10 )
DF_ANA (taxonomy =taxonomy,metadata2,analysisgroup,repea=0,label=label)  

##########################

label = paste("Sup Fig6c")
width = 2.5
height = 4
PCOA_ANA (taxonomy,metadata2,analysisgroup ,level,width=width, height=height,colors= colors,shapes = shapes,label=label)



#################Fungal_biomass

library(dplyr)
metadata=read.csv("R2R3_MO_biomass.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata[is.na(metadata)]<-"N"

colors <- c(
A="#3CE753",
D="#8EB9ED",
SynCom_A="#3CE753",
SynCom_D="#8EB9ED")
level = c("A","D","SynCom_A","SynCom_D")
widths = 2.6
heights = 4
analysisgroup = "group"
col_name = "Fungal_biomass"
ylabs = "Fungal biomass"

###########################Fig2e
label<-paste("Fig2e","_",col_name )
new_alpha2 <- metadata %>% filter( varieties == "NPB"  ) 

box_ANA (new_alpha2 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

###########################Fig2i
label<-paste("Fig2e","_",col_name )
new_alpha2 <- metadata %>% filter( varieties == "ZH11"  ) 

box_ANA (new_alpha2 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

###############################Fig2m
new_alpha2 <- metadata %>% filter( varieties == "9311"  ) 
label<-paste("Fig2m","_",col_name )

box_ANA (new_alpha2 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

##############################Relative lesion area

library(dplyr)
metadata=read.csv("R2R3_RL.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata[is.na(metadata)]<-"N"

colors <- c(
A="#3CE753",
D="#8EB9ED",
SynCom_A="#3CE753",
SynCom_D="#8EB9ED")
level = c("A","D","SynCom_A","SynCom_D")
widths = 2.6
heights = 4

analysisgroup = "group"
col_name = "Relative_lesion_area"
ylabs = "Relative lesion area (%)"


###############Fig2g
new_alpha2 <- metadata %>% filter( varieties == "NPB"  ) 
label<-paste("Fig2g","_",col_name )

box_ANA (new_alpha2 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

###############Fig2k
new_alpha2 <- metadata %>% filter( varieties == "ZH11"  ) 
label<-paste("Fig2k","_",col_name )

box_ANA (new_alpha2 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

###############Fig2o
new_alpha2 <- metadata %>% filter( varieties == "9311"  ) 
label<-paste("Fig2o","_",col_name )
box_ANA (new_alpha2 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

####################HCA

library(dplyr)
metadata=read.csv("R2R3_HCA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata[is.na(metadata)]<-"N"
analysisgroup = "group"
col_name = "HCA"
colors <- c(A="#3CE753", D="#8EB9ED")
level = c("A","D")
ylabs = "Concentration (mg/g)"
widths = 1.5
heights =3

###############################Fig2p
new_alpha1 <- metadata %>% filter( varieties == "NPB"  ) 
label<-paste("Fig2p","_",col_name )

box_ANA3 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

##############################Fig2q
new_alpha1 <- metadata %>% filter( varieties == "ZH11"  ) 
label<-paste("Fig2q","_",col_name )

box_ANA3 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

##############################Fig2r
new_alpha1 <- metadata %>% filter( varieties == "9311"  ) 
label<-paste("Fig2r","_",col_name )

box_ANA3 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

####################emzye

library(dplyr)
metadata=read.csv("R2R3_emzye.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata[is.na(metadata)]<-"N"


analysisgroup = "group"
col_name = "PAL_activity"
colors <- c(A="#3CE753", D="#8EB9ED")
level = c("A","D")
ylabs = "PAL activity (U/mg protein)"
widths = 1.5
heights =3

#####################################Fig2s
new_alpha1 <- metadata %>% filter( varieties == "NPB"  ) 
label<-paste("Fig2s","_",col_name )

box_ANA3 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

####################################Fig2t
new_alpha1 <- metadata %>% filter( varieties == "ZH11"  ) 
label<-paste("Fig2t","_",col_name )

box_ANA3 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

#####################################Fig2u
label<-paste("Fig2u","_",col_name )
new_alpha1 <- metadata %>% filter( varieties == "9311"  ) 
box_ANA3 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)

###################################Fig3i
new_alpha1 <- metadata %>% filter( varieties == "OsPAL1"  ) 
label<-paste("Fig3i","_",col_name )

widths = 1.5
heights =3
box_ANA3 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)






