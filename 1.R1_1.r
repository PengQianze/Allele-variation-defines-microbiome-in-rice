
################################################functions
################################################
################################################

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

STACK_ANA <- function (taxnonmy,metadata,analysisgroup, label = "ana_", level = level,width=4, height=3.5,
                       colorsF = c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                                   "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                                   "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39") ,
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
          panel.spacing = unit(0,"lines")
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
dfx = as.data.frame(melt(df,id=c("Group")))

sum<-plyr::ddply(dfx ,c("Group","variable"),summarise,
allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
n1=sum(!is.na(value)),se1=sd1/sqrt(n1))

sumorder_wide<- sum %>%  select(variable, Group, se1)  %>%
pivot_wider(names_from = Group, values_from = c(se1))

write.csv(as.data.frame(sumorder_wide),paste0(label,".",analysisgroup,"-wider-se.csv"))

orderstat.test<- dfx %>% group_by(variable)%>% t_test(value ~ Group)%>%adjust_pvalue(method = "fdr") %>%
add_significance(p.col = "p") # Benjamini & Yekutieli 
write.csv(as.data.frame(orderstat.test),paste0(label,".",analysisgroup,"-df-t-test.csv"))

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

#########################################
##########################################

Regression_sample1 <-function(data,ana_col,level,xlabs="x",ylabs="y", width=6.0, height=8.2,label="lab"
,colors=c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E")){

library(ggplot2)

sampFiletax2<-as.data.frame(data[, ana_col], row.names= row.names(data))
colnames(sampFiletax2)<-c("x","y","Group")
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
geom_point(aes(x,y,colour =Group),fill="transparent",alpha=0.8,size =1, shape = 21, stroke = 0.4)+

scale_colour_manual(values = colors)+
geom_smooth(method = "lm", formula = NULL,linewidth=0.2,se=FALSE,color="black",linetype="solid",aes(group=1))+
ggtitle(paste( " P = ",sprintf("%.4e",corr_p),   " R = ", sprintf("%.4f",corr_r), sep = ""))+
ylab(ylabs)+xlab(xlabs)+theme1

p21
ggsave(paste(label,"_",ana_col[1],"_",ana_col[2],"_","Line_pearson_R",".pdf"), p21, width=width, height=height, units="in")
ggsave(paste(label,"_",ana_col[1],"_",ana_col[2],"_","Line_pearson_R",".jpg"),p21, width=width, height=height, units="in")
}


########################################effect of microbiome

Effect_size_llm <-function(df,factors,random,label="llm_effect",width=6, height=6,cut_off_FC = 5, cut_off_FDR = 0.05,topn = 5){
bind_result <-data.frame()
for (col_name in colnames(df)) {
  if (col_name ==  {{factors}}){
    next
  }
    if (col_name ==  {{random}}){
    next
  }
    if (mean(df[[ {{col_name}} ]], na.rm = TRUE) < 0.0000001){
    next
  }
  
library(lme4)
library(lmerTest)
dfn2 <- df[c({{col_name}}, {{factors}},{{random}})]
names(dfn2) <- c("alpha", "factor","random")

fm1<-lmer(alpha~factor+(1|random),data=dfn2)
    presult<-car::Anova(fm1,type=2)
    coefs <- summary(fm1) $coefficients[2,1]  ## coefs
    SEs<-summary(fm1) $coefficients[2,2] ##standard errors
    t_vals<-summary(fm1) $coefficients[2,3]##t values
    p_vals <-summary(fm1) $coefficients[2,4]
    chisqP<-presult$`Pr(>Chisq)`
results <- data.frame(
  tax = col_name,
  factor = factors,
  Estimate = coefs,
  SE = SEs,
  t_value = t_vals,
  Pr_t = p_vals,
  ANOV_P = chisqP)
bind_result <-rbind(bind_result,results )
}
write.csv(bind_result, file = paste0(label,".","-LLMtest.csv"))


library(ggplot2)
library(ggrepel)

bind_result <- bind_result %>% filter(tax != "Other") %>%
mutate( p_adj = p.adjust(ANOV_P, method = "fdr"),
significance = case_when( p_adj < 0.05  ~ "sig", TRUE ~ "ns" ) )
bind_result$log10FDR <- log10(bind_result$ p_adj)

bind_result$Sig = ifelse(bind_result$p_adj < cut_off_FDR & abs(bind_result$Estimate) >= cut_off_FC, ifelse(bind_result$Estimate > cut_off_FC, 'Up', 'Down'), 'no')


write.csv(bind_result,paste0(label,".","LLM-test-plot.csv"))


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
  p1 <- ggplot(bind_result, aes(x = Estimate, y = -log10FDR)) +
    geom_point( aes(fill = Sig, size = abs(Estimate)), shape = 21,alpha = 0.8,colour = "black", stroke = 0.3) +
    scale_fill_manual(values = c(Up="#FF690B",Down="#4F95FF",no="#FFFFFF")) +
    geom_vline(xintercept = c(-cut_off_FC, cut_off_FC), lty = 4, col = "black", lwd = 0.25) +
    geom_hline(yintercept = -log10(cut_off_FDR), lty = 4, col = "black", lwd = 0.25) +
    labs(x = "Effect size", y = "-log10FDR") +
    theme1
 
  p1

  up <- bind_result[bind_result$Sig  == 'Up', ]
  top_5_up <- order(up$log10FDR, decreasing = FALSE)[1:topn]
  top_5_rows_up <- up[top_5_up, ]
  
  down <- bind_result[bind_result$Sig  == 'Down', ]
  top_5_down <- order(down$log10FDR, decreasing = FALSE)[1:topn]
  top_5_rows_down <- down[top_5_down, ]
  
  
  p2 <- p1 + geom_text_repel(
    data = top_5_rows_up,
    aes(label = tax), size = 1,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
  p2 <-p2 + geom_text_repel(
    data = top_5_rows_down,
    aes(label = tax), size = 1,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
  
  p2

  ggsave(paste0(label,"labled.","volcano.pdf"), p2, width=width, height=height, units="in")
  ggsave(paste0(label,"labled.","volcano.png"), p2, width=width, height=height, units="in")
  
}

#############lm

forest_plot <-function(data,ana_col,colors=c(Up="#FF690B",Down="#4F95FF"),label="order_top10_effect",width=5, height=8 ){

library(ggplot2)
library(dplyr)

bind_result<-as.data.frame(data[, ana_col], row.names= row.names(data))
colnames(bind_result)<-c("tax","Estimate","SE","p_adj","significance")
bind_result$tax <- factor(bind_result$tax, levels = unique(bind_result$tax))
bind_result$y <- reorder(bind_result$tax, -bind_result$Estimate)

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
      panel.spacing = unit(0,"lines"),legend.position = "top",plot.title = element_text(hjust = 0 , vjust = - 25 ,size=8),
      legend.background = element_blank(),legend.key = element_blank(),legend.title=element_blank() )

write.csv(bind_result,paste(label,".plot.csv") )
p1<- ggplot(bind_result, aes(x=Estimate, y=y,colour =significance)) + 
geom_errorbar(aes(xmin=(Estimate-SE), xmax=(Estimate+SE)), width =0.3,linewidth=0.8, position=position_dodge(0.05))+
geom_point(fill="white",alpha=1,size =2, shape = 21, stroke = 1)+
scale_colour_manual(values = colors)+
ylab("")+xlab("Effect size")+
geom_text(aes(x=max(Estimate)+max(Estimate)*0.15,label= sprintf("%.2e", p_adj)),size=3,
hjust = 0, vjust = 0)+
scale_x_continuous(expand = expansion(mult =c(0.05, 0.3)))+
theme1

ggsave(paste(label,"_","forest_plot",".pdf"), p1, width=width, height=height, units="in")
ggsave(paste(label,"_","forest_plot",".jpg"),p1, width=width, height=height, units="in")
p1

}



###############################

alpha_diversity <- function(x) {
  library(vegan)  
  observed_species <- specnumber(x, MARGIN = 1) #estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]#otu表需要为整数
  Shannon <- diversity(x, index = 'shannon',MARGIN = 1,base = 2)
  Simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  Pielou <- (Shannon) / log(observed_species)
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  Pielou <- sprintf("%0.4f", Pielou)
  result <- data.frame(observed_species, Chao1, Shannon, Simpson, Pielou)
  
  write.table(result, 'alpha_diversity.txt', quote = FALSE)
  return(result)
}

##########################################

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
  dfn <- dfn%>%mutate(Group=factor(Group,levels =level))
  stat.test<- dfn %>% t_test(Content ~ Group)%>%
  add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")
  readr::write_delim(as.data.frame(stat.test),paste0(label,"_",analysisgroup,col_name,"-df-t-test.txt"))

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
     stat_pvalue_manual(stat.test,label = "p = {p.adj}",hide.ns = F,tip.length = 0, label.size =3)+
     scale_size_continuous(range=c(1,3)) +  
     ylab(ylabs)+xlab("")+
     theme
p1
  
ggsave(paste0(label,"_",analysisgroup,col_name,".pdf"), p1, width=widths, height=heights, units="in")
ggsave(paste0(label,"_",analysisgroup,col_name,".png"), p1, width=widths, height=heights, units="in")

return(p1)
} 


##########################################
###############################
##############Sup Fig1
#################alpha_diversity

otu <- read.csv("R1_1_bacteria.S.A.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
otu[is.na(otu)]<-0
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

library(amplicon)

wild_samples <- rownames(metadata[metadata$group3 == "Wild", , drop = FALSE])
wild_otu_data <- otu[, colnames(otu) %in% wild_samples, drop = FALSE]

p1 = alpha_sample_rare(wild_otu_data, length=18, rep=30, count_cutoff=10)
ggsave(paste0("Sup Fig1a.pdf"), p1, width=6, height=4, units="in")
ggsave(paste0("Sup Fig1a.png"), p1, width=6, height=4, units="in")

JAP_samples <- rownames(metadata[metadata$group3 == "JAP", , drop = FALSE])
JAP_otu_data <- otu[, colnames(otu) %in% JAP_samples, drop = FALSE]

p2 = alpha_sample_rare(JAP_otu_data, length=18, rep=30,count_cutoff=10)
ggsave(paste0("Sup Fig1b.pdf"), p2, width=6, height=4, units="in")
ggsave(paste0("Sup Fig1b.png"), p2, width=6, height=4, units="in")


IND_samples <- rownames(metadata[metadata$group3 == "IND", , drop = FALSE])
IND_otu_data <- otu[, colnames(otu) %in% IND_samples, drop = FALSE]

p3 = alpha_sample_rare(IND_otu_data, length=18, rep=30, count_cutoff=10)
ggsave(paste0("Sup Fig1c.pdf"), p3, width=6, height=4, units="in")
ggsave(paste0("Sup Fig1c.png"), p3, width=6, height=4, units="in")


##################Extended Data Fig. 1
library(dplyr)
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
otu2 <- t(otu)
alpha <- alpha_diversity (otu2)

colors <- c(Wild="#8EB9ED",JAP ="#E2DB59", IND="#3CE753")
level = c("Wild","JAP","IND")
analysisgroup = "group3"
widths = 2.5
heights =4

metadata2 <-  metadata[,analysisgroup, drop = FALSE]
new_alpha = merge(alpha, metadata2, by= "row.names")
new_alpha <-tibble::column_to_rownames(new_alpha,var = "Row.names")  
new_alpha [is.na(new_alpha )] <- 0

write.csv(new_alpha,"Extended Data Fig. 1_new_alpha.csv")
new_alpha=read.csv("Extended Data Fig. 1_new_alpha.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

for (col_name in colnames(new_alpha)) {
  if (col_name == analysisgroup) {next}
  else {
col_name = col_name
ylabs = col_name
label<-paste("Extended Data Fig. 1","_",col_name )
box_ANA (new_alpha ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }


############Ext Fig1f

taxonomy=read.csv("R1_1_bacteria.O.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

taxonomy[is.na(taxonomy)]<-0
colors <- c(Wild="#8EB9ED",JAP ="#E2DB59", IND="#3CE753")
level = c("Wild","JAP","IND")
shapes<- c(Wild=21,JAP =22, IND=23)
label = "Ext Fig1f"
analysisgroup="group3"
width=5
height=5
PCOA_ANA (taxonomy,metadata,analysisgroup ,level,width=width, height=height,colors= colors,shapes = shapes,label=label)

################Fig1c
library(dplyr)
library(tidyverse)
library(tidyr)

genetic_distance<-read.csv("R1_1_genome_distance.csv", header=T, row.names=1, comment.char="")
otu_matrix<-read.table("Ext Fig1f_group3bray.txt", header=T, row.names=1, comment.char="")
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

df <- otu_matrix %>% as.data.frame() %>% tibble::rownames_to_column("row_name")
pair_data_tidy <- df %>%  pivot_longer( cols = -row_name, names_to = "col_name", values_to = "BCvalue" )     
clean_symmetric_df <- pair_data_tidy %>%mutate(new_1 = pmin(row_name, col_name), new_2 = pmax(row_name, col_name)) %>% 
distinct(new_1, new_2, .keep_all = TRUE) %>% select(-new_1, -new_2)

group_info <- metadata %>% select(tree_id, group3)  %>% tibble::rownames_to_column("row_name")
rownames(group_info) <- NULL  

group_info_row <-group_info %>%  rename(row_tree = 2, row_group = 3 )
group_info_col <-group_info %>%  rename(col_name = 1 ,col_tree = 2, col_group = 3 )

metadata4 <- merge( clean_symmetric_df,group_info_row, by = "row_name")
metadata5 <- merge( metadata4,group_info_col, by = "col_name")

result <- metadata5 %>%
rowwise() %>% mutate(distance = genetic_distance[row_tree, col_tree]) %>% ungroup()
result2 <-result %>% filter(distance > 0 ) %>% 
unite("group_mix", row_group , col_group, sep = "-", remove = FALSE)

group_mix <- result2 %>% select(group_mix)  %>%  distinct()
group_mix <- group_mix$group_mix
group_mix = c ("Wild-Wild", "JAP-JAP","IND-IND" ,"IND-Wild", "JAP-Wild","IND-JAP"  )

result3 <- result2 %>%
filter(group_mix %in% c("IND-Wild", "JAP-Wild","Wild-Wild"))%>%
 mutate(group_mix = case_when(
    group_mix == "IND-Wild" ~ "IND",
    group_mix == "JAP-Wild" ~ "JAP",
    group_mix == "Wild-Wild" ~ "Wild",TRUE ~ group_mix  ))

####均值
result4 <-result3 %>% unite("paired_mix", row_tree , col_tree, sep = "-", remove = FALSE)
result5 <- aggregate(BCvalue ~ paired_mix, result4 , mean, na.rm = TRUE)
result5 <- result5 %>%  rename(mean_BCvalue = 2 )
result4 <-merge(result4,result5,by="paired_mix")

group <- result4  %>% select(paired_mix,group_mix,distance) %>% distinct()
result6 <-merge(result5,group,by="paired_mix")
write.csv(result6,"Fig1c_meanBC_genome_distance.csv")

width=5
height=5
ana_col<-c("distance","mean_BCvalue","group_mix")
xlabs= "Genetic distance"
ylabs="Bray-curtis dissimilarity of phyllosphere microbiome"
colors <- c(Wild="#8EB9ED",JAP ="#E2DB59", IND="#3CE753")
level = c("Wild","JAP","IND")

label="Fig1c"
Regression_sample1 (result6 ,ana_col,level,xlabs=xlabs,ylabs=ylabs, width=width, height=height,label=label,colors=colors)

################Fig1d
taxonomy=read.csv("R1_1_bacteria.O.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)


taxonomy[is.na(taxonomy)]<-0
colorscat <-c("#5C8DC7","#EF7F29","#9A6731","#3EB1DB","#A09E21",
"#6ECCDD","#298C8C","#3E5E9A","#7B287C",
"#B4B4B4","#4B86C6","#AFD8E6",
"#AED8E5","#4B86C6","#FED700","#EC7722")

#colors <- c(Wild="#a4cfa9",JAP ="#B2FA3E", IND="#3DBEED")
colors <- c(Wild="#8EB9ED",JAP ="#E2DB59", IND="#3CE753")
level = c("Wild","JAP","IND")
shapes<- c(Wild=21,JAP =22, IND=23)
label = "Fig1d"
analysisgroup="group3"
width=4
height=4
STACK_ANA(taxnonmy,metadata,analysisgroup,label = label, level = level,width=width, height=height, colorsF = colorscat, topN = 10 )

################Ext Fig1i

library(dplyr)
taxonomy=read.csv("R1_1_bacteria.O.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
group_info <- metadata %>% select(repeats, Distance_toF393)  
df <-merge(group_info,t(taxonomy), by = "row.names")
df <- df %>% tibble::column_to_rownames(var = "Row.names")  

label_LMM="Ext Fig1i"
Effect_size_llm(df,factors="Distance_toF393",random="repeats",label=label_LMM,width=6, height=6,cut_off_FC = 2, cut_off_FDR = 0.05,topn = 10)

#####################Fig1e
data_ROW=read.csv(paste0(label=label_LMM,".","LLM-test-plot.csv"), header=T, row.names=1, comment.char="", stringsAsFactors=F)
data_ROW <- data_ROW %>% filter(Sig != "no" )

ana_col =c("tax","Estimate","SE","p_adj","Sig")
label_fores="Fig1e"
forest_plot (data_ROW ,ana_col,colors=c(Up="#FF690B",Down="#4F95FF"),label=label_fores,width=5, height=5 )

############Ext Fig1g

taxonomy=read.csv("R1_1_bacteria.C.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0
colorscat <-c("#5C8DC7","#EF7F29","#9A6731","#3EB1DB","#A09E21",
"#6ECCDD","#298C8C","#3E5E9A","#7B287C",
"#B4B4B4","#4B86C6","#AFD8E6",
"#AED8E5","#4B86C6","#FED700","#EC7722")

#colors <- c(Wild="#a4cfa9",JAP ="#B2FA3E", IND="#3DBEED")
colors <- c(Wild="#8EB9ED",JAP ="#E2DB59", IND="#3CE753")
level = c("Wild","JAP","IND")
shapes<- c(Wild=21,JAP =22, IND=23)
label = "Ext Fig1g"
analysisgroup="group3"
width=4
height=4
STACK_ANA(taxnonmy,metadata,analysisgroup,label = label, level = level,width=width, height=height, colorsF = colorscat, topN = 10 )

#######################Ext Fig1h

taxonomy=read.csv("R1_1_bacteria.G.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0
colorscat <-c("#5C8DC7","#EF7F29","#9A6731","#3EB1DB","#A09E21",
"#6ECCDD","#298C8C","#3E5E9A","#7B287C",
"#B4B4B4","#4B86C6","#AFD8E6",
"#AED8E5","#4B86C6","#FED700","#EC7722")

#colors <- c(Wild="#a4cfa9",JAP ="#B2FA3E", IND="#3DBEED")
colors <- c(Wild="#8EB9ED",JAP ="#E2DB59", IND="#3CE753")
level = c("Wild","JAP","IND")
shapes<- c(Wild=21,JAP =22, IND=23)
label = "Ext Fig1h"
analysisgroup="group3"
width=4
height=4
STACK_ANA(taxnonmy,metadata,analysisgroup,label = label, level = level,width=width, height=height, colorsF = colorscat, topN = 10 )

################Ext Fig1j

library(dplyr)

taxonomy3=read.csv("R1_1_bacteria.C.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy3[is.na(taxonomy3)] <-0
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

group_info <- metadata %>% select(repeats, Distance_toF393)  
df <-merge(group_info,t(taxonomy3), by = "row.names")
df <- df %>% tibble::column_to_rownames(var = "Row.names")  

label_LMM="Ext Fig1j"
Effect_size_llm(df,factors="Distance_toF393",random="repeats",label=label_LMM,width=6, height=6,cut_off_FC = 2, cut_off_FDR = 0.05,topn = 20)


################Ext Fig1k
library(dplyr)

taxonomy3=read.csv("R1_1_bacteria.G.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy3[is.na(taxonomy3)] <-0
metadata=read.csv("R1_1_metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

group_info <- metadata %>% select(repeats, Distance_toF393)  
df <-merge(group_info,t(taxonomy3), by = "row.names")
df <- df %>% tibble::column_to_rownames(var = "Row.names")  

label_LMM="Ext Fig1k"
Effect_size_llm(df,factors="Distance_toF393",random="repeats",label=label_LMM,width=6, height=6,cut_off_FC = 2, cut_off_FDR = 0.05,topn = 20)






