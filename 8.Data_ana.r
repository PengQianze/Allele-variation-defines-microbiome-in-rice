
###################################################################R2 ana
###################################################################R2 ana
###################################################################functions

PCOA_ANA2 <- function (taxonomy,metadata,analysisgroup = "TimeGrowth",width=8, height=8,
                      colors= c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                   "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                              "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39"),
                      shapes = c(21,22,23,24,25),label="ana_"
){
  
  library(vegan)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  idx = rownames(metadata) %in% colnames(taxonomy)
  metadata = metadata[idx, , drop = F]
  taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  
  otu_matrix <- vegdist(t(taxonomy),method = "bray")
  colnames(sampFile)[1] = "group"

  pcoa = cmdscale(otu_matrix, k = 3, eig = T)
  points = as.data.frame(pcoa$points)
  eig = pcoa$eig
  points = cbind(points, sampFile)
  colnames(points) = c("x", "y", "z", "group")
  #显著性
  samp_numbers<- nrow(metadata)
  pcoa.fit <- envfit(pcoa ~ group, data=sampFile, perm=100000)
  envfitR2 = paste(analysisgroup , "R2:", pcoa.fit$factors$r, "Pvalue:", pcoa.fit$factors$pvals,"Number:",samp_numbers,sep = ",")
  write.table(envfitR2,paste0(label, "ALL","-envfit.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

}
#####################################
##########################################

box_ANA2 <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,
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
  stat.test<- dfn %>% t_test(variable ~ Group) %>% add_significance(p.col = "p") 
 # readr::write_delim(as.data.frame(stat.test),paste0(label,"_",analysisgroup,col_name,"-df-t-test.txt"))

  test.p = paste(analysisgroup ,"G1:", stat.test$group1,"G2:",stat.test$group2,  "Pvalue:", stat.test$p,"N1:",stat.test$n1,"N2:",stat.test$n2,sep = ",")
  write.table(test.p,paste0(label,"_ALL","-t_test.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  
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
#  write.csv(as.data.frame(sumorder_wide),paste0(analysisgroup,".","-wider-means.csv"),append =T)
  means = paste("SNP:",sumorder_wide$SNP ,"G1:", sumorder_wide$G1,"G2:",sumorder_wide$G2,sep = ",")
  write.table(means, paste0(label,"_ALL","-mean.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  
  sumorder_wide_SE<- sumorder %>%  select( Group,SNP, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
  names(sumorder_wide_SE) <- c("SNP", "G1", "G2")
#  write.csv(as.data.frame(sumorder_wide_SE),paste0(analysisgroup,".","-wider-se.csv"),append =T)
  SE = paste("SNP:",sumorder_wide_SE$SNP ,"G1:", sumorder_wide_SE$G1,"G2:",sumorder_wide_SE$G2,sep = ",")
  write.table(SE, paste0(label,"_ALL","-se.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
} 

##########################################

box_ANA3 <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,
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
  stat.test<- dfn %>% t_test(variable ~ Group) %>% add_significance(p.col = "p") 

 # readr::write_delim(as.data.frame(stat.test),paste0(label,"_",analysisgroup,col_name,"-df-t-test.txt"))

  test.p = paste(col_name ,"G1:", stat.test$group1,"G2:",stat.test$group2,  "Pvalue:", stat.test$p,"N1:",stat.test$n1,"N2:",stat.test$n2,sep = ",")
  write.table(test.p,paste0(label,"ALL","-t_test.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  

  sumorder<-dfn %>%group_by(Group) %>%summarise(
    allmean = mean(variable, na.rm = TRUE),       
    sd1     = sd(variable, na.rm = TRUE),         
    SNP     = col_name,                      
    n1      = sum(!is.na(variable)),               
    se1     = sd1 / sqrt(n1),                   
    .groups = "drop")

#  write.csv(as.data.frame(sumorder),paste0(col_name,".","all_summary.csv"),append =T)

  sumorder_wide <- sumorder %>%  select( Group,SNP,allmean)  %>%
  pivot_wider(names_from = Group, values_from = c(allmean))
#  names(sumorder_wide) <- c("SNP", "G1", "G2")
  means = paste("tax:",sumorder_wide$SNP ,"G.A.type:", sumorder_wide$G,"T.D.type:",sumorder_wide$T,sep = ",")
  write.table(means, paste0(label,".ALL","-mean.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

  sumorder_wide_SE<- sumorder %>%  select( Group,SNP, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
#  names(sumorder_wide_SE) <- c("SNP", "G1", "G2")
  SE = paste("tax:",sumorder_wide_SE$SNP ,"G.A.type:", sumorder_wide_SE$G,"T.D.type:",sumorder_wide_SE$T,sep = ",")
  write.table(SE, paste0(label,".ALL","-se.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

} 


###################################################################R4 ana
#############################################

box_ANA4 <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,level=c(1,2),
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

  test.p = paste(analysisgroup ,"G1:", stat.test$group1,"G2:",stat.test$group2,  "Pvalue:", stat.test$p,"N1:",stat.test$n1,"N2:",stat.test$n2,sep = ",")
  write.table(test.p,paste0(label,"ALL","-t_test.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  

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
 # write.csv(as.data.frame(sumorder_wide),paste0(analysisgroup,".","-wider-means.csv"),append =T)

  means = paste("SNP:",sumorder_wide$SNP ,"G1:", sumorder_wide$G1,"G2:",sumorder_wide$G2,sep = ",")
  write.table(means, paste0(label,"ALL","-mean.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

  sumorder_wide_SE<- sumorder %>%  select( Group,SNP, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
  names(sumorder_wide_SE) <- c("SNP", "G1", "G2")
 # write.csv(as.data.frame(sumorder_wide_SE),paste0(analysisgroup,".","-wider-se.csv"),append =T)

  SE = paste("SNP:",sumorder_wide_SE$SNP ,"G1:", sumorder_wide_SE$G1,"G2:",sumorder_wide_SE$G2,sep = ",")
  write.table(SE, paste0(label,"ALL","-se.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
} 


################################################################################################

box_ANA5 <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,
                     colors = c("#6ED470","#9EA4E5","#00CD66","#7FFF00","#C0FF3E",
                                "#CAFF70","#8B864E","#CDCDB4","#CDCD00",
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

  test.p = paste(col_name ,"G1:", stat.test$group1,"G2:",stat.test$group2,  "Pvalue:", stat.test$p,"N1:",stat.test$n1,"N2:",stat.test$n2,sep = ",")
  write.table(test.p,paste0(label,".","ALL","-t_test.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  

  sumorder<-dfn %>%group_by(Group) %>%summarise(
    allmean = mean(variable, na.rm = TRUE),       
    sd1     = sd(variable, na.rm = TRUE),         
    SNP     = col_name,                      
    n1      = sum(!is.na(variable)),               
    se1     = sd1 / sqrt(n1),                   
    .groups = "drop")

 # write.csv(as.data.frame(sumorder),paste0(col_name,".","all_summary.csv"),append =T)

  sumorder_wide <- sumorder %>%  select( Group,SNP,allmean)  %>%
  pivot_wider(names_from = Group, values_from = c(allmean))
#  names(sumorder_wide) <- c("SNP", "G1", "G2")

  means = paste("tax:",sumorder_wide$SNP ,"HapA.type:", sumorder_wide$HapA,"HapD.type:",sumorder_wide$HapD,sep = ",")
  write.table(means, paste0(label,".",".ALL","-mean.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

  sumorder_wide_SE<- sumorder %>%  select( Group,SNP, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
#  names(sumorder_wide_SE) <- c("SNP", "G1", "G2")

  SE = paste("tax:",sumorder_wide_SE$SNP ,"HapA.type:", sumorder_wide_SE$HapA,"HapD.type:",sumorder_wide_SE$HapD,sep = ",")
  write.table(SE, paste0(label,".",".ALL","-se.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  
} 

##############################################

box_ANA6 <- function (metadata,analysisgroup = "TimeGrowth", col_name = "HCA",ylabs = "ylab",widths = 2.5,heights =4,
                     colors = c("#6ED470","#9EA4E5","#00CD66","#7FFF00","#C0FF3E",
                                "#CAFF70","#8B864E","#CDCDB4","#CDCD00",
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

  test.p = paste(col_name ,"G1:", stat.test$group1,"G2:",stat.test$group2,  "Pvalue:", stat.test$p,"N1:",stat.test$n1,"N2:",stat.test$n2,sep = ",")
  write.table(test.p,paste0(label,".","ALL","-t_test.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  

  sumorder<-dfn %>%group_by(Group) %>%summarise(
    allmean = mean(variable, na.rm = TRUE),       
    sd1     = sd(variable, na.rm = TRUE),         
    SNP     = col_name,                      
    n1      = sum(!is.na(variable)),               
    se1     = sd1 / sqrt(n1),                   
    .groups = "drop")

#  write.csv(as.data.frame(sumorder),paste0(col_name,".","all_summary.csv"),append =T)

  sumorder_wide <- sumorder %>%  select( Group,SNP,allmean)  %>%
  pivot_wider(names_from = Group, values_from = c(allmean))
#  names(sumorder_wide) <- c("SNP", "G1", "G2")

  means = paste("tax:",sumorder_wide$SNP ,"nivala:", sumorder_wide$nivala,"rufipogon:",sumorder_wide$rufipogon,sep = ",")
  write.table(means, paste0(label,".",".ALL","-mean.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)

  sumorder_wide_SE<- sumorder %>%  select( Group,SNP, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
#  names(sumorder_wide_SE) <- c("SNP", "G1", "G2")

  SE = paste("tax:",sumorder_wide_SE$SNP ,"nivala:", sumorder_wide_SE$nivala,"rufipogon:",sumorder_wide_SE$rufipogon,sep = ",")
  write.table(SE, paste0(label,".",".ALL","-se.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
 
} 

###############################################################
###############################################################Data
###############################################################
############order_level_all

library(dplyr)
taxonomy=read.csv("Data_bacteria.O.RA_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("Data_metadata4_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

for (col_name in colnames(metadata)) {
  if (col_name == "group3") {next}
  else {
metadata2 <- metadata %>% filter(.data[[col_name]] != "N"  )
analysisgroup = col_name
ylabs = col_name
label<-paste("all_locus_SNP_beta_effect","_")
colors <- c("#5C8DC7","#D8EF29", "#9A6731","#E2BA8F")
shapes<- c(21,22,23,24)
width=4
height=4
PCOA_ANA2 (taxonomy,metadata2,analysisgroup ,width=width, height=height,colors= colors,shapes = shapes,label=label)
  }
  }

#################alpha_diversity

library(dplyr)
metadata=read.csv("Data_metadata4_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
chao1_alpha=read.csv("Data_chao1_alpha.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

new_alpha = merge(chao1_alpha, metadata, by= "row.names")
rownames(new_alpha) <- new_alpha[, "Row.names"]
new_alpha <- new_alpha[, !colnames(new_alpha) %in% c("Row.names","group3")]
widths = 4
heights = 4

for (col_names in colnames(new_alpha)) {
  if (col_names == "Chao1") {next}
  else {
metadata2 <- new_alpha %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = "Chao1"
label<-paste("all_locus_","_",col_name )
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }
#################alpha_diversity

library(dplyr)
metadata=read.csv("Data_metadata4_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
chao1_alpha=read.csv("Data_shannon_alpha.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

new_alpha = merge(chao1_alpha, metadata, by= "row.names")
rownames(new_alpha) <- new_alpha[, "Row.names"]
new_alpha <- new_alpha[, !colnames(new_alpha) %in% c("Row.names","group3")]
widths = 4
heights = 4

for (col_names in colnames(new_alpha)) {
  if (col_names == "Shannon") {next}
  else {
metadata2 <- new_alpha %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = "Shannon"
label<-paste("all_locus_","_",col_name )
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

#################psedomonas_effect

library(dplyr)
taxonomy=read.csv("Data_bacteria.O.RA_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("Data_metadata4_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

taxonomy = t(taxonomy)
taxonomy2 <- taxonomy[, colnames(taxonomy) %in% c("Pseudomonadales","Enterobacterales")]

new_tax = merge(taxonomy2, metadata, by= "row.names")
rownames(new_tax) <- new_tax[, "Row.names"]
new_tax <- new_tax[, !colnames(new_tax) %in% c("Row.names","Enterobacterales","group3")]
widths = 4
heights = 4

for (col_names in colnames(new_tax)) {
  if (col_names == "Pseudomonadales") {next}
  else {
metadata2 <- new_tax %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = "Pseudomonadales"
label<-paste("all_locus_tax_effect","_",col_name )
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

#################Enterobacterales_effect

library(dplyr)
taxonomy=read.csv("Data_bacteria.O.RA_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("Data_metadata4_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

taxonomy = t(taxonomy)
taxonomy2 <- taxonomy[, colnames(taxonomy) %in% c("Pseudomonadales","Enterobacterales")]

new_tax = merge(taxonomy2, metadata, by= "row.names")
rownames(new_tax) <- new_tax[, "Row.names"]
new_tax <- new_tax[, !colnames(new_tax) %in% c("Row.names","Pseudomonadales","group3")]

widths = 4
heights = 4

for (col_names in colnames(new_tax)) {
  if (col_names == "Enterobacterales") {next}
  else {
metadata2 <- new_tax %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = "Enterobacterales"
label<-paste("all_locus_tax_effect","_",col_name )
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

########################A,D hapoltype

library(dplyr)
taxonomy=read.csv("Data_bacteria.O.RA_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("Data_metadata4_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

metadata2 <- metadata[, colnames(metadata) %in% c("chr.2.25894232","group3")]

new_tax = merge(metadata2, t(taxonomy), by= "row.names")
rownames(new_tax) <- new_tax[, "Row.names"]


new_tax <- new_tax[, !colnames(new_tax) %in% c("Row.names","group3")]
#col_name = "Chao1"
analysisgroup = "chr.2.25894232"
new_tax <- new_tax %>% filter( chr.2.25894232 != "N"  )

widths = 4
heights = 4

for (col_names in colnames(new_tax)) {
  if (col_names == "chr.2.25894232") {next}
  else {
metadata2 <- new_tax %>% filter( chr.2.25894232 != "N"  )
ylabs = col_names
col_name = col_names
label<-paste( "chr.2.25894232","_","effect order" )
box_ANA3 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

################# climate factor effect #######A,D hapoltype all rice 

library(dplyr)
trail=read.csv("Data_all_geo_climeate2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

trail[is.na(trail)]<-"N" 

new_tra <- trail %>% filter( HAP != "N"  ) %>% 
select(HAP, average_solar,temperature,precipitation,Longitude,Latitude)
analysisgroup = "HAP"
widths = 4
heights = 4

for (col_names in colnames(new_tra)) {
  if (col_names == "HAP") {next}
  if (col_names == "Latitude") {
metadata2 <- new_tra %>% filter(.data[[col_names]] != "N"  )
metadata2$Latitude <- abs(as.numeric(metadata2$Latitude))
ylabs = col_names
col_name = col_names
label<-paste("all_geo_climeate","_")
box_ANA5 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)  }
  else {
metadata2 <- new_tra %>% filter(.data[[col_names]] != "N"  )
ylabs = col_names
col_name = col_names
label<-paste("all_geo_climeate","_")
box_ANA5 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

########### climate factor effect ###### nivala_rufipogon

library(dplyr)
trail=read.csv("Data_all_geo_climeate2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

trail[is.na(trail)]<-"N" 

new_tra <- trail %>% filter( subgroup != "N"  ) %>%  filter( subgroup == "nivala" | subgroup == "rufipogon"  ) %>% 
select(subgroup, average_solar,temperature,precipitation,Longitude,Latitude)
analysisgroup = "subgroup"
widths = 4
heights = 4

for (col_names in colnames(new_tra)) {
  if (col_names == "subgroup") {next}
  if (col_names == "Latitude") {
metadata2 <- new_tra %>% filter(.data[[col_names]] != "N"  )
metadata2$Latitude <- abs(as.numeric(metadata2$Latitude))
ylabs = col_names
col_name = col_names
label<-paste("nivala_rufipogon_geo_climeate2","_")
box_ANA6 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)  }
  else {
metadata2 <- new_tra %>% filter(.data[[col_names]] != "N"  )
ylabs = col_names
col_name = col_names
label<-paste("nivala_rufipogon_geo_climeate2","_" )
box_ANA6 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

######################## A,D hapoltype

library(dplyr)
trail=read.csv("Data_rdp2_trails_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

trail[is.na(trail)]<-"N" 
new_tra <- trail %>% filter( HAP != "N"  )
analysisgroup = "HAP"
widths = 4
heights = 4

for (col_names in colnames(new_tra)) {
  if (col_names == "HAP") {next}
  else {
metadata2 <- new_tra %>% filter(.data[[col_names]] != "N"  )
ylabs = col_names
col_name = col_names
label<-paste("rdp2_rice_traits","_")
box_ANA5 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

########################traits _ solar  corrlations

library(dplyr)
library(linkET)
library(ggplot2)

env_trails <- read.csv("Data_cilmate_trails2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
env_trails[]<-lapply(env_trails, as.numeric)
env_trails <- na.omit(env_trails)

env <- env_trails[1:3]
trails<- env_trails[4:ncol(env_trails)]

library(psych)
ct <- corTest(env_trails, method = "pearson")

write.csv(ct$r,"cilmate_traits.pearson.r.csv")
write.csv(ct$p.adj,"cilmate_traits.pearson.p.adj.csv")
write.csv(ct$p,"cilmate_traits.pearson.p.csv")

################################################################### SNP_effect_trails
################################################################### SNP_effect_trails

library(dplyr)

trail=read.csv("Data_irri_snp_trial2_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
snp_info=read.csv("Data_SNP_cry2_trail.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
trail[is.na(trail)]<-"N"

trail2 <- trail[, colnames(trail) %in% c("cult_repro","liglt","hdg_80head","grwd","grwt100","mat")]
trail2[] <- lapply(trail2,as.numeric)
trail2[is.na(trail2)]<-"N"

##################################cult_repro

trail3<- trail2 %>% select(cult_repro) %>% filter(cult_repro != "N"  )
new_tax1 = merge(trail3, snp_info, by= "row.names")
rownames(new_tax1) <- new_tax1[, "Row.names"]
new_tax <- new_tax1[, !colnames(new_tax1) %in% c("Row.names")] 

widths = 4
heights = 4
effect_trial = "cult_repro"

for (col_names in colnames(new_tax)) {
  if (col_names == effect_trial) {next}
  else {
metadata2 <- new_tax %>% filter(.data[[col_names]] != "N"  )
#snp_n <- length(unique(metadata2[[col_names]]))
analysisgroup =  col_names
ylabs = col_names
col_name = effect_trial
label<-paste("trail_effect_",effect_trial,"_")
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

##################################hdg_80head

trail3<- trail2 %>% select(hdg_80head) %>% filter(hdg_80head != "N"  )

new_tax1 = merge(trail3, snp_info, by= "row.names")
rownames(new_tax1) <- new_tax1[, "Row.names"]
new_tax <- new_tax1[, !colnames(new_tax1) %in% c("Row.names")] 

widths = 4
heights = 4
effect_trial = "hdg_80head"

for (col_names in colnames(new_tax)) {
  if (col_names == effect_trial) {next}
  else {
metadata2 <- new_tax %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = effect_trial
label<-paste("trail_effect_",effect_trial,"_")
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

##################################liglt

trail3<- trail2 %>% select(liglt) %>% filter(liglt != "N"  )

new_tax1 = merge(trail3, snp_info, by= "row.names")
rownames(new_tax1) <- new_tax1[, "Row.names"]
new_tax <- new_tax1[, !colnames(new_tax1) %in% c("Row.names")] 

widths = 4
heights = 4
effect_trial = "liglt"

for (col_names in colnames(new_tax)) {
  if (col_names == effect_trial) {next}
  else {
metadata2 <- new_tax %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = effect_trial
label<-paste("trail_effect_",effect_trial,"_")
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

##################################grwd
trail3<- trail2 %>% select(grwd) %>% filter(grwd != "N"  )

new_tax1 = merge(trail3, snp_info, by= "row.names")
rownames(new_tax1) <- new_tax1[, "Row.names"]
new_tax <- new_tax1[, !colnames(new_tax1) %in% c("Row.names")] 

widths = 4
heights = 4
effect_trial = "grwd"

for (col_names in colnames(new_tax)) {
  if (col_names == effect_trial) {next}
  else {
metadata2 <- new_tax %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = effect_trial
label<-paste("trail_effect_",effect_trial,"_")
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

##################################grwt100
trail3<- trail2 %>% select(grwt100) %>% filter(grwt100 != "N"  )

new_tax1 = merge(trail3, snp_info, by= "row.names")
rownames(new_tax1) <- new_tax1[, "Row.names"]
new_tax <- new_tax1[, !colnames(new_tax1) %in% c("Row.names")] 

widths = 4
heights = 4
effect_trial = "grwt100"

for (col_names in colnames(new_tax)) {
  if (col_names == effect_trial) {next}
  else {
metadata2 <- new_tax %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = effect_trial
label<-paste("trail_effect_",effect_trial,"_" )
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

##################################"mat"
trail3<- trail2 %>% select(mat) %>% filter(mat != "N"  )

new_tax1 = merge(trail3, snp_info, by= "row.names")
rownames(new_tax1) <- new_tax1[, "Row.names"]
new_tax <- new_tax1[, !colnames(new_tax1) %in% c("Row.names")] 

widths = 4
heights = 4
effect_trial = "mat"

for (col_names in colnames(new_tax)) {
  if (col_names == effect_trial) {next}
  else {
metadata2 <- new_tax %>% filter(.data[[col_names]] != "N"  )
analysisgroup =  col_names
ylabs = col_names
col_name = effect_trial
label<-paste("trail_effect_",effect_trial,"_" )
box_ANA2 (metadata2 ,analysisgroup , col_name= col_name, label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

################################################################################################
################################################################################################link effect

library(dplyr)
metadata=read.csv("Data_trail_hap_sample.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata[is.na(metadata)]<-"N"

metadata2 <- metadata %>% filter( HAP != "N"  ) 
metadata2 <- metadata2[, colnames(metadata2) %in% c("HAP","cult_repro","liglt","hdg_80head","grwd","grwt100","mat")]
new_alpha1 <- metadata2 %>% filter( HAP %in% c("HapA_cryA", "HapD_cryB") ) 

colors <- c(
HapA_cryA="#3CE753",
HapA_cryB="#3CE7CD",
HapD_cryB="#8EB9ED",
HapD_cryA="#A1DFF0")

level = c("HapA_cryA","HapD_cryB")
analysisgroup = "HAP"
widths = 4
heights =4


for (col_name in colnames(new_alpha1)) {
  if (col_name == analysisgroup) {next}
  else {
col_name = col_name
ylabs = col_name
label<-paste("trial1","_","HapA_cryA" )
box_ANA4 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

##############################

new_alpha2 <- metadata2 %>% filter( HAP %in% c("HapA_cryB","HapD_cryA") ) 

colors <- c(
HapA_cryA="#3CE753",
HapA_cryB="#3CE7CD",
HapD_cryB="#8EB9ED",
HapD_cryA="#A1DFF0")

level = c("HapA_cryB","HapD_cryA")
analysisgroup = "HAP"
widths = 4
heights =4


for (col_name in colnames(new_alpha2)) {
  if (col_name == analysisgroup) {next}
  else {
col_name = col_name
ylabs = col_name
label<-paste("trial2","_","HapA_cryB" )
box_ANA4 (new_alpha2 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }


#############################

library(dplyr)
metadata=read.csv("Data_Micro.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata[is.na(metadata)]<-"N"


metadata2 <- metadata %>% filter( HAP != "N"  ) 
metadata2 <- metadata2[, colnames(metadata2) %in% c("HAP","Enterobacterales","Pseudomonadales")]
new_alpha1 <- metadata2 %>% filter( HAP %in% c("HapA_cryA", "HapD_cryB") ) 

colors <- c(
HapA_cryA="#3CE753",
HapA_cryB="#3CE7CD",
HapD_cryB="#8EB9ED",
HapD_cryA="#A1DFF0")

level = c("HapA_cryA","HapD_cryB")
analysisgroup = "HAP"
widths = 4
heights =4


for (col_name in colnames(new_alpha1)) {
  if (col_name == analysisgroup) {next}
  else {
col_name = col_name
ylabs = col_name
label<-paste("Micro1","_","HapA_cryA" )
box_ANA4 (new_alpha1 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }

##############################

new_alpha2 <- metadata2 %>% filter( HAP %in% c("HapA_cryB","HapD_cryA") ) 

colors <- c(
HapA_cryA="#3CE753",
HapA_cryB="#3CE7CD",
HapD_cryB="#8EB9ED",
HapD_cryA="#A1DFF0")

level = c("HapA_cryB","HapD_cryA")
analysisgroup = "HAP"
widths = 4
heights =4


for (col_name in colnames(new_alpha2)) {
  if (col_name == analysisgroup) {next}
  else {
col_name = col_name
ylabs = col_name
label<-paste("Micro2","_","HapA_cryB" )
box_ANA4 (new_alpha2 ,analysisgroup , col_name= col_name, colors = colors,level=level,  label=label, ylabs = ylabs,widths = widths,
heights =heights)
  }
  }
