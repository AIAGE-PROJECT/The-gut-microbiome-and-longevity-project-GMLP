################################################################################
#main codes for Figure 4 of gut microbiome and longevity studies
################################################################################
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggsci)
library(scales)
rm(list = ls())
data <- read.table("level-6.txt",header = T,row.names = 1,sep='\t',comment.char = '',check.names = F)

#Selection feature,Fig.4c
feature1 <- c("Klebsiella","Streptococcus","Enterobacter","Rhodococcus","Group")
df1 <- data %>% select(!!feature1)
#Log2 transformation
df1[,1:4]<- log2(df1[,1:4]+1)
#Delete 45-65 and 90-99
df1 <- subset(df1,Group !='45-65' & Group !='90-99')
df1$Group <- factor(df1$Group,levels = c("20-44","66-85","100-117"))
#pal_npg()
color <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF")
my_comparisons <- list(c("20-44", "66-85"),c("66-85", "100-117"))

#Plot
gc <- colnames(df1)
plist <- list()
for (i in 1:length(gc)){
  feature_group <- df1[,c(gc[i],"Group")]
  colnames(feature_group) <- c("abundance","group")
  pb1 <- ggboxplot(feature_group,
                 x = "group",
                 y = "abundance",
                 color = "group",
                 fill = NULL,
                 add = "jitter",alpha = 0.05,shape = 1,
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size = 1,
                 font.label = list(size = 12), 
                 palette = color)+
    theme(panel.background = element_blank())
  pb1 <- pb1+theme(axis.line = element_line(colour = "black"))+theme(axis.title.x = element_blank())
  pb1 <- pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 12,angle = 45,vjust = 1,hjust = 1))
  pb1 <- pb1+theme(axis.text.y = element_text(size = 12))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size=12,face="bold"))
  pb1 <- pb1+theme(legend.position = "NA",plot.title = element_text(hjust = 0.5,face = "italic"))
  pb1 <- pb1+stat_compare_means(method = "wilcox",hide.ns = F,comparisons = my_comparisons,label = "p.format")
  plist[[i]] <- pb1
} 
#Arrange multiple plots into a grid
pall1 <- plot_grid(plist[[1]],plist[[2]],plist[[3]],plist[[4]],ncol = 4)
ggsave(plot = pall1,'Fig4c_20-44vs66-85vs100-117.pdf',height = 4,width = 7,dpi = 600)

#################################################################################
#Selection feature,Fig.4e
feature2 <- c("Akkermansia", "Parabacteroides","Alistipes","Odoribacter","Group")
df2 <- data %>% select(!!feature2)
#Log2 transformation
df2[,1:4] <- log2(df2[,1:4]+1)
#Delete 45-65, 66-85 and 90-99
df2 <- subset(df2,Group !='45-65' & Group !='66-85' & Group !='90-99')
df2$Group <- factor(df2$Group,levels = c("20-44","100-117"))
my_comparisons <- list(c("20-44", "100-117"))

#Plot
gc <- colnames(df2)
plist <- list()
for (i in 1:length(gc)){
  feature_group <- df2[,c(gc[i],"Group")]
  colnames(feature_group) <- c("abundance","group")
  pb1 <- ggboxplot(feature_group,
                 x = "group",
                 y = "abundance",
                 color = "group",
                 fill = NULL,
                 add = "jitter",alpha = 0.05,shape = 1,
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size = 1,
                 font.label = list(size = 12), 
                 palette = color)+
    theme(panel.background = element_blank())
  pb1 <- pb1+theme(axis.line = element_line(colour = "black"))+theme(axis.title.x = element_blank())
  pb1 <- pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 12,angle = 0))
  pb1 <- pb1+theme(axis.text.y = element_text(size = 12))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))
  pb1 <- pb1+theme(legend.position = "NA",plot.title = element_text(hjust = 0.5,face = "italic"))
  pb1 <- pb1+stat_compare_means(method = "wilcox",hide.ns = F,comparisons = my_comparisons,label = "p.format")
  plist[[i]]<-pb1
} 
#Arrange multiple plots into a grid
pall2 <- plot_grid(plist[[1]],plist[[2]],plist[[3]],plist[[4]],ncol = 4)
ggsave(plot = pall2,'Fig4e_20-44vs100-117.pdf',height = 4,width = 7,dpi = 600)

################################################################################
#Supplementary Fig8
feature8 <- read.table("Supplementary Fig8.txt",header = T,row.names = 1,sep = '\t',comment.char = '',check.names = F)
#Log2 transformation
feature8[,1:8] <- log2(feature8[,1:8]+1)
feature8$group <- factor(feature8$group,levels = c("20-44","66-85"))
my_comparisons <- list(c("20-44", "66-85"))
color <- c("#E64B35FF","#4DBBD5FF")
#Plot
gc <- colnames(feature8)
plist <- list()
for (i in 1:length(gc)){
  feature_group <- feature8[,c(gc[i],"group")]
  colnames(feature_group) <- c("abundance","group")
  pb8 <- ggboxplot(feature_group,
                 x = "group",
                 y = "abundance",
                 color = "group",
                 fill = NULL,
                 add = "jitter",alpha = 0.05,shape = 1,
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size = 0.8,
                 font.label = list(size = 12), 
                 palette = color)+
    theme(panel.background = element_blank())
  pb8 <- pb8+theme(axis.line = element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb8 <- pb8+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 12,angle = 0))
  pb8 <- pb8+theme(axis.text.y = element_text(size = 12))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))
  pb8 <- pb8+theme(legend.position = "NA",plot.title = element_text(hjust = 0.5,face = "italic"))
  pb8 <- pb8+stat_compare_means(method = "wilcox",hide.ns = F,comparisons = my_comparisons,label = "p.format",paired = F)
  plist[[i]] <- pb8
} 
#Arrange multiple plots into a grid
pall8 <- plot_grid(plist[[1]],plist[[2]],plist[[3]],plist[[4]],
                 plist[[5]],plist[[6]],plist[[7]],plist[[8]],ncol = 4,nrow = 2)
ggsave(plot = pall8,'Supplementary Fig8_20-44vs66-85.pdf',height = 7,width = 8,dpi = 600)
