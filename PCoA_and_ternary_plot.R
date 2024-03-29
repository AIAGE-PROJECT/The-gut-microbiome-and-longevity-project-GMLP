################################################################################
#main codes for Figure 3 of gut microbiome and longevity studies
################################################################################

library(ggplot2)
library(ape)
library(tidyverse)
library(multcomp)
library(extrafont)
library(patchwork)
library(vegan)
library(ggtern)
library(reshape2)
library(ggsci)
library(scales)

rm(list = ls())
set.seed(1289)

#Import the Bray-Cutirs distance matrix, which generated by qiime diversity beta.
dist_BC <- read.table(as.matrix("beta_bray-curtis_distance-matrix.tsv"))
#Pcoa, Fig.3a
pcoa <- pcoa(dist_BC,correction = "none", rn = NULL) 
groups <- read.table("metadata_age_cohort1575.txt",row.names = 1,header = T,sep = '\t',comment.char = '',check.names = F)
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]
groups<-groups[match(rownames(pcoa$vectors),rownames(groups)),]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,groups$Group)
colnames(pcoadata) <- c("sample","PC1","PC2","group")

pcoadata$group <- factor(pcoadata$group,levels = c("20-44","45-65","66-85","90-99","100-117"))

#Adonis test
otu.adonis <- adonis(dist_BC~group,data = pcoadata,permutations = 9999)
write.table(otu.adonis$aov.tab,'beta_bray-curtis_adonis.tsv',sep = '\t',quote = F)
#-------------------------------------------------------------------------------
#Boxplot
df <- pcoadata
df1 <- df %>% group_by(group) %>% summarise(Max = max(PC1))
df2 <- df %>% group_by(group) %>% summarise(Max = max(PC2))
df1$Max <- df1$Max + max(df1$Max)*0.1
df2$Max <- df2$Max + max(df2$Max)*0.1

tuk1 <- aov(PC1~group,data = pcoadata) %>% 
  glht(alternative = 'two.sided',linfct = mcp(group = "Tukey")) %>% cld(alpah = 0.05)
#mod<-aov(PC1~group,data = pcoadata)
#glt<-glht(mod, alternative = 'two.sided', linfct = mcp(group = 'Tukey'))
#summary(glt)
tuk2 <- aov(PC2~group,data = pcoadata) %>% 
  glht(alternative = 'two.sided',linfct = mcp(group = "Tukey")) %>% cld(alpah = 0.05)

res <- data.frame(PC1 = tuk1$mcletters$Letters,PC2 = tuk2$mcletters$Letters,
                   df1 = df1$Max,df2 = df2$Max,group = df1$group)
res$group <- factor(res$group,levels = c("20-44","45-65","66-85","90-99","100-117"))
#-------------------------------------------------------------------------
p1 <- ggplot(pcoadata, aes(PC1, PC2,colour = group)) +
  geom_point(aes(colour = group),size = 1)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )")) +
  scale_colour_manual(values = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF")) +
  stat_ellipse(aes(PC1, PC2),level = 0.95)+
  theme(text = element_text(size = 15))+
  geom_vline(aes(xintercept = 0),linetype = "dotted")+
  geom_hline(aes(yintercept = 0),linetype = "dotted")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.text = element_text(colour = 'black',size = 12),
        legend.title = element_blank(),
        legend.key.height = unit(0.4,"cm"),
        legend.position = 'none',
        panel.grid = element_blank())
#-----------------------------------------------------------------------------------
p2 <- ggplot(pcoadata,aes(group,PC1)) +
  geom_boxplot(aes(fill = group,alpha = 0.1),outlier.colour = NA)+
  scale_fill_manual(values = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF"))+
  geom_text(data = res,aes(x = group,y = df1,label = PC1),
            size = 5,color = "black",fontface = "plain")+
  theme(panel.background = element_rect(fill='white',colour='black'),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = 'black',size = 12,face = "plain"),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())+coord_flip()
#--------------------------------------------------------------------------------------
p3 <- ggplot(pcoadata,aes(group,PC2)) +
  geom_boxplot(aes(fill = group,alpha = 0.1),outlier.colour = NA) +
  scale_fill_manual(values = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF"))+
  geom_text(data = res,aes(x = group,y = df2,label = PC2),
            size = 5,color = "black",fontface = "plain")+
  theme(panel.background = element_rect(fill = 'white',colour = 'black'),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 12,angle = 90,vjust = 0.5,hjust = 0.5,face = "plain"),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

p4 <- ggplot(pcoadata,aes(PC1, PC2))+
  geom_text(aes(x = -0.5,y = 0.6,
                label = paste("PERMANOVA:\ndf = ",otu.adonis$aov.tab$Df[1],"\nR2 = ",
                              round(otu.adonis$aov.tab$R2[1],4),"\np-value = ",
                              otu.adonis$aov.tab$`Pr(>F)`[1],
                              sep = "")),size = 2.8,family = "sans",fontface = 1)+
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

p <- p2+p4+p1+p3 + plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
ggsave(plot = p,"Fig3a_pcoa_of_bray-curtis_age1575.pdf",dpi = 400,width = 6,height = 5)

################################################################################

#Ternary plot, Fig.3b
otutab <- read.delim("feature-table_20-44_66-85_100-117.txt", header = T, row.names = 1,check.names = F)
design <- read.delim("metadata_20-44_66-85_100-117.txt", header = T, row.names = NULL)
design = design[,c("SampleID","Group")]

#Data processing function, 'threshold' indicates that at least one sample is larger than the threshold, 'times' is the scaling of the mean of the points.
data_clean <- function(otu, design, type = c("relative", "absolute"), threshold = 0.001, times = 100){
#Absolute abundance to relative abundance
  if (type == "absolute"){
    otu_relative <- apply(otutab, 2, function(x){x/sum(x)})
  }else {otu_relative <- otutab}
#Keep ASVs with at least one sample larger than the threshold
  idx <- rowSums(otu_relative > threshold) >= 1
  otu_threshold <- as.data.frame(otu_relative[idx, ])
  otu_threshold$OTUs <- row.names(otu_threshold)
#Pivot data from wide to long
  otu_longer <- pivot_longer(data = otu_threshold,
                             cols = -OTUs,
                             names_to = "SampleID",
                             values_to = "value")
#Adds Group in the metadata according to "SampleID"
  merge_data <- merge(otu_longer, design, by = "SampleID")
  otu <- subset(merge_data, select = c("Group","OTUs","value"))
#The mean based on OTUs and Group
  otu_mean <- otu %>% group_by(OTUs, Group) %>%
    summarise(value = mean(value))
#Pivot data from long to wide
  otu_tern <- otu_mean %>%
    group_by(Group, OTUs) %>%
    mutate(index = row_number()) %>%
    pivot_wider(names_from = Group,values_from = value)
  otu_tern <- otu_tern[,-2]
#Size of points
  otu_tern$size <- (apply(otu_tern[2:4], 1, mean))*times
  return(otu_tern)
}

otu_tern <- data_clean(otutab, design, type = "absolute", threshold = 0.05,times = 100)
#Adding taxonomy
taxa <- read.delim("taxonomy.txt", header = T,check.names = F)
otu_tern_taxa <- merge(otu_tern,taxa,by = "OTUs")
#Log2-transformed
otu_tern_taxa$`Log2(size)` <- log2(otu_tern_taxa[,5]+1)
#Plot
platte <- c("#3B4992FF","#EE0000FF","#008B45FF","#631879FF","#F39B7FFF","#4DBBD5FF","#808180FF")
otu_tern_taxa$Taxa <- factor(otu_tern_taxa$Taxa,levels = c("Actinobacteria","Bacteroidetes","Firmicutes","Fusobacteria",
                                                         "Proteobacteria","Verrucomicrobia","Others"))
p_ter <- ggtern(data = otu_tern_taxa, aes(x = `20-44`, y = `66-85`, z = `100-117`,color = Taxa,size = 4)) +
  scale_colour_manual(values = platte) +
  geom_point(aes(fill = Taxa,size = `Log2(size)`), alpha = 1,show.legend = T) +
  scale_size(range = c(0, 6)) + geom_mask() +theme_rgbw(base_size = 12 )+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        legend.position = c(1.1, 0.6),
        axis.text = element_text(colour = 'black',size = 12), axis.ticks = element_blank(),
        text = element_text(size = 15), 
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 15),
        legend.key = element_blank(),legend.text = element_text(size= 12 ),
        legend.box.spacing = unit(0,"line"),
        legend.key.height = unit(0.2,"line"),legend.key.width = unit(0.2,"line"),
        panel.spacing.y = unit(-1,"lines"),panel.spacing.x = unit(-1, "lines"))+
  guides(fill = guide_legend(override.aes = list(size=4)))
ggsave(plot = p_ter,'Fig3b_ternary plot_20-44_66-85_100-117.pdf',height = 4,width = 7,dpi = 600)
