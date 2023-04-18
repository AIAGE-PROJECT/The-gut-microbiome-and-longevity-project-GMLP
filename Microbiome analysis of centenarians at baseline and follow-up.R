################################################################################
#main codes for Figure 5 of gut microbiome and longevity studies
################################################################################

library(ggplot2)
library(tidyverse)
library(vegan)
library(ape)
library(ggpubr)
library(cowplot)
library(reshape)
library(ggsci)
library(scales)
rm(list=ls())
set.seed(7963)
#PCOA, Fig5a
data <- read.table("feature-table.txt", header = T,check.names = F,sep = "\t",row.names = 1) %>% t()
dist <- vegdist(data,method = "bray") %>% as.matrix()
pcoa <- pcoa(dist,correction = "none", rn = NULL)

groups <- read.table("metadata_follow-up.txt",row.names = 1,header = T,sep = '\t',comment.char = '',check.names = F)
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]

groups <- groups[match(rownames(pcoa$vectors),rownames(groups)),]
pcoadata <- data.frame(rownames(pcoa$vectors),PC1,PC2,groups$Group)
colnames(pcoadata) <- c("sample","PC1","PC2","group")

#Adonis test
otu.adonis <- adonis(dist~group,data = pcoadata,distance = "bray",permutations = 9999)

#Plot
p <- ggplot(pcoadata, aes(PC1, PC2,colour = group)) +
  geom_point(aes(colour = group),size = 2.5)+
  labs(title = paste0("PCoA of Bray-Curtis (P=",otu.adonis$aov.tab$`Pr(>F)`,",", "R2=",round(otu.adonis$aov.tab$R2,4),")"),
       x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PC1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PC2 ( ", ., "%", " )")) +
  scale_colour_manual(values = c("#ED0000FF","#42B540FF")) +
  theme(text = element_text(size = 10))+
  geom_vline(aes(xintercept = 0),linetype = "dotted")+
  geom_hline(aes(yintercept = 0),linetype = "dotted")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black',size = 12),
        legend.title = element_blank(),
        legend.key.height = unit(0.2,"cm"),legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = c(0.835, 0.915))+
  theme(panel.grid = element_blank())
p1<- p+stat_ellipse(aes(fill = group),geom = "polygon",
               level = 0.95,alpha = 0.05)
ggsave('Fig5a-PCoA_follow-up-bray.pdf', p1, width = 4, height = 4,dpi = 600)

################################################################################

#¦Á-diversity, Fig5e
alpha <- read.table('alpha-diversity.txt',header = T,row.names = 1,sep = '\t',check.names = F)
alpha$Group <- factor(alpha$Group,levels = c("Baseline","Follow-up"))
mycompare <- list(c('Baseline','Follow-up'))

pc <- ggplot(alpha,aes(x = Group,y = `chao1`,color = Group))+ 
  geom_boxplot(size = 1)+ 
  geom_jitter(aes(fill = Group),width = 0.2,shape = 1,size = 4,stroke = 1.5)+ 
  scale_color_manual(values = c("#ED0000FF","#42B540FF"))+ 
  stat_compare_means(comparisons = mycompare,label = "p.format",method = 'wilcox',
                     paired = T,label.y = 465)+ylim(60,500)+
  ggtitle("Chao1 index")+ 
  theme_bw()+ 
  theme(legend.position = "none", 
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("")+xlab("")
pp <- ggplot(alpha,aes(x = Group,y = `pielou_e`,color = Group))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(fill = Group),width = 0.2,shape = 1,size = 4,stroke = 1.5)+
  scale_color_manual(values = c("#ED0000FF","#42B540FF"))+
  stat_compare_means(comparisons = mycompare,label = "p.format",method = 'wilcox',
                     paired = T,label.y = 0.85)+ylim(0.45,0.89)+
  ggtitle("Pielou's index")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("")+xlab("")
ps <- ggplot(alpha,aes(x = Group,y = `shannon`,color = Group))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(fill = Group),width = 0.2,shape = 1,size = 4,stroke = 1.5)+
  scale_color_manual(values = c("#ED0000FF","#42B540FF"))+
  stat_compare_means(comparisons = mycompare,label = "p.format",method = 'wilcox',
                     paired = T,label.y = 7.25)+ylim(3.2,7.5)+
  ggtitle("Shannon index")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("")+xlab("")
p_alpha <- plot_grid(pc,pp,ps,ncol = 3)
ggsave(plot = p_alpha,'Fig5e_alpha.pdf',height = 4,width = 8,dpi = 600)
################################################################################

#¦Á-diversity by health status, Supplementary Fig7b

alpha_status <- alpha[!alpha$Status=='Frail',]
alpha_status <- alpha_status[,-c(4:6,8)]

alpha_status$Group_status <- factor(alpha_status$Group_status,levels = c("HB","HF","LHB","LHF"))
mycompare_status <- list(c('HB','HF'),c('LHB','LHF'),c('HB','LHB'),c('HF','LHF'))
mycolor <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")

#plot
gc <- colnames(alpha_status)
name <- list()
for (i in 1:(length(gc)-1)){
  feature_group <- alpha_status[,c(gc[i],"Group_status")]
  colnames(feature_group) <- c("abundance","group")
  name <-  as.character(colnames(alpha_status)[i])
  pb1 <- ggboxplot(feature_group,
                 x = "group",
                 y = "abundance",
                 color = "group",
                 fill = NULL,
                 add = "jitter",alpha = 0.01,shape = 1,
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size = 1,
                 palette = mycolor)+
    theme(panel.background = element_blank())
  pb1 <- pb1+theme(axis.line = element_line(colour = "black"))+theme(axis.title.x = element_blank())
  pb1 <- pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 12,hjust = 0.3,vjust = 0.5,angle = 0))
  pb1 <- pb1+theme(axis.text.y = element_text(size = 12))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))
  pb1 <- pb1+theme(legend.position = "NA",plot.title = element_text(hjust = 0.5,face = "plain"))
  pb1 <- pb1+stat_compare_means(method = "wilcox",comparisons = mycompare_status,label = "p.format",paired = F)
  ggsave(filename = paste(name, '_status.pdf', sep = ''),width = 3, height = 4, device = "pdf")
}

################################################################################

#Fig5d
dt <- dist
dt1 <- dt[seq(0,nrow(dt),2),]
dt1 <- dt1[,seq(1,ncol(dt1),2)]
df <- diag(as.matrix(dt1))
dt_paired <- data.frame(rownames(dt1),df)
rownames(dt_paired) <- dt_paired[,1]
colnames(dt_paired) <- c("sampleid","BC")
alpha_baseline <- subset(alpha,Group=='Baseline')
alpha_BC <- merge(alpha_baseline,dt_paired,by = "row.names", all = T)
#Plot
p1 <- ggplot(data = alpha_BC, aes(x = `chao1`, y = `BC`))+
  geom_point(color = "#E64B35FF",size = 3.5,alpha = 1.2,shape = 1,stroke = 1.5)+  
  stat_smooth(formula = y ~ x,method ="glm",se = TRUE,color = "#7E6148FF")+
  stat_cor(data = alpha_BC, method = "spearman",label.x.npc = 0.05,label.y.npc = 0.15,size = 5)+
  xlab("Baseline Chao1 index")+ylab("Bray-Curtis distance between paired samples")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size = 14),   
        axis.title.x = element_text(colour = 'black', size = 14),
        axis.title.y = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black',size = 14),
        legend.title = element_blank(),
        panel.grid = element_blank())

p2 <- ggplot(data = alpha_BC, aes(x = `pielou_e`, y = `BC`))+
  geom_point(color = "#E64B35FF",size = 3.5,alpha = 1.2,shape = 1,stroke = 1.5)+ 
  stat_smooth(formula = y ~ x,method = "glm",se = TRUE,color = "#7E6148FF")+
  stat_cor(data = alpha_BC, method = "spearman",label.x.npc = 0.05,label.y.npc = 0.15,size = 5)+
  xlab("Baseline Pielou's index")+ylab("")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size = 14), 
        axis.title.x = element_text(colour = 'black', size = 14),
        axis.title.y = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black',size = 14),
        legend.title = element_blank(),
        panel.grid = element_blank())

p3 <- ggplot(data = alpha_BC, aes(x = `shannon`, y = `BC`))+
  geom_point(color = "#E64B35FF",size = 3.5,alpha = 1.2,shape = 1,stroke = 1.5)+ 
  stat_smooth(formula = y ~ x,method = "glm",se = TRUE,color = "#7E6148FF")+
  stat_cor(data = alpha_BC, method = "spearman",label.x.npc = 0.05,label.y.npc = 0.15,size = 5)+
  xlab("Baseline Shannon index")+ylab("")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size = 14), 
        axis.title.x = element_text(colour = 'black', size = 14),
        axis.title.y = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black',size = 14),
        legend.title = element_blank(),
        panel.grid = element_blank())

p<-plot_grid(p1,p2,p3,ncol = 3)
ggsave(plot = p,'Fig5d-paied_BC_alpha.pdf',height = 4,width = 11,dpi = 600)

################################################################################

#Fig5b
#Bray distance of baseline
dt_b <- dt[seq(0,nrow(dt),2),]
dt_b <- dt_b[,seq(0,ncol(dt_b),2)]
df_b <- dt_b[upper.tri(dt_b, diag = F)]
group_b <- rep("Baseline",990)
bc_b <- data.frame(df_b,group_b)
colnames(bc_b) <- c("BC","group")

#Bray distance of folow up
dt_f <- dt[seq(1,nrow(dt),2),]
dt_f <- dt_f[,seq(1,ncol(dt_f),2)]
df_f <- dt_f[upper.tri(dt_f, diag = F)]
group_f <- rep("Follow-up",990)
bc_f <- data.frame(df_f,group_f)
colnames(bc_f) <- c("BC","group")

#Bray distance of paired samples
group_p <- rep("Paired",45)
bc_p <- data.frame(df,group_p)
colnames(bc_p) <- c("BC","group")

#Combine by rows
df_bc <- rbind(rbind(bc_b,bc_f),bc_p)

#Plot
df_bc$group <- factor(df_bc$group,levels = c("Baseline","Follow-up","Paired"))
mycompare <- list(c('Baseline','Follow-up'),c('Baseline','Paired'),c('Follow-up','Paired'))

p_bc <- ggplot(df_bc, aes(x = group, y = `BC`)) + 
  geom_violin(aes(fill = group),trim = T) + 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))+ 
  scale_fill_manual(values = c("#ED0000FF","#42B540FF","#00468BFF"))+
  stat_compare_means(comparisons = mycompare,label = "p.format",method = 'wilcox',paired = F,label.y = c(1.0,1.08,1.16))+
  ylim(0.45,1.22)+
  xlab("")+ylab("Bray-Curtis distance")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size = 18),
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.text = element_text(colour = 'black',size = 15),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')
ggsave(plot = p_bc,'Fig5b-BC_distance.pdf',height = 4,width = 3.5,dpi = 600)

#Fig5c
#Healthy baseline
name_hb <- c(rownames(subset(groups,Group=='Baseline' & Health_Status=='H')))
dt_hb <- dt[name_hb,]
dt_hb <- dt_hb[,name_hb]
df_hb <- dt_hb[upper.tri(dt_hb, diag = F)]
group_hb <- rep("HB",210)
bc_hb <- data.frame(df_hb,group_hb)
colnames(bc_hb) <- c("BC","group")

#Less healthy baseline
name_lhb <- c(rownames(subset(groups,Group=='Baseline' & Health_Status=='LH')))
dt_lhb <- dt[name_lhb,]
dt_lhb <- dt_lhb[,name_lhb]
df_lhb <- dt_lhb[upper.tri(dt_lhb, diag = F)]
group_lhb <- rep("LHB",171)
bc_lhb <- data.frame(df_lhb,group_lhb)
colnames(bc_lhb) <- c("BC","group")
 
#Healthy follow up
name_hf <- c(rownames(subset(groups,Group=='Follow-up' & Health_Status=='H')))
dt_hf <- dt[name_hf,]
dt_hf <- dt_hf[,name_hf]
df_hf <- dt_hf[upper.tri(dt_hf, diag = F)]
group_hf <- rep("HF",105)
bc_hf <- data.frame(df_hf,group_hf)
colnames(bc_hf) <- c("BC","group")

#Less healthy follow up
name_lhf <- c(rownames(subset(groups,Group=='Follow-up' & Health_Status=='LH')))
dt_lhf <- dt[name_lhf,]
dt_lhf <- dt_lhf[,name_lhf]
df_lhf <- dt_lhf[upper.tri(dt_lhf, diag = F)]
group_lhf <- rep("LHF",351)
bc_lhf <- data.frame(df_lhf,group_lhf)
colnames(bc_lhf) <- c("BC","group")

#Combine by rows
df_status <- rbind(rbind(rbind(bc_hb,bc_hf),bc_lhb),bc_lhf)
#Plot
df_status$group <- factor(df_status$group,levels = c("HB","HF","LHB","LHF"))
mycompare <- list(c('HB','HF'),c('LHB','LHF'),c('HB','LHB'),c('HF','LHF'))

p_status <- ggplot(df_status, aes(x = group, y = `BC`)) + 
  geom_violin(aes(fill = group),trim = T) + 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))+ 
  scale_fill_npg()+   
  stat_compare_means(comparisons = mycompare,label = "p.format",method = 'wilcox',paired = F,label.y = c(1.0,1.0,1.08,1.16))+
  ylim(0.5,1.22)+
  xlab("")+ylab("Bray-Curtis distance")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size = 15), 
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.text.x = element_text(hjust = 0.5,vjust = 0.5,angle = 0),
        axis.text = element_text(colour = 'black',size = 15),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')

ggsave(plot = p_status,'Fig5c_BC_distance_status.pdf',height = 4,width = 3.5,dpi = 600)
################################################################################

#Bray Curtis distance to 20-44,Fig5f
df_bf <- read.table("feature-table_20-44_baseline_follow-up.txt", header = T,sep = "\t",row.names = 1,check.names = F) %>% t()
dist_bf <- vegdist(df_bf,method = "bray") %>% as.matrix()
groups1 <- read.table("metadata-20-44_baseline_follow-up.txt",row.names = 1,header = T,sep='\t',comment.char = '',check.names = F)
name_y <- c(rownames(subset(groups1,Group=='20-44')))
name_b <- c(rownames(subset(groups1,Group=='Baseline')))
name_f <- c(rownames(subset(groups1,Group=='Follow-up')))
dt_bfy <- dist_bf[,name_y]

#Bray Curtis of baseline to 20-44 
dt_by <- dt_bfy[name_b,]
dt_by <- matrix(dt_by)
group_by <- rep("Baseline",nrow(dt_by))
bc_by <- data.frame(dt_by,group_by)
colnames(bc_by) <- c("BC","group")

#Bray Curtis of follow up to 20-44
dt_fy <- dt_bfy[name_f,]
dt_fy <- matrix(dt_fy)
group_fy <- rep("Follow-up",nrow(dt_fy))
bc_fy <- data.frame(dt_fy,group_fy)
colnames(bc_fy) <- c("BC","group")

#Combine by rows
df_bfy <- rbind(bc_by,bc_fy)
#Plot
df_bfy$group <- factor(df_bfy$group,levels = c("Baseline","Follow-up"))
mycompare <- list(c('Baseline','Follow-up'))

p_bfy <- ggplot(df_bfy, aes(x = group, y = `BC`)) + 
  geom_violin(aes(fill = group),trim = T) + 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))+
  scale_fill_npg()+
  stat_compare_means(comparisons=mycompare,label = "p.format",method = 'wilcox',
                     paired = F,label.y = 1.0)+ylim(0.55,1.05)+
  xlab("")+ylab("Bray-Curtis distance to 20-44")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size = 15), 
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.text.x = element_text(hjust = 0.5,vjust = 0.5,angle = 0),
        axis.text = element_text(colour = 'black',size = 15),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')
ggsave(plot = p_bfy,'Fig5f_BC_distance_to_20-44.pdf',height = 4,width = 3,dpi = 600)
################################################################################

#Bray Curtis distance to 20-44 by health status, Fig5g
name_hb <- c(rownames(subset(groups1,Group=='Baseline' & Health_status=='H')))
name_lhb <- c(rownames(subset(groups1,Group=='Baseline' & Health_status=='LH')))
name_hf <- c(rownames(subset(groups1,Group=='Follow-up' & Health_status=='H')))
name_lhf <- c(rownames(subset(groups1,Group=='Follow-up' & Health_status=='LH')))
#Bray Curtis of healthy baseline to 20-44 
dt_hby <- dt_bfy[name_hb,]
dt_hby <- matrix(dt_hby)
group_hby <- rep("HB",nrow(dt_hby))
bc_hby <- data.frame(dt_hby,group_hby)
colnames(bc_hby) <- c("BC","group")
#Bray Curtis of less healthy baseline to 20-44 
dt_lhby <- dt_bfy[name_lhb,]
dt_lhby <- matrix(dt_lhby)
group_lhby <- rep("LHB",nrow(dt_lhby))
bc_lhby <- data.frame(dt_lhby,group_lhby)
colnames(bc_lhby) <- c("BC","group")

#Bray Curtis of healthy follow up to 20-44
dt_hfy <- dt_bfy[name_hf,]
dt_hfy <- matrix(dt_hfy)
group_hfy <- rep("HF",nrow(dt_hfy))
bc_hfy <- data.frame(dt_hfy,group_hfy)
colnames(bc_hfy) <- c("BC","group")
#Bray Curtis of less healthy follow up to 20-44
dt_lhfy <- dt_bfy[name_lhf,]
dt_lhfy <- matrix(dt_lhfy)
group_lhfy <- rep("LHF",nrow(dt_lhfy))
bc_lhfy <- data.frame(dt_lhfy,group_lhfy)
colnames(bc_lhfy) <- c("BC","group")

#Combine by rows
df_status2 <- rbind(rbind(rbind(bc_hby,bc_lhby),bc_hfy),bc_lhfy)
#Plot
df_status2$group <- factor(df_status2$group,levels = c("HB","HF","LHB","LHF"))
mycompare <- list(c('HB','HF'),c('LHB','LHF'),c('HB','LHB'),c('HF','LHF'))
p_status2 <- ggplot(df_status2, aes(x = group, y = `BC`)) + 
  geom_violin(aes(fill = group),trim = T) + 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))+ 
  scale_fill_npg()+
  stat_compare_means(comparisons = mycompare,label = "p.format",method = 'wilcox',
                     paired = F,label.y = c(1.0,1.0,1.08,1.16))+ylim(0.55,1.22)+
  xlab("")+ylab("Bray-Curtis distance to 20-44")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size = 15), 
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.text.x = element_text(hjust = 0.5,vjust = 0.5,angle = 0),
        axis.text = element_text(colour = 'black',size = 15),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')
ggsave(plot = p_status2,'Fig5g_BC_distance_to_20-44_status.pdf',height = 4,width = 3.5,dpi = 600)
################################################################################

#Histogram, Fig5i
df_lv2 <- read.table("level-2_mean.txt", header = T, check.names = F,dec = ".", sep = "\t")
df_lv2$Phylum <- factor(df_lv2$Phylum,levels = rev(unique(df_lv2$Phylum)))
pdf("Fig5i_lv2-mean.pdf",height = 8,width = 7)
ggplot(melt(df_lv2), aes(x = variable, y = value, fill = Phylum)) + 
  geom_bar(stat = "identity", width = 0.5, col = 'black')  + 
  scale_fill_manual(values = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF",
                             "#20854EFF","#E18727FF","#0072B5FF","#BC3C29FF"))+
  labs(x = "",y = "Relative abundance(%)")+
  theme_bw()+ expand_limits(y = c(0,100))+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        panel.grid.minor = element_blank())+ 
  scale_y_continuous(expand = c(0,0))+
  geom_segment(data = df_lv2 %>% arrange(by = desc(Phylum)) %>% mutate(Baseline = cumsum(Baseline)) %>% mutate(`Follow-up` = cumsum(`Follow-up`)), aes(x = 1.25, xend = 1.75, y = Baseline, yend = `Follow-up`))
dev.off()
