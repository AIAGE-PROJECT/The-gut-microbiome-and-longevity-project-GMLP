################################################################################
#main codes for Figure 1 of gut microbiome and longevity studies
################################################################################
library(DirichletMultinomial)
library(cluster)
library(clusterSim)
library(vegan)
library(ade4)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsci)
library(scales)

rm(list = ls())
set.seed(8965)
# Load data
abundance <- read.table("level-6.txt", header=T,row.names = 1, dec=".", sep="\t",check.names =F )
#Delete metadata
df <- abundance[,-c(737:746)]
#Remove features with mean abundance less than 0.0001%
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  df_ab <- dataframe
  feature <- colMeans(df_ab) > percent 
  df_feature <- df_ab[,feature]
  print(percent)
  return(df_feature)
}
df <- noise.removal(df, percent = 0.0001)

#Fit the DMM model
count <- as.matrix(df)
fit <- lapply(1:10, dmn, count = count, verbose=TRUE)

#Judge the fitting effect and determine the optimal cluster number
laplace <- sapply(fit, laplace) 
write.table(laplace,'Sup_Fig3a_fit_DMM_laplace-0.0001.tsv',sep = '\t',quote=F)

#Supplementary Fig. 3a
pdf("Sup_Fig3a_fit_DMM_laplace-0.0001.pdf",height=4,width=5)
plot(laplace, type="b", xlab = "Number of Dirichlet Components", ylab="Laplace",lwd = 2)
dev.off()

##PAM clustering
#Based on bray-curtis distance
data.dist <- vegdist(df,method = "bray")
#Clustering algorithm
pam.clustering <- function(x,k) { # x is a distance matrix and k is the number of clusters
  require(cluster)
  cluster <- as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
data.cluster <- pam.clustering(data.dist, k = 4) #k=4 was used as the optimal cluster number for clustering
write.table(data.cluster,"cluster_pam-4.tsv",sep = '\t',quote=F)

#Used the function silhouette() for cluster validation
obs.silhouette <- mean(silhouette(data.cluster, data.dist)[,3])
cat(obs.silhouette)
#0.1181672

##Between-class analysis (BCA) was performed to support the clustering and identify the drivers for the enterotypes.
obs.pca <- dudi.pca(data.frame(df), scannf = F, nf = 4)
obs.bet <- bca(obs.pca, fac = as.factor(data.cluster), scannf = F, nf = 4) 
write.table(obs.bet$tab,"bca_table_drivers.txt",sep = '\t',quote=F) 

#Principal coordinates analysis 
obs.pcoa <- dudi.pco(data.dist, scannf = F, nf = 4)

#Plot with s.class(),Fig.1a
pdf("Fig1a_pcoa_bray_pam-4.pdf",height=5,width=5)
s.class(obs.pcoa$li, fac = as.factor(data.cluster), grid = F,
        col = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#EE4C97FF","#FFDC91FF"),
        #sub="Principal coordiante analysis",
        label = c("E1","E2","E3","E4"),
        xlim = c(-0.45,0.5),ylim = c(-0.5,0.5),
        clabel = 1, cpoint = 1,cstar = 1,cellipse = 1.5)
dev.off()

################################################################################

##Top drivers plot, Supplementary Fig. 3b-e
#Cluster1
cluster1 <- read.table("bca_table_drivers-cluster1.txt", header = T, row.names = 1, dec = ".", sep = "\t")
cluster1$Genus <- factor(cluster1$Genus,levels = rev(unique(cluster1$Genus)))
p1 <- ggplot(cluster1, aes(x = Genus, y = value1)) +
  geom_bar(stat = "identity",width = 0.8) +
  coord_flip() +theme_bw()+
  labs(title = paste("Top drivers: community type 1"))+
  xlab("Genus")+ylab("Value")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour='black', size = 15),
        axis.text.x = element_text(colour = 'black',size = 15),
        axis.text.y = element_text(colour = 'black',face="italic"),
        axis.text = element_text(colour = 'black',size = 15),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = p1,"Sup_Fig3b_bca_table_drivers-cluster1.pdf",height = 6,width =7 ,dpi = 600)

#Cluster2
cluster2 <- read.table("bca_table_drivers-cluster2.txt", header = T, row.names = 1, dec = ".", sep = "\t")
cluster2$Genus <- factor(cluster2$Genus,levels = rev(unique(cluster2$Genus)))
p2 <- ggplot(cluster2, aes(x = Genus, y = value2)) +
  geom_bar(stat = "identity",width = 0.8) +
  coord_flip() +theme_bw()+
  labs(title = paste("Top drivers: community type 2"))+
  xlab("Genus")+ylab("Value")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.text.x = element_text(colour='black',size = 15),
        axis.text.y = element_text(colour = 'black',face = "italic"),
        axis.text = element_text(colour = 'black',size = 15),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = p2,"Sup_Fig3c_bca_table_drivers-cluster2.pdf",height = 6,width = 7 ,dpi = 600)

#Cluster3
cluster3 <- read.table("bca_table_drivers-cluster3.txt", header = T, row.names = 1, dec = ".", sep = "\t")
cluster3$Genus <- factor(cluster3$Genus,levels = rev(unique(cluster3$Genus)))
p3 <- ggplot(cluster3, aes(x = Genus, y = value3)) +
  geom_bar(stat = "identity",width = 0.8) +
  coord_flip() +theme_bw()+
  labs(title = paste("Top drivers: community type 3"))+
  xlab("Genus")+ylab("Value")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour='black', size = 15),
        axis.text.x = element_text(colour = 'black',size = 15),
        axis.text.y = element_text(colour = 'black',face = "italic"),
        axis.text = element_text(colour = 'black',size = 15),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = p3,"Sup_Fig3d_bca_table_drivers-cluster3.pdf",height = 6,width = 7 ,dpi = 600)

#Cluster4
cluster4 <- read.table("bca_table_drivers-cluster4.txt", header = T, row.names = 1, dec = ".", sep = "\t")
cluster4$Genus <- factor(cluster4$Genus,levels = rev(unique(cluster4$Genus)))
p4 <- ggplot(cluster4, aes(x = Genus, y = value4)) +
  geom_bar(stat = "identity",width = 0.8) +
  coord_flip() +theme_bw()+
  labs(title = paste("Top drivers: community type 4"))+
  xlab("Genus")+ylab("Value")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.text.x = element_text(colour = 'black',size = 15),
        axis.text.y = element_text(colour = 'black',face = "italic"),
        axis.text = element_text(colour = 'black',size = 15),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = p4,"Sup_Fig3e_bca_table_drivers-cluster4.pdf",height = 6,width = 7 ,dpi = 600)

################################################################################

##Boxplot of the top contributors, Fig.1b
group <- matrix(paste('E',data.cluster,sep = ""))
colnames(group) <- "Enterotypes"
df_cbind <- cbind(abundance,group)
df_cbind$Enterotypes <- factor(df_cbind$Enterotypes,levels = c("E1","E2","E3","E4"))

pe1 <- ggplot(df_cbind, aes(x = Enterotypes, y = `Bacteroides`,fill = Enterotypes))+labs(title = "Bacteroides",x = "", y = "Relative abundance (%)")+
  stat_boxplot(geom = "errorbar",width = 0.3)+
  geom_boxplot(aes())+
  scale_fill_nejm()+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text.x = element_text(colour = 'black',size = 12,angle = 0,vjust = 0.5,hjust = 0.5,face = "plain"),
        axis.text = element_text(colour = 'black',size = 15),legend.title = element_blank(),
        legend.position = 'none',text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5,face = "italic"),panel.grid = element_blank())+
  stat_compare_means(comparisons=list(c('E1','E2'),c('E1','E3'),c('E1','E4')),
                     label = 'p.format',method = 'wilcox',label.y = c(85,100,115))+
  ylim(0,125)

pe2 <- ggplot(df_cbind, aes(x = Enterotypes, y = `Escherichia-Shigella`,fill = Enterotypes))+labs(title = "Escherichia-Shigella",x="", y="")+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot()+
  scale_fill_nejm()+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text = element_text(colour = 'black',size = 15),
        axis.text.x = element_text(colour = 'black',size = 12,angle = 0,vjust = 0.5,hjust = 0.5,face = "plain"),
        legend.title = element_blank(),legend.position = 'none',text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5,face = "italic"),panel.grid = element_blank())+
  stat_compare_means(comparisons=list(c('E2','E1'),c('E2','E3'),c('E2','E4')),
                     label = 'p.format',method = 'wilcox',label.y = c(95,110,125))+
  ylim(0,135)

pe3 <- ggplot(df_cbind, aes(x = Enterotypes, y = `Prevotella`,fill = Enterotypes))+labs(title = "Prevotella",x = "", y = "Relative abundance (%)")+
  stat_boxplot(geom = "errorbar",width = 0.3)+
  geom_boxplot()+
  scale_fill_nejm()+
  theme_set(theme_bw())+
  scale_x_discrete(labels = c("E1\nn=468","E2\nn=413","E3\nn=284","E4\nn=410"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text = element_text(colour = 'black',size = 15),
        axis.text.x = element_text(colour = 'black',size = 12,angle = 0,vjust = 0.5,hjust = 0.5,face = "plain"),
        legend.title = element_blank(),legend.position = 'none',text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5,face = "italic"),
        panel.grid = element_blank())+
  stat_compare_means(comparisons = list(c('E3','E1'),c('E2','E3'),c('E3','E4')),
                     label = 'p.format',method = 'wilcox',label.y = c(80,95,110))+
  ylim(0,120)

pe4 <- ggplot(df_cbind, aes(x = Enterotypes, y = `Blautia`,fill = Enterotypes))+labs(title = "Blautia",x="", y="")+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot()+
  scale_fill_nejm()+
  theme_set(theme_bw())+
  scale_x_discrete(labels = c("E1\nn=468","E2\nn=413","E3\nn=284","E4\nn=410"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text = element_text(colour ='black',size = 15),
        axis.text.x = element_text(colour = 'black',size = 12,angle = 0,vjust = 0.5,hjust = 0.5,face = "plain"),
        legend.title = element_blank(),legend.position = 'none',text = element_text(size=15),
        plot.title = element_text(hjust = 0.5,face = "italic"),panel.grid = element_blank())+
  stat_compare_means(comparisons = list(c('E4','E1'),c('E2','E4'),c('E3','E4')),
                     label = 'p.format',method = 'wilcox',label.y = c(95,110,125))+
  ylim(0,135)

pe1_4 <- pe1+pe2+pe3+pe4+plot_layout(ncol = 2,nrow = 2)
ggsave(plot = pe1_4,'Fig1b_E1-4-top-contributors.pdf',height = 6,width =6 ,dpi = 600)

################################################################################

##Bar chart of enterotypes,Fig.1c
#Fisher's test
fisher.test(table(df_cbind$Enterotypes,df_cbind$Group),simulate.p.value = T)
#Fisher's Exact Test for Count Data with simulated p-value (based on 2000 replicates)

#data:  table(df_cbind$Enterotypes, df_cbind$Group)
#p-value = 0.0004998
#alternative hypothesis: two.sided

enter_df <- data.frame(table(df_cbind$Enterotypes,df_cbind$Group))
colnames(enter_df) <- c('Enterotypes',"Group","value")
enter_df$Group <- factor(enter_df$Group,levels = c("20-44","45-65","66-85","90-99","100-117"))
p_bar <- ggplot(enter_df,mapping = aes(Group,value,fill = Enterotypes))+
  geom_bar(stat = 'identity',position = 'fill') +
  scale_fill_nejm()+
  labs(title = "P = 0.0005",x = '',y = 'Percentage') +
  theme(axis.title = element_text(size = 15),axis.text = element_text(size = 15, color = 'black'),text = element_text(size = 15))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust = 0.7,size = 18),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  theme(panel.grid = element_blank())
ggsave("Fig1c_enter-bar-pvalue.pdf", p_bar, width = 5, height = 6)

################################################################################

##Bar chart of enterotypes by health stats, Fig.1d
#Healthy
df_66_117H <- df_cbind[df_cbind$Status=='H',]
#Fisher test
fisher.test(table(df_66_117H$Enterotypes,df_66_117H$Group),workspace = 2e+7)
#Fisher's Exact Test for Count Data

#data:  table(df_66_117H$Enterotypes, df_66_117H$Group)
#p-value = 0.2471
#alternative hypothesis: two.sided

enter_66_100H <- data.frame(table(df_66_117H$Enterotypes,df_66_117H$Group))
colnames(enter_66_100H) <- c('Enterotypes',"Group","value")
enter_66_100H$Group <- factor(enter_66_100H$Group,levels = c("66-85","90-99","100-117"))
p_H <- ggplot(enter_66_100H,mapping = aes(Group,value,fill = Enterotypes))+
  geom_bar(stat = 'identity',position = 'fill') +
  scale_fill_nejm()+
  labs(title = "p-value = 0.25",x = '',y = 'Percentage') +
  theme(axis.title = element_text(size = 25),axis.text = element_text(size = 25, color = 'black'),text = element_text(size= 25 ))+
  theme(axis.text.x = element_text(angle = 0 , hjust = 0.5,vjust = 0.7),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  theme(panel.grid = element_blank())

#Less health
df_66_117LH <- df_cbind[df_cbind$Status=='LH',]
#Fisher test
fisher.test(table(df_66_117LH$Enterotypes,df_66_117LH$Group),workspace = 2e+7)
#Fisher's Exact Test for Count Data

#data:  table(df_66_117LH$Enterotypes, df_66_117LH$Group)
#p-value = 0.0003579
#alternative hypothesis: two.sided

enter_66_100LH <- data.frame(table(df_66_117LH$Enterotypes,df_66_117LH$Group))
colnames(enter_66_100LH) <- c('Enterotypes',"Group","value")
enter_66_100LH$Group <- factor(enter_66_100LH$Group,levels = c("66-85","90-99","100-117"))
p_LH <- ggplot(enter_66_100LH,mapping = aes(Group,value,fill = Enterotypes))+
  geom_bar(stat = 'identity',position = 'fill') +
  scale_fill_nejm()+
  labs(title = "p-value = 0.0004",x = '',y = 'Percentage') +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25, color = 'black'),
        text = element_text(size = 25))+
  theme(axis.text.y = element_blank(),   
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  theme(axis.text.x = element_text(angle = 0 , hjust = 0.5,vjust = 0.7),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  theme(panel.grid=element_blank())

#Frail
df_100F <- df_cbind[df_cbind$Status=='F',]
df_100F <- df_100F[df_100F$Group=='100-117',]

#fisher test
#fisher.test(table(df_66_117LH$Enterotypes,df_66_117LH$group),workspace = 2e+7)

enter_100F <- data.frame(table(df_100F$Enterotypes,df_100F$Group))
colnames(enter_100F) <- c('Enterotypes',"Group","value")
p_F <- ggplot(enter_100F,mapping = aes(Group,value,fill = Enterotypes))+
  geom_bar(stat = 'identity',position = 'fill') +
  scale_fill_nejm()+
  labs(title = "Frail",x = '',y = 'Percentage') +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25, color = 'black'),
        text = element_text(size=25))+
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  theme(axis.text.x = element_text(angle =0 , hjust = 0.5,vjust = 0.7),
        plot.title = element_text(hjust = 0.5))+
  theme(panel.grid = element_blank())

#Combine multiple plots
p_status <- p_H+ p_LH +p_F+ plot_layout(guides = 'collect',widths = c(3,3,1))
ggsave("Fig1d_enter-bar-66-117-pvalue-status.pdf", p_status, width =11, height = 6)

################################################################################

##Multivariate regression analysis of enterotypes
df_glm <- data.frame(df_cbind$Enterotypes,df_cbind$age,df_cbind$sex,df_cbind$BMI,
                     row.names = rownames(df_cbind))
colnames(df_glm) <- c("Enterotypes","Age","Sex","BMI")

#Transform classification variables into factor types
for (i in c("Enterotypes","Sex")){
  df_glm[,i] = as.factor(df_glm[,i]) }

str(df_glm)

#Glm
glm_enter <- glm(Enterotypes~Age+Sex+BMI,data = df_glm,family = binomial(link = 'logit'))
glmSum <- summary(glm_enter)

#Save the result, Supplementary Table 1
ResultMul <- c()
ResultMul <- rbind(ResultMul,glmSum$coef)
OR <- exp(glmSum$coef[,'Estimate'])
ResultMul <- cbind(ResultMul,confint(glm_enter),cbind(OR,exp(confint(glm_enter))))
write.csv(ResultMul,file = "gml_Mul_age_sex_BMI.csv")
