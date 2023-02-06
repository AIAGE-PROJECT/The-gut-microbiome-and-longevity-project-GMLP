################################################################################
#main codes for Figures 2,3 of gut microbiome and longevity studies
################################################################################
library(ggpubr)
library(ggplot2)
library(magrittr)
library(patchwork)
library(ggsci)
library(plyr)
library(scales)

rm(list = ls())
df<-read.table('alpha-diversity_phylum.txt',header = T,row.names = 1,sep = '\t',check.names=F)

#Differential analysis of alpha diversity and relative abundance, Fig.2a-c and Fig.3c-e
data<-df[,-c(9:11,13:18)]
data$Group<- factor(data$Group,levels = c("20-44","45-65","66-85","90-99","100-117"))
mycompare=list(c('45-65','66-85'),c('66-85','90-99'),c('20-44','66-85'),c('66-85','100-117'),c('20-44','100-117'))
mycolor<- c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF")

gene = list()
for (i in 1:(ncol(data)-1)) {
  gene <-  as.character(colnames(data)[i])
  p <- ggviolin(data,x = "Group", y = gene, fill = "Group",palette = mycolor, 
                font.label = list(size = 24, face = "bold", color ="red"),xlab = "",trim = T,
                add = "boxplot", add.params = list(fill = "white"),title = gene)
  p <- p+stat_compare_means(method = "wilcox",label = "p.format",comparisons = mycompare)+
    theme(axis.title.x=element_text(colour='black', size=15),
          axis.title.y=element_blank(),
          axis.text.x=element_text(colour='black',size=15,hjust = 0.5,vjust = 0.5,angle = 30),
          plot.title = element_text(hjust = 0.5,size=15),
          axis.text=element_text(colour='black',size=15),
          legend.position='none')
  ggsave(filename = paste(gene, '.pdf', sep = ''),width = 3.5, height = 4.5, device = "pdf")
}

################################################################################

#Differential analysis of alpha diversity and relative abundance by health status, Fig.2d and Fig.3f-g
set.seed(5579)
#100-117
data_cen<-df[df$Group=='100-117',]
F100_117<-data_cen[data_cen$sub_status =='F100-117',]
H100_117<-data_cen[data_cen$sub_status=='H100-117',]
LH100_117<-data_cen[data_cen$sub_status=='LH100-117',]

#Randomly selected 60 samples in LH of 100-117
index_select<-sample(nrow(LH100_117),60,replace=F)
LH100_117<-LH100_117[index_select,]
merge100_117<-merge(merge(F100_117,H100_117,all = TRUE),LH100_117,all = TRUE)

#66-85
data_old<-df[df$Group=='66-85',]
H66_85<-data_old[data_old$sub_status=='H66-85',]
LH66_85<-data_old[data_old$sub_status=='LH66-85',]

#Randomly selected 60 samples in H of 66-85
index_select<-sample(nrow(H66_85),60,replace=F)
H66_85<-H66_85[index_select,]
merge66_85<-merge(LH66_85,H66_85,all = TRUE)

#90-99
data_90<-df[df$Group=='90-99',]
H90_99<-data_90[data_90$sub_status=='H90-99',]
LH90_99<-data_90[data_90$sub_status=='LH90-99',]

#Randomly selected 60 samples in H of 90-99
index_select<-sample(nrow(H90_99),60,replace=F)
H90_99<-H90_99[index_select,]

#Randomly selected 60 samples in LH of 90-99
index_select<-sample(nrow(LH90_99),60,replace=F)
LH90_99<-LH90_99[index_select,]
merge90_99<-merge(H90_99,LH90_99,all = TRUE)

df_merge<-merge(merge(merge90_99,merge66_85,all = TRUE),merge100_117,all = TRUE)
df_merge<-df_merge[,-c(8:13,15:18)]

mycompare=list(c('H66-85','LH66-85'),c('H90-99','LH90-99'),c('H100-117','LH100-117'),c('LH100-117','F100-117'),
               c('H100-117','F100-117'),c('H66-85','H90-99'),c('H90-99','H100-117'),c('H66-85','H100-117'),
               c('LH66-85','LH90-99'),c('LH90-99','LH100-117'),c('LH66-85','LH100-117'))
df_merge$sub_status<- factor(df_merge$sub_status,levels = c("H66-85","LH66-85","H90-99","LH90-99",
                                                            "H100-117","LH100-117","F100-117"))
mycolor<- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")

#Boxplot
gc <- colnames(df_merge)
name<-list()
for (i in 1:(length(gc)-1)){
  feature_group<-df_merge[,c(gc[i],"sub_status")]
  colnames(feature_group)<-c("abundance","group")
  name <-  as.character(colnames(df_merge)[i])
  pb1<-ggboxplot(feature_group,
                 x="group",
                 y="abundance",
                 color="group",
                 fill=NULL,
                 add = "jitter",alpha=0.01,shape =1,
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size=0.6,
                 #font.label = list(size=15), 
                 palette = mycolor)+
    theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,hjust = 0.3,vjust = 0.5,angle = 30))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA",plot.title = element_text(hjust = 0.5,face = "italic"))
  pb1<-pb1+stat_compare_means(method="wilcox",comparisons =mycompare,label="p.format")
  ggsave(filename = paste(name, '-health-status.pdf', sep = ''),width = 4, height = 5, device = "pdf")
}

################################################################################

#Correlation analysis between pielou's and phylum, Fig.2e-f
df_corr <- read.table("level-2-7_100-117.txt",header = T,row.names = 1,sep = '\t',check.names=F)
p1 <- ggplot(data=df_corr, aes(x=`pielou_e`, y=`p_Bacteroidetes`))+
  geom_point(color="#E64B35FF",size =5,alpha=0.9,shape = 1,stroke = 1)+ 
  stat_smooth(formula = y ~ x,method="glm",se=TRUE,level = 0.95,color = "#7E6148FF")+
  stat_cor(data=df_corr, method = "spearman",label.x.npc = 0.02,label.y.npc = 0.85,size = 6,cor.coef.name = "R")+
  xlab("Pielou's")+ylab("Bacteroidetes relative abundance")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill='white', colour='black'),
        text = element_text(size = 20),
        axis.title.x=element_text(colour='black', size=20),
        axis.title.y=element_text(colour='black', size=20),
        axis.text=element_text(colour='black',size=20),
        legend.title=element_blank(),
        panel.grid=element_blank())
ggsave(plot = p1,'corr_Bacteroidetes_pielou.pdf',height = 5,width = 5,dpi = 600)

p2 <- ggplot(data=df_corr, aes(x=`pielou_e`, y=`p_Proteobacteria`))+
  geom_point(color="#E64B35FF",size =5,alpha=0.9,shape = 1,stroke = 1)+ 
  stat_smooth(formula = y ~ x,method="glm",se=TRUE,level = 0.95,color = "#7E6148FF")+
  stat_cor(data=df_corr, method = "spearman",label.x.npc = 0.02,label.y.npc = 0.2,size = 6,cor.coef.name = "R")+
  xlab("Pielou's")+ylab("Proteobacteria relative abundance")+
  theme_set(theme_bw())+
  theme(panel.background = element_rect(fill='white', colour='black'),
        text = element_text(size = 20),   
        axis.title.x=element_text(colour='black', size=20),
        axis.title.y=element_text(colour='black', size=20),
        axis.text=element_text(colour='black',size=20),
        legend.title=element_blank(),
        panel.grid=element_blank())
ggsave(plot = p2,'corr_Proteobacteria_pielou.pdf',height = 5,width = 5,dpi = 600)

#Associations between gut microbiome taxa and ¦Á-diversity, Fig.2g and Supplementary Table 2
#Transform classification variables into factor types
for (i in c("sex","hypertension")){
  df_corr[,i] = as.factor(df_corr[,i]) }

#Glm
for(n in c("chao1","pielou_e","shannon")) {
  glm_model<- function(x){
      FML<-as.formula(paste0(n,"~",x,"+sex+BMI+hypertension"))
      glm1<-glm(FML,data=df_corr,family = gaussian)
      glm2<-summary(glm1)
      OR<-round(exp(coef(glm1)),4)
      SE<-glm2$coefficients[,2]
      CI5<-round(exp(coef(glm1)-1.96*SE),4)
      CI95<-round(exp(coef(glm1)+1.96*SE),4)
      CI<-paste0(CI5,'-',CI95)
      P<-round(glm2$coefficients[,4],48)
      B<-round(glm2$coefficients[,1],48)
      
      glm_model <- data.frame('Characteristics'=c("Intercept",x,"sexmale","BMI","hypertensionyes"),
                              'OR' = OR,
                              'CI' = CI,
                              'Beta' = B,
                              'P' = P,
                              'SE' = SE)[-1,]
      return(glm_model)
    }
  
  #Selection of microbial features
  variable.names<- colnames(df_corr)[c(4:419)]
  Uni_glm <- lapply(variable.names, glm_model)
  Uni_glm<- ldply(Uni_glm,data.frame)
  write.table(Uni_glm,paste0(n,"_taxa_glm_Sex_BMI_Hypertension.txt"),sep = '\t',quote = FALSE,col.names = NA)
}

#Plot Fig2g
dt_phy<-read.table("glm_beta&p_phylum_adjust_p.txt",header = T,sep='\t',comment.char='',check.names=F)
dt_phy$Phylum<-factor(dt_phy$Phylum,levels = c("Fusobacteria","Actinobacteria","Synergistetes","Verrucomicrobia",
                                           "Proteobacteria","Bacteroidetes","Firmicutes"))
phy1<-ggplot(dt_phy,aes(x=chao1,y=Phylum))+
  geom_point(aes(color=`chao1_p`),size=6)+
  geom_segment(aes(x=0,xend=chao1,y=Phylum,yend=Phylum,color=`chao1_p`),size=2,linetype="solid")+
  scale_colour_manual(values=c("#0072B5FF","#BC3C29FF","#B09C85FF"))+
  labs(title = 'Chao1 index',y='Phylum',x='¦Â-coef.')+
  geom_vline(aes(xintercept = 0),linetype="dotted",lwd=1,col="black",alpha=0.5)+
  theme(legend.position = 'none')+
  guides(size=guide_legend(ncol = 4,order = 0,
                           label.position = 'bottom'),
         color=guide_legend(ncol = 2,order=1))+
  theme(panel.background = element_rect(fill='white', colour='black'),
        axis.title.x=element_text(colour='black', size=25),
        axis.title.y=element_text(colour='black', size=25),
        axis.text=element_text(colour='black',size=25),
        plot.title = element_text(hjust = 0.5,size = 25),panel.grid = element_blank()
  )

phy2<-ggplot(dt_phy,aes(x=pielou_e,y=Phylum))+
  geom_point(aes(color=`pielou_e_p`),size=6)+
  geom_segment(aes(x=0,xend=pielou_e,y=Phylum,yend=Phylum,color=`pielou_e_p`),size=2,linetype="solid")+
  scale_colour_manual(values=c("#6F99ADFF","#E18727FF","#BC3C29FF","#B09C85FF"))+
  labs(title = "Pielou's index",y='',x='¦Â-coef.')+
  geom_vline(aes(xintercept = 0),linetype="dotted",lwd=1,col="black",alpha=0.5)+
  theme(legend.position = 'none')+
  guides(size=guide_legend(ncol = 4,order = 0,
                           label.position = 'bottom'),
         color=guide_legend(ncol = 2,order=1))+
  theme(panel.background = element_rect(fill='white', colour='black'),
        axis.title.x=element_text(colour='black', size=25),
        axis.title.y=element_text(colour='black', size=25),
        axis.text.y = element_blank(),
        axis.text=element_text(colour='black',size=25),
        plot.title = element_text(hjust = 0.5,size = 25),panel.grid = element_blank()
  )

phy3<-ggplot(dt_phy,aes(x=shannon,y=Phylum))+
  geom_point(aes(color=`shannon_p`),size=6)+
  geom_segment(aes(x=0,xend=shannon,y=Phylum,yend=Phylum,color=`shannon_p`),size=2,linetype="solid")+
  scale_colour_manual(values=c("#BC3C29FF","#B09C85FF"))+
  labs(title = 'Shannon index',y='',x='¦Â-coef.')+
  geom_vline(aes(xintercept = 0),linetype="dotted",lwd=1,col="black",alpha=0.5)+
  theme(legend.position = 'none')+
  theme(panel.background = element_rect(fill='white', colour='black'),
        axis.title.x=element_text(colour='black', size=25),
        axis.title.y=element_text(colour='black', size=25),
        axis.text.y = element_blank(),
        axis.text=element_text(colour='black',size=25),
        plot.title = element_text(hjust = 0.5,size = 25),panel.grid = element_blank()
  )
phy=phy1+phy2+phy3+plot_layout(ncol = 3,nrow = 1)
ggsave(plot = phy,'Fig2g_corr_glm_beta_p-adjust.pdf.pdf',height = 6,width = 14,dpi = 600)

