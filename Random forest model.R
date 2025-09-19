################################################################################
#main codes for Figure 6 of gut microbiome and longevity studies
################################################################################
library(randomForest)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(caret)
library(pROC)
rm(list = ls())
set.seed(480)
feature <- read.table('level-6-20-44_66-85_100-117.txt',header = T,sep = "\t",check.names = T)
#The data is divided into the training set and the verification set by 3:1
inTrain <- createDataPartition(y = feature[,211],p = 0.25,list = F)
validation1 <- feature[inTrain,]
train1 <- feature[-inTrain,]
table(train1$Group)
table(validation1$Group)
write.table(validation1,'lv6-20-44_66-85_100-117-validation.txt',sep = '\t',quote = F)
write.table(train1,'lv6-20-44_66-85_100-117-train.txt',sep = '\t',quote = F)
train <- train1[,-1]
train$Group <- as.factor(train$Group)
#Optimize parameters ntree and mtry of RF model
n <- length(names(train))
set.seed(480)
min = 100
num = 0
for (i in 1:(n-1)){
  mtry_fit <- randomForest(Group~., data = train, mtry = i)
  err <- mean(mtry_fit$err.rate)
  print(err)
  if(err < min) {    
    min = err     
    num = i }
}
print(min)
print(num)
#[1] 74
ntree_fit <- randomForest(Group~.,data = train,mtry = 74,ntree = 1000)
pdf("ntree.pdf",height = 5,width = 5)
plot(ntree_fit)
dev.off()
write.table(ntree_fit$err.rate,'ntree.tsv',sep = '\t',quote = F)
#RF model
train.rf <- randomForest(Group ~ ., data = train, importance = TRUE, proximity = TRUE,ntree = 741,mtry = 74)
print(train.rf)

#Random Forest Cross-Valdidation for feature selection, the number of features is 29.
set.seed(290)
cross <- na.omit(train)
#mycross<-cbind(cross[1:210], matrix(runif(96 * nrow(cross)), nrow(cross), 96))
#Multiple cross verification
result <- replicate(5, rfcv(as.matrix(cross[,c(1:209)]), cross$Group,cv.fold = 10,step = 1.5), simplify = FALSE)
error.cv <- sapply(result, "[[", "error.cv")
write.table(error.cv,'error.cv.tsv',sep = '\t',quote = F)
pdf("CV_error.pdf",height = 5,width = 5)
matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type = "l",
        lwd = c(2, rep(1, ncol(error.cv))), col = 1, lty = 1, log = "x",
        xlab = "Number of variables", ylab = "CV Error")
dev.off()

#Variable Importance Top 29,which value of MeanDecreaseGini greater than 2, Fig6a
write.table(varImpPlot(train.rf),'varImpPlot_filter.tsv',sep = '\t',quote=F)
pdf("Fig6a_varImpPlot.pdf",height=7,width=10)
varImpPlot(train.rf, sort = TRUE, n.var = min(29, nrow(train.rf$importance)))
dev.off()

name_29 <- c("Escherichia.Shigella","Faecalibacterium","Clostridium.sensu.stricto.1","Parabacteroides",
           "Alistipes","Subdoligranulum","Eubacterium.coprostanoligenes.group","Fusicatenibacter",
           "Lachnospira","Erysipelotrichaceae.UCG.003","Desulfovibrio","Odoribacter","Ruminococcaceae.UCG.004",
           "Butyricimonas","Intestinibacter","Clostridium.innocuum.group","Eggerthella","Eisenbergiella",
           "Rhodococcus","Cloacibacillus","Staphylococcus","Acinetobacter","Cellulosilyticum","Epulopiscium",
           "Corynebacterium.1","Hydrogenoanaerobacterium","Pseudomonas","Sarcina","LogFB","Group")

#20-44 vs 100-117
set.seed(777) #20-44vs100-117
train_yc <- train[-c(223:507),][,name_29]
val_yc <- validation[-c(76:176),][,name_29]
train_yc$Group <- as.factor(train_yc$Group)
val_yc$Group <- as.factor(val_yc$Group)
#Cross-validate 10 fold
folds <- createFolds(y = train_yc[,1],k = 10)
fc1 <- as.numeric()
mod_pre1<-as.numeric()
for(i in 1:10){
  fold_test1 <- train_yc[folds[[i]],]
  fold_train1 <- train_yc[-folds[[i]],]
  model1 <- randomForest(Group ~ .,data = fold_train1,proximity = T,importance = T)
  model_pre1 <- predict(model1,newdata = fold_test1,type="prob")
  fc1 <- append(fc1,as.numeric(fold_test1$Group))
  mod_pre1 <- append(mod_pre1,model_pre1[,2])
}
df1 <- cbind(fc1,as.numeric(mod_pre1))
write.table(df1,"train_20-44vs100-117.tsv",sep = '\t',quote = F)
#ROC for train of 20-44vs100-117
pdf("Fig6b1_train_20-44vs100-117.pdf",height = 4,width = 4)
mycol <- c("slateblue","#4DBBD5FF","#000000")
a <- plot.roc(df1[,1],df1[,2],
              smooth = F,
              lwd = 2,
              ylim = c(0,1),
              xlim = c(1,0),
              legacy.axes = T,
              main = "Random Forest",
              ci = TRUE,
              col = mycol[2])
ciobj <- ci.se(a,specificities = seq(0, 1, 0.01))
plot(ciobj, type = "shape", col = "#4DBBD5FF")
legend.name <- c(paste("20-44vs100-117",sep = " "),paste("AUC","=",round(a$auc,3),sep = " "),
                 paste("95%CI =",round(a$ci[1],3),"-",round(a$ci[3],3),sep = " "))
legend("bottomright", 
       legend = legend.name,
       col = mycol[2:3],
       lwd = 2,
       bty = "n")
dev.off()
#Validation of 20-44vs100-117
fc2 <- as.numeric()
mod_pre2 <- as.numeric()
model_FR2 <- randomForest(Group ~ .,data = train_yc,proximity = T,importance =T )
model_FR_pre2 <- predict(model_FR2,newdata = val_yc,type = "prob")
fc2 <- append(fc2,as.numeric(val_yc$Group))
mod_pre2 <- append(mod_pre2,model_FR_pre2[,2])
df2 <- cbind(fc2,as.numeric(mod_pre2))
write.table(df2,"validation_20-44vs100-117.tsv",sep = '\t',quote = F)
#ROC for validation of 20-44vs100-117
pdf("Fig6c1_validation_20-44vs100-117.pdf",height = 4,width = 4)
b<- plot.roc(df2[,1],df2[,2],
             smooth = F,
             lwd = 2,
             ylim = c(0,1),
             xlim = c(1,0),
             legacy.axes = T,
             main = "Random Forest",
             ci = TRUE,
             col = mycol[2])
ciobj <- ci.se(b,specificities = seq(0, 1, 0.01))
plot(ciobj, type="shape", col = "#4DBBD5FF")
legend.name <- c(paste("66-85vs100-117",sep = " "),paste("AUC","=",round(b$auc,3),sep = " "),
                 paste("95%CI=",round(b$ci[1],3),"-",round(b$ci[3],3),sep = " "))
legend("bottomright", 
       legend = legend.name,
       col = mycol[2:3],
       lwd = 2,
       bty = "n")
dev.off()

################################################################################
#66-85vs100-117
set.seed(888)
train_oc <- train[-c(508:747),][,name_29]
val_oc <- validation[-c(177:250),][,name_29]
train_oc$Group <- as.factor(train_oc$Group)
val_oc$Group <- as.factor(val_oc$Group)
#cross-validate 10 fold
folds <- createFolds(y = train_oc[,1],k = 10)
fc3 <- as.numeric()
mod_pre3 <- as.numeric()
for(i in 1:10){
  fold_test3 <- train_oc[folds[[i]],]
  fold_train3 <- train_oc[-folds[[i]],]
  model3 <- randomForest(Group ~ .,data = fold_train3,proximity = T,importance = T)
  model_pre3 <- predict(model3,newdata = fold_test3,type = "prob")
  fc3 <- append(fc3,as.numeric(fold_test3$Group))
  mod_pre3 <- append(mod_pre3,model_pre3[,2])
}
df3 <- cbind(fc3,as.numeric(mod_pre3))
write.table(df3,"train_66-85vs100-117.tsv",sep = '\t',quote = F)
#ROC for train of 66-85vs100-117
pdf("Fig6b2_train_66-85vs100-117.pdf",height = 4,width = 4)
c <- plot.roc(df3[,1],df3[,2],
              smooth = F,
              lwd = 2,
              ylim = c(0,1),
              xlim = c(1,0),
              legacy.axes = T,
              main = "Random Forest",
              ci = TRUE,
              col = mycol[2])
ciobj <- ci.se(c,specificities=seq(0, 1, 0.01))
plot(ciobj, type = "shape", col="#4DBBD5FF")
legend.name <- c(paste("66-85vs100-117",sep = " "),paste("AUC","=",round(c$auc,3),sep = " "),
                 paste("95%CI =",round(c$ci[1],3),"-",round(c$ci[3],3),sep = " "))
legend("bottomright", 
       legend = legend.name,
       col = mycol[2:3],
       lwd = 2,
       bty = "n")
dev.off()
#Validation of 66-85vs100-117
fc4 <- as.numeric()
mod_pre4 <- as.numeric()
model_FR4 <- randomForest(Group ~ .,data = train_oc,proximity = T,importance = T)
model_FR_pre4 <- predict(model_FR4,newdata = val_oc,type = "prob")
fc4 <- append(fc4,as.numeric(val_oc$Group))
mod_pre4 <- append(mod_pre4,model_FR_pre4[,2])
df4 <- cbind(fc4,as.numeric(mod_pre4))
write.table(df4,"validation_66-85vs100-117.tsv",sep = '\t',quote = F)
#ROC for validation of 66-85vs100-117
pdf("Fig6c2_validation_66-85vs100-117.pdf",height = 4,width = 4)
d<- plot.roc(df4[,1],df4[,2],
             smooth = F,
             lwd = 2,
             ylim = c(0,1),
             xlim = c(1,0),
             legacy.axes = T,
             main = "Random Forest",
             ci = TRUE,
             col = mycol[2])
ciobj <- ci.se(d,specificities = seq(0, 1, 0.01))
plot(ciobj, type = "shape", col = "#4DBBD5FF")
legend.name <- c(paste("66-85vs100-117",sep=" "),paste("AUC","=",round(d$auc,3),sep = " "),
                 paste("95%CI =",round(d$ci[1],3),"-",round(d$ci[3],3),sep = " "))
legend("bottomright", 
       legend = legend.name,
       col = mycol[2:3],
       lwd = 2,
       bty = "n")
dev.off()

################################################################################

#Supplementary Fig9
df_29 <- feature[,name_29]
df_29[,1:28] <- df_29[,1:28]*100
df_29[,1:28] <- log2(df_29[,1:28]+1)
df_29$Group <- factor(df_29$Group,levels = c("20-85","100-117"))
color <- c("#00468BFF","#42B540FF","#925E9FFF")
my_comparisons <- list(c("20-85", "100-117"))

gc <- colnames(df_29)
plist <- list()
for (i in 1:length(gc)){
  feature_group <- df_29[,c(gc[i],"Group")]
  colnames(feature_group) <- c("abundance","group")
  pb29 <- ggboxplot(feature_group,
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
  pb29 <- pb29+theme(axis.line = element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb29 <- pb29+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 12))
  pb29 <- pb29+theme(axis.text.y = element_text(size = 12))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))
  pb29 <- pb29+theme(legend.position = "NA",plot.title = element_text(hjust = 0.5,face = "italic"))
  pb29 <- pb29+stat_compare_means(method = "wilcox",hide.ns = F,comparisons = my_comparisons,label = "p.format",paired = F)
  plist[[i]] <- pb29
} 
pall<-plot_grid(plist[[1]],plist[[3]],plist[[4]],plist[[5]],plist[[11]],
                plist[[12]],plist[[13]],plist[[14]],plist[[15]],plist[[16]],
                plist[[17]],plist[[18]],plist[[20]],plist[[21]],plist[[23]],
                plist[[24]],plist[[25]],plist[[26]],plist[[27]],plist[[28]],
                plist[[2]],plist[[6]],plist[[8]],plist[[9]],plist[[10]],
                plist[[19]],plist[[29]],plist[[7]],plist[[22]],ncol=5)
ggsave(plot = pall,'Supplementary Fig9_20-85vs100-117.pdf',height = 21,width = 12,dpi = 600)
