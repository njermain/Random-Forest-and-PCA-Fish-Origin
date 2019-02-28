library(shapeR)

setwd("C:/Users/w10007346/Dropbox/ATC shape analsis project (1)/ACM_ShapeAnalysis/Analysis")
load("partialrun18_Shapecoefs_minusItaly.RData")
head(getStdWavelet(shape))

# remove shape coefficients that covary significantly with length or year class
shape=stdCoefs(shape,classes="pop","length_cm", bonferroni=T) 
str(getMasterlist(shape)$AGE)
shape@master.list.org$AGE<-as.numeric(shape@master.list.org$age+1)
shape@master.list.org$YEARCLASS<-as.numeric(shape@master.list.org$year-shape@master.list.org$age)
getMasterlist(shape)$pop
which(getMasterlist(shape)$age>=0)
shape=stdCoefs(shape,classes="pop","YEARCLASS", bonferroni=T) 

# check for significant interactions for yearclass and age manually without function stdCoefs
head(str(getMasterlist(shape)))
coef.names<-colnames(getMasterlist(shape))[28:90]
# for age
ancovapvalues<-1
for (i in 28:90){
ancovapvalues[i]<-summary(aov(getMasterlist(shape)[,i]~getMasterlist(shape)$pop*getMasterlist(shape)$AGE))[[1]]["getMasterlist(shape)$pop:getMasterlist(shape)$AGE","Pr(>F)"]
}

length(which(ancovapvalues<.05))
rm.age.cols<-which(ancovapvalues<.05)  # which columns need to be removed for the analysis

# for year class
ancovapvalues.yc<-1
for (i in 28:90){
  ancovapvalues.yc[i]<-summary(aov(getMasterlist(shape)[,i]~getMasterlist(shape)$pop*getMasterlist(shape)$YEARCLASS))[[1]]["getMasterlist(shape)$pop:getMasterlist(shape)$YEARCLASS","Pr(>F)"]
}

length(which(ancovapvalues.yc<0.05))

rm.yc.cols<-which(ancovapvalues.yc<0.05)

# for length
ancovapvalues.length<-1
for (i in 28:90){
  ancovapvalues.length[i]<-summary(aov(getMasterlist(shape)[,i]~getMasterlist(shape)$pop*getMasterlist(shape)$length_cm))[[1]]["getMasterlist(shape)$pop:getMasterlist(shape)$length_cm","Pr(>F)"]
}

length(which(ancovapvalues.length<0.05))

rm.length.cols<-which(ancovapvalues.length<0.05)

# apply bonferonni corrections which is just adjust the critical value alpha to be 
# divided by the number of comparisons length(28:90), so .05/63
bonf.pvalue<-.05/length(28:90)

which(ancovapvalues.length<bonf.pvalue) # all good for length

which(ancovapvalues.yc<bonf.pvalue) # remove coefficient on column 66

which(ancovapvalues<bonf.pvalue)  # remove 29 38 40 51 52 75 76 77 78

# make stardardized coefficients list manually
library(dplyr)
select(getMasterlist(shape),Ws0c1)
usecols<-c(28,seq(30,37,1),39,seq(41,50,1),seq(53,74,1),seq(79,90,1))
standCfs<-getMasterlist(shape)[,usecols]

##########


# update master list
shape=enrich.master.list(shape)


# canonical analysis of princippal coordinates
library(vegan)
cap.res=capscale(standCfs ~ getMasterlist(shape)$pop, distance="euclidean")


anova(cap.res, by = "terms", step = 1000)

eig = eigenvals(cap.res,constrained = T)
eig.ratio = eig/sum(eig)



jpeg("PCoA_v2.jpg",width = 3000, height=2000,res = 400)
cluster.plot(scores(cap.res)$sites[,1:2],getMasterlist(shape)$pop, 
             
             xlim = c(-2,2),
             ylim = c(-2,2),
             xlab = paste("CAP1 (",round(eig[1]*100,1),"%)",sep = ""),
             ylab = paste("CAP2 (",round(eig[2]*100,1),"%)",sep = ""), 
             plotCI = F,conf.level = 0.95,
             las = 1, cex.lab=1.1, 
             cex.axis=1.1)

dev.off()



# Which Wavelet coefficients are most correlated with each axis
max(cor(scores(cap.res)$sites[,1],getStdWavelet(shape)))
max(cor(scores(cap.res)$sites[,2],getStdWavelet(shape)))

# simper analysis to determine which variables contribute most to significance of PERMANOVA
sim<-simper(getStdWavelet(shape), getMasterlist(shape)$pop)
sim$East_Gulf$overall # .1140949
sim$Gulf_West$overall # .1249572
sim$East_West$overall # .1291113

str(sim)
summary(sim)

# how many harmonics to include in ICC plot?
# Reconstruction of outlines; see how well they represent the real otoliths
est.list<-estimate.outline.reconstruction((shape))
save(est.list, file="estimate.outline.reconstruc.RData")
outline.reconstruction.plot(est.list,max.num.harmonics = 32, ref.w.level=0)


# Where on the otolith the variation is?

library(gplots)

jpeg("ICC.jpg",width = 3000, height=2000,res = 300)
plotWavelet(shape, level=5, class.name="pop", useStdcoef = T)
dev.off()

# plot overall shapes

jpeg("overlayshape.jpg",width = 3000, height=2100,res = 300)
plotWaveletShape(shape, "pop", show.angle=T)
dev.off()



summary(getMasterlist(shape)$pop)



# what do age compositions look like between east, west and gulf regions

dat<-read.csv("FISH.csv")

east<-subset(dat, dat$pop=="East")
west<-subset(dat, dat$pop=="West")
gulf<-subset(dat, dat$pop=="Gulf")
par(mfrow=c(1,3))
hist(east$age)
hist(west$age)
hist(gulf$age)


par(mfrow=c(1,3))
hist(east$month)
hist(west$month)
hist(gulf$month)



########### Random forest analysis ###########################

setwd("C:/Users/w10007346/Dropbox/ATC shape analsis project (1)/ACM_ShapeAnalysis/Analysis")
load("partialrun16_Shapecoefs_plusMauritania.RData")
head(getStdWavelet(shape))

library(randomForest)

class(getMasterlist(shape))
set.seed(100)
# divide into training and validation sets
train<-sample(nrow(getMasterlist(shape)),.7*nrow(getMasterlist(shape)), replace=F)
Trainset<-getMasterlist(shape)[train,]
Validset<-getMasterlist(shape)[-train,]
Trainset.coef<-standCfs[train,]
Validset.coef<-standCfs[-train,]
Train.df<-as.data.frame(Trainset.coef)
Train.df$pop<-Trainset$pop

Valid.df<-as.data.frame(Validset.coef)
Valid.df$pop<-Validset$pop


# random forest model, number of variables determined by for loop below
model1<-randomForest(pop~., data=Train.df, mtry=10, ntree=4000)
model1

predTrain<-predict(model1, Train.df, type="class")
table(predTrain, Train.df$pop)

# validation 
predValid<-predict(model1, Valid.df, type="class")
mean(predValid==Valid.df$pop) #90 percent classification accuracy
table(predValid, Valid.df$pop)


# evaluate the optimal number of variables to split at each node
a=c()

for (i in 3:length(standCfs)) {
  model3 <- randomForest(pop ~ ., data = Train.df, ntree = 4000, mtry = i, importance = TRUE)
  predValid <- predict(model3, Valid.df, type = "class")
  a[i-2] = mean(predValid == Valid.df$pop)
}

which(a==max(a))
a

plot(3:length(standCfs),a, type="l")



