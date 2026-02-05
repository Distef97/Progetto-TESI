setwd("~/GitHub/Tesis/PROF")
Rcpp::sourceCpp("rcpp1.cpp")
source("functions.R")
source("myfunctions.R")
set.seed(123)
library(knitr) #Per latex


####1
model<-gen_varma_noid(seed = 123,maxl = 5, part = "varma2")
model<-gen_varma_noid(seed = 123,maxl = 5, part = "varma3")
#model<-gen_varma_noid(seed = 123,maxl = 5, part = "var")
set.seed(123)
#model<-gen_varma_noid_grade(2, maxl = 5)
#####################à

mod1<-model$mod1
maxl=model$maxl
w<-model$w
s<-model$s
l<-model$l
irf1<-model$irf1
mod1
model$identified
saveRDS(mod1, file = "modellovarmainiziale.RDS")


n<-c(50, 100, 500, 1000, 2000, 3000,4000, 5000)
maxn<-5
n<-n[2:maxn]
simulation=1000*maxn

diff_irf<-c()
biases<-c()
VARIANZA<-c()
irf2<-list()
risultati_varma1 <- vector("list", simulation) # Pre-allocazione per efficienza

############################################################
########################FABLE
########################################################
for (k in 1:length(n[1:maxn])) {
  print(k)
  resih<-simulvarmanoidfable(NV = n[k],simulazioni = 1000)
  diff_irf[k]<-resih$diff_irf
  irf1c<-model$irf1[,,-1]
  resih$irf2<-resih$irf2[-(1:(w*w)),]
  irf2[[k]]<-resih$irf2
  risultati_varma1[(1+1000*(k-1)):(1000+1000*(k-1))]<-resih$risultati_varma
  bias<-matrix(0, length(irf1c), resih$simulazioni)
  for (i in 1:resih$simulazioni) {
    bias[,i]<-abs(irf1c-resih$irf2[,i])
  }
  biases[k]<-mean(rowMeans(bias))
  variance<-c()
  for (i in 1:resih$simulazioni) {
    variance[i]<-var(resih$irf2[,i])
  }
  VARIANZA[k]<-mean(variance)
}



#diff_irf
#biases
#VARIANZA
resvarma<-cbind(diff_irf, biases, VARIANZA)
saveRDS(resvarma, file = "tabellaIH.RDS")
kable(risultati_varma1[[c(1,2)]]$ar, format = "latex")
saveRDS(risultati_varma1, file = "risultatiIH.RDS")
saveRDS(irf2, file = "irf2IH.RDS")
risultati_varma1[[1000]]
dim(irf2[[1]])
varma1<-rvarma(3,1,1)
predict.varma(varma1, 10)
predict.varma(mod1, 10)

#####################à
############# KFAS
##########################



for (k in 1:length(n[1:4])) {
  print(k)
  resih<-simulvarmanoidkfas(NV = n[k],simulazioni = 1000)
  diff_irf[k]<-resih$diff_irf
  irf1c<-model$irf1[,,-1]
  resih$irf2<-resih$irf2[-(1:(w*w)),]
  irf2[[k]]<-resih$irf2
  risultati_varma1[(1+1000*(k-1)):(1000+1000*(k-1))]<-resih$risultati_varma
  
  bias<-matrix(0, length(irf1c), dim(resih$irf2)[2])
  for (i in 1:dim(resih$irf2)[2]) {
    bias[,i]<-abs(irf1c-resih$irf2[,i])
  }
  biases[k]<-mean(rowMeans(bias))
  variance<-c()
  for (i in 1:dim(resih$irf2)[2]) {
    variance[i]<-var(resih$irf2[,i])
  }
  VARIANZA[k]<-mean(variance)
}



#diff_irf
#biases
#VARIANZA
resvarma<-cbind(diff_irf, biases, VARIANZA)
saveRDS(resvarma, file = "tabellaKFAS.RDS")
kable(risultati_varma1[[c(1,2)]]$ar, format = "latex")
saveRDS(risultati_varma1, file = "risultatiKFAS.RDS")
saveRDS(irf2, file = "irf2KFAS.RDS")
risultati_varma1[[1000]]
dim(irf2[[1]])
varma1<-rvarma(3,1,1)
predict.varma(varma1, 10)
predict.varma(mod1, 10)

  }
  biases[k]<-mean(rowMeans(bias))
  variance<-c()
  for (i in 1:dim(resih$irf2)[2]) {
    variance[i]<-var(resih$irf2[,i])
  }
  VARIANZA[k]<-mean(variance)
}
