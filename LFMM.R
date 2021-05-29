library(lfmm)

load("gen.imp.rda")
e<-read.csv("e_data.csv")

setwd("/bulk/zhenbin/soybean_domestication/wj/G_E")
library(lfmm)
e<-e[,c("Latitude","Longitude","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","values")]
names(e)<-c("Latitude","Longitude","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","Altitude")

for(i in names(e)){
  mod.lfmm <- lfmm_ridge(Y = gen.imp, 
                        X = e$Latitude, 
                        K = 5)
                        
  pv <- lfmm_test(Y = gen.imp, 
                 X = e$Latitude, 
                 lfmm = mod.lfmm, 
                 calibrate = "gif")
  save(pv,file=paste0("pv.",i,".rda"))
}