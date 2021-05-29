library(geosphere)
library(ade4)
setwd("/bulk/zhenbin/soybean_domestication/wj")
pheno_geo<-read.csv("Geo_info_185_soja.csv",header=T)
names(pheno_geo)[1:5]<-c("Taxa","origin","groups","latitude","longitude")
pheno_geo<-pheno_geo[,1:5]
pheno_geo<-na.omit(pheno_geo)
ind_geo<-data.frame(pheno_geo[,c(5,4)])
ind_distance<-distm(ind_geo,fun = distGeo)
rownames(ind_distance)<-pheno_geo$Taxa
colnames(ind_distance)<-pheno_geo$Taxa
save(ind_distance,file="ind_distance.Rda")

# Genetic distance

library(SNPRelate)
vcf.fn <- "G.soja185.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "IBS.soja.gds", method="biallelic.only")
genofile <- snpgdsOpen("IBS.soja.gds")

ibs <- snpgdsIBS(genofile, num.thread=2)
ibs_ind<-ibs$ibs
row.names(ibs_ind)<-ibs$sample.id
colnames(ibs_ind)<-ibs$sample.id
ibs_ind<-(1-ibs_ind)/max(1-ibs_ind)
row.names(ibs_ind)<-gsub("NJAU.",'NJAU_',row.names(ibs_ind))
colnames(ibs_ind)<-gsub("NJAU.",'NJAU_',colnames(ibs_ind))

ibs_ind<-ibs_ind[which(row.names(ibs_ind)%in%row.names(ind_distance)),which(colnames(ibs_ind)%in%colnames(ind_distance))]
ind_distance<-ind_distance/1000
ind_distance[ind_distance<1]<-1
# ind_distance<-log10(ind_distance+1)
gen <- quasieuclid(as.dist(ibs_ind)) 
geo <- quasieuclid(as.dist(ind_distance))
r1 <- mantel.randtest(geo,log10(gen),nrepet = 9999)

# mantel test figure

df <- data.frame(geo<-c(geo),gen<-c(gen))
names(df)<-c("geo","gen")

## Use densCols() output to get density at each point
x <- densCols(geo,gen, colramp=colorRampPalette(c("blue","red")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(300)
df$col <- cols[df$dens]

png("mantel_geo_gen_scatter.png",width=8,height=4,units="in",type="cairo",res=600)
par(mfrow=c(1,2),mai=c(0.7,0.7,0.1,0.1),mgp=c(2,0.1,0),tck=-0.01,las=1)
plot(gen~geo, pch=20, col='gray90', cex=0.5,xlab="Geographic distance (Km)",ylab="Genetic distance",cex.lab=1.2)
abline(lm(gen~geo),lwd=2,col="red")

plot(r1,cex.lab=1.5,main="",ylim=c(0,3000),cex.axis=1.2,axes=F,cex.lab=1.2)
axis(1,at=seq(0,0.60,by=0.1),label=seq(0,0.60,by=0.1),las=1)
axis(2,at=seq(0,3000,by=500),label=seq(0,3000,by=500),las=1)
dev.off()
