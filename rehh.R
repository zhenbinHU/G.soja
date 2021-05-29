
vcftools<-"vcftools"

# subset the variants using vcftools
cmd<-paste0(vcftools, " --vcf G.soja185.vcf --chr 19 --from-bp 2000000 --to-bp 7000000 --maf 0.01 --recode --recode-INFO-all --out ","chr19")

system(cmd)

library(rehh)
hh <- data2haplohh(hap_file = "chr19.recode.vcf", polarize_vcf = FALSE,vcf_reader = "data.table",min_perc_geno.mrk = 50)
res <- calc_ehh(hh,mrk = "S19_04547654",include_nhaplo = TRUE,phased = T)

png("chr19.png",width=3,height=2,res=600,unit="in",type="cairo")
par(mai=c(0.3,0.3,0.1,0.1),mgp=c(1,0.1,0),las=1,cex=0.8)
plot(res,main="",tck=-0.01,cex.lab=0.8,cex.axis=0.8,xlim=c(4500000,4600000))
dev.off()
