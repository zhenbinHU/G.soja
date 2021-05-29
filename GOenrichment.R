library(topGO)

gene_ann<-read.delim("Gmax_275_Wm82.a2.v1.annotation_info.txt",header=T)
gene_go<-gene_ann[,c("locusName","GO")]
gene_go<-gene_go[grep("GO",gene_go$GO),]

GTOGO<-NULL
for (i in 1:dim(gene_go)[1]){
    tmp1<-gene_ann$locusName[i]
    tmp2<-as.character(gene_go$GO[i])
    go<-strsplit(tmp2,",")[[1]]
    gene_tmp<-rep(gene_go$locusName[i],length(go))
    gene_tmp<-data.frame(gene_tmp,go)
    GTOGO<-rbind(GTOGO,gene_tmp)
    if(i%%1000==0){ 
        cat(i," were done!\n")
    }
}

GTOGO<-as.data.frame(GTOGO)
row.names(GTOGO)<-NULL
names(GTOGO)<-c("gene_id","go_id")

geneID2GO <- by(GTOGO$go_id, GTOGO$gene_id, function(x) as.character(x))

all.genes <- sort(unique(as.character(GTOGO$gene_id)))

int.genes<-unique(gene_set_lfmm) # change the gene_set for each individual analyze
int.genes<-int.genes[which(int.genes%in%all.genes)]
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes

go.obj <- new("topGOdata", ontology='BP'
                 ,allGenes = int.genes
                 ,annot = annFUN.gene2GO
                 ,gene2GO = geneID2GO
                 )

results <- runTest(go.obj, algorithm = "weight01", statistic = "fisher")
results.tab <- GenTable(object = go.obj, elimFisher = results,topNodes = 60)

# GO graph

pdf("results_go_anno_gene.lfmm.pdf",width=6,height=6)
showSigOfNodes(go.obj, score(results), firstSigNode=10, useInfo ='all')
dev.off()

# figure the enrichment analysis
    
goEnrichment <- results.tab[as.numeric(results.tab$elimFisher)<0.05,]
    
goEnrichment <- goEnrichment[,c("GO.ID","Term","elimFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
goEnrichment$KS <- as.numeric(goEnrichment$elimFisher)

require(ggplot2)
goEnrichment$elimFisher<-as.numeric(goEnrichment$elimFisher)
p<-ggplot(goEnrichment, aes(x=Term, y=-log10(elimFisher))) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("Biological process") +
    ylab("Enrichment") +
    ggtitle("Title") +
    scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$elimFisher)), by = 2), 1)) +
    theme_bw(base_size=24) +
    theme(
        legend.position='none',
        legend.background=element_rect(),
        plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
        axis.text.x=element_text(angle=0, size=9, face="bold", hjust=1.10),
        axis.text.y=element_text(angle=0, size=9, face="bold", vjust=0.5),
        axis.title=element_text(size=12, face="bold"),
        legend.key=element_blank(),     #removes the border
        legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
        legend.text=element_text(size=18),  #Text size
        title=element_text(size=18)) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
png("goenrichment.gene.LFMM.png",width=6,height=6,type="cairo",res=600,unit="in")
p
dev.off()
write.csv(goEnrichment,'LFMM_goEnrichment.csv',row.names=F)
