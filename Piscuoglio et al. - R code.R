########## SUPPLEMENTARY MATERIAL #########
### Piscuoglio et al
### Genomic and transcriptomic heterogeneity of metaplastic breast carcinomas
###########################################


source("BACE.eset.97_CN.R")

library(RColorBrewer)
library(Biobase)
library(gplots)

########################## data preparation #####################
## expression
## read all the data in, annotate the probes
## the final object is metaplastic.mapped.lumi

metaplastic.lumi <- readBeads("20131101_beadstudio.txt", ann.library=NULL, flags.on=F, pval.thresh=0.01)

ensembl.61.ILMN <- read.csv("ensembl.61.ILMN.entrez.csv", sep=",", stringsAsFactors=F, na.strings=c("", " ", "NA", "#N/A"))
names(ensembl.61.ILMN) <- c("ensg", "description", "chrom", "start", "end", "cytoband", "status", "symbol", "ILMN", "entrez")
BS.file <- read.table("20131101_beadstudio.txt", sep="\t", as.is=T, header=T)
fData(metaplastic.lumi)$ILMN = BS.file$ILMN

metaplastic.probes <- ensembl.61.ILMN[match(fData(metaplastic.lumi)$ILMN, ensembl.61.ILMN$ILMN, nomatch=NA),]
metaplastic.probes.entrez <- ensembl.61.ILMN[match(fData(metaplastic.lumi)$entrez, ensembl.61.ILMN$entrez, nomatch=NA),]
metaplastic.probes[is.na(metaplastic.probes$ensg),] <- metaplastic.probes.entrez[is.na(metaplastic.probes$ensg),]

metaplastic.probes$chrom[metaplastic.probes$chrom == "X"] <- 23
metaplastic.probes$chrom[metaplastic.probes$chrom == "Y"] <- 24
metaplastic.probes$chrom <- as.integer(metaplastic.probes$chrom)

fvarLabels(metaplastic.lumi)[1] <- "PROBE_ID"
metaplastic.probes$description <- gsub("\\[\\S+\\ \\S+", "", metaplastic.probes$description)
metaplastic.probes$description <- gsub("\\ $", "", metaplastic.probes$description)

fData(metaplastic.lumi)$ILMN <- metaplastic.probes$ILMN
fData(metaplastic.lumi)$ensg <- metaplastic.probes$ensg
fData(metaplastic.lumi)$symbol <- metaplastic.probes$symbol
fData(metaplastic.lumi)$description <- metaplastic.probes$description
fData(metaplastic.lumi)$chrom <- metaplastic.probes$chrom
fData(metaplastic.lumi)$start <- metaplastic.probes$start
fData(metaplastic.lumi)$end <- metaplastic.probes$end
fData(metaplastic.lumi)$status <- metaplastic.probes$status
fData(metaplastic.lumi)$cytoband <- paste(metaplastic.probes$chrom, metaplastic.probes$cytoband, sep="")

sampleNames(metaplastic.lumi) <- gsub("\"", "", sampleNames(metaplastic.lumi))

metaplastic.lumi <- readPheno(metaplastic.lumi, "pheno.txt")
metaplastic.mapped.lumi <- removeUnmappedProbes(metaplastic.lumi)

save(metaplastic.mapped.lumi, file="20130519_metaplastic.mapped.lumi.RData")

## Copy number
## read all the data in
## remove chr Y
## final object is metaplasticCBS

system.time(metaplastic <- readBACE.cgh(exprsFile="Metaplastic.M.subset.txt", 
	featureDataFile="Metaplastic.fData.hg19.txt", stringsAsFactors=F, row.names=1))
metaplastic <- metaplastic[order(fData(metaplastic)$chrom, fData(metaplastic)$start),]

metaplastic <- readPheno(metaplastic, "pheno.txt")
system.time(metaplasticCBS <- cbsCGH(metaplastic))

SNP6mapping.noCNV <- read.delim("SNP6mapping.noCNV.toronto.txt", header=T, sep="\t", stringsAsFactors=F)
metaplasticCBS <- metaplasticCBS[is.element(featureNames(metaplasticCBS), SNP6mapping.noCNV$featureNames),]
metaplasticCBS <- callCGHStatesThreshold(metaplasticCBS, gainthresh=0.15, ampthresh=0.5, delthresh=-1, contig=3)

chrY <- which(fData(metaplasticCBS)$chrom %in% c(24,"Y"))
if (length(chrY)>0) {	
	metaplasticCBS <- metaplasticCBS[-chrY,]
}

library(gdata)
hd <- read.xls("allHD.ASCAT.ABSOLUTE.hg19.xlsx",1,stringsAsFactors=F) 

xx <- apply(assayDataElement(metaplasticCBS,"GL"),2,function(x){x[which(x== -2)] <- -1; x})
for (i in 1:nrow(hd)) {
	print(i)
	xx[which(fData(metaplasticCBS)$chrom==hd$chrom[i] & 
		((fData(metaplasticCBS)$start+fData(metaplasticCBS)$end)/2)>=hd$start[i] & 
		((fData(metaplasticCBS)$start+fData(metaplasticCBS)$end)/2)<=hd$end[i]),hd$Sample[i]] <- -2}

newCBS <- metaplasticCBS
assayDataElement(newCBS,"GL") <- xx

metaplasticCBS <- newCBS
rm(newCBS)
save(metaplasticCBS, file="20130519_metaplasticCBS.RData")


###################################### End of data preparation ###############################

### expression analysis
### hierarchical clustering and cluster stability assessment

metaplastic.heatmap <- removeProbesByMAD(metaplastic.mapped.lumi, mad.limit=1.2)
#  1566 probes with a median absolute deviation greater than 1.2 

metaplastic.heatmap <- centerGenes(metaplastic.heatmap)
metaplastic.heatmap <- removeReplicateProbesByMAD(metaplastic.heatmap)
# 1411 probes remaining

pheno_colours =brewer.pal(5,"Dark2")

gs <- rep("", length(fData(metaplastic.heatmap)$SYMBO))
gs[which(fData(metaplastic.heatmap)$symbol %in% c("CLDN3", "CLDN8", "CGN", "MYH11", "MYH14", "CDH1", "S100A1", "SOX10", "EPCAM"))] <- 
	fData(metaplastic.heatmap)$symbol[which(fData(metaplastic.heatmap)$symbol %in% 
	c("CLDN3", "CLDN8", "CGN", "MYH11", "MYH14", "CDH1", "S100A1", "SOX10", "EPCAM"))]

pdf("metaplastic.expression.heatmap.pdf",height=10,width=8)
heatmap.2(exprs(metaplastic.heatmap), trace="none",keysize=1,ColSideColors=pheno_colours[as.factor(pData(metaplastic.heatmap)$Cellular_type_simplified)], labRow=gs, scale="none", hclustfun=function(x){hclust(x, method="ward")}, distfun=function(x){as.dist(1-cor(t(x), method="pearson"))}, col=maPalette(high="darkorchid4", low="orange", mid="white", k=30))
dev.off()

library(pvclust)
set.seed(11237)
metaplastic.pv.exprs <- pvclust(exprs(metaplastic.heatmap), method.hclust="ward", method.dist="correlation")
pdf("pvclust_exprs.pdf")
plot(metaplastic.pv.exprs)
dev.off()
save(metaplastic.pv.exprs, file="20130611_pvclust_exprs_ward_correlation.RData")

### Differential expression analysis

metaplastic.de <- removeProbesByMAD(metaplastic.mapped.lumi, mad.limit=1)
#  2419 probes with a median absolute deviation greater than 1 

dietSAM2Class(metaplastic.de, metaplastic.de$chondroid.group, fData.ID="ILMN", 
	project="metaplastic.SAM.Chondroid.v.rest.q1", return.eset=F, q.value=1)
# 0 non-chondroid significant probes with Q value less than 1 %
# 48 Chondroid significant probes with Q value less than 1 %

dietSAM2Class(metaplastic.de, metaplastic.de$squamous.group, fData.ID="ILMN", 
	project="metaplastic.SAM.Squamous.v.rest.q1", return.eset=F, q.value=1)
# 19 Squamous significant probes with Q value less than 1 %
# 0 non-squamous significant probes with Q value less than 1 %

dietSAM2Class(metaplastic.de, metaplastic.de$spindle.group, fData.ID="ILMN", 
	project="metaplastic.SAM.Spindle.v.rest.q1", return.eset=F, q.value=1)
# 0 Spindle significant probes with Q value less than 1 %
# 190 non-spindle significant probes with Q value less than 1 %


## copy number analysis

plotFrequency(metaplasticCBS, device="PNG", project="Metaplastic.SNP6")

gl <- assayDataElement(metaplasticCBS,"GL")
band <- paste(fData(metaplasticCBS)$chrom, substr(fData(metaplasticCBS)$cytoband,1,1), sep="")
chrlen <- unlist(lapply(unique(band),function(cytoband) {
	max(fData(metaplasticCBS)$end[which(band== cytoband)])-min(fData(metaplasticCBS)$start[which(band== cytoband)])
}))
names(chrlen) <- unique(band)

newgl <- do.call("rbind", lapply(unique(band), function(arm) {
	smallgl <- gl[which(band==arm),]
	annot <- fData(metaplasticCBS)[which(band==arm),]
	for (i in 1:ncol(gl)) {
		rr <- rle(smallgl[,i])
		ends <- cumsum(rr$lengths)
		if(length(ends)==1) {
			starts=1
		} else {
			starts=c(1,ends[1:(length(ends)-1)]+1)
		}
		tab <- cbind(rr$values, starts, ends)
		for (n in which(tab[,1]==2)){
			if(annot$end[tab[n,3]]-annot $start[tab[n,2]] > chrlen[arm]/4) {
				smallgl[tab[n,2]:tab[n,3],i] <- 1
			}
		}
	}
	smallgl
}))

newCBS <- metaplasticCBS
assayDataElement(newCBS,"GL") <- newgl

plotFrequency(newCBS, device="PNG", project="Metaplastic.SNP6.focal")
listBreaksGL(newCBS,  contig=8, project="metaplastic.SNP6")
listBreaksGLAD(newCBS, gain.count=3, amp.count=2, project="metaplastic.SNP6", contig=8)

pdf("metaplastic.SNP6.heatmap.pdf", width=4, height=6)
cghHeatmap(metaplasticCBS, dist.method="euclidean", main="metaplastic SNP6", phenotypes=pData(metaplasticCBS)$Cellular_type_simplified, pheno.colours= pheno_colours)
dev.off()
gc()

png("metaplastic.SNP6.heatmap.png", width=1200, height=1800)
cghHeatmap(metaplasticCBS, dist.method="euclidean", main="metaplastic SNP6", phenotypes=pData(metaplasticCBS)$Cellular_type_simplified, pheno.colours= pheno_colours)
dev.off()
gc()


### copy number clusters

CNcluster1 = c("META41", "META42", "META52", "META39", "META30", "META55", "META31", "META59", "META64", "META32", "META49", "META50", "META36", "META47")
CNcluster2 = c("META53", "META62", "META37", "META40")

CNcluster_pheno <- rep(1, length(sampleNames(metaplasticCBS)))
CNcluster_pheno[which(sampleNames(metaplasticCBS) %in% CNcluster2)] <- 2

dietSAM2Class(metaplastic.mapped.lumi, CNcluster_pheno, fData.ID="ILMN", project="metaplastic.SAM.CNcluster1.v.CNcluster2", return.eset=F, q.value=1)
dietSAM2Class(metaplastic.mapped.lumi, CNcluster_pheno, fData.ID="ILMN", project="metaplastic.SAM.CNcluster1.v.CNcluster2", return.eset=F, q.value=5)
dietSAM2Class(metaplastic.mapped.lumi, CNcluster_pheno, fData.ID="ILMN", project="metaplastic.SAM.CNcluster1.v.CNcluster2", return.eset=F, q.value=10)
### NO DE GENES on at all q.value

metaplasticCBS.FE.CNcluster1.v.CNcluster2  <- fisherTestCGH(metaplasticCBS, pheno= CNcluster_pheno, 
	project="metaplastic_SNP6_Fisher_CNcluster1.v.CNcluster2")

listBreaksFisher(metaplasticCBS.FE.CNcluster1.v.CNcluster2, project="metaplasticCBS.FE.CNcluster1.v.CNcluster2", probeID=fData(metaplasticCBS.FE.CNcluster1.v.CNcluster2)$probeID)
latticePlotFishers(metaplasticCBS.FE.CNcluster1.v.CNcluster2, project="metaplasticCBS.FE.CNcluster1.v.CNcluster2")

save(metaplasticCBS.FE.CNcluster1.v.CNcluster2, file="20130520_metaplasticCBS.FE.CNcluster1.v.CNcluster2.RData")



### subtypes

metaplasticCBS.FE.chondroid.v.rest  <- fisherTestCGH(metaplasticCBS, 
	pheno= pData(metaplasticCBS)$chondroid.group, 
	project="metaplastic_SNP6_Fisher_chondoid.v.rest")
listBreaksFisher(metaplasticCBS.FE.chondroid.v.rest, project="metaplastic_SNP6_Fisher_chondoid.v.rest", 
	probeID=fData(metaplasticCBS.FE.chondroid.v.rest)$probeID)
latticePlotFishers(metaplasticCBS.FE.chondroid.v.rest, project="metaplastic_SNP6_Fisher_chondoid.v.rest")
save(metaplasticCBS.FE.chondroid.v.rest, file="20130519_metaplasticCBS.FE.chondroid.v.rest.RData")

metaplasticCBS.FE.spindle.v.rest  <- fisherTestCGH(metaplasticCBS, 
	pheno= pData(metaplasticCBS)$spindle.group, 
	project="metaplastic_SNP6_Fisher_spindle.v.rest")
listBreaksFisher(metaplasticCBS.FE.spindle.v.rest, project="metaplastic_SNP6_Fisher_spindle.v.rest", 
	probeID=fData(metaplasticCBS.FE.spindle.v.rest)$probeID)
latticePlotFishers(metaplasticCBS.FE.spindle.v.rest, project="metaplastic_SNP6_Fisher_spindle.v.rest")
save(metaplasticCBS.FE.spindle.v.rest, file="20130519_metaplasticCBS.FE.spindle.v.rest.RData")

metaplasticCBS.FE.squamous.v.rest  <- fisherTestCGH(metaplasticCBS, 
	pheno= pData(metaplasticCBS)$squamous.group, 
	project="metaplastic_SNP6_Fisher_squamous.v.rest")
listBreaksFisher(metaplasticCBS.FE.squamous.v.rest, project="metaplastic_SNP6_Fisher_squamous.v.rest", 
	probeID=fData(metaplasticCBS.FE.squamous.v.rest)$probeID)
latticePlotFishers(metaplasticCBS.FE.squamous.v.rest, project="metaplastic_SNP6_Fisher_squamous.v.rest")
save(metaplasticCBS.FE.squamous.v.rest, file="20130519_metaplasticCBS.FE.squamous.v.rest.RData")



#### overlay

metaplastic.exp.cgh.overlay <- medianCGH(metaplastic.mapped.lumi, metaplasticCBS)
metaplastic.exp.cgh.overlay <- expressionCGHCorrelation(metaplastic.exp.cgh.overlay)
metaplastic.exp.cgh.overlay <- wilcoxTestCGH(metaplastic.exp.cgh.overlay)
metaplastic.exp.cgh.overlay <- readPheno(metaplastic.exp.cgh.overlay, "pheno.txt")

save(metaplastic.exp.cgh.overlay, file="20130519_metaplastic.exp.cgh.overlay.RData")

correlated.probes <- fData(metaplastic.exp.cgh.overlay)[which(fData(metaplastic.exp.cgh.overlay)$pearson.adjp < 0.05),]
correlated.probes <- correlated.probes[order(correlated.probes$pearson.cor, decreasing=T),]
# 3071 significantly correlated genes
write.table(correlated.probes, "correlated.probes.xls", sep="\t", row.names=F, na="")

amp.overexpressed.probes <- fData(metaplastic.exp.cgh.overlay)[which(fData(metaplastic.exp.cgh.overlay)$Wilcox.adjp.amp < 0.05),]
amp.overexpressed.probes <- amp.overexpressed.probes[order(amp.overexpressed.probes$pearson.cor, decreasing=T),]
# 196 genes overexpressed when amplified
write.table(amp.overexpressed.probes, "amp.overexpressed.probes.xls", sep="\t", row.names=F, na="")

del.underexpressed.probes <- fData(metaplastic.exp.cgh.overlay)[which(fData(metaplastic.exp.cgh.overlay)$Wilcox.adjp.del < 0.05),]
del.underexpressed.probes <- del.underexpressed.probes[order(del.underexpressed.probes$pearson.cor, decreasing=T),]
# 0 genes under-expressed when deleted
write.table(del.underexpressed.probes, "del.underexpressed.probes.xls", sep="\t", row.names=F, na="")


metaplastic.exp.cgh.uni <- removeReplicateProbesByMAD(metaplastic.exp.cgh.overlay)
# 18931 probes remaining

splitHeatmap(metaplastic.exp.cgh.uni, chrom=8, start= 99413631, end= 104345094, project="8q22", device="PDF", heatmap.scale=1)
splitHeatmap(metaplastic.exp.cgh.uni, chrom=8, start= 117654369, end= 131029375, project="8q24.12-8q24.2", device="PDF", heatmap.scale=1)

