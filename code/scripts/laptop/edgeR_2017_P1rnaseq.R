#3/25/17
#script for my macbook
#R version 3.3.1
#EdgeR analysis of HTSEQcount output

setwd("/Users/Tara/GoogleDrive/P1_RNAseq/output/edgeR/")
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")

#install.packages("reshape2")
require("edgeR")
require("ggplot2")

#import data 
counts <- read.table("/Users/Tara/GoogleDrive/P1_RNAseq/data/processed/htseq/P1_htseq_all.txt")
#get rid of last five lines from htseq output
counts <- counts[1:5844,]


groups <- c("control", "control", "control", "la", "la", "la",
            "bdl", "bdl", "bdl")

#edgeR stores data in a simple list-based data object called a DGEList. 
#This type of object is easy to use because it can be manipulated like any list in R. 
#You can make this in R by specifying the counts and the groups in the function DGEList()
d <- DGEList(counts=counts,group=groups)


#First get rid of genes which did not occur frequently enough. 
#We can choose this cutoff by saying we must have at least 100 counts per million (calculated with cpm() in R) on any particular gene that we want to keep. 
#In this example, we're only keeping a gene if it has a cpm of 100 or greater for at least two samples.

cpm<-(cpm(d))
#returns counts per mILLION FROM A DGElist
#divides raw counts by library size and then multiplies by 1 mil

total_gene<-apply(d$counts, 2, sum) # total gene counts per sample

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep,]
dim(d)

#Recalculate library sizes of the DGEList object
d$samples$lib.size <- colSums(d$counts)
d$samples
names <- as.character(counts[,1])

#calcNormFactors" function normalizes for RNA composition
#it finds a set of scaling factors for the library sizes (i.e. number of reads in a sample)
#that minimize the log-fold changes between the samples for most genes 
d<- calcNormFactors(d)
#plot MDS - multi-dimensional scaling plot of the RNA samples in which distances correspond to leading log-fold-changes between each pair of RNA samples. 
tiff("P1_RNAseq.tiff")
plotMDS<-plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottom", as.character(unique(d$samples$group)), col=1:3, pch=20)
dev.off()

#create a design matrix using model.matrix function
fac <- paste(groups, sep=".")
fac <- factor(fac)
design <- model.matrix (~0+fac)
colnames(design) <- levels(fac)

#estimate dispersion
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
#mean variance plot
#grey = raw variance
#light blue = tagwise dispersion varaince
#dark blue = common dispersion
#black line = poission variance
tiff("Meanvarplot.tiff")
meanVarPlot <- plotMeanVar( d , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )
dev.off()

#BCV plot
#the BCV plot shows the common trended and genewise dispersions
#as a function of average logCPM
#plotBCV(y)

#fit genewise negative binomial general linear model
fit <- glmFit(d, design)

#the likelihood ratio stats are computed for the comparison of interest
#compare P1 in TH to P1 in TH + lactic acid 
##### positive logFC for genes upreg in la compared to control 
control.vs.la <- glmLRT(fit, contrast=c(0,-1,1))

#compare control to bdl 
####positive logFC means higher expressin in bdl
control.vs.bdl <- glmLRT(fit, contrast=c(1,-1,0))

##compare (pH 5 group - pH 7 group) to 0:
####positive logFC means higher expression in bdl 
bdl.vs.la <- glmLRT(fit, contrast=c(1,0,-1))
#make files with ALL genes
bdl.vs.la.tags <- topTags(bdl.vs.la, adjust="BH", n=Inf)
control.vs.bdl.tags<- topTags(control.vs.bdl, adjust="BH", n=Inf)
control.vs.la.tags <- topTags(control.vs.la, adjust="BH", n=Inf)

####FDR =.1 (10 in 100 sig genes wrong) AS CUT-OFF
### Write out stats to .csv and also filter by FDR (0.05)
cutoff=0.1

bdl.vs.la.dge <- decideTestsDGE(bdl.vs.la, p=cutoff, adjust="BH")
#obtain number of upreg, downreg genes
bdl.vs.la.dge.summary<-summary(bdl.vs.la.dge)
bdl.vs.la.dge.tags<-rownames(d)[as.logical(bdl.vs.la.dge)]
tiff("smearplot_bdl_la")
plotSmear(bdl.vs.la, de.tags = bdl.vs.la.dge.tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes Bdl compared to LA", cex.main=.8)
dev.off()

control.vs.la.dge <- decideTestsDGE(control.vs.la, p=cutoff, adjust="BH")
#obtain number of upreg, downreg genes
control.vs.la.dge.summary<-summary(control.vs.la.dge)
control.vs.la.dge.tags<-rownames(d)[as.logical(control.vs.la.dge)]
tiff("smearplot_control_la")
plotSmear(control.vs.la, de.tags = control.vs.la.dge.tags, xlab="LogCPM (counts per million)", 
          main= "Upregulated and downregulated genes in LA compared to control", cex.main=.8)
dev.off()

control.vs.bdl.dge <- decideTestsDGE(control.vs.bdl, p=cutoff, adjust="BH")
#obtain number of upreg, downreg genes
control.vs.bdl.dge.summary<-summary(control.vs.bdl.dge)
control.vs.bdl.dge.tags<-rownames(d)[as.logical(control.vs.bdl.dge)]
tiff("smearplot_control_bdl")
plotSmear(control.vs.bdl, de.tags = control.vs.bdl.dge.tags, xlab="LogCPM (counts per million)", 
          main= "Upreg and downreg genes in BDL compared to control", cex.main=.8)
dev.off()

#write out summaries from comparisons
summary_all <- cbind(control.vs.bdl.dge.summary, control.vs.la.dge.summary, bdl.vs.la.dge.summary)
colnames(summary_all) <- c("control_bdl", "control_la", "bdl_la")
write.table(summary_all, "summary_DGE_all.txt", sep="\t")

#export DGE (FDR <- 0.1)
out_control_bdl<- topTags(control.vs.bdl, n=Inf, adjust.method="BH")
keep_control_bdl <- out_control_bdl$table$FDR <= 0.1
de_control_bdl<-out_control_bdl[keep_control_bdl,]
write.table(x = de_control_bdl, file = "de_control_bdl_filter.txt", sep ="\t", quote = FALSE)

out_control_la<- topTags(control.vs.la, n=Inf, adjust.method="BH")
keep_control_la <- out_control_la$table$FDR <= 0.1
de_control_la<-out_control_la[keep_control_la,]
write.table(x = de_control_la, file = "de_control_la_filter.txt", sep ="\t", quote = FALSE)

out_bdl_la<- topTags(bdl.vs.la, n=Inf, adjust.method="BH")
keep_bdl_la<- out_bdl_la$table$FDR <= 0.1
de_bdl_la<-out_bdl_la[keep_bdl_la,]
write.table(x = de_bdl_la, file = "de_bdl_la_filter.txt", sep ="\t", quote = FALSE)

