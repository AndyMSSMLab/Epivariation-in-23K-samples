### Created by P. Garg
###     Program to permorm lumi and BMIQ normalization
###
###
### Example Usage: R CMD BATCH --no-save --no-restore '--args <input file name> <annotation file> <chr> <prefix> <remove_sample_file> runNormalization.r
### Input file must contain TargetID (ProbeID), Signal_A, Signal_B and Detection Pval for each sample
### <Annotation file> must contain ProbeID, location of targetted CpG (Chr_CpG, Start_CpG, End_CpG) and Infinium type (whether I or II in column Infinium_Type)
### <chr> can be auto for autosomal or chr1 for chromosome 1 probes etc
### <prefix> is required for output filename
### <remove_sample_file> remove sample file should contain list of samples to be removed. Can be left empty if no sample has to be removed
##########################################################################################################################################################

args=(commandArgs(TRUE))
if(length(args)==0)
{
    print("No arguments supplied. Please try again.")
    ##supply default values
   quit("no")
}else{
        args=c(args)
}

print(paste0("Input file: ", args[1]))
print(paste0("Annotation file: ", args[2]))
print(paste0("Chromosome: ", args[3]))
print(paste0("Prefix for Output file: ", args[4]))

input_file = args[1]
annotation_file = args[2]
chr = args[3]
prefix = args[4]
args
remove_samples = c()
if(length(args)==5){ remove_samples = unlist(read.table(args[5], as.is = T, sep = "\t"))}

###Load relevant libraries
.libPaths()
library(lumi)
library(methylumi)
library(FDb.InfiniumMethylation.hg19)
source("BMIQ_1.3.R")

print(paste0("Reading input file: ", input_file, " ..."))
data = lumiMethyR(input_file, lib="FDb.InfiniumMethylation.hg19", sep = "\t")
###print("Rows and Samples in input file: ",  paste(dim(data), collapse = ","))

print(paste0("Reading annotation file: ", input_file, " ..."))
anno = read.table(annotation_file, header = T, as.is = T, sep = "\t", check.names = F)
rownames(anno) = anno[,1]
if(chr == "auto"){
  probes = rownames(anno[grepl("chr[0-9]",anno$Chr_CpG),])  ###select autsomal probes
}else{
  probes = rownames(anno[anno$Chr_CpG == chr,])  ### probes on selected chr
}

data2 = data[intersect(probes, featureNames(data)),setdiff(sampleNames(data),remove_samples)]
rm(data)

### Lumi Normalizaition

data2.c.adj <- lumiMethyC(data2, verbose = TRUE); ##Color correction
data2.b.adj <- lumiMethyB(data2.c.adj); rm(data2.c.adj); ## Background correction
data2.c.q <- lumiMethyN(data2.b.adj, method='quantile'); rm(data2.b.adj) ## Normalization across samples

###make before and after normalization plots
pdf(paste0(prefix,"_norm_plot.pdf"));
  plotColorBias1D(data2, main = "ColorBias1D before normalization");
  plotColorBias1D(data2.c.q, main = "ColorBias1D after normalization");
  density(estimateIntensity(data2), xlab="log2(CpG-site Intensity)", addLegend = F, main = "Density: Intensity before normalization");
  density(estimateIntensity(data2.c.q), xlab="log2(CpG-site Intensity)", addLegend = F,  main = "Density: Intensity after normalization");
  density(data2, main="Density: M-value before quantile normalization", addLegend = F);
  density(data2.c.q, main="Density: M-value after quantile normalization", addLegend = F);
  
### BMIQ normalization to correct for distribution of beta values for Infinium I and Infinium II

beta = estimateBeta(data2.c.q, returnType="matrix", offset = 100); rm(data2); rm(data2.c.q)
type = anno[rownames(beta),]$Infinium_Type
design = rep(1,length(type))
design[type=="II"] = 2
beta.norm = sapply(colnames(beta), function(t){BMIQ(unlist(beta[,t]), design, sampleID=t, plots=FALSE, nfit=3000)$nbeta})
names(design) = rownames(beta)

### Make plot of beta values before and after BMIQ normalization
plot(0,0, col = "white", ylim = c(0,6), xlim = c(-0.1,1), ylab = "Density", xlab = "Beta-values", main = "Infinium I and II probe normalization")
for(i in 1:ncol(beta.norm)){
  points(density(beta[names(design[design==1]),i]), type = "l")
}
for(i in 1:ncol(beta.norm)){
  points(density(beta[names(design[design==2]),i]), type = "l", col = "red")
}
for(i in 1:ncol(beta.norm)){
  points(density(beta.norm[names(design[design==2]),i]), type = "l", col = "green")
}
legend("topright", legend = c("Infinium I probes", "Infinium II probes", "Infinium II probes BMIQ"), text.col = c("black","red","green"))
dev.off()

print("Done with normalization")

beta.anno = cbind(anno[rownames(beta.norm),], round(beta.norm,4))

write.table(beta.anno,paste0(prefix,"_norm.txt"), row.names = F, quote = F, sep = "\t")
