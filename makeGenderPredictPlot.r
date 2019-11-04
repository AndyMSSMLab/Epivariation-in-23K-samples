##########################################################################################################################################################
### Program to check Gender using mean Beta on chrX and missingness on chrY
###
###
### Meant to be run in BATCH providing a filename in --args that corresponds to chromosome. Can be run on each chromosome in this way.
###
###
### Example Usage: R CMD BATCH --no-save --no-restore '--args <input file name> <annotation file> <prefix> <Sample Info> <remove_sample_file> checkGender.r
### <chr> can be auto for autosomal or chr1 for chromosome 1 probes etc
### <Sample Info> sample sheet with column names sampleID with same name as samples in input file. Can be left empty
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
print(paste0("Prefix for Output file: ", args[3]))

##/sc/orga/projects/sharpa01a/Paras/Dataset/illumina_450K_annotation

input_file = args[1]
annotation_file = args[2]
prefix = args[3]
args
sample_info = NA
remove_samples = c()
if(length(args)>=4){ 
  if(file.exists(args[4])){
    sample_info = read.table(args[4], as.is = T, sep = "\t", header = T, check.names = F, quote = "\"")
    rownames(sample_info) <- sample_info[,1]
  }
}
if(length(args)==5){ remove_samples = unlist(read.table(args[5], as.is = T, sep = "\t"))}

###Load relevant libraries
.libPaths()
library(lumi)
library(methylumi)
library(FDb.InfiniumMethylation.hg19)
library(ggplot2)

print(paste0("Reading input file: ", input_file, " ..."))
data = lumiMethyR(input_file, lib="FDb.InfiniumMethylation.hg19", sep = "\t")

anno = read.table(annotation_file, header = T, as.is = T, sep = "\t", check.names = F)
rownames(anno) = anno[,1]
probes = rownames(anno[anno$Mapping=="UM"&is.na(anno$SNP_in_5bp)&anno$Chr_CpG == "chrX",])  ### probes on chrX
dataX = data[intersect(probes, featureNames(data)),setdiff(sampleNames(data),remove_samples)]
probes = rownames(anno[anno$Mapping=="UM"&is.na(anno$SNP_in_5bp)&anno$Chr_CpG == "chrY",])  ### probes on chrY
dataY = data[intersect(probes, featureNames(data)),setdiff(sampleNames(data),remove_samples)]
rm(data)

###Get mean methylation for chrX and ratio of failed detection P value on chrY
methX <- apply(betas(dataX),2, mean, na.rm = T)
detPY <- apply(pvals(dataY),2, function(t){ round(sum(t>0.01)/length(t),4)})

out <- data.frame(sampleID = names(methX), meanBetaOnChrX = round(methX,4), fractionfailedDetPOnChrY = detPY, stringsAsFactors = F)
write.table(out,paste0(prefix,"_genderCheck.txt"), row.names = F, quote = F, sep = "\t")

if(is.data.frame(sample_info)){
  out = cbind(out, sample_info[out$sampleID,])
  head(out)
}

jpeg(paste0(prefix,"_genderCheck_plot.jpeg"), height = 2400, width = 3000, res = 300, pointsize = 12)
g = ggplot(out, aes(x = meanBetaOnChrX, y = fractionfailedDetPOnChrY))

if(any(colnames(out)=="gender")){
  g = g + geom_point(size = 2, aes(color = gender))
}else{
  g = g+geom_point(size = 2)
}

g = g +  xlab("Mean Beta values for probes on chrX") +
  ylab("Fraction of probes with failed detection P values on chrY") +
  ggtitle(paste0("Gender Check for ", prefix)) +
  theme_bw()
print(g)
dev.off()
