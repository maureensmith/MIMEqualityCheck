#!/usr/local/bin Rscript --vanilla --slave


##############################################################################
# sam2counts Sanity Check:
# plot  mutation frequencies for each file for 1d, 2d, 3d (if present)
# plot coverage (distribution) for each position
##############################################################################

# clear the environment variables
#rm(list = ls())

# install needed packages
dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")"))))
    return(TRUE)
  
  install.packages(package)
  return (eval(parse(text=paste("require(",package,")"))))
}

#"ggformula"
for(p in c("ggplot2", "seqinr", "stringr" , "gridExtra", "grid")) {
  dynamic_require(p)
}


#Read in arguments and normalize paths
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  cat("\nCall the script with at least 3 arguments: countDirectory referenceFile resultDirectory (sampleSheetFile)\n
The countDirectory contains the count files in subdirectories (1d, 2d, 3d).\n
The reference file contains the reference sequence in a fasta file.\n
The plots are saved in the resultDirectory.\n
An optional fourth argument sampleSheetFile can be given, in order to change the sample names for the barcodes of the count filenames if it contains the columns 'Sample' and 'Encoding' \n\n")
  
  #terminate without saving workspace
  quit("no")
}

#cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
countDir <- normalizePath(args[1])
referenceFile <- normalizePath(args[2])
resultDir <- file.path(args[3])
sample_sheet_file=""
if(length(args) > 3)
  sample_sheet_file=normalizePath(args[4])

if(!dir.exists(resultDir))
  dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)
resultDir <- normalizePath(resultDir)

cat(c("Arguments: ", countDir, referenceFile, resultDir, sample_sheet_file, "\n"), sep = "\n")





# set working directory to call other R Scripts
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

this.dir <- getSrcDirectory(function(x) {x})

if (rstudioapi::isAvailable()) {
  this.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
}else {
  this.dir <- getScriptPath()
}
setwd(this.dir)

source("plotRoutinesQualityCheck.R")
source("readCountData.R")


nucl.df = getNucleotide() 
refSeq.df <- getRefSeq(referenceFile, nucl.df)

coverage1d_all.df <- readCountData_1d(countDir, refSeq.df, nucl.df)
coverage2d_all.df <- readCountData_2d(countDir, refSeq.df, nucl.df)

# if available, use sample names instead of encoding for the plots
if(file.exists(sample_sheet_file)) {
  "samle"
  sample_sheet=read.table(sample_sheet_file,header=T, sep=";")
  if("Sample" %in% colnames(sample_sheet) & "Encoding" %in% colnames(sample_sheet)) {
    coverage1d_all.df["sample"] = sample_sheet$Sample[unlist(sapply(coverage1d_all.df$sample, function(x) which(sample_sheet$Encoding==as.numeric(x))))]
    coverage2d_all.df["sample"] = sample_sheet$Sample[unlist(sapply(coverage2d_all.df$sample, function(x) which(sample_sheet$Encoding==as.numeric(x))))]
  }
}

##### plots 1d ##########
samples <- unique(coverage1d_all.df$sample)
ps_mutCounts <- lapply(samples,function(s) plotMutCountsForEachSample_1d(coverage1d_all.df[coverage1d_all.df$sample==s,],
                                                                    refSeq.df = refSeq.df, nucl.df= nucl.df))
#p_mutCounts_perSample<-multiplot(plotlist = ps_mutCounts, cols = floor(sqrt(length(samples))))
p_mutCounts_perSample<-grid.arrange(grobs=ps_mutCounts,ncol=floor(sqrt(length(samples))))
outputFile = paste0(resultDir, "/mut_counts_perSample_1d.pdf")
ggsave(file=outputFile, plot = p_mutCounts_perSample, device = "pdf",width=40, heigh=10*floor(sqrt(length(samples))), units = "cm")


ps_mutFreqs <- lapply(samples,function(s) plotMutFreqPerWtForEachSample_1d(coverage1d_all.df[coverage1d_all.df$sample==s,],
                                                                         refSeq.df = refSeq.df, nucl.df= nucl.df))
#p_mutCounts_perSample<-multiplot(plotlist = ps_mutCounts, cols = floor(sqrt(length(samples))))
p_mutFreqs_perSample<-grid.arrange(grobs=ps_mutFreqs,ncol=floor(sqrt(length(samples))))
outputFile = paste0(resultDir, "/mut_freqs_perSample_1d.pdf")
ggsave(file=outputFile, plot = p_mutFreqs_perSample, device = "pdf",width=40, heigh=10*floor(sqrt(length(samples))), units = "cm")


# plot mutation frequencey per position and wild type
ps_mutFreqPerWt<- lapply(samples, function(s) {
  ps<-plotMutFreqPerPosforEachWt(coverage1d_all.df[coverage1d_all.df$sample==s,], nucl.df)
  p=grid.arrange(grobs=ps,ncol=4,top = textGrob(s,gp=gpar(fontsize=20)) )
  outputFile = paste0(resultDir, "/mut_freq_perPos_perWt_1d_",s,".pdf")
  ggsave(filename = outputFile, plot = p, device = "pdf",width=40, height=10, units = "cm")
  return(p)})
# p_mutFreqPerWt_perSample<-grid.arrange(grobs=unlist(ps_mutFreqPerWt),ncol=4)
# outputFile = paste0(resultDir, "/mut_freq_perPos_perWt_1d.pdf")
# ggsave(filename = outputFile, plot = p_mutFreqPerWt_perSample, device = "pdf",width=40, height=10*length(samples), units = "cm")

# plot coverage
p_cov_all <- plotCoverageAll_1d(coverage1d_all.df) 
outputFile = paste0(resultDir, "/all_coverage_1d.pdf")
ggsave(file=outputFile, plot = p_cov_all, device = "pdf",width=40, heigh=10, units = "cm")

p_cov_all_log  <- plotCoverageAll_1d(coverage1d_all.df, logCov=T) 
outputFile = paste0(resultDir, "/all_log_coverage_1d.pdf")
ggsave(file=outputFile, plot = p_cov_all_log, device = "pdf",width=40, heigh=10, units = "cm")

#plot mutFreq
p_mut_all_log<- plotLogMutFreqAll_1d(coverage1d_all.df)
outputFile = paste0(resultDir, "/all_log_mutRate_1d.pdf")
ggsave(file=outputFile, plot = p_mut_all_log, device = "pdf",width=30, heigh=10, units = "cm")


#plot entropy 
p_entropy<- plotEntropy_1d(coverage1d_all.df)
outputFile = paste0(resultDir, "/all_entropy_perPos_1d.pdf")
ggsave(file=outputFile, plot = p_entropy, device = "pdf",width=30, heigh=10, units = "cm")

###### plots 2d ##########

p_cov_all_2d<-plotCoverageAll_2d(coverage2d_all.df)
outputFile = paste0(resultDir, "/all_coverage_2d.pdf")
ggsave(file=outputFile, plot = p_cov_all_2d, device = "pdf",width=40, heigh=10, units = "cm")

p_cov_all_log_2d <- plotCoverageAll_2d(coverage2d_all.df, logCov=T)
outputFile = paste0(resultDir, "/all_log_coverage_2d.pdf")
ggsave(file=outputFile, plot = p_cov_all_log_2d, device = "pdf",width=40, heigh=10, units = "cm")



