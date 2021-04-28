##############################################################################
# sam2counts Sanity Check:
# plot  mutation frequencies for each file for 1d, 2d, 3d (if present)
# plot coverage (distribution) for each position
##############################################################################

require(ggplot2)
require(scales)

#### 1d plots #######
plotMutCountsForEachSample_1d <- function(coverage1d.df, refSeq.df, nucl.df) {
  # has been filtered before, to only contain one sample  
  s <- coverage1d.df$sample[1]
  L <- nrow(refSeq.df)
  #       #transform mut1 2 and 3 to the respective nucleotide
    mut_1 =sapply(refSeq.df$asInt, function(x) nucl.df$asInt[-x][1])
    mut_2 = sapply(refSeq.df$asInt, function(x) nucl.df$asInt[-x][2])
    mut_3 = sapply(refSeq.df$asInt, function(x) nucl.df$asInt[-x][3])
    
    mut_counts.df=data.frame(wt=rep(nucl.df$asChar[coverage1d.df$wt],3),
                             counts=c(coverage1d.df$mutCount1, 
                                      coverage1d.df$mutCount2, 
                                      coverage1d.df$mutCount3),
                             mutNo=c(rep("mut1", L), rep("mut2", L),rep("mut3", L)),
                             mut=nucl.df$asChar[c(mut_1, mut_2, mut_3)])
  
    p_mutCounts <- ggplot(data = mut_counts.df)+
      geom_boxplot(aes(x=wt, y=log10(counts), fill=mut), alpha=0.8)+
      labs(x="wt", y=expression(log[10](counts)))+
      ggtitle(s)+
    theme(axis.text = element_text(size=10),
          axis.title = element_text(size=10,face="bold"),
          title = element_text(size=10),
          legend.text = element_text(size=8), 
          legend.title = element_text(size=8, face="bold"),
          legend.position="bottom", 
          legend.box="vertical", 
          legend.margin=margin())

  return(p_mutCounts)
}

plotMutFreqPerWtForEachSample_1d <- function(coverage1d.df, refSeq.df, nucl.df) {
  # has been filtered before, to only contain one sample  
  s <- coverage1d.df$sample[1]
  L <- nrow(refSeq.df)
  #       #transform mut1 2 and 3 to the respective nucleotide
  mut_1 =sapply(refSeq.df$asInt, function(x) nucl.df$asInt[-x][1])
  mut_2 = sapply(refSeq.df$asInt, function(x) nucl.df$asInt[-x][2])
  mut_3 = sapply(refSeq.df$asInt, function(x) nucl.df$asInt[-x][3])
  mut_counts.df=data.frame(wt=rep(nucl.df$asChar[coverage1d.df$wt],3),
                           counts=c(coverage1d.df$mutCount1/coverage1d.df$coverage, 
                                    coverage1d.df$mutCount2/coverage1d.df$coverage, 
                                    coverage1d.df$mutCount3/coverage1d.df$coverage),
                           mutNo=c(rep("mut1", L), rep("mut2", L),rep("mut3", L)),
                           mut=nucl.df$asChar[c(mut_1, mut_2, mut_3)])
  
  p_mutCounts <- ggplot(data = mut_counts.df)+
    geom_boxplot(aes(x=wt, y=log10(counts), fill=mut), alpha=0.8)+
    labs(x="wt", y=expression(log[10](frequency)))+
    ggtitle(s)+
    theme(axis.text = element_text(size=10),
          axis.title = element_text(size=10,face="bold"),
          title = element_text(size=10),
          legend.text = element_text(size=8), 
          legend.title = element_text(size=8, face="bold"),
          legend.position="bottom", 
          legend.box="vertical", 
          legend.margin=margin())
  
  return(p_mutCounts)
}


plotMaxMut_1d <- function(refSeqChar, max_mut_nucl) {
  p_maxMut<- ggplot()+
  geom_bar(aes(x=refSeqChar, fill=max_mut_nucl))+
  labs(x="wt", y="max mut count per pos.", fill="")+
  ggtitle(f_name)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))
  return(p_maxMut)
}


# plot mutation frequencey per position and wild type
plotMutFreqPerPosPerWt <- function(i, coverage1d_all.df, nucl.df, mut=0) {
  wtIdx1 <- which(coverage1d_all.df$wt == i & coverage1d_all.df$mutRate1>0)
  wtIdx2 <- which(coverage1d_all.df$wt == i & coverage1d_all.df$mutRate2>0)
  wtIdx3 <- which(coverage1d_all.df$wt == i & coverage1d_all.df$mutRate3>0)
    # wtIdx1 <- which(ref.df$asInt == i & coverage1d_all.df$mutRate1>0)
    # wtIdx2 <- which(ref.df$asInt == i & coverage1d_all.df$mutRate2>0)
    # wtIdx3 <- which(ref.df$asInt == i & coverage1d_all.df$mutRate3>0)
    p_mutFreq_perPos_perWt_1d <- ggplot()
    if(mut==0) {
      p_mutFreq_perPos_perWt_1d<-p_mutFreq_perPos_perWt_1d + geom_line(data=coverage1d_all.df[wtIdx1,], 
                                            aes(x=position, y=log10(mutRate1),col=nucl.df$asChar[-i][1], group=sample), size=1, alpha=0.7,na.rm = T)+
        geom_line(data=coverage1d_all.df[wtIdx2,], aes(x=position, y=log10(mutRate2), 
                                                       col=nucl.df$asChar[-i][2], group=sample), size=1, alpha=0.7,na.rm = T)+
        geom_line(data=coverage1d_all.df[wtIdx3,], aes(x=position, y=log10(mutRate3), 
                                                       col=nucl.df$asChar[-i][3], group=sample), size=1, alpha=0.7,na.rm = T)+
        ggtitle(paste0("wild type ", nucl.df$asChar[i]))+
        scale_color_manual(values=nucl.df$color[-i])
    } else{
      p_mutFreq_perPos_perWt_1d<-p_mutFreq_perPos_perWt_1d + 
        ggtitle(paste0("wild type ", nucl.df$asChar[i], " -> ", nucl.df$asChar[-i][mut]))
    } 
    if(mut==1)  
      p_mutFreq_perPos_perWt_1d <- p_mutFreq_perPos_perWt_1d + geom_line(data=coverage1d_all.df[wtIdx1,], aes(x=position, y=log10(mutRate1),col=sample, group=sample), size=1, alpha=0.7,na.rm = T)
    if( mut==2)  
      p_mutFreq_perPos_perWt_1d <- p_mutFreq_perPos_perWt_1d + geom_line(data=coverage1d_all.df[wtIdx2,], aes(x=position, y=log10(mutRate2), col=sample, group=sample), size=1, alpha=0.7,na.rm = T)
    if(mut==3)  
      p_mutFreq_perPos_perWt_1d <- p_mutFreq_perPos_perWt_1d + geom_line(data=coverage1d_all.df[wtIdx3,], aes(x=position, y=log10(mutRate3), col=sample, group=sample), size=1, alpha=0.7,na.rm = T)
    p_mutFreq_perPos_perWt_1d <- p_mutFreq_perPos_perWt_1d + labs(x="position", y=expression(log[10] ("mut. freq.")), col="mutation")+
      #ggtitle(paste0("Sample ", f_name, ", wild type ", nucl.df$asChar[i]))+
      ylim(-6,0)+
      theme(axis.text = element_text(size=8),
            axis.title = element_text(size=8),
            legend.text = element_text(size=4),
            legend.title = element_text(size=4),
            title = element_text(size=8),
            legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
            guides(color=F)
    
    return(p_mutFreq_perPos_perWt_1d)
}

plotMutFreqPerSamplePerWtBoxplot <- function(i, coverage1d_all.df, nucl.df) {
  wtIdx1 <- which(coverage1d_all.df$wt == i & coverage1d_all.df$mutRate1>0)
  wtIdx2 <- which(coverage1d_all.df$wt == i & coverage1d_all.df$mutRate2>0)
  wtIdx3 <- which(coverage1d_all.df$wt == i & coverage1d_all.df$mutRate3>0)
  # wtIdx1 <- which(ref.df$asInt == i & coverage1d_all.df$mutRate1>0)
  # wtIdx2 <- which(ref.df$asInt == i & coverage1d_all.df$mutRate2>0)
  # wtIdx3 <- which(ref.df$asInt == i & coverage1d_all.df$mutRate3>0)
  mutFreq.df = data.frame(sample=c(coverage1d_all.df$sample[wtIdx1], coverage1d_all.df$sample[wtIdx2], coverage1d_all.df$sample[wtIdx3]),
                          mutRate=c(coverage1d_all.df$mutRate1[wtIdx1], coverage1d_all.df$mutRate2[wtIdx2], coverage1d_all.df$mutRate3[wtIdx3]),
                          mutation=c(rep(nucl.df$asChar[-i][1], length(wtIdx1)), 
                                rep(nucl.df$asChar[-i][2], length(wtIdx2)),
                                rep(nucl.df$asChar[-i][3], length(wtIdx3))))
  p_mutFreq_perSamples_perWt_1d <- ggplot(data=mutFreq.df) +
    geom_boxplot(aes(x=sample, y=log10(mutRate),fill=mutation), alpha=0.8,na.rm = T)+
    ggtitle(paste0("wild type ", nucl.df$asChar[i]))+
    scale_fill_manual(values=nucl.df$color[-i])+
    labs(x="position", y=expression(log[10] ("mut. freq.")), fill="mutation")+
    #ggtitle(paste0("Sample ", f_name, ", wild type ", nucl.df$asChar[i]))+
    ylim(-6,0)+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1, size=12),
          axis.title = element_text(size=14),
          legend.text = element_text(size=14),
          legend.title = element_text(size=14))
  
  return(p_mutFreq_perSamples_perWt_1d)
}

plotMutFreqPerPosforEachWt <- function(coverage1d_all.df, nucl.df) {
  ps_mutFreqPerWt <- lapply(seq(4), plotMutFreqPerPosPerWt, coverage1d_all.df=coverage1d_all.df, nucl.df=nucl.df)
    
  #outputFile = paste0(resultDir, "/",f_name,"_mutFreq_perPos_wt",nucleotides[i],"_1d.pdf")
  #ggsave(filename = outputFile, plot = p_mutFreq_perPos_perWt_1d, device = "pdf",width=16, height=8, units = "cm")
  return(ps_mutFreqPerWt)
}


plotCoverageAll_1d <- function(coverage1d_all.df, logCov=F) {
  p_cov_all <- ggplot(data=coverage1d_all.df) +
    scale_y_continuous(labels = scales::comma, limits = c(0,NA))+
    theme(axis.text = element_text(size=10),
          axis.title = element_text(size=12,face="bold"),
          legend.text = element_text(size=10), 
          legend.title = element_text(size=10, face="bold"))
    if(logCov)
      p_cov_all <- p_cov_all + geom_line(aes(x=position, y= log10(coverage), color=sample), alpha=0.8, size=1)
    else
      p_cov_all <- p_cov_all + geom_line(aes(x=position, y= (coverage), color=sample), alpha=0.8, size=1)
    
    
  # only change scale of y axis to log10
  # if(logCov) 
  #   p_cov_all <- p_cov_all + scale_y_continuous(labels = scales::comma, limits = c(0,NA), trans=pseudo_log_trans(base=10))
  #     
  return(p_cov_all)
}

plotLogMutFreqAll_1d <- function(coverage1d_all.df) {
  p_mut_all_log<-
    ggplot(data=coverage1d_all.df, aes(x=factor(sample), y= log10(mutRate), fill=sample)) +
    geom_boxplot(alpha=0.8, show.legend = F, na.rm = T)+
    labs(x="sample", y=expression(log[10](mutRate)))+
    ylim(-6, 0) +
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          legend.text = element_text(size=12), 
          legend.title = element_text(size=12, face="bold"))
  return(p_mut_all_log)
}

plotEntropy_1d <- function(coverage1d_all.df) {
  p_entropy<-
    ggplot(data=coverage1d_all.df[coverage1d_all.df$coverage>0,]) +
    geom_line(aes(x=position, y= entropy, col=sample), alpha=0.8)+
    labs(x="sample", y=expression(entropy))+
    #xlim(0,300)+
    #ylim(0,0.3)
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          legend.text = element_text(size=12), 
          legend.title = element_text(size=12, face="bold"))
  return(p_entropy)
}

#### 2d plots #######

plotCoverageAll_2d <- function(coverage2d_all.df, logCov=F) {
  p_cov_all_2d = ggplot(coverage2d_all.df) +
    geom_ribbon(aes(x=position, ymin = (cov75), ymax = (maxCov), fill=sample), alpha=.2) +
    geom_line(aes(x=position, y=(maxCov), col=sample), size=1, alpha=0.8) + 
    ylab("coverage (max to 75th percentile)") +
    scale_y_continuous(labels = scales::comma, limits = c(0,NA))+
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=12), 
          legend.title = element_text(size=12, face="bold"))
  # only change scale of y axis to log10
  if(logCov) 
    p_cov_all_2d <- p_cov_all_2d + scale_y_continuous(trans=pseudo_log_trans(base=10), breaks =pretty_breaks())
  
  return(p_cov_all_2d)
}


########## Multiplot from Cookbook for R
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

