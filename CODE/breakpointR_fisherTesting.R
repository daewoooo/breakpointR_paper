## Load required libraries
library(ggplot2)
library(breakpointR)

## Parameter settings
## Calculate strand-state distribution for all possible combination of 1000 directional reads
c.reads <- 0:1000
w.reads <- 1000:0 
roi.reads <- c.reads + w.reads
test.backG <- c(0.01, 0.025, 0.05, 0.1, 0.2)

plots <- list()
comp <- list()
for (backG in test.backG) {
  bestFit.fisher <- list()
  bestFit.binom <- list()
  for (i1 in 1:length(c.reads)) {
      roi.reads <- c.reads[i1] + w.reads[i1]
      fisher.test <- genotype.fisher(cReads = c.reads[i1], wReads = w.reads[i1], roiReads = roi.reads, background = backG)
      bestFit.fisher[[i1]] <- data.frame(pvalue=fisher.test$pval, state=fisher.test$bestFit)
      binom.test <- countBinProb(minusCounts = w.reads[i1], plusCounts = c.reads[i1], alpha = backG, log = FALSE)
      best.fit <- which.max(binom.test)
      if (best.fit == 1) {
        best.fit.df <- data.frame(pvalue=binom.test[best.fit], state='ww')
      } else if (best.fit == 2) {
        best.fit.df <- data.frame(pvalue=binom.test[best.fit], state='cc')
      } else {
        best.fit.df <- data.frame(pvalue=binom.test[best.fit], state='wc')
      }
      bestFit.binom[[i1]] <- best.fit.df
  }
  plt.df1 <- do.call(rbind, bestFit.fisher)
  plt.df2 <- do.call(rbind, bestFit.binom)
  ## Compare strand assignment between fisher and binom test
  agree <-  plt.df1$state == plt.df2$state
  comp.df <- data.frame(c=0:1000, w=1000:0, back=backG, fisher=plt.df1$state, binom=plt.df2$state, agree=agree)
  comp[[length(plots) + 1]] <- comp.df
  
  plt1 <- ggplot(plt.df1) + 
    geom_linerange(aes(x=0:1000, ymin=1, ymax=pvalue, color=state), size=3) + 
    scale_color_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown")) + 
    annotate(geom = "text", x = c(0,250,500,750,1000), y = -0.1, label = c('1000','750','500','250','0'), size=4) + 
    xlab("WC reads counts") + 
    ylab("pValue") + 
    ggtitle(paste('Fisher: Background level', backG)) +
    theme_bw()
  
  plt2 <- ggplot(plt.df2) + 
    geom_linerange(aes(x=0:1000, ymin=1, ymax=pvalue, color=state), size=3) + 
    scale_color_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown")) + 
    annotate(geom = "text", x = c(0,250,500,750,1000), y = -0.1, label = c('1000','750','500','250','0'), size=4) + 
    xlab("WC reads counts") + 
    ylab("pValue") + 
    ggtitle(paste('Binom: Background level', backG)) +
    theme_bw()
  
  plt <- plot_grid(plt1, plt2, nrow = 1)
  plots[[length(plots) + 1]] <- plt
}  
## Calculate level of agreement between Fisher and binom test
comp.all <- do.call(rbind, comp)
comp.all[comp.all$agree == FALSE,]
nrow(comp.all[comp.all$agree == FALSE,]) / nrow(comp.all)

final.plt <- plot_grid(plotlist = plots, ncol = 1)
destination <- "/home/porubsky/WORK/BreakpointR_paper/Revisions/fisherVSbinom_test.pdf"
ggsave(filename = destination, plot = final.plt, device = 'pdf', width = 10, height = 12)
destination <- "/home/porubsky/WORK/BreakpointR_paper/Revisions/fisherVSbinom_test.png"
ggsave(filename = destination, plot = final.plt, width = 10, height = 12, units = 'in', limitsize = FALSE, dpi = 300)

#==================================================================================================================================
## breakpointR genotyping function
genotype.fisher <- function(cReads, wReads, roiReads, background=0.02, minReads=10) {
  ## FISHER EXACT TEST
  result <- list(bestFit=NA, pval=NA)
  if (length(roiReads)==0) {
    return(result)
  }
  if (is.na(roiReads)) {
    return(result)
  }
  if ( roiReads >= minReads ) {
    m <- matrix(c(cReads, wReads, round(roiReads*(1-background)), round(roiReads*background)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Cs','Ws')))
    CCpVal <- stats::fisher.test(m, alternative="greater")[[1]]
    m <- matrix(c(cReads, wReads, round(roiReads*0.5), round(roiReads*0.5)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Cs','Ws')))
    WCpVal <- 1 - stats::fisher.test(m, alternative="two.sided")[[1]]
    m <- matrix(c(wReads, cReads, round(roiReads*(1-background)), round(roiReads*background)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Ws','Cs')))
    WWpVal <- stats::fisher.test(m, alternative="greater")[[1]]
    
    pVal <- c(wc=WCpVal, cc=CCpVal, ww=WWpVal)
    result <- list(bestFit=names(pVal)[which.min(pVal)], pval=min(pVal))
    return(result)
    
  } else { 
    return(result)
  }
}

countBinProb <- function(minusCounts, plusCounts, alpha=0.1, log=FALSE) {
  
  sumCounts <- minusCounts + plusCounts
  #calculate that given region is WW
  prob.ww <- stats::dbinom(minusCounts, size = sumCounts, prob = 1-alpha, log = log)
  
  #calculate that given region is CC
  prob.cc <- stats::dbinom(minusCounts, size = sumCounts, prob = alpha, log = log)
  
  #calculate that given region is WC
  prob.mix <- stats::dbinom(minusCounts, size = sumCounts, prob = 0.5, log = log)
  
  prob.m <- cbind(prob.ww, prob.cc, prob.mix)
  
  return(prob.m)
}
