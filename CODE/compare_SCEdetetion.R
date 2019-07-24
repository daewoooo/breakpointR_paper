## Load required libraries
library(breakpointR)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(scales)

## Chromosomes to analyze
chroms <- paste0('chr', c(1:22, 'X'))

## Run breakpointR for 200kb window ##
######################################
ptm <- proc.time()
breakpointr(inputfolder = "/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac01_back01/", 
            outputfolder = "BreakpointR_results_200kb", windowsize = 200000, binMethod = 'size', 
            pairedEndReads = FALSE, chromosomes = chroms, min.mapq = 10, filtAlt = TRUE, background = 0.05, minReads = 100, conf = 0.95)
proc.time() - ptm

ptm <- proc.time()
breakpointr(inputfolder = "/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac025_back01/", 
            outputfolder = "BreakpointR_results_200kb", windowsize = 200000, binMethod = 'size', 
            pairedEndReads = FALSE, chromosomes = chroms, min.mapq = 10, filtAlt = TRUE, background = 0.05, minReads = 100, conf = 0.95)
proc.time() - ptm

ptm <- proc.time()
breakpointr(inputfolder = "/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac05_back01/", 
            outputfolder = "BreakpointR_results_200kb", windowsize = 200000, binMethod = 'size', 
            pairedEndReads = FALSE, chromosomes = chroms, min.mapq = 10, filtAlt = TRUE, background = 0.05, minReads = 100, conf = 0.95)
proc.time() - ptm

ptm <- proc.time()
breakpointr(inputfolder = "/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac075_back01/", 
            outputfolder = "BreakpointR_results_200kb", windowsize = 200000, binMethod = 'size', 
            pairedEndReads = FALSE, chromosomes = chroms, min.mapq = 10, filtAlt = TRUE, background = 0.05, minReads = 100, conf = 0.95)
proc.time() - ptm

## Calculate stat for 200kb window in BAIT ##
#############################################
setwd("/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac01_back01/")
COV01_BACK01 <- compareResults(simul.breaks = "COV01_BACK01_SCE_list.txt", 
                  bait.sce = "Summary_files_BAIT/BAIT_TOTAL_LIBRARY_SCE.BED", 
                  breakpointr.sce = "BreakpointR_results_200kb/breakpoints/breakPointSummary.txt", 
                  index = "COV01_BACK01")

setwd("/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac025_back01/")
COV025_BACK01 <- compareResults(simul.breaks = "COV025_BACK01_SCE_list.txt", 
                               bait.sce = "Summary_files_BAIT/BAIT_TOTAL_LIBRARY_SCE.BED", 
                               breakpointr.sce = "BreakpointR_results_200kb/breakpoints/breakPointSummary.txt", 
                               index = "COV025_BACK01")

setwd("/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac05_back01/")
COV05_BACK01 <- compareResults(simul.breaks = "COV05_BACK01_SCE_list.txt", 
                                bait.sce = "Summary_files_BAIT/BAIT_TOTAL_LIBRARY_SCE.BED", 
                                breakpointr.sce = "BreakpointR_results_200kb/breakpoints/breakPointSummary.txt", 
                                index = "COV05_BACK01")

setwd("/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac075_back01/")
COV075_BACK01 <- compareResults(simul.breaks = "COV075_BACK01_SCE_list.txt", 
                               bait.sce = "Summary_files_BAIT/BAIT_TOTAL_LIBRARY_SCE.BED", 
                               breakpointr.sce = "BreakpointR_results_200kb/breakpoints/breakPointSummary.txt", 
                               index = "COV075_BACK01")

## Calculate summary statistics (specificity) for 200kb window ##
#################################################################
all.libs.stat <- rbind(COV01_BACK01$specificity, 
                       COV025_BACK01$specificity, 
                       COV05_BACK01$specificity, 
                       COV075_BACK01$specificity)

all.libs.stat.long <- melt(all.libs.stat, id.vars = c('index','ID'))
plt.df <- all.libs.stat.long %>% group_by(variable, ID) %>% summarise(count=sum(value))
plt.df$perc <- (plt.df$count / max(plt.df$count)) * 100
## Save plotted table
destination <- "/home/porubsky/WORK/BreakpointR_paper/Revisions/summary_BAIT_breakpointR_comparison_200kb.csv"
write.table(plt.df, file = destination, quote = FALSE, row.names = FALSE, sep=",")

## Get ROC curve counts ##
##########################
roc.df <- plt.df[plt.df$ID %in% c('BAIT', 'breakPR', 'breakPR.CI'),]
roc.df <- roc.df[roc.df$variable != 'total.sce',]
total.counts <- plt.df[plt.df$variable == 'total.sce',]
roc.df$TP <- roc.df$count
roc.df$FN <- total.counts$count[total.counts$ID == 'simul'] - roc.df$TP
roc.df$FP <- 0
roc.df$FP[roc.df$ID == 'BAIT'] <- total.counts$count[total.counts$ID == 'BAIT'] - roc.df$count[roc.df$ID == 'BAIT']
roc.df$FP[roc.df$ID == 'breakPR'] <- total.counts$count[total.counts$ID == 'breakPR'] - roc.df$count[roc.df$ID == 'breakPR'] 
roc.df$FP[roc.df$ID == 'breakPR.CI'] <- total.counts$count[total.counts$ID == 'breakPR.CI'] - roc.df$count[roc.df$ID == 'breakPR.CI'] 
roc.df$TN <- 0
roc.df$TN[roc.df$ID == 'BAIT'] <- total.counts$count[total.counts$ID == 'simul'] - total.counts$count[total.counts$ID == 'BAIT']
roc.df$TN[roc.df$ID == 'breakPR'] <- total.counts$count[total.counts$ID == 'simul'] - total.counts$count[total.counts$ID == 'breakPR'] 
roc.df$TN[roc.df$ID == 'breakPR.CI'] <- total.counts$count[total.counts$ID == 'simul'] - total.counts$count[total.counts$ID == 'breakPR.CI']

roc.df$TPR <- roc.df$TP / (roc.df$TP + roc.df$FN)
roc.df$FPR <- 1 - (roc.df$FP / (roc.df$TN + roc.df$FP))
roc.df$recall <- roc.df$TP / (roc.df$TP + roc.df$FN)
roc.df$precision <- roc.df$TP / (roc.df$TP + roc.df$FP)

roc.200kb <- ggplot(roc.df) + geom_line(aes(x=recall, y=precision, group=ID, color=ID), size=1) + 
  geom_point(aes(x=recall, y=precision, group=variable, shape=variable, color=ID, size=0.5)) + 
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_shape_manual(values = c(15:18)) + 
  scale_color_manual(values = brewer.pal(n = 4, name = 'Set2')) +
  theme_bw() + ggtitle("200kb window size")

plt1 <- ggplot() + geom_col(data=plt.df, aes(x=ID, y=perc, fill=ID)) +
  facet_wrap(variable ~ ., nrow = 1) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) +
  scale_fill_manual(values = brewer.pal(n = 4, name = 'Set2'), guide='none') +
  xlab("") +
  ylab("# of detected SCEs")

## Calculate summary statistics (resolution)
BAIT.dists <- c(COV01_BACK01$BAIT.dists, COV025_BACK01$BAIT.dists, COV05_BACK01$BAIT.dists, COV075_BACK01$BAIT.dists)
BreakPR.dists <- c(COV01_BACK01$BreakPR.dists, COV025_BACK01$BreakPR.dists, COV05_BACK01$BreakPR.dists, COV075_BACK01$BreakPR.dists)
BreakPR.CI.dists <- c(COV01_BACK01$BreakPR.CI.dists, COV025_BACK01$BreakPR.CI.dists, COV05_BACK01$BreakPR.CI.dists, COV075_BACK01$BreakPR.CI.dists)
BAIT.dists.mean <- mean(BAIT.dists)
BreakPR.dists.mean <- mean(BreakPR.dists)
BreakPR.CI.dists.mean <- mean(BreakPR.CI.dists)
BAIT.widths <- c(COV01_BACK01$BAIT.width, COV025_BACK01$BAIT.width, COV05_BACK01$BAIT.width, COV075_BACK01$BAIT.width)
BreakPR.widths <- c(COV01_BACK01$BreakPR.width, COV025_BACK01$BreakPR.width, COV05_BACK01$BreakPR.width, COV075_BACK01$BreakPR.width)
BreakPR.CI.widths <- c(COV01_BACK01$BreakPR.CI.width, COV025_BACK01$BreakPR.CI.width, COV05_BACK01$BreakPR.CI.width, COV075_BACK01$BreakPR.CI.width)
BAIT.widths.mean <- mean(BAIT.widths)
BreakPR.widths.mean <- mean(BreakPR.widths)
BreakPR.CI.widths.mean <- mean(BreakPR.CI.widths)

plt.df <- data.frame(BAIT.dists.mean=BAIT.dists.mean, BreakPR.dists.mean=BreakPR.dists.mean, BreakPR.CI.dists.mean=BreakPR.CI.dists.mean)
plt.df <- melt(plt.df)
plt.df$ID <- c('BAIT', 'BreakPR', 'BreakPR.CI')
plt2 <- ggplot() + geom_col(data=plt.df, aes(x=ID, y=value, fill=ID)) +
  scale_fill_manual(values = brewer.pal(n = 6, name = 'Set2'), guide='none') +
  scale_y_continuous(labels = comma) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust = 1)) +
  ylab("Distance to simulated breakpoint (bp)")

BAIT.widths <- c(median(COV01_BACK01$BAIT.width), median(COV025_BACK01$BAIT.width), median(COV05_BACK01$BAIT.width), median(COV075_BACK01$BAIT.width))
BreakPR.widths <- c(median(COV01_BACK01$BreakPR.width), median(COV025_BACK01$BreakPR.width), median(COV05_BACK01$BreakPR.width), median(COV075_BACK01$BreakPR.width))
BreakPR.CI.widths <- c(median(COV01_BACK01$BreakPR.CI.width), median(COV025_BACK01$BreakPR.CI.width), median(COV05_BACK01$BreakPR.CI.width), median(COV075_BACK01$BreakPR.CI.width))
plt.df <- data.frame(BAIT.widths=BAIT.widths, BreakPR.widths=BreakPR.widths, BreakPR.CI.widths=BreakPR.CI.widths)
plt.df <- melt(plt.df)
plt.df$ID <- factor(rep(c('BAIT', 'BreakPR', 'BreakPR.CI'), each=4))
plt.df$categ <- rep(c('COV01','COV025','COV05','COV075'), 3)

plt3 <- ggplot() + geom_col(data=plt.df, aes(x=categ, y=value, fill=ID), position=position_dodge()) +
  #theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) +
  scale_fill_manual(values = brewer.pal(n = 6, name = 'Set2')) +
  scale_y_continuous(labels = comma) +
  xlab("") +
  theme_bw() +
  ylab("Median breakpoint width (bp)")
plt.final.200kb <- plot_grid(plt2, plt3, nrow = 1, rel_widths = c(1,3), align = 'h', axis = 'b')


## Run breakpointR for 1Mb window ##
####################################
ptm <- proc.time()
breakpointr(inputfolder = "/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac01_back01/", 
            outputfolder = "BreakpointR_results_1Mb", windowsize = 1000000, binMethod = 'size', 
            pairedEndReads = FALSE, chromosomes = chroms, min.mapq = 10, filtAlt = TRUE, background = 0.05, minReads = 100, conf = 0.95)
proc.time() - ptm

ptm <- proc.time()
breakpointr(inputfolder = "/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac025_back01/", 
            outputfolder = "BreakpointR_results_1Mb", windowsize = 1000000, binMethod = 'size', 
            pairedEndReads = FALSE, chromosomes = chroms, min.mapq = 10, filtAlt = TRUE, background = 0.05, minReads = 100, conf = 0.95)
proc.time() - ptm

ptm <- proc.time()
breakpointr(inputfolder = "/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac05_back01/", 
            outputfolder = "BreakpointR_results_1Mb", windowsize = 1000000, binMethod = 'size', 
            pairedEndReads = FALSE, chromosomes = chroms, min.mapq = 10, filtAlt = TRUE, background = 0.05, minReads = 100, conf = 0.95)
proc.time() - ptm

ptm <- proc.time()
breakpointr(inputfolder = "/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac075_back01/", 
            outputfolder = "BreakpointR_results_1Mb", windowsize = 1000000, binMethod = 'size', 
            pairedEndReads = FALSE, chromosomes = chroms, min.mapq = 10, filtAlt = TRUE, background = 0.05, minReads = 100, conf = 0.95)
proc.time() - ptm

## Calculate stat for 1Mb window in BAIT ##
###########################################
setwd("/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac01_back01/")
COV01_BACK01 <- compareResults(simul.breaks = "COV01_BACK01_SCE_list.txt", 
                               bait.sce = "BAIT_1Mb/Summary_files_BAIT/BAIT_TOTAL_LIBRARY_SCE.BED", 
                               breakpointr.sce = "BreakpointR_results_1Mb/breakpoints/breakPointSummary.txt", 
                               index = "COV01_BACK01")

setwd("/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac025_back01/")
COV025_BACK01 <- compareResults(simul.breaks = "COV025_BACK01_SCE_list.txt", 
                                bait.sce = "BAIT_1Mb/Summary_files_BAIT/BAIT_TOTAL_LIBRARY_SCE.BED", 
                                breakpointr.sce = "BreakpointR_results_1Mb/breakpoints/breakPointSummary.txt", 
                                index = "COV025_BACK01")

setwd("/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac05_back01/")
COV05_BACK01 <- compareResults(simul.breaks = "COV05_BACK01_SCE_list.txt", 
                               bait.sce = "BAIT_1Mb/Summary_files_BAIT/BAIT_TOTAL_LIBRARY_SCE.BED", 
                               breakpointr.sce = "BreakpointR_results_1Mb/breakpoints/breakPointSummary.txt", 
                               index = "COV05_BACK01")

setwd("/home/porubsky/WORK/BreakpointR_paper/Revisions/simulations_GRCh38/BAMs_frac075_back01/")
COV075_BACK01 <- compareResults(simul.breaks = "COV075_BACK01_SCE_list.txt", 
                                bait.sce = "BAIT_1Mb/Summary_files_BAIT/BAIT_TOTAL_LIBRARY_SCE.BED", 
                                breakpointr.sce = "BreakpointR_results_1Mb/breakpoints/breakPointSummary.txt", 
                                index = "COV075_BACK01")

## Calculate summary statistics (specificity)
all.libs.stat <- rbind(COV01_BACK01$specificity, 
                       COV025_BACK01$specificity, 
                       COV05_BACK01$specificity, 
                       COV075_BACK01$specificity)

all.libs.stat.long <- melt(all.libs.stat, id.vars = c('index','ID'))
plt.df <- all.libs.stat.long %>% group_by(variable, ID) %>% summarise(count=sum(value))
plt.df$perc <- (plt.df$count / max(plt.df$count)) * 100
## Save plotted table
destination <- "/home/porubsky/WORK/BreakpointR_paper/Revisions/summary_BAIT_breakpointR_comparison_1Mb.csv"
write.table(plt.df, file = destination, quote = FALSE, row.names = FALSE, sep = ",")

## Get ROC curve counts ##
##########################
roc.df <- plt.df[plt.df$ID %in% c('BAIT', 'breakPR', 'breakPR.CI'),]
roc.df <- roc.df[roc.df$variable != 'total.sce',]
total.counts <- plt.df[plt.df$variable == 'total.sce',]
roc.df$TP <- roc.df$count
roc.df$FN <- total.counts$count[total.counts$ID == 'simul'] - roc.df$TP
roc.df$FP <- 0
roc.df$FP[roc.df$ID == 'BAIT'] <- total.counts$count[total.counts$ID == 'BAIT'] - roc.df$count[roc.df$ID == 'BAIT']
roc.df$FP[roc.df$ID == 'breakPR'] <- total.counts$count[total.counts$ID == 'breakPR'] - roc.df$count[roc.df$ID == 'breakPR'] 
roc.df$FP[roc.df$ID == 'breakPR.CI'] <- total.counts$count[total.counts$ID == 'breakPR.CI'] - roc.df$count[roc.df$ID == 'breakPR.CI'] 
roc.df$TN <- 0
roc.df$TN[roc.df$ID == 'BAIT'] <- total.counts$count[total.counts$ID == 'simul'] - total.counts$count[total.counts$ID == 'BAIT']
roc.df$TN[roc.df$ID == 'breakPR'] <- total.counts$count[total.counts$ID == 'simul'] - total.counts$count[total.counts$ID == 'breakPR'] 
roc.df$TN[roc.df$ID == 'breakPR.CI'] <- total.counts$count[total.counts$ID == 'simul'] - total.counts$count[total.counts$ID == 'breakPR.CI']

roc.df$TPR <- roc.df$TP / (roc.df$TP + roc.df$FN)
roc.df$FPR <- 1 - (roc.df$FP / (roc.df$TN + roc.df$FP))
roc.df$recall <- roc.df$TP / (roc.df$TP + roc.df$FN)
roc.df$precision <- roc.df$TP / (roc.df$TP + roc.df$FP)

roc.1Mb <- ggplot(roc.df) + geom_line(aes(x=recall, y=precision, group=ID, color=ID), size=1) + 
  geom_point(aes(x=recall, y=precision, group=variable, shape=variable, color=ID, size=0.5)) + 
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_shape_manual(values = c(15:18)) + 
  scale_color_manual(values = brewer.pal(n = 4, name = 'Set2')) +
  theme_bw() + ggtitle("1Mb window size")

plt1 <- ggplot() + geom_col(data=plt.df, aes(x=ID, y=perc, fill=ID)) +
  facet_wrap(variable ~ ., nrow = 1) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) +
  scale_fill_manual(values = brewer.pal(n = 4, name = 'Set2'), guide='none') +
  xlab("") +
  ylab("# of detected SCEs")

plt1 <- ggplot() + geom_col(data=plt.df, aes(x=ID, y=perc, fill=ID)) +
  facet_wrap(variable ~ ., nrow = 1) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) +
  scale_fill_manual(values = brewer.pal(n = 4, name = 'Set2'), guide='none') +
  xlab("") +
  ylab("# of detected SCEs")

## Calculate summary statistics (resolution)
BAIT.dists <- c(COV01_BACK01$BAIT.dists, COV025_BACK01$BAIT.dists, COV05_BACK01$BAIT.dists, COV075_BACK01$BAIT.dists)
BreakPR.dists <- c(COV01_BACK01$BreakPR.dists, COV025_BACK01$BreakPR.dists, COV05_BACK01$BreakPR.dists, COV075_BACK01$BreakPR.dists)
BreakPR.CI.dists <- c(COV01_BACK01$BreakPR.CI.dists, COV025_BACK01$BreakPR.CI.dists, COV05_BACK01$BreakPR.CI.dists, COV075_BACK01$BreakPR.CI.dists)
BAIT.dists.mean <- mean(BAIT.dists)
BreakPR.dists.mean <- mean(BreakPR.dists)
BreakPR.CI.dists.mean <- mean(BreakPR.CI.dists)
BAIT.widths <- c(COV01_BACK01$BAIT.width, COV025_BACK01$BAIT.width, COV05_BACK01$BAIT.width, COV075_BACK01$BAIT.width)
BreakPR.widths <- c(COV01_BACK01$BreakPR.width, COV025_BACK01$BreakPR.width, COV05_BACK01$BreakPR.width, COV075_BACK01$BreakPR.width)
BreakPR.CI.widths <- c(COV01_BACK01$BreakPR.CI.width, COV025_BACK01$BreakPR.CI.width, COV05_BACK01$BreakPR.CI.width, COV075_BACK01$BreakPR.CI.width)
BAIT.widths.mean <- mean(BAIT.widths)
BreakPR.widths.mean <- mean(BreakPR.widths)
BreakPR.CI.widths.mean <- mean(BreakPR.CI.widths)

plt.df <- data.frame(BAIT.dists.mean=BAIT.dists.mean, BreakPR.dists.mean=BreakPR.dists.mean, BreakPR.CI.dists.mean=BreakPR.CI.dists.mean)
plt.df <- melt(plt.df)
plt.df$ID <- c('BAIT', 'BreakPR', 'BreakPR.CI')
plt2 <- ggplot() + geom_col(data=plt.df, aes(x=ID, y=value, fill=ID)) +
  scale_fill_manual(values = brewer.pal(n = 6, name = 'Set2'), guide='none') +
  scale_y_continuous(labels = comma) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust = 1)) +
  ylab("Distance to simulated breakpoint (bp)")

BAIT.widths <- c(median(COV01_BACK01$BAIT.width), median(COV025_BACK01$BAIT.width), median(COV05_BACK01$BAIT.width), median(COV075_BACK01$BAIT.width))
BreakPR.widths <- c(median(COV01_BACK01$BreakPR.width), median(COV025_BACK01$BreakPR.width), median(COV05_BACK01$BreakPR.width), median(COV075_BACK01$BreakPR.width))
BreakPR.CI.widths <- c(median(COV01_BACK01$BreakPR.CI.width), median(COV025_BACK01$BreakPR.CI.width), median(COV05_BACK01$BreakPR.CI.width), median(COV075_BACK01$BreakPR.CI.width))
plt.df <- data.frame(BAIT.widths=BAIT.widths, BreakPR.widths=BreakPR.widths, BreakPR.CI.widths=BreakPR.CI.widths)
plt.df <- melt(plt.df)
plt.df$ID <- factor(rep(c('BAIT', 'BreakPR', 'BreakPR.CI'), each=4))
plt.df$categ <- rep(c('COV01','COV025','COV05','COV075'), 3)

plt3 <- ggplot() + geom_col(data=plt.df, aes(x=categ, y=value, fill=ID), position=position_dodge()) +
  #theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) +
  scale_fill_manual(values = brewer.pal(n = 6, name = 'Set2')) +
  scale_y_continuous(labels = comma) +
  xlab("") +
  theme_bw() +
  ylab("Median breakpoint width (bp)")
plt.final.1Mb <- plot_grid(plt2, plt3, nrow = 1, rel_widths = c(1,3), align = 'h', axis = 'b')

## Create final summary plot
#final.plt <- plot_grid(plt1, plt2, ncol = 1, rel_heights = c(2,1))
roc.plts <- plot_grid(roc.200kb, roc.1Mb, nrow = 1)
final.plt <- plot_grid(plt.final.200kb, plt.final.1Mb, roc.plts, ncol = 1, rel_heights = c(2,2,2.5))
destination <- "/home/porubsky/WORK/BreakpointR_paper/Revisions/summary_BAIT_breakpointR_comparison.pdf"
ggsave(filename = destination, plot = final.plt, device = 'pdf', width = 10, height = 10, useDingbats=FALSE)
destination <- "/home/porubsky/WORK/BreakpointR_paper/Revisions/summary_BAIT_breakpointR_comparison.png"
ggsave(filename = destination, plot = final.plt, width = 10, height = 10, units = 'in', limitsize = FALSE, dpi = 300)

#================================================================================================================================
## Helper function

#' This fuction evaluates breakpoint detection performance between BAIT and breakpointR
#'
#' @param simul.breaks A path to simulted set of breakpoints (SCEs).
#' @param bait.sce A path to detected set of breakpoints (SCEs) by BAIT.
#' @param breakpointr.sce A path to detected set of breakpoints (SCEs) by breakpoiintR.
#' @param index An user defined ID.
#' @return A \code(list} object.
#' @author David Porubsky

compareResults <- function(simul.breaks=NULL, bait.sce=NULL, breakpointr.sce=NULL, index=NULL) {
  ## Read in simulated breakpoints
  simul.breaks <- read.table(simul.breaks, header = TRUE)
  simul.breaks.gr <- GRanges(seqnames=simul.breaks$Chr, ranges=IRanges(start=simul.breaks$Pos, end=simul.breaks$Pos+1), ID=simul.breaks$Filename, method='simul', level=1)

  ## Read in BAIT results
  bait.sce <- read.table(bait.sce, header = FALSE)
  #ID <- paste0('COV01_BACK01_', bait.sce$V4)
  ID <- gsub(bait.sce$V4, pattern = '_srt_dedup', replacement = '')
  bait.sce.gr <- GRanges(seqnames=bait.sce$V1, ranges=IRanges(start=bait.sce$V2, end=bait.sce$V3), ID=ID, method='BAIT', level=2)

  ## Read breakpointR results
  breakpointr.sce <- read.table(breakpointr.sce, header = TRUE)
  breakpointr.sce$filenames <- gsub(breakpointr.sce$filenames, pattern = '_srt_dedup.bam.RData.\\d', replacement = '')
  breakpointr.sce.gr <- GRanges(seqnames=breakpointr.sce$seqnames, ranges=IRanges(start = breakpointr.sce$start, end=breakpointr.sce$end), ID=breakpointr.sce$filenames, method='BreakPR', level=3)
  breakpointr.sce.CI.gr <- GRanges(seqnames=breakpointr.sce$seqnames, ranges=IRanges(start = breakpointr.sce$CI.start, end=breakpointr.sce$CI.end), ID=breakpointr.sce$filenames, method='BreakPR.CI', level=4)

  ## Merge all ranges
  ranges.all <- c(simul.breaks.gr, bait.sce.gr, breakpointr.sce.gr, breakpointr.sce.CI.gr)

  ranges.all.df <- as.data.frame(ranges.all)
  plt <- ggplot() + geom_rect(data=ranges.all.df, aes(xmin=start, xmax=end, ymin=0+level, ymax=1+level, fill=method, color=method)) +
    facet_wrap(ID ~ seqnames, scales='free')

  ## Calculate distance of predicted breakpoints to a simulated breakpoint
  BAIT.dists <- list()
  BreakPR.dists <- list()
  BreakPR.CI.dists <- list()
  ranges.all.grl <- split(ranges.all, ranges.all$ID)
  for (i in seq_along(ranges.all.grl)) {
    chr.SCEs <- split(ranges.all.grl[[i]], seqnames(ranges.all.grl[[i]]))
    for (j in seq_along(chr.SCEs)) {
      chr.break <- chr.SCEs[[j]]
      if ('BAIT' %in% chr.break$method) {
        dist.BAIT <-  distanceToNearest(chr.break[chr.break$method == 'BAIT'], chr.break[chr.break$method == 'simul'])
        BAIT.dists[[length(BAIT.dists) + 1]] <- dist.BAIT@elementMetadata$distance
      }

      if ('BreakPR' %in% chr.break$method) {
        dist.BreakPR <-  distanceToNearest(chr.break[chr.break$method == 'BreakPR'], chr.break[chr.break$method == 'simul'])
        BreakPR.dists[[length(BreakPR.dists) + 1]] <- dist.BreakPR@elementMetadata$distance
      }
      
      if ('BreakPR.CI' %in% chr.break$method) {
        dist.BreakPR <-  distanceToNearest(chr.break[chr.break$method == 'BreakPR.CI'], chr.break[chr.break$method == 'simul'])
        BreakPR.CI.dists[[length(BreakPR.CI.dists) + 1]] <- dist.BreakPR@elementMetadata$distance
      }
    }
  }
  BAIT.dists <- unlist(BAIT.dists, use.names = FALSE)
  BreakPR.dists <- unlist(BreakPR.dists, use.names = FALSE)
  BreakPR.CI.dists <- unlist(BreakPR.CI.dists, use.names = FALSE)
  
  BAIT.medDist <- median(BAIT.dists)
  BreakPR.medDist <- median(BreakPR.dists)
  BreakPR.CI.medDist <- median(BreakPR.CI.dists)
  BAIT.meanDist <- mean(BAIT.dists)
  BreakPR.meanDist <- mean(BreakPR.dists)
  BreakPR.CI.meanDist <- mean(BreakPR.CI.dists)
  BAIT.medWidth <- median(width(bait.sce.gr))
  BreakPR.medWidth <- median(width(breakpointr.sce.gr))
  BreakPR.CI.medWidth <- median(width(breakpointr.sce.CI.gr))

  breakpointr.sce.found <- length(subsetByOverlaps(breakpointr.sce.gr, simul.breaks.gr))
  breakpointr.sceCI.found <- length(subsetByOverlaps(breakpointr.sce.CI.gr, simul.breaks.gr))
  bait.sce.found <- length(subsetByOverlaps(bait.sce.gr, simul.breaks.gr))
  
  breakpointr.sce.found.10kb <- length(subsetByOverlaps(breakpointr.sce.gr, simul.breaks.gr, maxgap = 10000))
  breakpointr.sceCI.found.10kb <- length(subsetByOverlaps(breakpointr.sce.CI.gr, simul.breaks.gr, maxgap = 10000))
  bait.sce.found.10kb <- length(subsetByOverlaps(bait.sce.gr, simul.breaks.gr, maxgap = 10000))
  
  breakpointr.sce.found.100kb <- length(subsetByOverlaps(breakpointr.sce.gr, simul.breaks.gr, maxgap = 100000))
  breakpointr.sceCI.found.100kb <- length(subsetByOverlaps(breakpointr.sce.CI.gr, simul.breaks.gr, maxgap = 100000))
  bait.sce.found.100kb <- length(subsetByOverlaps(bait.sce.gr, simul.breaks.gr, maxgap = 100000))
  
  breakpointr.sce.found.500kb <- length(subsetByOverlaps(breakpointr.sce.gr, simul.breaks.gr, maxgap = 500000))
  breakpointr.sceCI.found.500kb <- length(subsetByOverlaps(breakpointr.sce.CI.gr, simul.breaks.gr, maxgap = 500000))
  bait.sce.found.500kb <- length(subsetByOverlaps(bait.sce.gr, simul.breaks.gr, maxgap = 500000))
  
  ## Compile summary tables for export
  resolution.df <- data.frame(median.breakpoint.width=c(BAIT.medWidth, BreakPR.medWidth, BreakPR.CI.medWidth, 0),
                              median.distance2simulSCE=c(BAIT.medDist, BreakPR.medDist, BreakPR.CI.medDist, 0), 
                              mean.distance2simulSCE=c(BAIT.meanDist, BreakPR.meanDist, BreakPR.CI.meanDist, 0),
                              ID=c('BAIT', 'breakPR', 'breakPR.CI', 'simul'), index=index)
  
  specificity.df <- data.frame(overlap.direct=c(bait.sce.found, breakpointr.sce.found, breakpointr.sceCI.found, length(simul.breaks.gr)),
                          overlap.10kb=c(bait.sce.found.10kb, breakpointr.sce.found.10kb, breakpointr.sceCI.found.10kb, length(simul.breaks.gr)),
                          overlap.100kb=c(bait.sce.found.100kb, breakpointr.sce.found.100kb, breakpointr.sceCI.found.100kb, length(simul.breaks.gr)),
                          overlap.500kb=c(bait.sce.found.500kb, breakpointr.sce.found.500kb, breakpointr.sceCI.found.500kb, length(simul.breaks.gr)),
                          total.sce=c(length(BAIT.dists), length(BreakPR.dists), length(BreakPR.CI.dists), length(simul.breaks.gr)),
                          ID=c('BAIT', 'breakPR', 'breakPR.CI', 'simul'), index=index)
  
  ## Return results
  return(list(brk.plot=plt, resolution=resolution.df, specificity=specificity.df, 
              BAIT.dists = BAIT.dists, BreakPR.dists = BreakPR.dists, BreakPR.CI.dists = BreakPR.CI.dists,
              BAIT.width = width(bait.sce.gr), BreakPR.width = width(breakpointr.sce.gr), BreakPR.CI.width = width(breakpointr.sce.CI.gr)))
}

