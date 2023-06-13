#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/REVAMP/outdir/PATH/blast_results/ASV_blastn_nt.btab" #POINT TO FILE
#args[2]<-"/REVAMP/outdir/PATH/blast_results/" #Point to outdir
#args[3]<-326 #blast length cutoff
########################################
library(dplyr)
setwd(as.character(args[2]))

all_results <- read.delim(as.character(args[1]), header=FALSE, stringsAsFactors=FALSE)
colnames(all_results) <- c("ASV", "percent", "length", "taxid", "accession")
all_results_trimLength <- all_results[all_results$length>=(as.numeric(args[3])), ]
all_results_trimLength$correction <- all_results_trimLength$percent/100 *all_results_trimLength$length
all_results_trimLength_filter <- all_results_trimLength %>% group_by(ASV) %>% filter(percent == max(percent))
all_results_bestHit <- all_results_trimLength_filter %>% group_by(ASV)%>%summarise_all(list(~toString(unique(.))))
write.table(all_results_bestHit, file='ASV_blastn_nt_formatted.txt', sep='\t', quote=FALSE, row.names=FALSE)
