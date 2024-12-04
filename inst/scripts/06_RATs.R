library(data.table)
library(DRIMSeq)
library(here)
library(rats)

d <- readRDS("tutorials/06_differential_transcript_usage/www/dmDSdata.rds")

samples(d)
counts(d)

countA <- data.table(counts(d)[,c(2,which(samples(d)$group=="G02")+2)])
countB <- data.table(counts(d)[,c(2,which(samples(d)$group=="G03")+2)])

tx2gene <- counts(d)[,2:1]
names(tx2gene) <- c("target_id","parent_id")

st_count <- system.time({
  rats_count <- call_DTU(annot=tx2gene, count_data_A=countA, 
                         count_data_B=countB, qboot=FALSE)
})

plot_gene(rats_count,"ENSG00000117682",style="bycondition")

plot_gene(rats_count,"ENSG00000117682",style="byisoform")
