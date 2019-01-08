#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library("rtracklayer")
library("dplyr")



args = commandArgs(trailingOnly=TRUE)


grupt_frqs.phial <-args[1]
frq_shifts.phial <- args[2]

#freqCompare.bg <- import.bedGraph(grupt_frqs.phial)
freqCompare.bg <- import.bedGraph(grupt_frqs.phial)
names(mcols(freqCompare.bg)) <-  c("score", "name", "blup", "ref", "alt", "parent1_count", "parent1_alt_af", "parent2_count", "parent2_alt_af","offspring_count", "offspring_alt_af")

freqCompare.bg$parent2_introg_deltaF <- (freqCompare.bg$offspring_alt_af - freqCompare.bg$parent2_alt_af) * sign(freqCompare.bg$parent2_alt_af-freqCompare.bg$parent1_alt_af)
freqCompare.bg$parent1_depletion_deltaF <- (freqCompare.bg$offspring_alt_af - freqCompare.bg$parent1_alt_af) * sign(freqCompare.bg$parent2_alt_af-freqCompare.bg$parent1_alt_af)




freqCompare.df <- data.frame(seqnames=seqnames(freqCompare.bg),
								   starts=start(freqCompare.bg)-1,
								   ends=end(freqCompare.bg),
								   names=c(rep(".", length(freqCompare.bg))),
								   score=c(rep(".", length(freqCompare.bg))),
								   strand=strand(freqCompare.bg)  )

freqCompare.df <- cbind(freqCompare.df, mcols(freqCompare.bg) %>% as.data.frame() %>% select(c("ref","alt","parent1_count","parent1_alt_af","parent2_count","parent2_alt_af","offspring_count","offspring_alt_af","parent2_introg_deltaF","parent1_depletion_deltaF")))


write.table(freqCompare.df, file=frq_shifts.phial, quote=F, sep="\t", row.names=F)#, col.names=F)



