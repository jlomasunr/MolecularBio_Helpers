library(reticulate)
library(dplyr)
library(ggplot2)
library(data.table)

setwd("/Users/jslomas/Box/Cushman Lab/SynCAM/SynCAM_redesign/MolecularBio_Helpers")

use_condaenv(condaenv = 'dnachisel')
source_python("codonUsageCalculations.py")

translateSeqToRSCU("ATGTGTGAGAGCTTCAAGAGAGAATACCAGCTCTGCGAGGAGATCGGTCGAGGCCGCTTCGGCGTCGTCTACCGCTGCTACAACCCCTCCTCCACCGAAGATTCCGACACTCCTTTGGCTGTCAAGTCCATTGATAAGCGCCTCCTCCTCGACGACGAGACCGACCGTGAGTGCCTTGACAAGGAACCCAAAATCCTCCACCTTCTCTCTCCTCACCCTAATATTCTGCAGATCCATAACCTCTTTGATTCTGATACCCATCTCCTCATTGTTACCGATCTCTGCCAAGAGGAGACCCTTTACGAACGCATTATCTCTAACGGCCCCTTCTCTGAGCCTGACGCTGCGGCTATTTTCTGTCAGCTCGCCGAGGCCCTCGCTCACTGCCACCGCAACTACGTCGCTCACCGTGATATCAAACCAGACAACATACTATTCGACTCAAGGAACAGGCTGAAGCTCTGCGACTTCGGCTCGGCCGAGTGGTTCGGAGCAGGAGACAGAGAGATGCGCGGTGTGGTAGGAACACCCTACTACGTGGCGCCGGAGGTGCTTTCCGGAAAAGACTACAATGAGAAGGCGGATGTGTGGAGTGCTGGTGTTATTCTCTACATCATGCTCGGTGGAGTTCCACCCTTCTATGGTGAGACTGTTGAGGAGACTTTTGAGGCTGTTCTTAGGGGAAATCTTCGCTTTCCGGCCAGGATTTTCCGAAACGTTTCTACTCAGGCAAGGGATTTGTTGAGGAAGATGATGTGCAAGGATGTTTCCAGAAGGTTTTCTGCTGAACAAGTCTTAAGGCATCCTTGGGTAACCAGTGGAGGATTGGCCAACATGTAA", "../Athaliana_167_TAIR10.cds_primaryTranscriptOnly.rscu")

data = read.csv("/Users/jslomas/Box/Cushman Lab/SynCAM/SynCAM_redesign/Actin-Ath_Os.csv")
data <- data %>% select(Ath, Os)

cor(as.vector(data$Ath), as.vector(data$Os), method = "pearson")
cor.test(data$Ath, data$Os, method=c("pearson"))

#data$Os = data$Os*-1
# data <- data.frame(a, o)
data <- cbind(Location = rownames(data), data)
rownames(data) <- 1:nrow(data)
setDT(data)
data <- melt(data = data,
                  id.vars = c("Location"),
                  variable.name = "Sequence",
                  value.name = "RSCU")
data$Location <- as.numeric(data$Location)

ggplot(data, aes(x=Location, y=RSCU, fill=Sequence, alpha=0.5)) +
  geom_smooth() +
  geom_point() +
  #geom_bar(stat = "identity", width=1) +
  #geom_line(group = "Sequence")
  facet_wrap(vars(Sequence))


