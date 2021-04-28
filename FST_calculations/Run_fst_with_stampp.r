library("StAMPP")
library("adegenet")



CFRL18.genotype.df <- read.table("CFRL18_genotype.txt", header = TRUE, row.names = 1, check.names = FALSE)

CFRL18.binary.genotype.df <- data.frame(lapply(CFRL18.genotype.df, function(x){as.numeric(factor(x, levels = names(sort(-table(x)))))-1}), check.names =  FALSE)

row.names(CFRL18.binary.genotype.df) <- row.names(CFRL18.genotype.df)

CFRL18.LungGl <- new("genlight", as.matrix(CFRL18.binary.genotype.df))

CFRL18.pops <- read.table("CFRL18_sample_pops.txt", row.names = 1, stringsAsFactors = TRUE)

CFRL18.LungGl$pop <- CFRL18.pops$V2

CFRL18.FST <- stamppFst(CFRL18.LungGl)



CFRL18.fst_output <- data.frame(cbind("Pop1" = CFRL18.FST$Bootstraps$Population1, "pop2" = CFRL18.FST$Bootstraps$Population2, "Lower_CI" = CFRL18.FST$Bootstraps$`Lower bound CI limit`, "Upper_CI" = CFRL18.FST$Bootstraps$`Upper bound CI limit`, "P-value" = CFRL18.FST$Bootstraps$`p-value`, "Fst" = CFRL18.FST$Bootstraps$Fst))

write.csv(CFRL18.fst_output, "CFRL18_Fst_results.csv", row.names = FALSE)


