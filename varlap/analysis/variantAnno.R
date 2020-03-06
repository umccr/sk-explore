library(VariantAnnotation)

#read inputs
bcbio <- readVcf("/Users/kanwals/Documents/UMCCR/data/projects/varlap/0000.vcf", "hg38")
dragen <- readVcf("/Users/kanwals/Documents/UMCCR/data/projects/varlap/0001.vcf", "hg38")
both <- readVcf("/Users/kanwals/Documents/UMCCR/data/projects/varlap/0002.vcf", "hg38")

#names(info(bcbio))
#introduce a new variable in the INFO field
info(bcbio)$CLASS <- "bcbio"
info(dragen)$CLASS <- "dragen"
info(both)$CLASS <- "both"

#update header with the new variable description
info(header(bcbio)) <- rbind(info(header(bcbio)), data.frame(Number=1, Type="String", 
                                                             Description="Variant class - bcbio, dragen or both", row.names = "CLASS"))
info(header(dragen)) <- rbind(info(header(dragen)), data.frame(Number=1, Type="String", 
                                                         Description="Variant class - bcbio, dragen or both", row.names = "CLASS"))
info(header(both)) <- rbind(info(header(both)), data.frame(Number=1, Type="String", 
                                                     Description="Variant class - bcbio, dragen or both", row.names = "CLASS"))


#write output to file
writeVcf(bcbio, "/Users/kanwals/Documents/UMCCR/data/projects/varlap/bcbio.vcf")
write_vcf(dragen, "/Users/kanwals/Documents/UMCCR/data/projects/varlap/dragen.vcf")
write_vcf(both, "/Users/kanwals/Documents/UMCCR/data/projects/varlap/both.vcf")
