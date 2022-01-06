library("data.table")

#read hotspots tsv file
hotspots <- as.data.frame(read.table("/Users/kanwals/Documents/UMCCR/data/mosdepth/hotspots/Hotspot.tsv.gz", 
                                        header = FALSE, sep = "\t", stringsAsFactors = FALSE))
#assign column names
colnames(hotspots) <- c("chrom", "start", "ref", "alt")

#create a vector that stores the end posiion value for a hotspot - using the start position and number of characters of ref
pos <- vector()
pos <- apply(hotspots[,c('start', 'ref')], 1, function(x){
  if(nchar(x['ref']) > 1){
    numchar <- nchar(x['ref'])
    val <- as.integer(x['start']) + numchar
  } else {
    val = as.integer(x['start']) + 1
  }
})

#create a new column using the above created vector
hotspots['end'] <- pos

#change columns order in hotspots df
hotspots <- hotspots[,c('chrom', 'start', 'end', 'ref', 'alt')]

#write final output df to a file
fwrite(hotspots, 
       file = "~/Documents/UMCCR/data/mosdepth/hotspots/hotspots.bed",
       sep = "\t",
       col.names = FALSE)
