
#input directory to search for sample files
d <- "/home/gunnar/cluster/becgsc_016/analysis"

l <- list.files(d, pattern = "sample.table.csv", recursive = T)

for (i in 1:length(l))
{
  tmp <- as.character(l[i])
  sample.file <- paste(d,tmp,sep="/")
  path <- 
  if (exists("st", inherits = F)) {
  
  st <- read.csv(paste(d,sf,sep="/"))
}
