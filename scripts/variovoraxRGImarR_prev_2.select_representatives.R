library(ohchibi)

set.seed(130816)

mdf_chosen <- read.table(file = "./cleandata/df_genomes_bymarkers.tsv",header = F,sep = "\t")
strains_chosen <- mdf_chosen$V1 %>% unique

#Read metadata
Map <- read.table(file = "./rawdata/metadata_refseq.txt",
                  header = T,sep = "\t",comment.char = "",quote = "")

Map$Assembly <- Map$local_filename %>% gsub(pattern = "\\.fna\\.gz",replacement = "") %>%
  gsub(pattern = "\\.\\/",replacement = "")

Map$Genus <- Map$organism_name %>% gsub(pattern = " .*",replacement = "")

Map <- which(Map$Assembly %in% strains_chosen) %>%
  Map[.,] %>% droplevels 


df_genus <- Map$Genus %>% table %>%
  data.frame
colnames(df_genus)[1] <- "Genus"


#Sample randomly each genus to select 1 representative
mgenus <- df_genus$Genus %>% as.character %>% unique

#Load the marker table
df_markers <- read.table(file = "./cleandata/df_genomes_bymarkers.tsv",
                         header = F,sep = "\t",comment.char = "",quote = "")

df_markers$Count <- 1
df_sum_markers <- aggregate(Count~V1,df_markers,sum)

Map$CountMarkers <- match(Map$Assembly,df_sum_markers$V1) %>%
  df_sum_markers$Count[.]

Map <- with(Map,order(-CountMarkers)) %>%
  Map[.,]

reps <- NULL
for(genus in mgenus){
  Map_sub <- Map %>% subset(Genus == genus)  %>% droplevels
  mnum <- Map_sub$CountMarkers %>% max
  mc <- Map_sub %>% subset(CountMarkers == mnum) %$%Assembly %>%
    sample(size = 1,replace = F) 
  reps <- data.frame(Genus = genus,Rep = mc,NumMarkers = mnum) %>%
    rbind(reps,.)
}


write(x = reps$Rep,file = "./cleandata/reps_224.txt",append = F)


