library(ohchibi)


Tab <- read.table(file = "./rawdata/markers_count_table.tsv",
                  header = T,sep = "\t",row.names = 1,quote = "",check.names = F)
Tab[Tab>1] <- 0

colSums(Tab)

Tab <- colnames(Tab) %>% grep(pattern = "rpsK",invert = T) %>%
  Tab[,.]
colSums(Tab)


df_freq <- rowSums(Tab) %>% 
  data.frame(Assembly = names(.),NumMarkers = .,row.names = NULL)

#Read metadata
Map <- read.table(file = "./rawdata/metadata_refseq.txt",
                  header = T,sep = "\t",comment.char = "",quote = "")
Map$Assembly <- Map$local_filename %>% gsub(pattern = "\\.fna\\.gz",replacement = "") %>%
  gsub(pattern = "\\.\\/",replacement = "")
Map$Genus <- Map$organism_name %>% gsub(pattern = " .*",replacement = "")

df_freq$Genus <- match(df_freq$Assembly,Map$Assembly) %>%
  Map$Genus[.]

#224 Genus
df_freq$Genus %>% unique

#Check distirbution
df_freq$NumMarkers %>% table
df_touse <- df_freq %>% subset(NumMarkers > 24) %>% droplevels

df_freq %>% subset(NumMarkers <=24) %$% Genus %>% table

Tab_sub  <-match(df_touse$Assembly,rownames(Tab)) %>%
  Tab[.,]
colSums(Tab_sub)
rowSums(Tab_sub) %>% min

dim(Tab_sub)
dim(df_touse)

colSums(Tab_sub)
rowSums(Tab_sub) %>% min


df_touse$Genus %>% unique

#Write a file with markers and genomes to use
genomes <- rownames(Tab)
markers <- colnames(Tab)

melted <- Tab_sub %>%  as.matrix %>% melt
melted <- melted %>% subset(value == 1)
melted <- melted[,-3]

write.table(x = melted,file = "./cleandata/df_genomes_bymarkers.tsv",quote = F,sep = "\t",
            append = F,row.names = F,col.names = F)

