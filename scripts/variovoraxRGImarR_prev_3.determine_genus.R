library(ohchibi)
set.seed(130816)

df <- read.table(file = "./rawdata/all_df_summary_iaadegradation.tsv",
                 header = F,sep = "\t")
df_freq <- df$V7 %>% table %>%
  data.frame
colnames(df_freq)[1] <- "Pattern"

colnames(df) <- c("Assembly","Hotspot","Scaffold","Border","Position","Length","Pattern","COGs")

#Read metadata
Map <- read.table(file = "./rawdata/metadata_refseq.txt",
                  header = T,sep = "\t",comment.char = "",quote = "")
Map$Assembly <- Map$local_filename %>% gsub(pattern = "\\.fna\\.gz",replacement = "") %>%
  gsub(pattern = "\\.\\/",replacement = "")


df <- merge(df,Map[,c("Assembly","taxid","species_taxid","organism_name","infraspecific_name","refseq_category")])
df$Genus <- df$organism_name %>% gsub(pattern = " .*",replacement = "")

#Create core patterns for iac operon
df$COG4638 <- 0
df$COG4638[df$COGs %>% grep(pattern = "COG4638")] <- 1

df$COG5517 <- 0
df$COG5517[df$COGs %>% grep(pattern = "COG5517")] <- 1

df$iacB <- 0
df$iacB[df$COGs %>% grep(pattern = "iacB")] <- 1

df$iacI <- 0
df$iacI[df$COGs %>% grep(pattern = "iacI")] <- 1

#Create patterns for Variovovrax operon
df$COG4231 <- 0
df$COG4231[df$COGs %>% grep(pattern = "COG4231")] <- 1

df$COG1014 <- 0
df$COG1014[df$COGs %>% grep(pattern = "COG1014")] <- 1


df$MinimalPattern <- paste0(df$COG4638,df$COG5517,df$iacB,df$iacI,df$COG4231,df$COG1014) 

df$MinimalPattern %>% table

chosen <- c("111100","110011","110111","111101","111110","111111")

df_sub <- which(df$MinimalPattern %in% chosen) %>%
  df[.,] %>% droplevels 

### Define all Genus ###
mgenus <- df_sub$Genus %>% as.character %>%
  unique


#Subset the strain to scan the from
Map$Genus <- Map$organism_name %>%
  gsub(pattern = " .*",replacement = "")
Map_sub <- which(Map$Genus %in% mgenus) %>%
  Map[.,] %>% droplevels 

write(x = Map_sub$Assembly,
      file = "./cleandata/ids_buildtree.txt",ncolumns = 1,append = F)


#Determine the total 
df_tot <- df[,c("Assembly","Genus")] %>% unique %$% Genus %>% table %>%
  sort(decreasing = T) %>% as.data.frame

df_found <- df_sub[,c("Assembly","Genus")] %>% unique %$% Genus %>% table %>%
  sort(decreasing = T) %>% as.data.frame
merged <- merge(df_tot,df_found, by = ".")
colnames(merged) <- c("Genus","TotalRefSeqGenomes","PositiveRefSeqGenomes")
merged$Perc <- (merged$PositiveRefSeqGenomes/merged$TotalRefSeqGenomes)*100
merged <- with(merged,order(-Perc)) %>%
  merged[.,]
merged


#Subset Genus with higher than 5
merged %>% subset(TotalRefSeqGenomes > 5)

mgenus <- c("Acinetobacter","Pseudomonas","Paraburkholderia",
            "Burkholderia","Bradyrhizobium","Variovorax","Enterobacter",
            "Ruegeria","Alcaligenes","Achromobacter","Leclercia","Lelliottia","Raoultella","Paenibacillus")

Map_sub <- which(Map$Genus %in% mgenus) %>%
  Map[.,] %>% droplevels
ids_subphylogeny <- Map_sub$Assembly %>% as.character %>% unique
write(x = ids_subphylogeny,file = "./cleandata/ids_genusbuildtree.txt",append = F)
