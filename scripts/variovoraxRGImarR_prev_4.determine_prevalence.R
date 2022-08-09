library(ohchibi)
library(ggtree)
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

df_sub$Type <- c("Hybrid")
df_sub$Type[df_sub$MinimalPattern %>% grep(pattern = "111100")] <- "iacLike"
df_sub$Type[df_sub$MinimalPattern %>% grep(pattern = "110011")] <- "varioLike"
df_sub$Type[df_sub$MinimalPattern %>% grep(pattern = "110111")] <- "varioLike"


#Determine the repated
df_sub$Count <- 1
Tab <- acast(data = df_sub,formula = Assembly~Type,fun.aggregate = sum,value.var = "Count")



df_repeated <- apply(X = Tab,MARGIN = 1,FUN = function(x)which(x>0) %>% length) %>%
  data.frame(Assembly = names(.),Num = .,row.names = NULL) %>%
  subset(Num > 1)

#Los repeated tienen mas de un operon
df_melt <- Tab %>% melt
df_melt <- df_melt %>% subset(value != 0) %>% droplevels

msum <- aggregate(value~Var1,df_melt,sum)

df_melt$NumHSGenome <- match(df_melt$Var1,msum$Var1) %>%
  msum$value[.]


colnames(df_melt) <- c("Assembly","Type","Freq","NumHsGenome")
df_melt$Perc <- df_melt$Freq/df_melt$NumHsGenome
df_melt$Genus <- match(df_melt$Assembly,df_sub$Assembly) %>%
  df_sub$Genus[.]

df_melt %>% subset(Genus == "Variovorax") 



#Count based on Type
df_total <- dcast(data = df_melt,formula = Genus~Type,fun.aggregate = sum,value.var = "Perc")
df_total$TotalPositive <- df_total$Hybrid + df_total$iacLike + df_total$varioLike


### Define all Genus ###
mgenus <- df_sub$Genus %>% as.character %>%
  unique

#Subset the strain to scan the from
Map$Genus <- Map$organism_name %>%
  gsub(pattern = " .*",replacement = "")
Map_sub <- which(Map$Genus %in% mgenus) %>%
  Map[.,] %>% droplevels 

#Determine the total 
df_tot <- df[,c("Assembly","Genus")] %>% unique %$% Genus %>% table %>%
  sort(decreasing = T) %>% as.data.frame

df_total$TotalGenus <- match(df_total$Genus,df_tot$.) %>%
  df_tot$Freq[.]

df_total$PercPositive <- df_total$TotalPositive/df_total$TotalGenus

#Convert it to percentage
df_total$PercHybrid <- (df_total$Hybrid/df_total$TotalGenus)
df_total$PerciacLike <- (df_total$iacLike/df_total$TotalGenus)
df_total$PercvarioLike <- (df_total$varioLike/df_total$TotalGenus)
df_total$PercMissing <- 1-df_total$PercPositive

#Create a classification scheme based on percentage
#Classify each genus based on majority rule
melted <- melt(data = df_total,id.vars = c("Genus"),
               measure.vars = c("PercHybrid","PerciacLike","PercvarioLike"))
mgenus <- melted$Genus %>% unique
Res_classification <- NULL
for(mg in mgenus){
  melted_temp <- melted %>% subset(Genus == mg) %>% droplevels
  melted_temp <- with(melted_temp,order(-value)) %>%
    melted_temp[.,]
  mval <- melted_temp$value[1]
  melted_temp <- melted_temp %>% subset(value == mval)
  if(nrow(melted_temp) == 1){
    Res_classification <- data.frame(Genus = mg,Majority = melted_temp$variable[1] %>% as.character) %>%
      rbind(Res_classification,.)
  }else{
    cat(mg,"\t",nrow(melted_temp),"\n")
    print(melted_temp)
  }
  
}
Res_classification <- data.frame(Genus = c("Celeribacter","Cupriavidus","Ideonella","Labrenzia","Zoogloea"),
                                 Majority = c("PerciacLike","PerciacLikePercvarioLike","PerciacLikePercvarioLike","PerciacLikePercvarioLike","PerciacLike")) %>%
  rbind(Res_classification,.)

df_total$ClassificationMajority <- match(df_total$Genus,Res_classification$Genus) %>%
  Res_classification$Majority[.]
df_total$ClassificationMajority <- df_total$ClassificationMajority %>% gsub(pattern = "Perc",replacement = "")

#Append the taxonomic width information
tree <- read.tree(file = "./cleandata/concatenated.aln.fasttree.wag_gamma.newick")

## Extract the information about diversity per genus
Map_tree <- Map
Map_tree <- match(tree$tip.label,Map_tree$Assembly) %>%
  Map_tree[.,] %>% droplevels
mgenus <- Map_tree$Genus %>% unique

#Loop over the genus to compute total span diversity by positive hits and total diversity in branch
Res_phylo <- NULL
counter <- 0
for(mg in mgenus){
  counter <- counter+1
  cat("Working on ",mg,"\t",counter,"\n")
  #Add the Acidovorax caveat
  if(mg == "Acidovorax"){
    next
  }
  Map_temp <- Map_sub %>% subset(Genus == mg) %>% droplevels
  ids <- Map_temp$Assembly %>% as.character
  todrop <- which(!(tree$tip.label %in% ids)) %>%
    tree$tip.label[.]
  tree_temp <- ape::drop.tip(phy = tree,tip = todrop)
  if(length(tree_temp$tip.label) ==1){
    av_val <- tree_temp$edge.length
    mean_val <- tree_temp$edge.length
    max_val <- tree_temp$edge.length
  }else{
    melted_cop <- cophenetic.phylo(tree_temp) %>% melt_dist
    av_val <- melted_cop$dist %>% sample(x = .,replace = TRUE,size = 1000) %>% mean
    mean_val <- melted_cop$dist %>% mean
    max_val <- melted_cop$dist %>% max
  }
  
  #Select the only positive ones
  mids <- df_sub %>% subset(Genus == mg) %$% Assembly %>% as.character
  todrop <- which(!(tree$tip.label %in% mids)) %>%
    tree$tip.label[.]
  tree_temp <- ape::drop.tip(phy = tree,tip = todrop)
  if(length(tree_temp$tip.label) ==1){
    av_val_positive <- tree_temp$edge.length
    mean_val_positive <- tree_temp$edge.length
    max_val_positive <- tree_temp$edge.length
  }else{
    melted_cop <- cophenetic.phylo(tree_temp) %>% melt_dist
    av_val_positive <- melted_cop$dist %>% sample(x = .,replace = TRUE,size = 1000) %>% mean
    mean_val_positive <- melted_cop$dist %>% mean
    max_val_positive <- melted_cop$dist %>% max
  }
  temp <- data.frame(Genus = mg,MeanDistPhyloTotal = mean_val,MaxDistPhyloTotal = max_val,SampledDistPhyloTotal = av_val,
                     MeanDistPhyloPositive = mean_val_positive,MaxDistPhyloPositive = max_val_positive,SampledDistPhyloPositive = av_val_positive)
  Res_phylo <- rbind(Res_phylo,temp)
}

#Appendn the Acidovorax information
mg <- "Acidovorax"
Map_temp <- Map_sub %>% subset(Genus == mg) %>% droplevels
ids <- Map_temp$Assembly %>% as.character
todrop <- which(!(tree$tip.label %in% ids)) %>%
  tree$tip.label[.]
tree_temp <- ape::drop.tip(phy = tree,tip = todrop)
if(length(tree_temp$tip.label) ==1){
  av_val <- tree_temp$edge.length
  mean_val <- tree_temp$edge.length
  max_val <- tree_temp$edge.length
}else{
  melted_cop <- cophenetic.phylo(tree_temp) %>% melt_dist
  av_val <- melted_cop$dist %>% sample(x = .,replace = TRUE,size = 1000) %>% mean
  mean_val <- melted_cop$dist %>% mean
  max_val <- melted_cop$dist %>% max
}
#Addd value
av_val_positive <- 0.01
mean_val_positive <- 0.01
max_val_positive <-0.01
temp <- data.frame(Genus = mg,MeanDistPhyloTotal = mean_val,MaxDistPhyloTotal = max_val,SampledDistPhyloTotal = av_val,
                   MeanDistPhyloPositive = mean_val_positive,MaxDistPhyloPositive = max_val_positive,SampledDistPhyloPositive = av_val_positive)
Res_phylo <- rbind(Res_phylo,temp)


#Calculate the ratios
Res_phylo$RatioSampledMeanDistPhylo <- Res_phylo$SampledDistPhyloPositive/Res_phylo$SampledDistPhyloTotal

#Crete the distribution of chosen genus
plots_genus <- NULL
paleta <- c("black","red")
names(paleta) <- c("No","Yes")
mgenus <- c("Pseudomonas","Acinetobacter","Enterobacter","Paraburkholderia","Burkholderia",
            "Ruegeria","Variovorax","Bradyrhizobium","Achromobacter","Alcaligenes",
            "Escherichia","Leclercia","Acetobacter")

for(mg in mgenus){
  cat("Plotting ",mg,"\n")
  Map_temp <- Map_sub %>% subset(Genus == mg) %>% droplevels
  ids <- Map_temp$Assembly %>% as.character
  todrop <- which(!(tree$tip.label %in% ids)) %>%
    tree$tip.label[.]
  tree_temp <- ape::drop.tip(phy = tree,tip = todrop)
  Map_temp <- Map %>% subset(Genus == mg) %>% droplevels
  Map_temp <- match(tree_temp$tip.label,Map_temp$Assembly) %>%
    Map_temp[.,]
  rownames(Map_temp) <- Map_temp$Assembly
  Map_temp$IAADeg <- "No"
  Map_temp$IAADeg[which(Map_temp$Assembly %in% df_sub$Assembly)] <- "Yes"
  Map_temp$IAADeg <- Map_temp$IAADeg %>% factor
  Map_temp <- Map_temp[,c(24,1:23,25,26)]
  
  p <- ggtree(tr = tree_temp) + scale_y_reverse()
  p <- p %<+% Map_temp + geom_tippoint(aes(color=IAADeg), size=1)  +
    scale_color_manual(values = paleta)
  plots_genus[[mg]] <- p
  
}

#Create merged structure
df_total <- merge(df_total,Res_phylo,by = "Genus",all.x = TRUE)

#Append the representative information
ids_rep <- read.table(file = "./cleandata/reps_224.txt") %$% V1

df_rep <- data.frame(RepresentativeId = ids_rep, Genus = match(ids_rep,Map$Assembly) %>%
                       Map$Genus[.])

df_total$RepId <- match(df_total$Genus,df_rep$Genus) %>%
  df_rep$RepresentativeId[.]


#Save structures
mlist <- list(
  df_positive_raw = df_sub,
  df_iaadeg_genus = df_total,
  plots_genus = plots_genus,
  Map_sub = Map_sub
)
saveRDS(object = mlist,file = "./cleandata/dat_prev_genus.RDS")
rm(list=ls())
