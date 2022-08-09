library(ohchibi)
library(ggtree)

set.seed(130816)

Res_reg <- readRDS(file = "./cleandata/res_reg_cluster.RDS")

#Detine specific clsuters with unique ids
muids_vario <- Res_reg %>%
  subset(V5 == "BorderContigNo") %>%
  subset(V10 == "3021") %>%
  subset(Type == "varioLike") %>%
  subset(EndCount == 1 & COG == "COG1846" ) %$% UId %>% as.character %>% unique

vario_ids_two <- which(Res_reg$UId %in% muids_vario) %>%
  Res_reg[.,] %>%
  subset(EndCount == 2 & COG == "COG1846" ) %$% UId %>% as.character %>% unique

ids_iac_a <- Res_reg %>%
  subset(V10 == "3303") %>%
  subset(V5 == "BorderContigNo") %>%
  subset(Type == "iacLike") %>%
  subset(EndCount == 1 & COG == "COG1846" ) %$% UId %>% as.character %>% unique


ids_iac_temp <- Res_reg %>%
  subset(V10 == "3303")  %>% 
  subset(V5 == "BorderContigNo") %>%
  subset(Type == "iacLike") %>%
  subset(EndCount == 1 & COG != "COG0583"  & ClusterId == "C4") %$% UId %>% as.character %>% unique

ids_iac_b <- which(Res_reg$UId %in% ids_iac_temp) %>%
  Res_reg[.,] %>%
  subset(V10 == "3303")  %>% 
  subset(V5 == "BorderContigNo") %>%
  subset(EndCount == 2 & COG == "COG1846" ) %$% UId %>% as.character %>%
  unique

chosen_uids <- c(vario_ids_two,ids_iac_a,ids_iac_b) %>% unique

Res_chosen <- which(Res_reg$UId %in% chosen_uids) %>%
  Res_reg[.,] %>%
  subset(EndCount == 1 | EndCount == 2) %>%
  subset(COG == "COG1846") %>% droplevels

Res_a <- Res_chosen %>% subset(Type == "varioLike" & EndCount == 1) 
Res_b <- Res_chosen %>% subset(Type == "varioLike" & EndCount == 2) 
Res_a$Specific <- "varioLikemarR73"
Res_b$Specific <- "varioLikemarR50"
Res_c <- Res_chosen %>% subset(Type == "iacLike") 
Res_c$Specific <- "iacLike"
Res_chosen <- rbind(Res_a,Res_b,Res_c)

id_genes <- Res_chosen$V4 %>% unique

#Select markers from other cog
#This is a similar marR
control <- Res_reg %>%
  subset(V5 == "BorderContigNo") %>%
  subset(COGId == "COG5631") %$% V4 %>% unique

#MAke sure the marRs from control are included
one <- Res_reg$V1 %>%
  grep(pattern = "Eso") %>%
  Res_reg[.,] %>%
  subset(COG == "COG1846") %$% V4

two <- Res_reg$V1 %>%
  grep(pattern = "ASM1136v1") %>%
  Res_reg[.,] %>%
  subset(COG == "COG1846") %$% V4

three <- Res_reg$V1 %>%
  grep(pattern = "ASM1196v2") %>%
  Res_reg[.,] %>%
  subset(COG == "COG1846") %$% V4

four <- Res_reg$V1 %>%
  grep(pattern = "ASM47302v1") %>%
  Res_reg[.,] %>%
  subset(COG == "COG1846") %$% V4

five <- Res_reg$V1 %>%
  grep(pattern = "ASM73714v1") %>%
  Res_reg[.,] %>%
  subset(COG == "COG1846") %$% V4

six <- Res_reg$V1 %>%
  grep(pattern = "ASM164267v1") %>%
  Res_reg[.,] %>%
  subset(COG == "COG1846") %$% V4

seven <- Res_reg$V1 %>%
  grep(pattern = "ASM190831v1") %>%
  Res_reg[.,] %>%
  subset(COG == "COG1846")  %$% V4

eight <- Res_reg$V1 %>%
  grep(pattern = "ASM508068v1") %>%
  Res_reg[.,] %>%
  subset(COG == "COG1846")  %$% V4


chosen_ids <- c(id_genes ,control,one,two,three,four,five,six,seven,eight) %>% unique

#Read the chosen ids
chosen_ids <- read.table(file = "./rawdata/chosen_ids_marR_tree.txt") %$% V1

##
Res_chosen <- match(chosen_ids,Res_reg$V4) %>%
  Res_reg[.,] %>% droplevels

df <- Res_chosen[,c("V1","V4","Type","EndCount")] 
df$Specific <- "iacLike"
df$Specific[which( (df$Type == "varioLike") & (df$EndCount == 1) )] <- "varioLikeMarR73"
df$Specific[which( (df$Type == "varioLike") & (df$EndCount == 2) )] <- "varioLikeMarR50"

df_gt <- read.table(file = "./rawdata/gtdbtk.bac120.summary.tsv",
                    header = T,sep = "\t",quote = "",comment.char = "")
df_gt$user_genome <- df_gt$user_genome %>% gsub(pattern = "\\.fna",replacement = "")
df$Classification <- match(df$V1,df_gt$user_genome) %>%
  df_gt$classification[.]
tobind <- df$Classification %>% strsplit(split = "\\;") %>%
  unlist %>%
  gsub(pattern = "[dpcofgs]__",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T)  %>%
  as.data.frame
colnames(tobind) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
df <- cbind(df,tobind)
df$ShortGenus <- df$Genus %>% gsub(pattern = "_[A-Z]+$",replacement = "")
df$PhylumShortGenus <- paste0(df$Phylum,"|",df$ShortGenus)

#Append the representative column
df$Representative <- NA
colnames(df)[1:2] <- c("Id","GeneId")

df$Representative[df$Id %>% grep(pattern = "Eso")] <- "Yes"
df$Representative[df$Id %>% grep(pattern = "ASM1136v1")] <- "Yes"
df$Representative[df$Id %>% grep(pattern = "ASM1196v2")] <- "Yes"
#df$Representative[df$Id %>% grep(pattern = "ASM47302v1")] <- "Yes"
df$Representative[df$Id %>% grep(pattern = "ASM73714v1")] <- "Yes"
#df$Representative[df$Id %>% grep(pattern = "ASM164267v1")] <- "Yes"
#df$Representative[df$Id %>% grep(pattern = "ASM190831v1")] <- "Yes"
df$Representative[df$Id %>% grep(pattern = "ASM508068v1")] <- "Yes"
df$Representative[df$Id %>% grep(pattern = "^2")] <- "Yes"


df_freq_acineto <- df %>% subset(Genus == "Acinetobacter") %$% Species %>% table %>%
  sort(decreasing = T) %>% data.frame


df_other <- df %>% subset(Genus != "Acinetobacter") %>%
  droplevels


x <- df %>% subset(Species == "Acinetobacter baumannii" & Representative == "Yes") %$% GeneId
pa  <- df %>% subset(Species == "Acinetobacter baumannii") %$% GeneId  %>%
  sample(size = 10) %>% unique

y <- df %>% subset(Species == "Acinetobacter seifertii") %$% GeneId  %>%
  sample(size = 10) %>% unique
z <- df %>% subset(Species == "Acinetobacter nosocomialis") %$% GeneId  %>%
  sample(size = 10) %>% unique
p  <- df %>% subset(Species == "Acinetobacter pittii") %$% GeneId  %>%
  sample(size = 10) %>% unique

chosen_acineto <- c(x,pa,y,z,p) %>% unique
df_chose_acineto <- which(df$GeneId %in% chosen_acineto) %>%
  df[.,] %>% droplevels
df_end <- rbind(df_other,df_chose_acineto)



###What do I need to add here for this phylogeny
tree <- read.tree(file = "./rawdata/chosen_ids_marR_tree.mafft.fasttree.wag_gamma.rooted.newick")

todrop <- which(!(tree$tip.label %in% df_end$GeneId)) %>%
  tree$tip.label[.]

tree <- ape::drop.tip(phy = tree,tip = todrop)


df_end <- match(tree$tip.label,df_end$GeneId) %>%
  df_end[.,] %>% droplevels

rownames(df_end) <- df_end$GeneId
df_end <- df_end[,c(2,1,3:16)]

df_freq_genus <- df_end$PhylumShortGenus %>% table %>%
  sort(decreasing = T) %>% data.frame
chosen_Genus <- c(df_freq_genus$.[1:7] %>% as.character,"Proteobacteria|Ruegeria","Proteobacteria|Enterobacter")

df_end$ColorGenus <- "Other"
mindices <- which(df_end$PhylumShortGenus %in% chosen_Genus)
df_end$ColorGenus[mindices] <- df_end$ShortGenus[mindices]
df_end$ColorGenus <- df_end$ColorGenus %>% factor()

df_end$Representative[which(is.na(df_end$Representative))] <- "No"
df_end$Representative <- df_end$Representative  %>% factor(levels = c("No","Yes") %>% rev)

#Identify contradictory
mchange <- df_end %>% subset(Type == "varioLike" & Specific == "iacLike") %$% GeneId
mids <- df_end$Id[which(df_end$GeneId %in% mchange)]
df_t <- df_end[which(df_end$Id %in% mids),] %>% droplevels
df_nt <-  df_end[which(!(df_end$Id %in% mids)),] %>% droplevels

df_t$Specific <- df_t$Specific %>% gsub(pattern = "varioLikeMarR50",replacement = "varioLikeMarR73") %>%
  gsub(pattern = "iacLike",replacement = "varioLikeMarR50")

df_end <- rbind(df_nt,df_t)

df_end$Label <- NA
mindices <- which(df_end$Representative == "Yes")
df_end$Label[mindices] <- df_end$ShortGenus[mindices]


paleta_type <- c("#BEA2EB","#8AEB9F","#F28F9B")
names(paleta_type) <- c("iacLike","varioLikeMarR73","varioLikeMarR50")


p <- ggtree(tr = tree) +
  scale_y_reverse()
p_end <- p %<+% df_end +
  geom_tippoint(aes(color = Specific,shape = Representative),size = 1)+
  scale_shape_manual(values = c(16,15) %>% rev) +
  geom_tiplab(aes(label = Label)) +
  scale_color_manual(values = paleta_type)

write.tree(phy = tree,file = "./cleandata/fig5b_tree.newick",append = FALSE)

Mapa <- df_end[,c("GeneId","Specific","Representative","Label")]
Mapa$Specific <- Mapa$Specific %>% 
  gsub(pattern = "varioLike",replacement = "iad")
rownames(Mapa) <- NULL
colnames(Mapa)[2] <- "Type"
write.table(x = Mapa,file = "./cleandata/fig5b_metadata.tsv",
            append = FALSE,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


oh.save.pdf(p = p_end,outname = "phylogeny_marR.pdf",outdir = "./figures/",width = 8,height = 14)
rm(list=ls())
gc()
