library(ohchibi)
library(ggtree)
library(ggrepel)

set.seed(130816)

#Prepare data for ITOL tree
#Read list of prevalence
Dat <- readRDS(file = "./cleandata/dat_prev_genus.RDS")
Map <- Dat$df_iaadeg_genus

Map$Size <- 0.2
Map$Size[which(Map$TotalGenus>6 & Map$TotalGenus <= 30)] <- 1
Map$Size[which(Map$TotalGenus>31 & Map$TotalGenus <= 100)] <- 5
Map$Size[which(Map$TotalGenus>101 & Map$TotalGenus <= 500)] <- 10
Map$Size[which(Map$TotalGenus>501 & Map$TotalGenus <= 5000)] <- 20
Map$Size[which(Map$TotalGenus>5000)] <- 30

#Fix the NaN
Map$RatioSampledMeanDistPhylo %>% grep(pattern = "NaN") %>%
  Map[.,]
Map$RatioSampledMeanDistPhylo[Map$RatioSampledMeanDistPhylo %>% grep(pattern = "NaN")] <- 1

#Remove marine and bacterium from the phylogenetic tree
mnoms <- which(Map$Genus %in% c("bacterium","marine","Proteobacteria")) %>%
  Map$RepId[.]

Map <- which(!(Map$Genus %in% mnoms)) %>%
  Map[.,] %>% droplevels

#By definition if there is Prevalence of 1 there should be Phylogeny of 1
Map$RatioSampledMeanDistPhylo[which(Map$PercPositive ==1)] <- 1

Map$RatioSampledMeanDistPhylo[which(Map$RatioSampledMeanDistPhylo > 1.05)] <- 1.05

## Tree colors for labels
paleta_type <- c("#BEA2EB","#8AEB9F","#EBB573","#8AE7FA")
names(paleta_type) <- c("iacLike","varioLike","Hybrid","iacLikevarioLike")

Map$Label <- NA
mgenus <- c("Pseudomonas","Acinetobacter","Enterobacter","Paraburkholderia","Burkholderia",
            "Ruegeria","Variovorax","Bradyrhizobium","Achromobacter","Alcaligenes",
            "Leclercia","Lelliottia","Raoultella","Escherichia","Acetobacter")

Map$Label[which(Map$Genus %in% mgenus)] <- Map$Genus[which(Map$Genus %in% mgenus)]



p1 <- ggplot(data = Map,aes(PercPositive,RatioSampledMeanDistPhylo)) +
  geom_hline(yintercept = 0.8,size =1,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0.4,size =1,color = "#D9D9D9",linetype = "longdash")+
  geom_vline(xintercept = 0.8,size =1,color = "#D9D9D9",linetype = "longdash")+
  geom_vline(xintercept = 0.4,size =1,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 1.0,size =1,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(size = Size,color =ClassificationMajority)) +
  geom_label_repel(aes(label = Label,fill = ClassificationMajority),
                   size = 5,alpha = 0.8,segment.alpha = 1,label.size = 0.5) +
  theme_ohchibi() +
  xlab(label = "Prevalence in Genus") +
  ylab(label = "Ratio phylogenetic width") +
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1.06)) +
  scale_x_continuous(breaks = seq(0,1,0.2),limits = c(0,1)) +
  scale_size(name = "Number of genomes",breaks = c(0.2,1,5,10,20,30),
             labels = c("1-5","6-30","31-100","101-500","501-5000",">5000")) +
  scale_color_manual(values = paleta_type)+
  scale_fill_manual(values = paleta_type)



oh.save.pdf(p = p1,outname = "scatter_prevalence.pdf",outdir = "./figures/",width = 16,height = 10)

#Save process table
df_end <- Map[,c("Genus","Hybrid","iacLike","varioLike","TotalGenus",
       "PercHybrid","PerciacLike","PercvarioLike","ClassificationMajority",
       "RatioSampledMeanDistPhylo")]

colnames(df_end)[4] <- "iadLike"
colnames(df_end)[8] <- "PropiadLike"
colnames(df_end)[6] <- "PropHybrid"
colnames(df_end)[7] <- "PropiacLike"

df_end$ClassificationMajority <- df_end$ClassificationMajority %>%
  gsub(pattern = "varioLike",replacement = "iadLike") %>% 
  gsub(pattern = "iacLikeiadLike",replacement = "Hybrid") 

#wirte table
write.table(x = df_end,file = "./cleandata/fig3bc_data.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
rm(list=ls())
gc()
