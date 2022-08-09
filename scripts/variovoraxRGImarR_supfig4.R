library(ohchibi)
library(gggenes)

set.seed(130816)

#Read the structure
merged <- read.table(file = "./rawdata/df_map.iaa_degradation_padbwithreps.tsv",
                     header = T,sep = "\t",quote = "",comment.char = "")

#First analysis using the representatives
merged <- merged %>% subset(Representative == "Yes") %>% droplevels


merged_sub <- merged$Classification %in% c("iac putatively conserved","Vario putatively conserved") %>%
  merged[.,] %>% droplevels() %>% subset(Representative == "Yes") %>% droplevels

aggregate(Hotspot ~taxon_oid,merged_sub,FUN = function(x)table(x) %>% length)

##Check representatives
merged_sub %>% subset(taxon_oid == "GCF_005080685.1")
merged_sub %>% subset(taxon_oid == "GCF_000020125.1")
merged_sub %>% subset(taxon_oid == "GCF_000737145.1")
merged_sub %>% subset(taxon_oid == "GCF_001654845.1")

#Variovorax
merged_sub %>% subset(taxon_oid == "2643221508")



#Bradyrhizobium
merged_sub %>% subset(taxon_oid == "GCF_000807315.1")


#Ruegeria
merged_sub %>% subset(taxon_oid == "GCF_000011965.2")



#Create matrix
Tab <- acast(data = merged_sub,formula = Hotspot~COG,fun.aggregate = sum,value.var = "Count")
Tab <- colnames(Tab) %>% grep(pattern = "NA",invert = T) %>%
  Tab[,.]
#Binarize
Tab[Tab>1] <- 1
melted <- Tab %>% melt
colnames(melted) <- c("Hotspot","COG","value")

colnames(merged_sub)
Map <- merged_sub[,c("taxon_oid","Hotspot","GeneralGenus","Species")] %>%
  unique


## Count the number of enzymes 
df_variofreq <- aggregate(varioFlavor~Hotspot,data = merged_sub,FUN = function(x)table(x) %>% length)
df_variofreq <- merge(Map,df_variofreq,by = "Hotspot")


df_variofreq <- with(df_variofreq,order(taxon_oid)) %>%
  df_variofreq[.,]

#Check the frequency to gather the one we keep
df_variofreq

#Define selectable operons to draw structure tih
mids <- c("HS3278PADB","HS17Reps","HS24Reps","HS57Reps","HS67Reps","HS78Reps","HS90Reps","HS97Reps","HS97Reps","HS106Reps")
merged_sub <- which(merged_sub$Hotspot %in% mids) %>%
  merged_sub[.,]

#Read the gff structure and prepare structure for plotting
df_gff <- read.table(file = "./rawdata/df_gff.tsv")
colnames(df_gff) <- c("taxon_oid","contig","gene_oid","start","end","strand")
df_gff$taxon_oid <- df_gff$taxon_oid %>% gsub(pattern = "\\.gbk\\.gff",replacement = "")



df_gff <- merge(df_gff,merged_sub, by = c("taxon_oid","gene_oid"))
df_gff$direction <- ifelse(df_gff$strand == "+", 1, -1)


#Adjust level of species
df_gff$Species <- df_gff$Species %>%
  factor(levels = c("Pseudomonas putida 1290",
                    "Enterobacter soli ATCC BAA-2102","Acinetobacter baumannii ATCC 19606",
                    "Paraburkholderia phytofirmans PsJN",
                    "Variovorax paradoxus CL14","Ruegeria pomeroyi DSS-3","Bradyrhizobium diazoefficiens USDA 110","Bradyrhizobium japonicum E109"))

df_gff <- which(!(is.na(df_gff$Species))) %>%
  df_gff[.,] %>% droplevels

##Create colors for the genes
df_core <- data.frame(
  COG = c("COG4638","COG5517","COG3347","COG1018"),
  Name = c("iacC","iacD","iacE","iacF"),
  Letter = c("c","d","e","f"),
  Color  = c("#E31A1C","#E31A1C","#E31A1C","#CAB2D6")
)

df_iac <- data.frame(
  COG = c("COG1960","iacB","COG1853","COG0154","iacI"),
  Name = c("iacA","iacB","iacG","iacH","iacI"),
  Letter = c("a","b","g","h","i"),
  Color = c("#6A3D9A","#6A3D9A","#CAB2D6","#CAB2D6","#6A3D9A")
)


df_vario <- data.frame(
  COG =c("COG1014","COG4231","COG1878"),
  Name =c("iorA","iorB","COG1878"),
  Letter = c("j","k","l"),
  Color = c("#33A02C","#33A02C","#33A02C")
)

df_regulation <- data.frame(
  COG = c("COG1846","COG0583"),
  Name = c("MarR","LysR"),
  Letter = c("t","u"),
  Color = c("#FFFF99","#B15928")
)

df_extra <- data.frame(
  COG = c("COG0318","COG0411","COG0183","COG0683","COG0179","COG0559","COG4177"),
  Name = c("COG0318","COG0411","COG0183","COG0683","COG0179","COG0559","COG4177"),
  Letter= c("n","o","m","r","p","s","q"),
  Color = c("#B2DF8A","#1F78B4","#33A02C","#A6CEE3","#1F78B4","#A6CEE3","#1F78B4")
)



df_noms <- rbind(df_core,df_iac,df_vario,df_regulation,df_extra)


#Create a variable name
df_gff$Name <- NA


indices <- which(df_gff$COG %in% df_noms$COG)
df_gff$Name[indices] <- match(df_gff$COG[indices],df_noms$COG) %>% df_noms$Letter[.]


paleta_genes <- df_noms$Color
names(paleta_genes) <- df_noms$Letter

# ggplot(df_gff, aes(xmin = start, xmax = end, y = Species,fill = Name)) +
#   geom_gene_arrow(aes(forward = direction)) +
#   facet_wrap(~ Species, scales = "free", ncol = 1) +
#   theme_genes() +
#   scale_fill_manual(values = paleta_genes)


dummies <- make_alignment_dummies(
  df_gff,
  aes(xmin = start, xmax = end, y = Species, id = Name),
  on = "COG5517"
)

p <- ggplot(df_gff, aes(xmin = start, xmax = end, y = Species,fill = Name,label = Name)) +
  #geom_gene_arrow(aes(forward = direction)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),aes(forward = direction)) +
  geom_gene_label(align = "left") +
  geom_blank(data = dummies) +
  facet_wrap(~ Species, scales = "free", ncol = 1) +
  theme_genes()  +
  scale_fill_manual(values = paleta_genes)

oh.save.pdf(p = p,outname = "genemap_reps.pdf",outdir = "./figures/",width = 16,height = 12)  


#
write.table(x = df_gff,file = "./cleandata/supfig4_data.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

rm(list=ls())
gc()
dev.off()

