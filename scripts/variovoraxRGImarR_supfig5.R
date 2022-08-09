library(ohchibi)

set.seed(130816)

df_trans <- data.frame(Id = c("GCF_003293615.1_ASM329361v1_genomic","GCF_000011365.1_ASM1136v1_genomic","GCA_017831945.1_ASM1783194v1_genomic",
                              "2643221508","GCF_001654845.1_EsoBAA2102_DRAFTv1_genomic","GCF_005080685.1_ASM508068v1_genomic","GCF_000020125.1_ASM2012v1_genomic"),
                       Translation = c("Achromobacter_xylosoxidans","Bradyrhizobium_japonicum_USDA500","Bradyrhizobium_USDA500",
                                       "Variovorax_paradoxus_CL14","Enterobacter_soli","Pseudomonas_putida","Paraburkholderia_phytofirmans"))


order_genomes <- c("Variovorax_paradoxus_CL14","Achromobacter_xylosoxidans","Bradyrhizobium_USDA110","Bradyrhizobium_japonicum_USDA500",
                   "Enterobacter_soli","Pseudomonas_putida","Paraburkholderia_phytofirmans")
### 
mfile <- "./rawdata/COG4638.fasta.mafft.sident.pure"
Tab <- read.table(file = mfile,header = F,sep = "\t",row.names = 1) %>% as.matrix
Tab <- Tab[,1:7]
rownames(Tab)  <- rownames(Tab) %>% gsub(pattern = "\\|.*",replacement = "")
rownames(Tab) <- match(rownames(Tab),df_trans$Id) %>%
  df_trans$Translation[.]
colnames(Tab) <- rownames(Tab)

melted <- Tab %>% melt

df_a <- Tab %>% melt_dist() %>% dplyr::rename(.data = .,Strain1 = iso1,Strain2 = iso2,PercIdentity = dist) %>%
  dplyr::mutate(.data = .,COG = "COG4637_iacC_iadD")

melted$Var1 <- melted$Var1 %>% factor(levels = order_genomes)
melted$Var2 <- melted$Var2 %>% factor(levels = order_genomes)
melted$value <- melted$value %>% gsub(pattern = "^0$",replacement = "100") %>% as.numeric

p1 <- melted %>%
  subset(Var1 != "Bradyrhizobium_USDA500" & Var2!= "Bradyrhizobium_USDA500") %>%
  ggplot(data = .,aes(Var1,Var2)) +
  geom_tile(aes(fill = value),color = "black") +
  geom_text(aes(label = round(value,2))) +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",limits = c(0,100),name = "Alignment\nidentity") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))  +
  ggtitle(label = "COG4638: iacC/iadC")

melted$Type <- "COG4638"
melted_COG4638 <- melted

### 
mfile <- "./rawdata/COG5517.fasta.mafft.sident.pure"
Tab <- read.table(file = mfile,header = F,sep = "\t",row.names = 1) %>% as.matrix
Tab <- Tab[,1:7]
rownames(Tab)  <- rownames(Tab) %>% gsub(pattern = "\\|.*",replacement = "")
rownames(Tab) <- match(rownames(Tab),df_trans$Id) %>%
  df_trans$Translation[.]
colnames(Tab) <- rownames(Tab)

melted <- Tab %>% melt

df_b <- Tab %>% melt_dist() %>% dplyr::rename(.data = .,Strain1 = iso1,Strain2 = iso2,PercIdentity = dist) %>%
  dplyr::mutate(.data = .,COG = "COG5517_iacD_iadE")

melted$Var1 <- melted$Var1 %>% factor(levels = order_genomes)
melted$Var2 <- melted$Var2 %>% factor(levels = order_genomes)
melted$value <- melted$value %>% gsub(pattern = "^0$",replacement = "100") %>% as.numeric

p2 <- melted %>%
  subset(Var1 != "Bradyrhizobium_USDA500" & Var2!= "Bradyrhizobium_USDA500") %>%
  ggplot(data = .,aes(Var1,Var2)) +
  geom_tile(aes(fill = value),color = "black") +
  geom_text(aes(label = round(value,2))) +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",limits = c(0,100),name = "Alignment\nidentity") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))  +
  ggtitle(label = "COG5517: iacD/iadD")

melted$Type <- "COG5517"
melted_COG5517 <- melted

#Check core components
mfile <- "./rawdata/COG1014.fasta.mafft.sident.pure"
Tab <- read.table(file = mfile,header = F,sep = "\t",row.names = 1) %>% as.matrix
Tab <- Tab[,1:4]
rownames(Tab)  <- rownames(Tab) %>% gsub(pattern = "\\|.*",replacement = "")
rownames(Tab) <- match(rownames(Tab),df_trans$Id) %>%
  df_trans$Translation[.]
colnames(Tab) <- rownames(Tab)

melted <- Tab %>% melt


df_c <- Tab %>% melt_dist() %>% dplyr::rename(.data = .,Strain1 = iso1,Strain2 = iso2,PercIdentity = dist) %>%
  dplyr::mutate(.data = .,COG = "COG1014_iorB")

melted$Var1 <- melted$Var1 %>% factor(levels = order_genomes)
melted$Var2 <- melted$Var2 %>% factor(levels = order_genomes)
melted$value <- melted$value %>% gsub(pattern = "^0$",replacement = "100") %>% as.numeric

p3 <- melted %>%
  subset(Var1 != "Bradyrhizobium_USDA500" & Var2!= "Bradyrhizobium_USDA500") %>%
  ggplot(data = .,aes(Var1,Var2)) +
  geom_tile(aes(fill = value),color = "black") +
  geom_text(aes(label = round(value,2))) +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",limits = c(0,100),name = "Alignment\nidentity") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))  +
  ggtitle(label = "COG1014: iorB")

melted$Type <- "COG1014"
melted_COG1014 <- melted

mfile <- "./rawdata/COG4231.fasta.mafft.sident.pure"
Tab <- read.table(file = mfile,header = F,sep = "\t",row.names = 1) %>% as.matrix
Tab <- Tab[,1:4]
rownames(Tab)  <- rownames(Tab) %>% gsub(pattern = "\\|.*",replacement = "")
rownames(Tab) <- match(rownames(Tab),df_trans$Id) %>%
  df_trans$Translation[.]
colnames(Tab) <- rownames(Tab)

melted <- Tab %>% melt

df_d <- Tab %>% melt_dist() %>% dplyr::rename(.data = .,Strain1 = iso1,Strain2 = iso2,PercIdentity = dist) %>%
  dplyr::mutate(.data = .,COG = "COG4231_iotA")

melted$Var1 <- melted$Var1 %>% factor(levels = order_genomes)
melted$Var2 <- melted$Var2 %>% factor(levels = order_genomes)
melted$value <- melted$value %>% gsub(pattern = "^0$",replacement = "100") %>% as.numeric

p4 <- melted %>%
  subset(Var1 != "Bradyrhizobium_USDA500" & Var2!= "Bradyrhizobium_USDA500") %>%
  ggplot(data = .,aes(Var1,Var2)) +
  geom_tile(aes(fill = value),color = "black") +
  geom_text(aes(label = round(value,2))) +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",limits = c(0,100),name = "Alignment\nidentity") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))  +
  ggtitle(label = "COG4231: iotA")

df_end <- rbind(df_a,df_b) %>%
  rbind(df_c) %>%
  rbind(df_d)

write.table(x = df_end,
            file = "./cleandata/supfig5_data.tsv",append = F,
            quote = F,sep = "\t",row.names = F,col.names = T)

rm(list=ls())
gc()

#Save figures
#oh.save.pdf(p = p1,outname = "iacc.idn.pdf",outdir = "./",width = 9,height = 8)
#oh.save.pdf(p = p2,outname = "iacd.idn.pdf",outdir = "./",width = 9,height = 8)
#oh.save.pdf(p = p3,outname = "iorB.idn.pdf",outdir = "./",width = 9,height = 8)
#oh.save.pdf(p = p4,outname = "iotA.idn.pdf",outdir = "./",width = 9,height = 8)
