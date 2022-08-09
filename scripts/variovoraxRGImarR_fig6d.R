library(ohchibi)

set.seed(130816)

Dat <- readRDS(file = "./cleandata/dat_amplicon_sc175_withdropout.RDS")

#Raw spiking
Dat_sub <- Dat$Dat_copies

melted <- Dat_sub$Tab %>% melt
colnames(melted) <- c("Id","TLMiSeq","value")

melted <- merge(melted,Dat_sub$Map, by = "TLMiSeq")
melted <- melted$Id %>% grep(pattern = "^ASV",invert = T) %>%
  melted[.,] %>% droplevels

melted$Bacteria <- melted$Bacteria %>% 
  factor(levels = c("PelletSC175","Alone","10Vario","CL14","MarRS28AR46A","Delta73","Delta68-69",
                    "B.japonicum","P.putida1290","Root562","PsJN","E.soli","CL71","CL69"))

melted$Id <- melted$Id %>% factor(levels = c("CL14","VarioASV1","VarioASV2",
                                             "B.japonicum","P.putida1290","Root562","PsJN","E.soli","CL71","CL69"))
mids <- melted$Id %>% levels
paleta <- paletteer_d("ggthemes::Tableau_20")[1:10]
names(paleta) <- mids

melted_end <- melted[which(!(melted$Bacteria %in% c("Alone","PelletSC175"))),] %>%
  droplevels  

#Put colors as illustrator
melted_end$Lvalue <- log10(melted_end$value)
melted_end$Lvalue[is.infinite(melted_end$Lvalue)] <- NA

#Establish color palette
paleta <- c("#39A048","#C2B59B","#D9D9D9","#B4D88C","#E11F28","#F6999B","#6B4099","#CAB4D7","#B25A29","#FFDE17")
noms <- c("CL14","VarioASV1","VarioASV2","B.japonicum","P.putida1290","Root562","PsJN","E.soli","CL71","CL69")
names(paleta) <- noms

#Keep only consistent pairs
df_a <- melted_end %>%
  subset(Treatment == "175SynCom+_Delta68-69") %>% droplevels
df_a$Lvalue[which((df_a$Id == "CL14") & (is.na(df_a$Lvalue)))] <- 0


df_b <- melted_end %>%
  subset(Treatment == "175SynCom+_Delta73") %>% droplevels
df_b$Lvalue[which((df_b$Id == "CL14") & (is.na(df_b$Lvalue)))] <- 0

df_c <- melted_end %>%
  subset(Treatment == "175SynCom+_B.japonicum") %>% droplevels
df_c$Lvalue[which((df_c$Id == "B.japonicum") & (is.na(df_c$Lvalue)))] <- 0

df_d <- melted_end %>%
  subset(Treatment == "175SynCom+_CL69") %>% droplevels
df_d$Lvalue[which((df_d$Id == "CL69") & (is.na(df_d$Lvalue)))] <- 0


df_e <- melted_end %>%
  subset(Treatment == "175SynCom+_10vario") %>% droplevels
df_e$Lvalue[which((df_e$Id == "CL14") & (is.na(df_e$Lvalue)))] <- 0
df_e$Lvalue[which((df_e$Id == "VarioASV1") & (is.na(df_e$Lvalue)))] <- 0
df_e$Lvalue[which((df_e$Id == "VarioASV2") & (is.na(df_e$Lvalue)))] <- 0

df_f <- melted_end %>%
  subset(Treatment == "175SynCom+_Root562") %>% droplevels
df_f$Lvalue[which((df_f$Id == "Root562") & (is.na(df_f$Lvalue)))] <- 0


df_g <- melted_end %>%
  subset(Treatment == "175SynCom+_CL14") %>% droplevels
df_g$Lvalue[which((df_g$Id == "CL14") & (is.na(df_g$Lvalue)))] <- 0

df_h <- melted_end %>%
  subset(Treatment == "175SynCom+_CL71") %>% droplevels
df_h$Lvalue[which((df_h$Id == "CL71") & (is.na(df_h$Lvalue)))] <- 0

df_i <- melted_end %>%
  subset(Treatment == "175SynCom+_P.putida1290") %>% droplevels
df_i$Lvalue[which((df_i$Id == "P.putida1290") & (is.na(df_i$Lvalue)))] <- 0

df_j <- melted_end %>%
  subset(Treatment == "175SynCom+_MarRS28AR46A") %>% droplevels
df_j$Lvalue[which((df_j$Id == "CL14") & (is.na(df_j$Lvalue)))] <- 0


df_k <- melted_end %>%
  subset(Treatment == "175SynCom+_PsJN") %>% droplevels
df_k$Lvalue[which((df_k$Id == "PsJN") & (is.na(df_k$Lvalue)))] <- 0

df_l <- melted_end %>%
  subset(Treatment == "175SynCom+_E.soli") %>% droplevels
df_l$Lvalue[which((df_l$Id == "E.soli") & (is.na(df_l$Lvalue)))] <- 0

melted_end <- rbind(df_a,df_b) %>%
  rbind(df_c) %>%
  rbind(df_d) %>%
  rbind(df_e) %>%
  rbind(df_f) %>%
  rbind(df_g) %>%
  rbind(df_h) %>%
  rbind(df_i) %>%
  rbind(df_j) %>%
  rbind(df_k) %>%
  rbind(df_l) 

melted_end$Bacteria <- melted_end$Bacteria %>% 
  factor(levels = c("10Vario","CL14","Delta68-69","Delta73","MarRS28AR46A","B.japonicum","P.putida1290","Root562","PsJN","E.soli","CL71","CL69"))

melted_end <- melted_end %>% subset(Id != "VarioASV2") %>% droplevels

#Save structure

#Save table
df_end <- melted_end[,c("Id","Bacteria","Lvalue")]%>% 
  na.omit %>%
  dplyr::rename(.data = .,Background = Id,Treatment = Bacteria, LogCFU = Lvalue)


write.table(x = df_end,file = "./cleandata/fig6d_dropin_copies.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

p <- ggplot(data = melted_end,aes(Id,Lvalue)) +
  geom_point(aes(color = Id),alpha = 0.5,size = 3,shape = 21) +
  stat_summary(fun  = median,geom = "point",aes(fill = Id),size = 5,shape =23,stroke = 0.3)  +
  facet_grid(.~Bacteria,space = "free") +
  theme_ohchibi(size_panel_border = 0.1) +
  theme(
    axis.text.x = element_text(size = 5,angle = 45,vjust = 1,hjust = 1),
    strip.background.x = element_blank()
  ) +
  ylab("Number of 16s copies") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7),labels = c(0,10,100,1000,10000,100000,1000000,10000000)) +
  scale_x_discrete(expand = c(0,1.5)) +
  theme(
    panel.grid.major.y = element_line(color = "grey89",size = 0.25),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),axis.title.x = element_blank()
  ) +
  scale_colour_manual(values = paleta) +
  scale_fill_manual(values = paleta)

oh.save.pdf(p = p,outname = "16s_end.revision.pdf",outdir = "./figures/",width = 10,height = 10)


dev.off()
rm(list=ls())
gc()

