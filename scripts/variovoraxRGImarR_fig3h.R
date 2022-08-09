library(ohchibi)
library(DESeq2)
set.seed(130816)


###Added block for refferee request
Dat_all <- readRDS(file = "./cleandata/dat_deseq2_variovoraxmarRRNASeq.RDS")
df_lrt <- Dat_all$dds_lrt   %>%
  results %>% as.data.frame %>% 
  tibble::rownames_to_column()  %>%
  dplyr::rename(.data = .,Gene = rowname )

df_exp <- read.table(file = "./cleandata/df_avexpression_treatments_variovorax_rnaseq.csv",sep = ",",
                     comment.char = "",quote = "",header = TRUE,check.names = FALSE)


df_exp <- match(df_lrt$Gene,df_exp$Gene) %>%
  df_exp[.,]
merged <- cbind(df_lrt,df_exp[,-1])

chosen <- seq(2643613649,2643613677,1) 

merged$Figure2H <- "No"
merged$Figure2H[which(merged$Gene %in% chosen)] <- "Yes"

write.table(x = merged,file = "./cleandata/fig2h_data.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)


####
Dat_all <- readRDS(file = "./cleandata/dat_deseq2_variovoraxmarRRNASeq.RDS")
melted_av <- Dat_all$melted_av  

#Create a matrix to share
df_all <- dcast(data = melted_av,formula = Gene~Treatment,value.var = "value")
Tab <- acast(data = melted_av,formula = Gene~Treatment,value.var = "value")
df_var <- apply(X = Tab,MARGIN = 1,FUN = var) %>%
  data.frame(Gene = names(.),Variance = .,row.names = NULL)
df_all <- df_all[,c("Gene","CL14wt_minus","CL14wt_IAA",
                    "d6869_minus","d6869_IAA",
                    "d73_minus","d73_IAA",
                    "73SR_A_minus","73SR_A_IAA",
                    "d50_minus","d50_IAA"
)]
df_all <- merge(df_all,df_var, by = "Gene")
df_all$Gene <- df_all$Gene %>% factor
df_all <- with(df_all,order(Gene)) %>%
  df_all[.,]

#write.table(x = df_all,file = "../cleandata/df_avexpression_treatments_variovorax_rnaseq.csv",
#            append = F,quote = F,sep = ",",row.names = F,col.names = T)

chosen <- seq(2643613649,2643613677,1)

melted <- which(melted_av$Gene %in%chosen) %>%
  melted_av[.,] %>% droplevels

melted$Gene <- melted$Gene %>% factor
melted$Bacteria <- melted$Bacteria %>% 
  factor(levels = c("CL14wt","d6869","d73","73SR_A","d50"))
melted$value %>% sort %>% plot
melted$Label <- melted$value %>% round(1)

melted <- melted %>%
  subset(Bacteria != "73SR_A") %>% droplevels

p <- ggplot(data = melted,aes(IAA,Gene)) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = Label)) +
  facet_grid(.~Bacteria,space = "free",scales = "free") +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2,2), oob = squish)+
  theme_ohchibi(font_family = "Helvetica") +
  theme(
    strip.background.x = element_blank(),
    # axis.title.x = element_blank()
  ) +
  xlab("IAA supplementation")

oh.save.pdf(p = p,outname = "heatmaphs33_rnaseq.pdf",
            outdir = "./figures/",width = 10,height = 14)

rm(list=ls())
gc()
