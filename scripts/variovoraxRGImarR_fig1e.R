library(ohchibi)
library(xlsx)
library(emmeans)

set.seed(130816)


df <- read.xlsx(file = "./rawdata/Root Elongation quantification IAA and CL28 with CL14 68-69 mutants.xlsx",sheetIndex = 1)

df$condition %>% table


melted <- df[,c(3,5:14)] %>%
  melt(data = .,id.variables = "condition") %>%
  na.omit


melted$Facet <- "NB"
melted$Facet[melted$condition %>% grep(pattern = "^CL28")] <- "CL28"
melted$Facet[melted$condition %>% grep(pattern = "IAA")] <- "IAA"

melted$Facet <- melted$Facet %>%
  factor(levels = c("NB","IAA","CL28"))


melted$BacteriaLong <- melted$condition <- melted$condition %>%
  gsub(pattern = ".*\\+",replacement = "") %>%
  gsub(pattern = " ",replacement = "") %>%
  gsub(pattern = "IAA",replacement = "NB") %>%
  gsub(pattern = "^CL28$",replacement = "NB") %>%
  factor(levels = c("NB","CL14","Δ68-69","Δ68-70","ΔHS33",
                    "pBBR::68-69","pBBR::68-70","ΔHS33pBBR1:68-69","ΔHS33pBBR1:68-70"))


melted$BacteriaLong <- melted$BacteriaLong %>% 
  gsub(pattern = "^pBBR::68-69",replacement = "ΔHS33pBBR1:68-69") %>%
  gsub(pattern = "pBBR::68-70",replacement = "ΔHS33pBBR1:68-70") %>%
  gsub(pattern = "ΔHS33p",replacement = "ΔHS33 p") %>%
  factor(levels = c("NB","CL14","ΔHS33","Δ68-69","Δ68-70",
                    "ΔHS33 pBBR1:68-69","ΔHS33 pBBR1:68-70"))


colnames(melted)[3] <- "Length_cm"

paleta <- c("#A6CEE3","#33A02C","#F5B452","#F5823D","#A81E0C","#6122F5","#00F5E1")
names(paleta) <- levels(melted$BacteriaLong)
#names(paleta) <- c("NB","CL14","ΔHS33","Δ68-69")

#a <- paletteer_d("awtools::spalette",n = 6,direction = -1)
#a <- a[c(1,3,6)]
#paleta <- c(paleta,a)
#noms <- melted$BacteriaLong %>% levels
#names(paleta)[5:7] <- noms[5:7]

df_cuenta <- melted
df_cuenta$Count <- 1
aggregate(Count~BacteriaLong+Facet,df_cuenta,sum)

## Replot using boxplot
p <- ggplot(data = melted,aes(BacteriaLong,Length_cm)) +
  geom_boxplot(outlier.colour = NA,outlier.size = NA,aes(color = BacteriaLong)) +
  geom_sina(aes(color = BacteriaLong),size = 1.5,stroke = 0.04,alpha = 0.3) +
  facet_grid(.~Facet,space = "free",scales = "free") +
  #geom_text(data = Res_em,aes(label = Letters,x = BacteriaLong, y = 5),inherit.aes = F) +
  theme_ohchibi(font_family = "Helvetica",size_panel_border = 0.75,
                size_axis_text.x = 12,size_axis_text.y = 12,size_axis_title.x = 15,size_axis_title.y = 15) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 15),
    panel.grid.major.y = element_line(color = "grey89",size = 0.25),
    panel.grid.major.x = element_line(color = "grey89",size = 0.25),
    panel.grid.minor.x = element_line(color = "grey89",size = 0.25)
  ) +
  ylab(label = "Primary root elongation (cm)") +
  xlab(label = "Bacterial treatment")  +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8),limits = c(0,8) ,oob = rescale_none) +
  scale_color_manual(values = paleta,name = "Bacterial\ntreatment")

oh.save.pdf(p = p,outname = "rgi_cl14_conditions.boxplot.pdf",outdir = "./figures/",width = 12,height = 8)



### Perform testing
df <- melted
df$Background <- df$Facet
mback <- df$Background %>% unique %>% as.character %>%
  grep(pattern = "NB",invert = T,value = T)
Res_em <- NULL
for(mb in mback){
  df_sub <- df %>% subset(Background == mb) %>% droplevels
  mtemp <- df_sub %>% 
    aov((Length_cm) ~ BacteriaLong,data = .) %>%
    emmeans(specs = "BacteriaLong") %>% CLD %>% 
    as.data.frame
  mtemp$Background <- mb
  Res_em <- rbind(Res_em,mtemp)
}

Res_em %>% subset(Background == "CL28")

melted$condition %>% table
melted$condition %>% unique %>% na.omit


#Save tables
df <- melted[,c("BacteriaLong","Facet","Length_cm")] %>%
  dplyr::rename(.data = ., Treatment = BacteriaLong , Background = Facet) %>%
  droplevels

write.table(x = df,file = "./cleandata/fig1e_data.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)  

rm(list=ls())
gc()
