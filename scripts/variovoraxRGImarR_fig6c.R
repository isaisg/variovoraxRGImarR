library(ohchibi)
library(emmeans)

set.seed(130816)

#######################
### Rep 1 #############
df <- read.table(file = "./rawdata/data_dropin_rep1.csv",header = T,sep = ",")
nfactor <- 2920
df$Length_cm <- (df$Length_Pixel*12)/(nfactor)


df$Background <- df$Treatment %>% gsub(pattern = "_.*",replacement = "") %>%
  gsub(pattern = "Mono",replacement = "Non RGI treatment") %>%
  factor(levels = c("Non RGI treatment","IAA","SC33","SC175"))
df$Bacteria <- df$Treatment %>% gsub(pattern = ".*_",replacement = "") %>%
  gsub(pattern = "CL14Delta",replacement = "CL14Delta6869") %>%
  factor(levels =c("NB","10Vario","CL14CL69","CL14","CL14Delta6869","Brady",
                   "Pse1290","Pse562","Paraburk","Esoli","AcineCL69","AcineCL71"))

df$Rep <- "Rep1"
df_rep1 <- df

### Test inside each background
mback <- df$Background %>% unique %>% as.character
Res_em <- NULL
for(mb in mback){
  df_sub <- df %>% subset(Background == mb) %>% droplevels
  mtemp <- df_sub %>% 
    aov((Length_cm) ~ Bacteria,data = .) %>%
    emmeans(specs = "Bacteria") %>% multcomp::cld() %>% 
    as.data.frame
  mtemp$Background <- mb
  Res_em <- rbind(Res_em,mtemp)
}

#Transform letters so it is concistent
Res_em$.group <- Res_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em)[which(colnames(Res_em) == ".group")] <- "group"
#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>% 
  strsplit(split = "") %>% 
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop 
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>% 
    unlist %>% toupper  %>% paste0(collapse = "")
}

Res_em$Background <- Res_em$Background %>% factor(levels = df$Background %>% levels)
Res_em$Bacteria <- Res_em$Bacteria %>% factor(levels = df$Bacteria %>% levels)

## Adjust the names so they are already meaningful
df$BacteriaLong <- df$Bacteria %>%
  gsub(pattern = "10Vario",replacement = "10 Variovorax") %>%
  gsub(pattern = "^CL14$",replacement = "Variovorax CL14") %>%
  gsub(pattern = "CL14CL69",replacement = "CL14 + CL69") %>%
  gsub(pattern = "CL14Delta",replacement = "Variovorax CL14 Delta") %>%
  gsub(pattern = "Pse1290",replacement = "Pseudomonas 1290") %>%
  gsub(pattern = "Pse562",replacement = "Pseudomonas 562") %>%
  gsub(pattern = "Brady",replacement = "Bradyrhizobium japonicum") %>%
  gsub(pattern = "Esoli",replacement = "Enterobacter soli") %>%
  gsub(pattern = "Paraburk",replacement = "Parabukholderia PsJN") %>%
  gsub(pattern = "AcineCL69",replacement = "Acinetobacter CL69") %>%
  gsub(pattern = "AcineCL71",replacement = "Acinetobacter CL71")

order_long <- with(df,order(Bacteria)) %>%
  df$BacteriaLong[.] %>% as.character %>% unique
df$BacteriaLong <- df$BacteriaLong %>% factor(levels = order_long)

Res_em$BacteriaLong <- Res_em$Bacteria %>%
  gsub(pattern = "10Vario",replacement = "10 Variovorax") %>%
  gsub(pattern = "^CL14$",replacement = "Variovorax CL14") %>%
  gsub(pattern = "CL14CL69",replacement = "CL14 + CL69") %>%
  gsub(pattern = "CL14Delta",replacement = "Variovorax CL14 Delta") %>%
  gsub(pattern = "Pse1290",replacement = "Pseudomonas 1290") %>%
  gsub(pattern = "Pse562",replacement = "Pseudomonas 562") %>%
  gsub(pattern = "Brady",replacement = "Bradyrhizobium japonicum") %>%
  gsub(pattern = "Esoli",replacement = "Enterobacter soli") %>%
  gsub(pattern = "Paraburk",replacement = "Parabukholderia PsJN") %>%
  gsub(pattern = "AcineCL69",replacement = "Acinetobacter CL69") %>%
  gsub(pattern = "AcineCL71",replacement = "Acinetobacter CL71")
Res_em$BacteriaLong <- Res_em$BacteriaLong %>% factor(levels = order_long)

## Create a palette
paleta <- c("#A6CEE3","black","#B8B50F","#33A02C","#F0A928",
            "#B2DF8A",
            "#E31A1C","#FB9A99","#6A3D9A","#CAB2D6","#FFFF99","#B15928")
names(paleta) <- df$BacteriaLong %>% levels


###
p <- ggplot(data = df,aes(BacteriaLong,Length_cm)) +
  geom_sina(shape = 21,stroke = 0.04,size = 2) +
  stat_summary(fun.data = "mean_cl_normal",geom = "pointrange",aes(color = BacteriaLong),size = 0.75) +
  facet_grid(.~Background,space = "free",scales = "free") +
  geom_text(data = Res_em,aes(label = Letters,x = BacteriaLong, y = 6),inherit.aes = F) +
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
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6),limits = c(0,6)) +
  scale_color_manual(values = paleta,name = "Bacterial treatment")

df_uno <- df %>% subset(Background == "Non RGI treatment" | Background == "IAA") %>% droplevels

pr1 <- ggplot(data = df %>% subset(Background == "Non RGI treatment" | Background == "IAA") %>% droplevels,aes(BacteriaLong,Length_cm)) +
  geom_sina(shape = 21,stroke = 0.04,size = 2) +
  stat_summary(fun.data = "mean_cl_normal",geom = "pointrange",aes(color = BacteriaLong),size = 0.75) +
  facet_grid(.~Background,space = "free",scales = "free") +
  geom_text(data = Res_em %>% subset(Background == "Non RGI treatment" | Background == "IAA") %>%droplevels,aes(label = Letters,x = BacteriaLong, y = 6),inherit.aes = F) +
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
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6),limits = c(0,6)) +
  scale_color_manual(values = paleta,name = "Bacterial treatment")

p1 <- ggplot(data = df %>% subset(Background == "SC175"),aes(BacteriaLong,Length_cm)) +
  geom_sina(shape = 21,stroke = 0.04,size = 2) +
  stat_summary(fun.data = "mean_cl_normal",geom = "pointrange",aes(color = BacteriaLong),size = 0.75) +
  facet_grid(.~Background,space = "free",scales = "free") +
  geom_text(data = Res_em %>% subset(Background == "SC175"),aes(label = Letters,x = BacteriaLong, y = 5),inherit.aes = F) +
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
  scale_y_continuous(breaks = c(0,1,2,3,4,5),limits = c(0,5.1)) +
  scale_color_manual(values = paleta,name = "Bacterial treatment")


#########################
####################Rep2
df <- read.table(file = "./rawdata/data_dropin_rep2.csv",header = T,sep = ",")
nfactor <- 2920
df$Length_cm <- (df$Length*12)/(nfactor)

df$Background <- df$Background %>%
  gsub(pattern = "NB",replacement = "Non RGI treatment")
df$Rep <- "Rep2"
df_rep2 <- df

df <- df %>% subset(Background == "SC175") %>% droplevels
df$Bacteria <- df$Bacteria %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels =c("NB","10Vario","CL14","CL14Delta6869","CL14Delta73","CL14Delta73Comp","Brady",
                   "Pse1290","Pse562","Paraburk","Esoli","AcineCL69","AcineCL71"))

### Test inside each background
mback <- df$Background %>% unique %>% as.character
Res_em <- NULL
for(mb in mback){
  df_sub <- df %>% subset(Background == mb) %>% droplevels
  mtemp <- df_sub %>% 
    aov((Length_cm) ~ Bacteria,data = .) %>%
    emmeans(specs = "Bacteria") %>% multcomp::cld() %>% 
    as.data.frame
  mtemp$Background <- mb
  Res_em <- rbind(Res_em,mtemp)
}

#Transform letters so it is concistent
Res_em$.group <- Res_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em)[which(colnames(Res_em) == ".group")] <- "group"
#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>% 
  strsplit(split = "") %>% 
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop 
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>% 
    unlist %>% toupper  %>% paste0(collapse = "")
}

Res_em$Bacteria <- Res_em$Bacteria %>% factor(levels = df$Bacteria %>% levels)

## Adjust the names so they are already meaningful
df$BacteriaLong <- df$Bacteria %>%
  gsub(pattern = "10Vario",replacement = "10 Variovorax") %>%
  gsub(pattern = "^CL14$",replacement = "Variovorax CL14") %>%
  gsub(pattern = "CL14Delta6869",replacement = "Variovorax CL14 Delta6869") %>%
  gsub(pattern = "CL14Delta73",replacement = "Variovorax CL14 Delta 73") %>%
  gsub(pattern = "CL14Delta73Comp",replacement = "Variovorax CL14 Delta 73 Comp") %>%
  gsub(pattern = "Pse1290",replacement = "Pseudomonas 1290") %>%
  gsub(pattern = "Pse562",replacement = "Pseudomonas 562") %>%
  gsub(pattern = "Brady",replacement = "Bradyrhizobium japonicum") %>%
  gsub(pattern = "Esoli",replacement = "Enterobacter soli") %>%
  gsub(pattern = "Paraburk",replacement = "Parabukholderia PsJN") %>%
  gsub(pattern = "AcineCL69",replacement = "Acinetobacter CL69") %>%
  gsub(pattern = "AcineCL71",replacement = "Acinetobacter CL71")


order_long <- with(df,order(Bacteria)) %>%
  df$BacteriaLong[.] %>% as.character %>% unique
df$BacteriaLong <- df$BacteriaLong %>% factor(levels = order_long)

Res_em$BacteriaLong <- Res_em$Bacteria %>%
  gsub(pattern = "10Vario",replacement = "10 Variovorax") %>%
  gsub(pattern = "^CL14$",replacement = "Variovorax CL14") %>%
  gsub(pattern = "CL14Delta6869",replacement = "Variovorax CL14 Delta6869") %>%
  gsub(pattern = "CL14Delta73",replacement = "Variovorax CL14 Delta 73") %>%
  gsub(pattern = "CL14Delta73Comp",replacement = "Variovorax CL14 Delta 73 Comp") %>%
  gsub(pattern = "Pse1290",replacement = "Pseudomonas 1290") %>%
  gsub(pattern = "Pse562",replacement = "Pseudomonas 562") %>%
  gsub(pattern = "Brady",replacement = "Bradyrhizobium japonicum") %>%
  gsub(pattern = "Esoli",replacement = "Enterobacter soli") %>%
  gsub(pattern = "Paraburk",replacement = "Parabukholderia PsJN") %>%
  gsub(pattern = "AcineCL69",replacement = "Acinetobacter CL69") %>%
  gsub(pattern = "AcineCL71",replacement = "Acinetobacter CL71")
Res_em$BacteriaLong <- Res_em$BacteriaLong %>% factor(levels = order_long)

df$BacteriaLong %>% levels
names(paleta)[3] <- "Variovorax CL14 Delta 73"
z <- "#73722C"
names(z) <- "Variovorax CL14 Delta 73Comp"
paleta <- c(paleta,z)

###
p2 <- ggplot(data = df,aes(BacteriaLong,Length_cm)) +
  geom_sina(shape = 21,stroke = 0.04,size = 2) +
  stat_summary(fun.data = "mean_cl_normal",geom = "pointrange",aes(color = BacteriaLong),size = 0.75) +
  facet_grid(.~Background,space = "free",scales = "free") +
  geom_text(data = Res_em,aes(label = Letters,x = BacteriaLong, y = 5),inherit.aes = F) +
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
  scale_y_continuous(breaks = c(0,1,2,3,4,5),limits = c(0,5.1)) +
  scale_color_manual(values = paleta,name = "Bacterial treatment")

composition <- egg::ggarrange(p1+ theme(legend.position = "none"),p2,nrow = 1)


###### Merged analysis #########
colnames(df_rep1)
colnames(df_rep2)

##Count samples ###

df_cuenta <- rbind(df_rep1[,c("Bacteria","Background","Length_cm","Rep")],
                   df_rep2[,c("Bacteria","Background","Length_cm","Rep")])
df_cuenta$Count <- 1  
aggregate(Count~Background+Bacteria,df_cuenta,sum)



df <- rbind(df_rep1[,c("Bacteria","Background","Length_cm","Rep")] %>% subset(Background == "SC175"),
            df_rep2[,c("Bacteria","Background","Length_cm","Rep")] %>% subset(Background == "SC175")
) %>% droplevels

df$Bacteria %>% unique


chosen_bacteria <- df[,c("Bacteria","Rep")] %>% unique %$% Bacteria %>% table %>%
  data.frame %>% subset(Freq == 2) %$% . %>% as.character

df <- which(df$Bacteria %in% chosen_bacteria) %>%
  df[.,] %>% droplevels

#Order the factors
df$Bacteria <- df$Bacteria %>%
  factor(levels =c("NB","10Vario","CL14","CL14Delta6869","CL14Delta73","CL14Delta73Comp","Brady",
                   "Pse1290","Pse562","Paraburk","Esoli","AcineCL69","AcineCL71")) %>%
  droplevels


Res_em <- aov(formula = Length_cm~Bacteria+Rep,data = df) %>%
  emmeans(specs = "Bacteria") %>% multcomp::cld() %>% 
  as.data.frame

#Transform letters so it is concistent
Res_em$.group <- Res_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em)[which(colnames(Res_em) == ".group")] <- "group"
#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>% 
  strsplit(split = "") %>% 
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop 
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>% 
    unlist %>% toupper  %>% paste0(collapse = "")
}


Res_em$Bacteria <- Res_em$Bacteria %>% factor(levels = df$Bacteria %>% levels)

## Adjust the names so they are already meaningful
df$BacteriaLong <- df$Bacteria %>%
  gsub(pattern = "10Vario",replacement = "10 Variovorax") %>%
  gsub(pattern = "^CL14$",replacement = "Variovorax CL14") %>%
  gsub(pattern = "CL14Delta6869",replacement = "Variovorax CL14 Delta6869") %>%
  gsub(pattern = "CL14Delta73",replacement = "Variovorax CL14 Delta 73") %>%
  gsub(pattern = "CL14Delta73Comp",replacement = "Variovorax CL14 Delta 73 Comp") %>%
  gsub(pattern = "Pse1290",replacement = "Pseudomonas 1290") %>%
  gsub(pattern = "Pse562",replacement = "Pseudomonas 562") %>%
  gsub(pattern = "Brady",replacement = "Bradyrhizobium japonicum") %>%
  gsub(pattern = "Esoli",replacement = "Enterobacter soli") %>%
  gsub(pattern = "Paraburk",replacement = "Parabukholderia PsJN") %>%
  gsub(pattern = "AcineCL69",replacement = "Acinetobacter CL69") %>%
  gsub(pattern = "AcineCL71",replacement = "Acinetobacter CL71")

order_long <- with(df,order(Bacteria)) %>%
  df$BacteriaLong[.] %>% as.character %>% unique
df$BacteriaLong <- df$BacteriaLong %>% factor(levels = order_long)

Res_em$BacteriaLong <- Res_em$Bacteria %>%
  gsub(pattern = "10Vario",replacement = "10 Variovorax") %>%
  gsub(pattern = "^CL14$",replacement = "Variovorax CL14") %>%
  gsub(pattern = "CL14Delta6869",replacement = "Variovorax CL14 Delta6869") %>%
  gsub(pattern = "CL14Delta73",replacement = "Variovorax CL14 Delta 73") %>%
  gsub(pattern = "CL14Delta73Comp",replacement = "Variovorax CL14 Delta 73 Comp") %>%
  gsub(pattern = "Pse1290",replacement = "Pseudomonas 1290") %>%
  gsub(pattern = "Pse562",replacement = "Pseudomonas 562") %>%
  gsub(pattern = "Brady",replacement = "Bradyrhizobium japonicum") %>%
  gsub(pattern = "Esoli",replacement = "Enterobacter soli") %>%
  gsub(pattern = "Paraburk",replacement = "Parabukholderia PsJN") %>%
  gsub(pattern = "AcineCL69",replacement = "Acinetobacter CL69") %>%
  gsub(pattern = "AcineCL71",replacement = "Acinetobacter CL71")

Res_em$BacteriaLong <- Res_em$BacteriaLong %>% factor(levels = order_long)


p <- ggplot(data = df,aes(BacteriaLong,Length_cm)) +
  geom_jitter(aes(shape = Rep,color = BacteriaLong),size = 1.5,stroke = 0.04) +
  scale_shape_manual(values = c(21,24)) +
  stat_summary(fun.data = "mean_cl_normal",geom = "pointrange",aes(color = BacteriaLong),size = 0.75) +
  facet_grid(.~Background,space = "free",scales = "free") +
  geom_text(data = Res_em,aes(label = Letters,x = BacteriaLong, y = 6),inherit.aes = F) +
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
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6),limits = c(0,6)) +
  scale_color_manual(values = paleta,name = "Bacterial treatment")

df_dos <- df

#Save table
df_end <- rbind(match(colnames(df_dos),colnames(df_uno)) %>%
  df_uno[,.],
  df_dos)

df_end <- df_end[,-1] %>% 
  dplyr::rename(.data = .,Bacteria = BacteriaLong) %>%
  dplyr::select(.data = .,c("Background","Bacteria","Rep","Length_cm"))

write.table(x = df_end,file = "./cleandata/fig6c_data.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

merged <- egg::ggarrange(pr1 + theme(legend.position = "none"),p + theme(legend.position = "none",axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()),nrow =1,widths = c(1,0.5))


oh.save.pdf(p = merged,outname = "rgi_iaa_operons_merged.final.pdf",outdir = "./figures/",
            width = 12,height = 7)


rm(list=ls())
gc()
dev.off()
