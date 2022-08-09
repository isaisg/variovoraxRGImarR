library(ohchibi)

set.seed(130816)

Dat_binary <- readRDS(file = "./cleandata/dat_variovoraxRGI_marR_pco.RDS")


distfun <- function(x,method) vegan::vegdist(x = x, method = "jaccard")

mpco <- oh.pco(Tab = Dat_binary$Tab %>% t,
               Map = Dat_binary$Map,
               ndim = 3,eig = T,
               distfun = distfun,
               id_var = "UId")



p <- ggplot(data = mpco$Map_pco ,aes(PCo1,PCo2)) + 
  geom_vline(xintercept = 0,linetype = "longdash",color = "#D9D9D9",size = 1) +
  geom_hline(yintercept = 0,linetype = "longdash",color = "#D9D9D9",size = 1) +
  geom_point(data = mpco$Map_pco %>% subset(Representative == "No"),
             aes(color = ColorGenus),alpha = 0.75,stroke = 0,size = 2) +
  geom_point(data = mpco$Map_pco %>% subset(Representative == "Yes"),
             aes(color = ColorGenus),alpha = 1,stroke = 0,size = 4,shape = 15) + 
  scale_color_paletteer_d("ggthemes::Tableau_20",name = "Taxonomic\nclassification") +
  theme_ohchibi(size_panel_border = 0.75) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )  +
  xlab(label = "PCo1 (11.05%)") +
  ylab(label = "PCo2 (3.96%)")

oh.save.pdf(p = p,outname = "pco_end.pdf",outdir = "./figures/",
            width = 10,height = 6)

#Save data for reproducibility
Tab <- Dat_binary$Tab
Map <- mpco$Map_pco



#Save files
write.table(x = Tab,file = "./cleandata/fig3a_matrix.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(x = Map,file = "./cleandata/fig3a_metadata_projection.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)

rm(list=ls())
gc()
