library(ohchibi)
library(ggtree)

set.seed(130816)

#Prepare data for ITOL tree
#Read list of prevalence
Dat <- readRDS(file = "./cleandata/dat_prev_genus.RDS")
Map <- Dat$df_iaadeg_genus

Map$Size <- 0.5
Map$Size[which(Map$TotalGenus>6 & Map$TotalGenus <= 30)] <- 3
Map$Size[which(Map$TotalGenus>31 & Map$TotalGenus <= 100)] <- 5
Map$Size[which(Map$TotalGenus>101 & Map$TotalGenus <= 500)] <- 7
Map$Size[which(Map$TotalGenus>501 & Map$TotalGenus <= 5000)] <- 9
Map$Size[which(Map$TotalGenus>5000)] <- 11

Map$Genus %>% table %>% names
#Remove marine and bacterium from the phylogenetic tree
mnoms <- which(Map$Genus %in% c("bacterium","marine")) %>%
  Map$RepId[.]
tree <- read.tree(file = "./cleandata/reps_224.concatenated.aln.fasttree.wag_gamma.newick")
tree <- ape::drop.tip(phy = tree,tip = mnoms)
write.tree(phy = tree,file = "./cleandata/reps_222.concatenated.aln.fasttree.wag_gamma.newick")
write.tree(phy = tree,file = "./cleandata/fig2b_tree.newick")


Map <- which(!(Map$RepId %in% mnoms)) %>%
  Map[.,] %>% droplevels

#Map$itolName <- paste0(Map$Genus," (n=",Map$TotalGenus,", prev=",round(Map$PercPositive,2),", phyl=",
#       round(Map$RatioSampledMeanDistPhylo,2),")")
Map$itolName <- paste0(Map$Genus," (n=",Map$TotalGenus,")")

#LABELS
outfile <- "./cleandata//itol_reps_taxonoid2names.txt"
mline <- c("LABELS","Separator TAB","DATA")
write(x = mline,file = outfile,append = F,sep = "")
df <- data.frame(taxon_oid =Map$RepId,Names = Map$itolName)
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)


## Tree colors for labels
paleta_type <- c("#BEA2EB","#8AEB9F","#EBB573","#8AE7FA")
names(paleta_type) <- c("iacLike","varioLike","Hybrid","iacLikevarioLike")
Map$ColorLabel <- match(Map$ClassificationMajority,names(paleta_type)) %>%
  paleta_type[.] %>% as.character

outfile <- "./cleandata//itol_colors_labels.txt"
mline <- c("TREE_COLORS","Separator TAB","DATA")
write(x = mline,file = outfile,append = F,sep = "")
df <- data.frame(taxon_oid = Map$RepId,Type = "label_background",Color = Map$ColorLabel)
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)

#Prepare the bar charts one dffor prevalence
outfile <- "./cleandata/itol_bar_prevalence.txt"
mline <- c("DATASET_MULTIBAR","Separator TAB","DATASET_LABEL\tPrev","COLOR\tblack","FIELD_COLORS\t#BEA2EB\t#8AEB9F\t#EBB573\twhite","FIELD_LABELS\tiac\tvario\thybrid\tblank",
           "LEGEND_TITLE\tPrev","LEGEND_SHAPES\t1\t1\t1\t1","LEGEND_COLORS\t#BEA2EB\t#8AEB9F\t#EBB573\twhite","LEGEND_LABELS\tiac\tvario\thybrid\tblank",
           "BORDER_WIDTH\t2","WIDTH\t125","HEIGHT_FACTOR\t1.5","MARGIN\t10",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")
df <- data.frame(taxon_oid =Map$RepId,iac=Map$PerciacLike,vario=Map$PercvarioLike,hybrid=Map$PercHybrid,blank=Map$PercMissing)
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)



#Second barchart for phylogenetic breadth
outfile <- "./cleandata/itol_bar_phyl.txt"
mline <- c("DATASET_MULTIBAR","Separator TAB","DATASET_LABEL\tPhyl","COLOR\tred","FIELD_COLORS\t#D9D9D9","FIELD_LABELS\tRatio",
           "DATASET_SCALE\t0-0.0-black-1-1-5\t0.5-0.5-black-1-1-5\t1-1.0-black-1-1-5",
           "LEGEND_TITLE\tPhyl","LEGEND_SHAPES\t1","LEGEND_COLORS\t#D9D9D9","LEGEND_LABELS\tRatio",
           "BORDER_WIDTH\t2","WIDTH\t125","HEIGHT_FACTOR\t1.5","MARGIN\t10",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")
df <- data.frame(taxon_oid =Map$RepId,Ratio=Map$RatioSampledMeanDistPhylo)
df$Ratio[which(df$Ratio>1.0)] <- 1.1
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)


#Create  a legend for the dataset symbol based on size of the node
outfile <- "./cleandata/itol_shape_tip.txt"
mline <- c("DATASET_SYMBOL","Separator TAB","DATASET_LABEL\tTip","COLOR\tblue",
           "LEGEND_TITLE\tTip","LEGEND_SHAPES\t2","LEGEND_COLORS\t#D9D9D9","LEGEND_LABELS\tSize",
           "MAXIMUM_SIZE\t20",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")
df <- data.frame(taxon_oid =Map$RepId,Shape = 2,Size=Map$Size,Color="rgba(0,0,0,0.2)",Fill =1,Position = 1)
#df$Ratio[which(df$Ratio>1.0)] <- 1.05
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)


#Create a colorstrip dataset for the size so it is easier to compare
paleta_size <- c("#ffffcc","#ffeda0","#feb24c","#fc4e2a","#e31a1c","#800026")
names(paleta_size) <- c("0.5","3","5","7","9","11")
outfile <- "./cleandata/itol_strip_size.txt"
mline <- c("DATASET_COLORSTRIP","Separator TAB","DATASET_LABEL\tStripSize","COLOR\tpurple",
           "LEGEND_TITLE\tSize","LEGEND_SHAPES\t1\t1\t1\t1\t1\t1","LEGEND_COLORS\t#F6F5E2\t#ffeda0\t#feb24c\t#fc4e2a\t#e31a1c\t#800026",
           "LEGEND_LABELS\t0-5\t6-30\t31-100\t101-500\t501-5000\t>5000",
           "BORDER_WIDTH\t2",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")
Map$ColorLabel <- match(Map$Size,names(paleta_size)) %>%
  paleta_size[.] %>% as.character
df <- data.frame(taxon_oid =Map$RepId,Color=Map$ColorLabel)
#df$Ratio[which(df$Ratio>1.0)] <- 1.05
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)

rm(list=ls())
gc()
