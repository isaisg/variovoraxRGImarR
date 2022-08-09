library(ohchibi)
library(DESeq2)

set.seed(130816)

df <- read.table(file = "./rawdata/resfc.0.tsv",header = T,
                 sep = "\t",row.names = 1,check.names = F,quote = "",comment.char = "")

colnames(df)  <- colnames(df) %>% 
  gsub(pattern = ".*bowtie2\\.",replacement = "") %>%
  gsub(pattern = "\\.sam",replacement = "")

#Create Map
Map <- data.frame(RawId = colnames(df) )
Map$Treatment <- Map$RawId %>% gsub(pattern = "_S[0-9].*",replacement = "") %>%
  gsub(pattern = "^[0-9]_[0-9]_",replacement = "") %>%
  gsub(pattern = "_[ACTG][ACTG][ACTG].*",replacement = "")
Map$Bacteria <- Map$Treatment %>% gsub(pattern = "_minus|_IAA",replacement = "")
Map$IAA <- "No"
Map$IAA[Map$Treatment %>% grep(pattern = "IAA")] <- "Yes"
Map$SampleId <- Map$RawId %>% gsub(pattern = ".*_S",replacement = "") %>%
  gsub(pattern = "_.*",replacement = "") %>%
  paste0("S",.)
colnames(df) <- match(colnames(df),Map$RawId) %>%
  Map$SampleId[.]
rownames(Map) <- Map$SampleId
Map <- Map[,c("SampleId","Bacteria","IAA","Treatment","RawId")]

Tab <- df %>% as.matrix
Tab <- Tab[which(rowSums(Tab)!=0),]

Dat <- create_dataset(Tab = Tab,Map = Map)
Dat <- subset.Dataset(x = Dat,subset = Treatment != "NTC",drop = T,clean = T)
Dat$Map$Bacteria <- Dat$Map$Bacteria %>% factor
Dat$Map$IAA <- Dat$Map$IAA %>% factor
Dat$Map$Treatment <- Dat$Map$Treatment %>% factor

### Perform Deseq model
dds <- DESeqDataSetFromMatrix(countData = Dat$Tab,
                              colData = Dat$Map,
                              design = ~Treatment)
dds <- DESeq(object = dds,fitType = "local")

#Create object to plot
mvst <- vst(object = dds,blind = F,fitType = "local")

mat <- assay(mvst)

Tab_z <- mat %>% t %>% scale (center = T,scale = F)

Dat_z <- create_dataset(Tab = Tab_z %>% t,Map = Dat$Map)

#Create summary object for plotting
melted <- Dat_z$Tab %>% melt
colnames(melted) <- c("Gene","SampleId","value")
melted <- merge(melted,Dat_z$Map, by = "SampleId")
Tab_mean <- acast(data = melted,formula = Gene~Treatment,
                  fun.aggregate = mean,value.var = "value") 
melted_av <- Tab_mean %>% melt
colnames(melted_av) <- c("Gene","Treatment","value")
melted_av <- merge(melted_av,Dat_z$Map[,c("Treatment","Bacteria","IAA")] %>% unique, 
                   by = "Treatment")

###LRT model
### Perform Deseq model
dds_lrt <- DESeqDataSetFromMatrix(countData = Dat$Tab,
                                  colData = Dat$Map,
                                  design = ~Treatment)
dds_lrt <- DESeq(object = dds_lrt,
                 fitType = "local",test = "LRT",reduced = ~1)



mlist <- list(
  Dat_raw = Dat,
  Dat_z = Dat_z,
  dds = dds,
  melted_av = melted_av,
  dds_lrt = dds_lrt
)
saveRDS(object = mlist,
        file = "./cleandata/dat_deseq2_variovoraxmarRRNASeq.RDS")
rm(list=ls())
gc()
