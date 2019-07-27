rm(list = objects()); ls()

library(tidyverse)
library(tidylog)
library(data.table)
library(janitor)

amKeyFile<- fread("./GenoDatabase/AMPanel_gbs.csv")
bpKeyFile<- fread("./GenoDatabase/BP_KeyFiles/PYN_2019.txt")

glumeKeyFile<- amKeyFile %>% 
  bind_rows(bpKeyFile) %>% 
  mutate(FullSampleName = str_replace(FullSampleName,
                                      "Smoky Hill","smokyhill")) %>% 
  mutate(FullSampleName = str_replace_all(FullSampleName,
                                      " ","")) %>% 
  mutate(FullSampleName = str_replace_all(FullSampleName,
                                      "^","-")) %>% 
  mutate(FullSampleName = tolower(FullSampleName))

write.table(glumeKeyFile, "./GenoDatabase/GlumeMaster.txt", quote = F,
            sep = "\t", row.names = F, col.names = T)

###############################################################################

snpChip <- read_delim(
  "./GenoDatabase/Glume_imputed.hmp.txt", 
  "\t", escape_double = FALSE, trim_ws = TRUE)
snpChip<- snpChip %>% 
  clean_names()
snpChip<- separate(snpChip,alleles,c("allele_a","allele_b"), sep = "/")

missAmbiguous = c('0', '+', '-')
hetCodes = c('R','Y','S','W','K','M','B','D','H','V')
hapgeno=as.matrix(snpChip[,13:ncol(snpChip)])
hapgeno[hapgeno %in% missAmbiguous]=NA
hapgeno[hapgeno=='N']=NA
hapgeno[hapgeno %in% hetCodes]='H'
snpChip=cbind(snpChip[,1:12], hapgeno)
rm(hapgeno)

write.table(snpChip, file="./GenoDatabase/Glume_imputedBeagle.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./GenoDatabase/Glume_imputedBeagle.txt", 
                header=TRUE, check.names=F, sep = "\t")

snpChip[snpChip == snpChip$allele_a] = -1
snpChip[snpChip == snpChip$allele_b] = 1
snpChip[snpChip == "H"] = 0
snpChip[snpChip == "C"] = NA
snpChip[snpChip == "A"] = NA
snpChip[snpChip == "T"] = NA
snpChip[snpChip == "G"] = NA
snpChip[snpChip == "-"] = NA
snpChip[snpChip == "."] = NA

snpChip<- snpChip[ ,c(1,4,5,13:1210)]

write.table(snpChip, file="./GenoDatabase/Glume_imputedNumeric.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./GenoDatabase/Glume_imputedNumeric.txt", 
                header=TRUE, check.names=F, sep = "\t")

chrSum<- plyr::count(snpChip, vars = "chrom")
snpChip<- snpChip %>% 
  filter(chrom != "UN")
snpMatrix<- t(snpChip[ , c(-1, -2, -3)])

pcaMethods::checkData(snpMatrix)  #Check PCA assumptions

pcaAM<- pcaMethods::pca(snpMatrix, nPcs = 15) #SVD PCA

sumPCA<- as.data.frame(summary(pcaAM))

Scores<- as.data.frame(pcaMethods::scores(pcaAM)) 
Scores<- setDT(Scores, keep.rownames = TRUE)

Scores<- Scores %>% 
  mutate(Year = str_sub(rn,3,4),
         Year = as.numeric(Year),
         Year = replace_na(Year,"Unknown"),
         Year = str_replace(Year,"44","Unknown"),
         dh = str_detect(rn,"dh"))

Scores %>% 
  ggplot(aes(x = PC1, y = PC2, colour = factor(Year), shape = dh)) +
  geom_point(alpha = 1) +
  scale_color_manual(name = "Year", 
                     values = c('#1b9e77','#d95f02','#7570b3',
                                '#e7298a','#66a61e','#e6ab02',
                                '#a6761d','#666666','#000000')) +
  scale_shape_manual(name = "DoubleHaploid", 
                     values = c(17,19)) +
  theme_bw() +
  theme(axis.title = element_text(colour = "black", size = 16),
        axis.text = element_text(colour = "black", size = 14),
        aspect.ratio = 1:1,
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid.major = element_line(colour = "#b7b7b7",linetype = 2),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(colour = "black", size = 14)) +
  labs(title = "PCA of Glume Colour Lines",
       x = expression(paste("PC1 ",R^2, " = 4.7%")),
       y =  expression(paste("PC2 ",R^2, " = 3.8%")))