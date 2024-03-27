##########
# load packages
library(stats)
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(AICcmodavg)
library(grid)
library(forestploter)
library(lm.beta)

##############
# preprocess #
##############
TotScoreBySample <- readRDS('data/TotScoreBySample.rds')
TotScoreBySample$projid <- as.numeric(TotScoreBySample$projid)
phenotypes <- read.csv(file = 'data/raw/dataset_978_basic_01-13-2022.csv', header = TRUE)

# remove white space
phenotypes$study <- substring(phenotypes$study, 1,3)

# table with every sample's scores & phenotype
scoreWithPheno <- merge(phenotypes, TotScoreBySample, by='projid')
scoreWithPheno <- scoreWithPheno[which(scoreWithPheno$study.x==scoreWithPheno$study.y),]
scoreWithPheno$study.y <-NULL
colnames(scoreWithPheno)[2] <- 'study'

######
# PC #
######
pc_egval <- read.csv(file = 'data/rosmap_pca.eigenval', header = FALSE)
pc_egvec <- read.table(file = 'data/rosmap_pca.eigenvec', sep = '\t', header = FALSE)

# coverting AD to ROSMAP id
ROSmaster <- readRDS("data/ROSmaster.rds")
ROSmaster_keylist <- ROSmaster[c('projid', 'study', 'IID')]
rm(ROSmaster)

ROSmaster_keylist$study <- as.character(ROSmaster_keylist$study)
ROSmaster_keylist$projid <- as.integer(ROSmaster_keylist$projid)
ROSmaster_keylist$study <- gsub('\\s+', '', ROSmaster_keylist$study)

pc_egvec$V1 <- NULL
AD <- pc_egvec[which(grepl('AD', pc_egvec$V2)),]

temp <- merge(AD, ROSmaster_keylist, by.x = "V2", by.y = "IID")
temp$V2<-NULL
temp <- cbind(temp[,c(11,12)], temp[,c(1:10)])

temp2 <- pc_egvec[-which(grepl('AD', pc_egvec$V2)), ]
study <- substring(temp2$V2, 1, 3)
projid <- as.numeric(substring(temp2$V2, 4))
temp2 <- cbind(projid, study, temp2[,c(2:11)])

# missing 3
pc_egvec <- rbind(temp, temp2)
colnames(pc_egvec) <- c('projid', 'study', 'egvec1', 'egvec2', 'egvec3',
                        'egvec4', 'egvec5','egvec6', 'egvec7', 'egvec8',
                        'egvec9', 'egvec10')

scoreWithPheno <- merge(scoreWithPheno, pc_egvec, by=c("projid", "study"))

rm(temp, temp2, study, AD, pc_egvec, pc_egval)

#########################
# phenotype_long <- read.csv(file = 'data/dataset_978_long_01-13-2022.csv', header = TRUE)
# all.last <- lapply(unique(phenotype_long$projid), function(id){
#   ss <- subset(phenotype_long, projid==id & !is.na(cogn_global))
#   max.fu <- max(ss$fu_year)
#   return(ss[which(ss$fu_year == max.fu),])
# })
#
# all.last <- do.call(rbind,all.last)
#
# saveRDS(all.last, 'pheno_long_lastobs.rds')

pheno_long <- readRDS('data/pheno_long_lastobs.rds')
pheno_long$study <- substring(pheno_long$study, 1, 3) # remove white space
scoreWithPheno <- merge(scoreWithPheno, pheno_long, by=c('projid'))

## remove that 1 outlier
outlier <- which(grepl(max(scoreWithPheno$loeuf_dup), scoreWithPheno$loeuf_dup))
scoreWithPheno <- scoreWithPheno[-outlier, ]

# make a new column named cogdx_binom: cogdx = {1,2} -> 0; cogdx = {4} -> 1;
scoreWithPheno$cogdx_binom <- scoreWithPheno$cogdx
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==1)] <- 0
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==2)] <- 0
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==3)] <- NA
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==4)] <- 1
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==5)] <- NA

saveRDS(scoreWithPheno, 'data/scoreWithPheno.rds')
