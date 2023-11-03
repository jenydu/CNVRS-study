library(ggplot2)
library(dplyr)
library(tidyr)
library(geomtextpath)
library(ggpubr)
library(gridExtra)

#######################################################
CNVcalls <- read.csv(file = 'data/raw/CNVcalls.tsv',
                       sep = '\t', header = TRUE)

# del_count <- colSums(CNVcalls[which(CNVcalls$SVTYPE=='DEL'),]!= 2)[-(1:6)]
# dup_count <- colSums(CNVcalls[which(CNVcalls$SVTYPE=='DUP'),]!= 2)[-(1:6)]

del <- 2 - CNVcalls[which(CNVcalls$SVTYPE=='DEL'),-(1:6)]
dup <- CNVcalls[which(CNVcalls$SVTYPE=='DUP'),-(1:6)] -2
del_count <- colSums(del)
dup_count <- colSums(dup)

outlier <- which(grepl(max(dup_count), dup_count))
dup_count <- dup_count[-outlier]
del_count <- del_count[-outlier]
plot(del_count, dup_count, pch=20, title('ROSMAP (with copy number as multiplier)'))
Hmisc::rcorr(del_count, dup_count)

#### load annotations (VEP) ####
VEP1 <- read.csv(file = 'data/raw/01_04_2023_VEP_output_filtered.txt',
                   sep = '\t', header = TRUE)

VEP1 <- VEP1[c(2, 3, 6)]
VEP<-VEP1[!duplicated(VEP1),]
colnames(VEP) <- c('location', 'SVTYPE', 'gene')
VEP <- VEP[!grepl("-",VEP$gene),]

VEP2<-tidyr::separate(data = VEP, col = location, into = c("CHROM", "loc"), sep = "\\:")
VEP2<-tidyr::separate(data = VEP2, col = loc, into = c("START", "END"), sep = "\\-")
VEP2[VEP2 == 'deletion'] <- 'DEL'
VEP2[VEP2 == 'duplication'] <- 'DUP'
VEP2[, c(1:3)] <- sapply(VEP2[, c(1:3)], as.integer)

VEP <- VEP2

rm(VEP1, VEP2)


####

haplo_triplo <- readRDS("data/pHaplo_pTriplo_data.rds")
colnames(haplo_triplo) <- c('gene', 'pHaplo', 'pTriplo')

pli_loeuf <- readRDS("data/pLI_LOEUF_data.rds")
pli_loeuf <- pli_loeuf[, c(1, 21, 30)]
pli_loeuf$oe_lof_upper <- pli_loeuf$oe_lof_upper ** (-1)

pHaploLst<-haplo_triplo
pHaploLst$pTriplo = NULL

pTriploLst<-haplo_triplo
pTriploLst$pHaplo = NULL

rm(haplo_triplo)

########## thresholding ###################
# pHaplo >= 0.86
thresh1 <- 0.86
pHaplo_thresh <- pHaploLst
pHaplo_thresh[which(pHaplo_thresh$pHaplo < thresh1),][2] <- 0
pHaplo_thresh[which(pHaplo_thresh$pHaplo >= thresh1),][2] <- 1

# pTriplo >= 0.94
thresh2 <- 0.94
pTriplo_thresh <- pTriploLst
pTriplo_thresh[which(pTriplo_thresh$pTriplo < thresh2),][2] <- 0
pTriplo_thresh[which(pTriplo_thresh$pTriplo >= thresh2),][2] <- 1

rm(thresh1, thresh2)

############## gene count ################

geneContentTable <- crossing(CNVcalls[,c(1,2,3,5)], VEP, .name_repair = "minimal")
colnames(geneContentTable) <- c('CHROM.x', 'START.x', 'END.x', 'SVTYPE.x',
                                'CHROM.y', 'START.y', 'END.y', 'SVTYPE.y', 'gene')
geneContentTable <- filter(geneContentTable,
                           geneContentTable$SVTYPE.x==geneContentTable$SVTYPE.y)
geneContentTable <- filter(geneContentTable,
                           geneContentTable$CHROM.x==geneContentTable$CHROM.y)
geneContentTable <- filter(geneContentTable,
                           geneContentTable$START.x-1<=geneContentTable$START.y)
geneContentTable <- filter(geneContentTable,
                           geneContentTable$END.x+1>=geneContentTable$END.y)

geneContentTable$SVTYPE.y = NULL
geneContentTable$CHROM.y = NULL
geneContentTable$START.y = NULL
geneContentTable$END.y = NULL
colnames(geneContentTable) <- c('CHROM', 'START', 'END', 'SVTYPE', 'gene')

# remove duplicate rows
geneContentTable <- unique(geneContentTable)

####### PLOT CNV count per chrom (total & genic) ##########

CNV_pos <- unique(geneContentTable[1:3])

VEPinput <- CNVcalls[,c(1,2,3,5)]
VEPinput <- VEPinput[which(VEPinput$SVTYPE!="mCNV"),]
VEPinput$size <- VEPinput$END-VEPinput$START

p2 <- ggplot() +
  geom_histogram(VEPinput[which(VEPinput$SVTYPE=='DEL'),],
                 mapping=aes(x = size,fill="a"),
                 color="black") +
  scale_x_log10() +
  geom_histogram(VEPinput[which(VEPinput$SVTYPE=='DUP'),],
                 mapping=aes(x = size,fill="b"),
                 color="black") +
  theme_bw() +
  labs(tag = "A")+
  ggtitle('Distribution of CNV Sizes') +
  xlab("Size (# base pairs; log scale)") +
  ylab("Number of CNVs") +
  scale_fill_manual(name = NULL,
                    values =c('a'="#4DBBD5CC",'b'="#F39B7FCC"),
                    labels = c('DEL CNVs','DUP CNVs')) +
  theme(legend.position = 'right')
p2

p1 <- ggplot() +
  geom_bar(VEPinput, mapping=aes(x = CHROM, fill="a"),color="black") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 22, by = 1)) +
  labs(tag = "B")+
  ggtitle('Distribution of CNVs Across Chromosomes') +
  xlab("Chromosome Number") +
  ylab("Number of CNVs (log scale)") +
  geom_bar(data=CNV_pos, mapping=aes(x = CHROM, fill="b"), color="black", position = 'dodge') +
  scale_y_log10() +
  scale_fill_manual(name = NULL, values = c('a'="gold",'b'="forestgreen"),
                    labels = c('All CNVs','Genic CNVs')) +
  theme(legend.position = 'right')
p1
grid.arrange(p2, p1)

rm(CNV_pos,VEPinput, p1, p2)
##############################################

# load ID conversion files
ROSmaster <- readRDS("data/ROSmaster.rds")
ROSmaster_keylist <- ROSmaster[c('projid', 'study', 'IID')]
rm(ROSmaster)

ROSmaster_keylist$study <- as.character(ROSmaster_keylist$study)
ROSmaster_keylist$projid <- as.integer(ROSmaster_keylist$projid)
ROSmaster_keylist$study <- gsub('\\s+', '', ROSmaster_keylist$study)

## sm -> MAP/AD ##
WGS_keyfile <- read.table(file = 'data/raw/WGS_keyfile.txt',
                         sep = ' ', header = TRUE)

trnsfm <- function(count) {
    # edit this function for adjusting how the counts are multiplied (e.g. take the log of counts)
    # currently it is left as linear multiplication
    return (count)
}

########## sum up scores per patient ###############
TotScoreBySample <- as.data.frame(matrix(nrow = ncol(CNVcalls), ncol = 14))

colnames(TotScoreBySample) <- c('study','projid', 'size_del', 'size_dup',
                                'gcount_del', 'gcount_dup',
                                'pli_del', 'pli_dup',
                                'loeuf_del', 'loeuf_dup',
                                'pHI', 'pTS', 'pHI_thresh', 'pTS_thresh')

for(i in 7:ncol(CNVcalls)) {       # for-loop over columns
  # i <- 8    # for debugging purpose
  lstCNV <- CNVcalls[,c(1:3, 5, i)]
  lstCNV <- lstCNV[which(lstCNV$SVTYPE!="mCNV"),]

  smid <- colnames(lstCNV)[5]
  colnames(lstCNV)[5] <- 'count'

  # convert smid -> iid -> projid #
  study <- projid <- NULL
  if (substring(smid, 1, 3) == 'ROS' || substring(smid, 1, 3) == "MAP") {
    study <- substring(smid, 1, 3)
    projid <- as.numeric(substr(smid, 4, nchar(smid)))
  } else {
    smid <- gsub("\\.", "-", smid)
    iid <- WGS_keyfile[which(WGS_keyfile$V2 == smid),]$V1
    if (length(iid) > 0) {
      if (substr(iid, 3, 4)!='AD') {
        study <- substr(iid, 1, 3)
        projid <- as.numeric(substr(iid, 4, nchar(iid)))
      } else {
        study <- ROSmaster_keylist[(which(ROSmaster_keylist$IID == iid)),][2]
        projid <- as.numeric(ROSmaster_keylist[
          (which(ROSmaster_keylist$IID == iid)),][1])
      }
    }
  }

  # change in copy number
  # 2 = diploid normal
  lstCNV[which(lstCNV$SVTYPE == 'DUP'),5] <- lstCNV[which(lstCNV$SVTYPE == 'DUP'),5] - 2
  lstCNV[which(lstCNV$SVTYPE == 'DEL'),5] <- 2 - lstCNV[which(lstCNV$SVTYPE == 'DEL'),5]

  # size
  CNVsizes <- lstCNV$END - lstCNV$START
  lstCNV <- cbind(lstCNV, CNVsizes)

  lstCNV$count <- trnsfm(lstCNV$count)

  # METHOD 1: include count (i.e. score of single CNV * count of CNV copies)
  ######### CNV size ####################
  sum_size_del <- lstCNV[which(lstCNV$SVTYPE=='DEL'), c(5, 6)]
  sum_size_del <- sum(sum_size_del$count * sum_size_del$CNVsizes)

  sum_size_dup <- lstCNV[which(lstCNV$SVTYPE=='DUP'), c(5, 6)]
  sum_size_dup <- sum(sum_size_dup$count * sum_size_dup$CNVsizes)

  ######### gene count ######################
  genes <- merge(lstCNV, geneContentTable)

  sum_gcount_del <- sum(genes[which(genes$SVTYPE=='DEL'), 5])
  sum_gcount_dup <- sum(genes[which(genes$SVTYPE=='DUP'), 5])

  # pli/loeuf ###################
  pli <- unique(merge(genes, pli_loeuf, by='gene'))
  pli <- cbind(pli, pli$count * pli$pLI)
  sum_pli_del <- sum(pli[which(pli$SVTYPE=='DEL'), 10], na.rm = TRUE)
  sum_pli_dup <- sum(pli[which(pli$SVTYPE=='DUP'), 10], na.rm = TRUE)

  # LOEUF ############
  pli <- cbind(pli, pli$count * pli$oe_lof_upper)
  sum_loeuf_del <- sum(pli[which(pli$SVTYPE=='DEL'), 11], na.rm = TRUE)
  sum_loeuf_dup <- sum(pli[which(pli$SVTYPE=='DUP'), 11], na.rm = TRUE)

  ### pHI/pTS
  phi <- unique(merge(genes, pHaploLst, by = 'gene'))
  phi <- phi[which(phi$SVTYPE=='DEL'),]

  pts <- unique(merge(genes, pTriploLst, by = 'gene'))
  pts <- pts[which(pts$SVTYPE=='DUP'),]

  sum_phi <- sum(phi$count * phi$pHaplo)
  sum_pts <- sum(pts$count * pts$pTriplo)

  ### pHI/pTS with threshold
  phi_thresh <- unique(merge(genes, pHaplo_thresh, by = 'gene'))
  pts_thresh <- unique(merge(genes, pTriplo_thresh, by = 'gene'))

  phi_thresh <- phi_thresh[which(phi_thresh$SVTYPE=='DEL'),]
  sum_phi_thresh <- sum(phi_thresh$count * phi_thresh$pHaplo)
  pts_thresh <- pts_thresh[which(pts_thresh$SVTYPE=='DUP'),]
  sum_pts_thresh <- sum(pts_thresh$count * pts_thresh$pTriplo)

  # METHOD 2: do not consider count, duplication/deletion exists or doesn't exist

#   ######### CNV size ####################
#   sum_size_del <- lstCNV[which(lstCNV$SVTYPE=='DEL'), c(5, 6)]
#   sum_size_del <- sum_size_del[(which(sum_size_del$count != 0)), ]
#   sum_size_del <- sum(sum_size_del$CNVsizes)
#
#   sum_size_dup <- lstCNV[which(lstCNV$SVTYPE=='DUP'), c(5, 6)]
#   sum_size_dup <- sum_size_dup[(which(sum_size_dup$count != 0)), ]
#   sum_size_dup <- sum(sum_size_dup$CNVsizes)
#
#   ######### gene count ######################
#
#   genes <- merge(lstCNV, geneContentTable)
#
#   sum_gcount_del <- nrow(genes[which(genes$SVTYPE=='DEL' & genes$count!=0),])
#   sum_gcount_dup <- nrow(genes[which(genes$SVTYPE=='DUP' & genes$count!=0),])
#
#   ######### pli/loeuf ###################
#   pli <- unique(merge(genes, pli_loeuf, by='gene'))
#
#   pli <- cbind(pli, pli$count * pli$pLI)
#   sum_pli_del <- sum(pli[which(pli$SVTYPE=='DEL'), 10], na.rm = TRUE)
#   sum_pli_dup <- sum(pli[which(pli$SVTYPE=='DUP'), 10], na.rm = TRUE)
#
#   ######### LOEUF ############
#   pli <- cbind(pli, pli$count * pli$oe_lof_upper)
#
#   pli_del <- pli[which(pli$count!=0 & pli$SVTYPE== 'DEL'),]
#   sum_pli_del <- sum(pli_del$pLI,na.rm=TRUE)
#   sum_loeuf_del <- sum(pli_del$oe_lof_upper, na.rm = TRUE)
#
#   pli_dup <- pli[which(pli$count!=0 & pli$SVTYPE== 'DUP'),]
#   sum_pli_dup <- sum(pli_dup$pLI, na.rm = TRUE)
#   sum_loeuf_dup <- sum(pli_dup$oe_lof_upper, na.rm = TRUE)
#
#   ####### pHI/pTS ########
#   phi <- unique(merge(genes, pHaploLst, by = 'gene'))
#   phi <- phi[which(phi$SVTYPE=='DEL'),]
#
#   pts <- unique(merge(genes, pTriploLst, by = 'gene'))
#   pts <- pts[which(pts$SVTYPE=='DUP'),]
#
#   sum_phi <- phi[which(phi$count!=0),]
#   sum_phi <- sum(sum_phi$pHaplo)
#   sum_pts <- pts[which(pts$count!=0),]
#   sum_pts <- sum(sum_pts$pTriplo)
#
#   ### pHI/pTS with threshold
#   phi_thresh <- unique(merge(genes, pHaplo_thresh, by = 'gene'))
#   pts_thresh <- unique(merge(genes, pTriplo_thresh, by = 'gene'))
#
#   phi_thresh <- phi_thresh[which(phi_thresh$SVTYPE=='DEL'),]
#   pts_thresh <- pts_thresh[which(pts_thresh$SVTYPE=='DUP'),]
#
#   sum_phi_thresh <- phi_thresh[which(phi_thresh$count!=0),]
#   sum_phi_thresh <- sum(sum_phi_thresh$pHaplo)
#   sum_pts_thresh <- pts_thresh[which(pts_thresh$count!=0),]
#   sum_pts_thresh <- sum(sum_pts_thresh$pTriplo)

  ########
  if (!is.null(study)) {
    TotScoreBySample[i,] <- c(study, projid, sum_size_del, sum_size_dup,
                              sum_gcount_del, sum_gcount_dup, sum_pli_del, sum_pli_dup,
                              sum_loeuf_del, sum_loeuf_dup, sum_phi, sum_pts, sum_phi_thresh, sum_pts_thresh)
  }
  rm(study, projid, genes)
}

TotScoreBySample <- TotScoreBySample[-c(1:6),]
TotScoreBySample[, 3:14] <- sapply(TotScoreBySample[, 3:14], as.double)

# size_del <- ggplot(TotScoreBySample, aes(x=size_del)) +
#   geom_histogram(color="black", fill="lightblue") + theme_bw()
# size_del
#
# size_dup <- ggplot(TotScoreBySample, aes(x=size_dup)) +
#   geom_histogram(color="black", fill="lightblue") + theme_bw()
# size_dup
#
# gcount_del <- ggplot(TotScoreBySample, aes(x=gcount_del)) +
#   geom_histogram(color="black", fill="lightblue") + theme_bw()
# gcount_del
#
# gcount_dup <- ggplot(TotScoreBySample, aes(x=gcount_dup)) +
#   geom_histogram(color="black", fill="lightblue") + theme_bw()
# gcount_dup

pli_del <- ggplot(TotScoreBySample, aes(x=pli_del)) +
  geom_histogram(color="black", fill="skyblue") + theme_bw() +
  xlab("pLI Score (DEL)") +
  ylab("# CNVs") + scale_x_log10()

pli_dup <- ggplot(TotScoreBySample, aes(x=pli_dup)) +
  geom_histogram(color="black", fill="salmon") + theme_bw()+
  ylab(" ")+
  xlab("pLI Score (DUP)") + scale_x_log10()

loeuf_del <- ggplot(TotScoreBySample, aes(x=loeuf_del)) +
  geom_histogram(color="black", fill="skyblue") + theme_bw()+
  xlab("LOEUF Score (DEL)") +
  ylab("# CNVs")+scale_x_log10()

loeuf_dup <- ggplot(TotScoreBySample, aes(x=loeuf_dup)) +
  geom_histogram(color="black", fill="salmon") + theme_bw()+
  xlab("LOEUF Score (DUP)") +
  ylab(" ")+
  scale_x_log10()

pHI <- ggplot(TotScoreBySample, aes(x=pHI)) +
  geom_histogram(color="black", fill="skyblue") + theme_bw()+
  xlab("pHI Score") +
  ylab("# CNVs")+ scale_x_log10()

pTS <- ggplot(TotScoreBySample, aes(x=pTS)) +
  geom_histogram(color="black", fill="salmon") + theme_bw()+
  ylab(" ")+
  xlab("pTS Score") +scale_x_log10()

pHI_thresh <- ggplot(TotScoreBySample, aes(x=pHI_thresh)) +
  geom_histogram(color="black", fill="skyblue") + theme_bw()+
  xlab("pHI Score (binarized)") +
  ylab("# CNVs")+scale_x_log10()

pTS_thresh <- ggplot(TotScoreBySample, aes(x=pTS_thresh)) +
  geom_histogram(color="black", fill="salmon") + theme_bw()+
  ylab(" ")+
  xlab("pTS Score (binarized)") +scale_x_log10()

grid.arrange(pli_del, pli_dup, loeuf_del, loeuf_dup, pHI, pTS, pHI_thresh, pTS_thresh, ncol=2)

saveRDS(TotScoreBySample[complete.cases(TotScoreBySample), ], "data/TotScoreBySample.rds")


# # a small overlap test
# a<-CNVcalls[which(CNVcalls$SVTYPE!="mCNV"),]
#
# for (i in 1:nrow(a)) {
#   chr <- a[i,1]
#   start <-a[i,2]
#   end <-a[i,3]
#   type <- a[i,5]
#   smaller <- a[which(a$CHROM==chr &
#                        a$START>=start &
#                        a$END<=end &
#                        a$SVTYPE==type &
#                        !(a$START==start &
#                            a$END==end)),]
#   if (nrow(smaller)>0) {
#     a[i,] <- NA
#     print(i)
#   }
#
# }
