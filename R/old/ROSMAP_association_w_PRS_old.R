##########
# load packages
library(dplyr)
library(AICcmodavg)
library(ggplot2)
#####

scoreWithPheno <- readRDS('data/scoreWithPheno.rds')

# PRS <- read.csv('data/raw/polygenic_risk_score_Bellenguez.tsv', sep = '\t', header = TRUE)
# colnames(PRS) <- c('IID', 'PRS')
PRS <- read.csv('data/raw/raw_score_kunkle_full.tsv', sep = '\t', header = TRUE)

# PRS_subset <- read.csv('data/raw/raw_score_kunkle_subset.tsv', sep = '\t', header = TRUE)

##############################

# coverting AD to projid
ROSmaster <- readRDS("data/ROSmaster.rds")
ROSmaster_keylist <- ROSmaster[c('projid', 'study', 'IID')]
rm(ROSmaster)

ROSmaster_keylist$study <- as.character(ROSmaster_keylist$study)
ROSmaster_keylist$projid <- as.integer(ROSmaster_keylist$projid)
ROSmaster_keylist$study <- gsub('\\s+', '', ROSmaster_keylist$study)

AD <- PRS[which(grepl('AD', PRS$IID)),]

temp <- merge(AD, ROSmaster_keylist, by.x = "IID", by.y = "IID")
temp$IID<-NULL

ROS <- PRS[-which(grepl('AD', PRS$IID)),]
study <- substring(ROS$IID, 1, 3)
projid <- as.numeric(substring(ROS$IID, 4))
temp2 <- data.frame(ROS$PRS, projid, study)
colnames(temp2) <- c('PRS','projid','study')
PRS <- rbind(temp, temp2)

scoreWithPheno <- merge(scoreWithPheno, PRS, by.x=c('projid','study.x'), by.y=c('projid','study'))
rm(temp, temp2, AD, ROS, ROSmaster_keylist)

############

# list of 10 eigenvectors
str_egvec <- '+ egvec1 + egvec2 + egvec3 + egvec4 + egvec5 + egvec6 + egvec7 + egvec8 + egvec9 + egvec10'

# other covariates (neurological)
str_covar_autop <- '+ msex + age_death + pmi'

# other covariates (cognitive)
str_covar_cogn <- '+ msex + educ + age_at_visit'

# 4 models, 2 scores per model
str_scores <- c('~ pli_del + pli_dup ',
                '~ loeuf_del + loeuf_dup ',
                '~ pHI + pTS ',
                '~ pHI_thresh + pTS_thresh ')

# continuous variables
lst_pheno_autop <- c('sqrt(tangles)', 'sqrt(amyloid)', 'arteriol_scler', 'cvda_4gp2')
lst_pheno_cog <- c('cogn_global')

lst_risk_scores <- c('pli_del','pli_dup', 'loeuf_del','loeuf_dup',
                     'pHI','pTS', 'pHI_thresh','pTS_thresh')

scitf_note <- function(num, digit) {
  # converts to scientific notation
  return (formatC(num, format = "E", digits = digit))
}

p_t_val_calculation <- function(scoreWithPheno, lst_risk_scores, str_scores,
                                lst_pheno_autop, lst_pheno_cog) {
  lst_pheno <- c(lst_pheno_autop, lst_pheno_cog)
  p_t_val <- matrix(nrow = length(lst_pheno) * length(lst_risk_scores), ncol = 13)
  colnames(p_t_val) <- c('phenotype', 'risk_score', 't','p', 'adj.r','t_w_prs','p_w_prs','adj.r_prs','anova_p',
                         'AIC', 'AIC_w_prs', 'prs_ad_t','prs_ad_p')

  for (i in 1:length(lst_pheno)) {
    for (j in 1:length(str_scores)) {
      if (i <= length(lst_pheno_autop)) {
        formula <- paste0(lst_pheno[i], str_scores[j], str_covar_autop, str_egvec)
        formula_prs <- paste0(lst_pheno[i], str_scores[j], str_covar_autop, str_egvec, '+ PRS')
      } else {
        formula <- paste0(lst_pheno[i], str_scores[j], str_covar_cogn, str_egvec)
        formula_prs <- paste0(lst_pheno[i], str_scores[j], str_covar_cogn, str_egvec, '+ PRS')
      }
      fit <- lm(formula, data=scoreWithPheno)
      fit_prs <- lm(formula_prs, data=scoreWithPheno)

      anova <- anova(fit, fit_prs)
      #define list of models
      models <- list(fit, fit_prs)

      #specify model names
      mod.names <- c('fit','fit_w_prs')

      #calculate AIC of each model
      aic <- aictab(cand.set = models, modnames = mod.names)

      fit_coeff <- scitf_note(summary(fit)$coefficients, 2)
      fit_prs_coeff <- scitf_note(summary(fit_prs)$coefficients, 2)
      aic$AICc <- scitf_note(aic$AICc, 2)

      a <-(i-1)*8+2*j-1
      p_t_val[a,] <- c(lst_pheno[i],
                       lst_risk_scores[2*j-1],
                       fit_coeff[,3][2],
                       fit_coeff[,4][2],
                       scitf_note(summary(fit)$adj.r.squared, 2),
                       fit_prs_coeff[,3][2],
                       fit_prs_coeff[,4][2],
                       scitf_note(summary(fit_prs)$adj.r.squared, 2),
                       scitf_note(anova$`Pr(>F)`[2], 2),
                       aic$AICc[which(aic$Modnames=='fit')],
                       aic$AICc[which(aic$Modnames=='fit_w_prs')],
                       fit_prs_coeff[,3][17],
                       fit_prs_coeff[,4][17]
                       )

      p_t_val[a+1,] <- c(lst_pheno[i],
                         lst_risk_scores[2*j],
                         fit_coeff[,3][3],
                         fit_coeff[,4][3],
                         NA,
                         fit_prs_coeff[,3][3],
                         fit_prs_coeff[,4][3],
                         NA,
                         NA,
                         aic$AICc[which(aic$Modnames=='fit')],
                         aic$AICc[which(aic$Modnames=='fit_w_prs')],
                         NA,
                         NA
                         )
    }
  }
  ### categorical
  # cogdx
  for (j in 1:length(str_scores)) {
    formula <- paste0('cogdx_binom', str_scores[j], str_covar_cogn, str_egvec)
    formula_prs <- paste0('cogdx_binom', str_scores[j], str_covar_cogn, str_egvec, '+ PRS')

    fit <- glm(formula, data=scoreWithPheno, family='binomial')
    fit_prs <- glm(formula_prs, data=scoreWithPheno, family='binomial')
    anova <- anova(fit, fit_prs)

    #define list of models
    models <- list(fit, fit_prs)

    #specify model names
    mod.names <- c('fit','fit_w_prs')

    #calculate AIC of each model
    aic <- aictab(cand.set = models, modnames = mod.names)

    new_rows <- matrix(nrow=2, ncol=13)
    colnames(new_rows) <- c('phenotype', 'risk_score', 't','p', 'adj.r','t_w_prs','p_w_prs','adj.r_prs','anova_p',
                            'AIC', 'AIC_w_prs', 'prs_ad_t','prs_ad_p')

    fit_coeff <- scitf_note(summary(fit)$coefficients, 2)
    fit_prs_coeff <- scitf_note(summary(fit_prs)$coefficients, 2)
    aic$AICc <- scitf_note(aic$AICc, 2)

    new_rows[1, ] <- c('cogdx',
                       lst_risk_scores[2*j-1],
                       fit_coeff[,3][2],
                       fit_coeff[,4][2],
                       NA,
                       fit_prs_coeff[,3][2],
                       fit_prs_coeff[,4][2],
                       NA,
                       NA,
                       aic$AICc[which(aic$Modnames=='fit')],
                       aic$AICc[which(aic$Modnames=='fit_w_prs')],
                       fit_prs_coeff[,3][17],
                       fit_prs_coeff[,4][17]
                       )

    new_rows[2, ] <- c('cogdx',
                       lst_risk_scores[2*j],
                       fit_coeff[,3][3],
                       fit_coeff[,4][3],
                       NA,
                       fit_prs_coeff[,3][3],
                       fit_prs_coeff[,4][3],
                       NA,
                       NA,
                       aic$AICc[which(aic$Modnames=='fit')],
                       aic$AICc[which(aic$Modnames=='fit_w_prs')],
                       fit_prs_coeff[,3][18],
                       fit_prs_coeff[,4][18]
                       )
    p_t_val <- rbind(p_t_val, new_rows)
  }

  p_t_val <- as.data.frame(p_t_val)
  p_t_val[p_t_val == "pTS_thresh"] <- "pTS (binarized)"
  p_t_val[p_t_val == "pHI_thresh"] <- "pHI (binarized)"
  p_t_val[p_t_val == "pli_del"] <- "pLI (DEL)"
  p_t_val[p_t_val == "pli_dup"] <- "pLI (DUP)"
  p_t_val[p_t_val == "loeuf_del"] <- "LOEUF (DEL)"
  p_t_val[p_t_val == "loeuf_dup"] <- "LOEUF (DUP)"

  return(p_t_val)
}

######
p_t_val <- p_t_val_calculation(scoreWithPheno, lst_risk_scores, str_scores,
                               lst_pheno_autop, lst_pheno_cog)

## multiple testing correction
p_vals <- p_t_val[, c('p_w_prs', 'prs_ad_p', 'p')]
combined_vector <- c(p_vals$p_w_prs, as.numeric(p_vals$prs_ad_p), p_vals$p)

fdr_p <- scitf_note(p.adjust(combined_vector, method = 'fdr'), 2)
bonferroni_p <- scitf_note(p.adjust(combined_vector, method = "bonferroni"), 2)

p_t_val$fdr_RS_p <- data.frame(fdr_p)[97:144,]
p_t_val$fdr_RS_wprs_p <- data.frame(fdr_p)[1:48,]
p_t_val$fdr_prs_ad_p <- data.frame(fdr_p)[49:96,]
p_t_val$bonferroni_RS_p <- data.frame(bonferroni_p)[97:144,]
p_t_val$bonferroni_RS_wprs_p <- data.frame(bonferroni_p)[1:48,]
p_t_val$bonferroni_prs_ad_p <- data.frame(bonferroni_p)[49:96,]

write.csv(p_t_val, 'output/CNVRS_PRS_association_kunkle.csv', row.names = FALSE)
