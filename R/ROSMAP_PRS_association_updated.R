# Note: instead of checking the ANOVA between pheno ~ RS and pheno ~ RS + ADPRS,
# here it checks ANOVA between pheno ~ ADPRS and pheno ~ ADPRS + RS

##########
# load packages
library(dplyr)
library(AICcmodavg)
library(ggplot2)
library(pscl)
library(lmtest)
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
str_scores <- c('+ pli_del + pli_dup ',
                '+ loeuf_del + loeuf_dup ',
                '+ pHI + pTS ',
                '+ pHI_thresh + pTS_thresh ')

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
  p_t_val <- matrix(nrow = length(lst_pheno) * length(lst_risk_scores), ncol = 9)
  colnames(p_t_val) <- c('phenotype', 'risk_score',
                         'PRS_t', 'PRS_p', 'adj.r', 'adj.r_rs',
                         'anova_p', 'AIC', 'AIC_w_RS')
  # i=1
  # j=1
  for (i in 1:length(lst_pheno)) {
    for (j in 1:length(str_scores)) {
      if (i <= length(lst_pheno_autop)) {
        formula <- paste0(lst_pheno[i], '~ PRS', str_covar_autop, str_egvec)
        formula_rs <- paste0(lst_pheno[i], '~ PRS', str_scores[j], str_covar_autop, str_egvec)
      } else {
        formula <- paste0(lst_pheno[i], '~ PRS', str_covar_cogn, str_egvec)
        formula_rs <- paste0(lst_pheno[i], '~ PRS', str_scores[j], str_covar_cogn, str_egvec)
      }
      fit <- lm(formula, data=scoreWithPheno)
      fit_rs <- lm(formula_rs, data=scoreWithPheno)

      anova <- anova(fit, fit_rs)
      #define list of models
      models <- list(fit, fit_rs)

      #specify model names
      mod.names <- c('fit','fit_rs')

      #calculate AIC of each model
      aic <- aictab(cand.set = models, modnames = mod.names)
      # aic$AICc <- scitf_note(aic$AICc, 2)

      fit_coeff <- scitf_note(summary(fit)$coefficients, 2)
      fit_rs_coeff <- scitf_note(summary(fit_rs)$coefficients, 2)

      a <-(i-1)*8+2*j-1
      p_t_val[a,] <- c(lst_pheno[i],
                       lst_risk_scores[2*j-1],
                       fit_coeff[,3][2],
                       fit_coeff[,4][2],
                       scitf_note(summary(fit)$adj.r.squared, 2),
                       scitf_note(summary(fit_rs)$adj.r.squared, 2),
                       scitf_note(anova$`Pr(>F)`[2], 2),
                       aic$AICc[which(aic$Modnames=='fit')],
                       aic$AICc[which(aic$Modnames=='fit_rs')]
      )

      p_t_val[a+1,] <- c(lst_pheno[i],
                         lst_risk_scores[2*j],
                         NA,
                         NA,
                         NA,
                         NA,
                         NA,
                         NA,
                         NA
      )
    }
  }

  for (j in 1:length(str_scores)) {
    formula <- paste0('cogdx_binom', '~ PRS', str_covar_cogn, str_egvec)
    formula_rs <- paste0('cogdx_binom', '~ PRS', str_scores[j], str_covar_cogn, str_egvec)

    fit <- glm(formula, data=scoreWithPheno, family='binomial')
    fit_rs <- glm(formula_rs, data=scoreWithPheno, family='binomial')

    #define list of models
    models <- list(fit, fit_rs)

    #specify model names
    mod.names <- c('fit','fit_rs')

    #calculate AIC of each model
    # aic <- aictab(cand.set = models, modnames = mod.names)
    # aic$AICc <- scitf_note(aic$AICc, 2)
    aic_fit <- AIC(fit)
    aic_fit_rs <- AIC(fit_rs)
    aic_data <- data.frame(Model = c('fit', 'fit_rs'), AIC = c(aic_fit, aic_fit_rs))

    fit_coeff <- scitf_note(summary(fit)$coefficients, 2)
    fit_rs_coeff <- scitf_note(summary(fit_rs)$coefficients, 2)

    pseudo_rsq_model1 <- pR2(fit, type = "mcfadden")
    pseudo_rsq_model2 <- pR2(fit_rs, type = "mcfadden")
    model_comparison <- lrtest(fit, fit_rs)

    new_rows <- matrix(nrow=2, ncol=9)
    colnames(new_rows) <- c('phenotype', 'risk_score',
                            'PRS_t', 'PRS_p', 'adj.r', 'adj.r_rs',
                            'anova_p', 'AIC', 'AIC_w_RS')

    new_rows[1, ] <- c('cogdx',
                       lst_risk_scores[2*j-1],
                       fit_coeff[,3][2],
                       fit_coeff[,4][2],
                       scitf_note(pseudo_rsq_model1['r2ML'], 2),
                       scitf_note(pseudo_rsq_model2['r2ML'], 2),
                       scitf_note(model_comparison$`Pr(>Chisq)`[2], 20),
                       aic_data$AIC[1],
                       aic_data$AIC[2]
    )

    new_rows[2, ] <- c('cogdx',
                       lst_risk_scores[2*j],
                       NA,
                       NA,
                       NA,
                       NA,
                       NA,
                       NA,
                       NA
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
p_t_val$PRS_p_fdr <- p.adjust(p_t_val$PRS_p)
write.csv(p_t_val, 'output/CNVRS_ANOVA_updated.csv', row.names = FALSE)
