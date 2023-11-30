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
# all the code above are used to generate the table with Scores & Phenotypes.
scoreWithPheno <- readRDS('data/scoreWithPheno.rds')
scoreWithPheno_nocopynum <- readRDS('data/scoreWithPheno_nocopynumber.rds')

# list of 10 eigenvectors
str_egvec <- '+ egvec1 + egvec2 + egvec3 + egvec4 + egvec5 + egvec6 + egvec7 + egvec8 + egvec9 + egvec10'
str_egvec_stdz <- '+ scale(egvec1) + scale(egvec2) + scale(egvec3) + scale(egvec4) + scale(egvec5) + scale(egvec6) + scale(egvec7) + scale(egvec8) + scale(egvec9) + scale(egvec10)'

# other covariates (neurological)
str_covar_autop <- '+ msex + age_death + pmi'
str_covar_autop_stdz <- '+ scale(msex) + scale(age_death) + scale(pmi)'

# other covariates (cognitive)
str_covar_cogn <- '+ msex + educ + age_at_visit'
str_covar_cogn_stdz <- '+ scale(msex) + scale(educ) + scale(age_at_visit)'

# 4 models, 2 scores per model
str_scores <- c('~ pli_del + pli_dup ',
                '~ loeuf_del + loeuf_dup ',
                '~ pHI + pTS ',
                '~ pHI_thresh + pTS_thresh ')
str_scores_stdz <- c('~ scale(pli_del) + scale(pli_dup)',
                     '~ scale(loeuf_del) + scale(loeuf_dup)',
                     '~ scale(pHI) + scale(pTS)',
                     '~ scale(pHI_thresh) + scale(pTS_thresh)')

# continuous variables
lst_pheno_autop <- c('sqrt(tangles)', 'sqrt(amyloid)', 'arteriol_scler', 'cvda_4gp2')
lst_pheno_cog <- c('cogn_global')

lst_risk_scores <- c('pli_del','pli_dup','loeuf_del','loeuf_dup',
                     'pHI','pTS', 'pHI_thresh','pTS_thresh')

#######################
p_t_val_calculation <- function(scoreWithPheno, lst_risk_scores, str_scores,
                                lst_pheno_autop, lst_pheno_cog) {
  lst_pheno <- c(lst_pheno_autop, lst_pheno_cog)
  p_t_val <- matrix(nrow = length(lst_pheno) * length(lst_risk_scores), ncol = 7)
  colnames(p_t_val) <- c('phenotype', 'risk_score', 'standardized_beta',
                         'p_val', 'ci_lw', 'ci_up', 't_val')

  for (i in 1:length(lst_pheno)) {
    for (j in 1:length(str_scores)) {
      if (i <= length(lst_pheno_autop)) {
        formula <- paste0(lst_pheno[i], str_scores[j], str_covar_autop, str_egvec)
        stdz_formula <- paste0('scale(', lst_pheno[i], ')', str_scores_stdz[j],
                               str_covar_autop_stdz, str_egvec_stdz)
      } else {
        formula <- paste0(lst_pheno[i], str_scores[j],
                          str_covar_cogn, str_egvec)
        stdz_formula <- paste0('scale(', lst_pheno[i], ')', str_scores_stdz[j],
                               str_covar_cogn_stdz, str_egvec_stdz)
      }

      fit <- lm(formula, data=scoreWithPheno)
      fit_stdz <- lm(stdz_formula, data=scoreWithPheno)

      confint1 <- confint(fit_stdz, paste0('scale(',lst_risk_scores[2*j-1],')'),level=0.95)
      confint2 <- confint(fit_stdz, paste0('scale(',lst_risk_scores[2*j],')'), level=0.95)
      standardized_betas <- lm.beta(fit_stdz)

      curr_row <-(i-1)*8+2*j-1
      p_t_val[curr_row, ] <- c(lst_pheno[i], lst_risk_scores[2*j-1],
                               unname(standardized_betas$coefficients[2]),
                               summary(fit)$coefficients[,4][2],
                               confint1[1], confint1[2],
                               summary(fit)$coefficients[,3][2])
      p_t_val[curr_row+1, ] <- c(lst_pheno[i], lst_risk_scores[2*j],
                                 unname(standardized_betas$coefficients[3]),
                                 summary(fit)$coefficients[,4][3],
                                 confint2[1], confint2[2],
                                 summary(fit)$coefficients[,3][3])
    }
  }

  ### categorical
  # cogdx
  for (j in 1:length(str_scores)) {
    formula <- paste0('cogdx_binom', str_scores[j], str_covar_cogn, str_egvec)
    stdz_formula <- paste0('scale(cogdx_binom)', str_scores_stdz[j],
                           str_covar_cogn_stdz, str_egvec_stdz)

    fit <- glm(formula, data=scoreWithPheno, family='binomial')
    fit_stdz <- lm(stdz_formula, data=scoreWithPheno)

    confint1 <- confint(fit_stdz, paste0('scale(',lst_risk_scores[2*j-1],')'),level=0.95)
    confint2 <- confint(fit_stdz, paste0('scale(',lst_risk_scores[2*j],')'), level=0.95)
    standardized_betas <- lm.beta(fit_stdz)

    new_rows <- matrix(nrow=2, ncol=7)
    new_rows[1, ] <- c('cogdx',
                       lst_risk_scores[2*j-1],
                       unname(standardized_betas$coefficients[2]),
                       summary(fit)$coefficients[,4][2],
                       confint1[1],
                       confint1[2],
                       summary(fit)$coefficients[,3][2])
    new_rows[2, ] <- c('cogdx',
                       lst_risk_scores[2*j],
                       unname(standardized_betas$coefficients[3]),
                       summary(fit)$coefficients[,4][3],
                       confint2[1],
                       confint2[2],
                       summary(fit)$coefficients[,3][3])

    p_t_val <- rbind(p_t_val, new_rows)
  }

  p_t_val <- as.data.frame(p_t_val)
  p_t_val$standardized_beta <- as.double(p_t_val$standardized_beta)
  p_t_val$p_val <- as.double(p_t_val$p_val)
  p_t_val$t_val <- as.double(p_t_val$t_val)
  p_t_val$ci_lw <- as.double(p_t_val$ci_lw)
  p_t_val$ci_up <- as.double(p_t_val$ci_up)

  # multiple testing correction
  p_t_val$p_val_adj <- p.adjust(p_t_val$p_val, 'fdr')
  return(p_t_val)
}

p_t_val <- p_t_val_calculation(scoreWithPheno, lst_risk_scores, str_scores,
                                   lst_pheno_autop, lst_pheno_cog)
p_t_val_ncn <- p_t_val_calculation(scoreWithPheno_nocopynum, lst_risk_scores, str_scores,
                               lst_pheno_autop, lst_pheno_cog)

# renaming things to make the plots look prettier
p_t_val[p_t_val == "pTS_thresh"] <- "pTS (binarized)"
p_t_val[p_t_val == "pHI_thresh"] <- "pHI (binarized)"
p_t_val[p_t_val == "pli_del"] <- "pLI (DEL)"
p_t_val[p_t_val == "pli_dup"] <- "pLI (DUP)"
p_t_val[p_t_val == "loeuf_del"] <- "LOEUF (DEL)"
p_t_val[p_t_val == "loeuf_dup"] <- "LOEUF (DUP)"

p_t_val[p_t_val == "sqrt(tangles)"] <- 'Tangle Density\n(n=1001)'
p_t_val[p_t_val == "sqrt(amyloid)"] <- 'Overall Amyloid Level\n(n=1005)'
p_t_val[p_t_val == "cvda_4gp2"] <- "Cerebral Atherosclerosis\nRating\n(n=1004)"
p_t_val[p_t_val == "cogn_global"] <- "Global Cognitive Function\n(n=1011)"
p_t_val[p_t_val == "cogdx"] <- 'Final Consensus\nCognitive Diagnosis\n(n=946)'
p_t_val[p_t_val == "arteriol_scler"] <- 'Arteriolosclerosis\n(n=1003)'


# Renaming for p_t_val_ncn
p_t_val_ncn[p_t_val_ncn == "pTS_thresh"] <- "pTS (binarized)"
p_t_val_ncn[p_t_val_ncn == "pHI_thresh"] <- "pHI (binarized)"
p_t_val_ncn[p_t_val_ncn == "pli_del"] <- "pLI (DEL)"
p_t_val_ncn[p_t_val_ncn == "pli_dup"] <- "pLI (DUP)"
p_t_val_ncn[p_t_val_ncn == "loeuf_del"] <- "LOEUF (DEL)"
p_t_val_ncn[p_t_val_ncn == "loeuf_dup"] <- "LOEUF (DUP)"

p_t_val_ncn[p_t_val_ncn == "sqrt(tangles)"] <- 'Tangle Density\n(n=1001)'
p_t_val_ncn[p_t_val_ncn == "sqrt(amyloid)"] <- 'Overall Amyloid Level\n(n=1005)'
p_t_val_ncn[p_t_val_ncn == "cvda_4gp2"] <- "Cerebral Atherosclerosis\nRating\n(n=1004)"
p_t_val_ncn[p_t_val_ncn == "cogn_global"] <- "Global Cognitive Function\n(n=1011)"
p_t_val_ncn[p_t_val_ncn == "cogdx"] <- 'Final Consensus\nCognitive Diagnosis\n(n=946)'
p_t_val_ncn[p_t_val_ncn == "arteriol_scler"] <- 'Arteriolosclerosis\n(n=1003)'
# heatmap
rs_order <- c('pLI (DEL)', 'pLI (DUP)',
              'LOEUF (DEL)', 'LOEUF (DUP)',
              'pHI', 'pTS',
              'pHI (binarized)', 'pTS (binarized)')



#  p< 0.05 (*)
# 0.01 (**)
# 0.001 (***)

# ... (Your existing code)

# Create heatmap for scoreWithPheno
plot_scoreWithPheno <- ggplot(data = p_t_val, aes(x=factor(risk_score, level = rs_order),
                                                  y=phenotype, fill=t_val)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-6,6), space = "Lab", name="t-statistic") +
  theme_minimal() +
  xlab('Risk Score') +
  ylab('Phenotype') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val < 0.05 & p_val_adj >= 0.05 , as.character(sprintf("%.3f", t_val)), ' ')),
            color = "grey30", size = 3.5) + theme(axis.text.x = element_text(size = 8)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.05 & p_val_adj >= 0.01, paste0(as.character(sprintf("%.2f", t_val)), "*"), ' '), fontface = 'bold'),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.01 & p_val_adj >= 0.001, paste0(as.character(sprintf("%.2f", t_val)), "**"), ' '), fontface = 'bold'),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.001, paste0(as.character(sprintf("%.3f", t_val)), "***"), ' '), fontface = 'bold'),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8))

# Create heatmap for scoreWithPheno_nocopynum
plot_scoreWithPheno_nocopynum <- ggplot(data = p_t_val_ncn, aes(x=factor(risk_score, level = rs_order),
                                                            y=phenotype, fill=t_val)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-6,6), space = "Lab", name="t-statistic") +
  theme_minimal() +
  xlab('Risk Score') +
  ylab('Phenotype') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val < 0.05 & p_val_adj >= 0.05 , as.character(sprintf("%.3f", t_val)), ' ')),
            color = "grey30", size = 3.5) + theme(axis.text.x = element_text(size = 8)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.05 & p_val_adj >= 0.01, paste0(as.character(sprintf("%.2f", t_val)), "*"), ' '), fontface = 'bold'),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.01 & p_val_adj >= 0.001, paste0(as.character(sprintf("%.2f", t_val)), "**"), ' '), fontface = 'bold'),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.001, paste0(as.character(sprintf("%.3f", t_val)), "***"), ' '), fontface = 'bold'),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8))

# Arrange the plots side by side
library(ggpubr)
ggarrange(plot_scoreWithPheno, plot_scoreWithPheno_nocopynum, ncol = 2, labels = c("* weights", "no weights"))


# forest plot
p_t_val$phenotype = factor (p_t_val$phenotype)
p_t_val <- p_t_val %>% mutate(Annotation = case_when(
  p_val < 0.05 & p_val_adj >= 0.05 ~ '*',
  p_val < 0.05 & p_val_adj < 0.05 ~ '**',
  TRUE ~ ''))
forest_plot <- ggplot(data=p_t_val, aes(y=phenotype, x=standardized_beta, xmin=ci_lw, xmax=ci_up, col=risk_score,fill=risk_score)) +
  geom_linerange(size=0.7,position=position_dodge(width = 0.5), alpha=ifelse(p_t_val$p_val >= 0.05, 0.5, 1)) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(size=3, shape=20, stroke = 0.5,position=position_dodge(width = 0.5), alpha=ifelse(p_t_val$p_val >= 0.05, 0.3, 1)) +
  scale_x_continuous(limits = c(-0.15, 0.25), breaks = c(-0.2,-0.1,0, 0.1, 0.2), name = 'Standardized Effect (Beta Coefficient)') +
  scale_y_discrete(name="Outcomes") +
  theme_minimal() +
  scale_color_discrete(name = "CNV-RS") +  # Rename the legend title
  scale_fill_discrete(name = "CNV-RS") +   # Rename the legend title
  theme(legend.position = "bottom")

forest_plot <- forest_plot + geom_text(aes(x = ci_up, label = Annotation), hjust = -1, vjust=0.7, position = position_dodge(width = 0.5), size = 5, color='black')

# grid.arrange(forest_plot, plot_bothsex)
# g <- arrangeGrob(forest_plot, plot_bothsex, heights = 3:2)
ggsave("output/main_fig.png", forest_plot, width = 20, height = 25, units = "cm")
ggsave(file="output/main_fig.svg", plot=forest_plot, width = 20, height = 25, units = "cm")

a <- p_t_val
a1 <- plot_bothsex
a2 <- forest_plot
