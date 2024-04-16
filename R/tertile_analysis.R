# load packages
library(stats)
library(dplyr)
library(scales)
library(lm.beta)

scoreWithPheno <- readRDS('data/scoreWithPheno.rds')

lst_risk_scores <- c('pli_del', 'pli_dup', 'loeuf_del', 'loeuf_dup',
                     'pHI', 'pTS', 'pHI_thresh', 'pTS_thresh')

output_df <- data.frame(matrix(ncol = 5, nrow = 8))

colnames(output_df) <- c("CNV_RS",
                         "Fisher_p",
                         "Odds_ratio",
                         "Lower_CI",
                         "Upper_CI")

for (i in 1:length(lst_risk_scores)) {
  score <- lst_risk_scores[i]
  print(score)
  quantiles <- quantile(scoreWithPheno[[score]], probs = c(0.33, 0.66))
  scoreWithPheno$CNV_RS_group <- cut(scoreWithPheno[[score]],
                                     breaks = c(-Inf, quantiles, Inf),
                                     labels = c("Low", "Medium", "High"),
                                     include.lowest = TRUE)
  filtered <- scoreWithPheno %>%
    filter(cvda_4gp2 == 0 | cvda_4gp2 == 3)

  contingency_table <- table(filtered$CNV_RS_group,
                             filtered$cvda_4gp2)

  contingency_table <- contingency_table[rownames(contingency_table) != "Medium", ]

  print(contingency_table)

  fisher_test_result <- fisher.test(contingency_table, conf.int = TRUE)
  odds_ratio <- fisher_test_result$estimate
  conf_int <- fisher_test_result$conf.int

  output_df[i, ] <- c(score,
                      fisher_test_result$p.value,
                      odds_ratio,
                      conf_int[1],
                      conf_int[2])
}

print(output_df)
