
fulltable <- read.csv('CNVRS_PRS_association_kunkle.csv', header = TRUE)

colnames(fulltable) <- gsub("\\.+", "_", colnames(fulltable))

fulltable_subset <- fulltable[, c("phenotype",
                                  "adj_R_squared_PRS_LOAD_only_", "adj_R_squared_PRS_LOAD_RS_",
                                  "AIC_PRS_LOAD_only_", "AIC_PRS_LOAD_RS_",
                                  "ANOVA_p_val")]
colnames(fulltable_subset) <- c("Phenotype",
                                "PRS-LOAD only", "PRS-LOAD + RS1 + RS2",
                                "AIC_PRS_LOAD_only_", "AIC_PRS_LOAD_RS_",
                                "ANOVA_p_val")
df <- aggregate(. ~ Phenotype, data = fulltable_subset, FUN = mean, na.action = na.omit)
library(tidyr)
library(ggplot2)
library(reshape2)
scitf_note <- function(num, digit) {
  # converts to scientific notation
  return (formatC(num, format = "e", digits = digit))
}

# Melt the data frame
dfm <- melt(df[, c('Phenotype', 'PRS-LOAD only', 'PRS-LOAD + RS1 + RS2')], id.vars = 'Phenotype')

dfm[dfm == "sqrt(tangles)"] <- 'Tangle Density'
dfm[dfm == "sqrt(amyloid)"] <- 'Overall Amyloid Level'
dfm[dfm == "cvda_4gp2"] <- "Cerebral Atherosclerosis\nRating"
dfm[dfm == "cogn_global"] <- "Global Cognitive Function"
dfm[dfm == "cogdx"] <- 'Final Consensus\nCognitive Diagnosis'
dfm[dfm == "arteriol_scler"] <- 'Arteriolosclerosis'


labels <- df$ANOVA_p_val
labels <- c(rep(NA, 5), labels)
# Create the grouped bar plot
g<-ggplot(dfm, aes(x = Phenotype, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
  geom_text(aes(label = if_else(labels < 0.05, paste0('p = ', scitf_note(labels, 2), '***'), ' '), group = Phenotype, fontface = 'bold'),
            position = position_dodge(width = 0.9), vjust = -1, size = 3.5, color = 'black') +
  geom_text(aes(label = if_else(labels >= 0.05, paste0('p = ', as.character(sprintf("%.3f", labels))), ' '), group = Phenotype),
            position = position_dodge(width = 0.9), vjust = -1, size = 3.5, color = "grey30") +

  labs(fill = "Category") +
  scale_fill_manual(values = c('PRS-LOAD only' = 'darkgoldenrod1', 'PRS-LOAD + RS1 + RS2' = 'skyblue')) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.85),
        legend.background = element_rect(linetype="solid", colour ="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1)) +
  guides(fill = guide_legend(title = 'Models')) +
  labs(x = "Outcomes") +
  scale_y_continuous(limits = c(0, 0.4), name = "Adjusted R-Squared (Average of 4 Models)")



g




##########

cogdx <- read.csv('./output/CNVRS_ANOVA_updated.csv', header = TRUE)
cogdx <- cogdx[cogdx$phenotype=='cogdx',]

cogdx_avg <- data.frame(
  Phenotype = c('Final Consensus\nCogitive Diagnosis', 'Final Consensus\nCogitive Diagnosis'),
  variable = c('PRS-LOAD only', 'PRS-LOAD + RS1 + RS2'),
  value = c(mean(cogdx$adj.r, na.rm=TRUE), mean(cogdx$adj.r_rs, na.rm=TRUE))
)


labels2 <- mean(cogdx$anova_p, na.rm=TRUE)
labels2 <- c(rep(NA, 1), labels2)
# Create the grouped bar plot
# Reorder the levels of the "Category" variable
cogdx_avg$variable <- factor(cogdx_avg$variable, levels = c('PRS-LOAD only', 'PRS-LOAD + RS1 + RS2'))

# Create the grouped bar plot
g2 <- ggplot(cogdx_avg, aes(x = Phenotype, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
  geom_text(aes(label = if_else(labels2 < 0.05, paste0('p = ', scitf_note(labels2, 2), '***'), ' '), group = Phenotype, fontface = 'bold'),
            position = position_dodge(width = 0.9), vjust = -1, size = 3.5, color = 'black') +
  geom_text(aes(label = if_else(labels2 >= 0.05, paste0('p = ', as.character(sprintf("%.3f", labels2))), ' '), group = Phenotype),
            position = position_dodge(width = 0.9), vjust = -1, size = 3.5, color = "grey30") +
  labs(fill = "Category") +
  scale_fill_manual(values = c('PRS-LOAD only' = 'darkgoldenrod1', 'PRS-LOAD + RS1 + RS2' = 'skyblue')) +
  theme_bw() +
  theme(legend.position = 'none',  # Remove the legend
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1)) +
  guides(fill = FALSE) +  # Remove the legend
  labs(x = "Outcome (Binomial)") +
  scale_y_continuous(limits = c(0, 0.4), name = "Adjusted McFadden Pseudo R-Squared (Average of 4 Models)")

g2
library(patchwork)

# Set the same y-axis limits for both plots
y_limits <- c(0, 0.4)

# Combine g and g2 using patchwork
combined_plot <- g + g2 +
  plot_layout(ncol = 2, byrow = TRUE, widths = c(5, 1)) +
  plot_annotation(tag_levels = "A")

# Print the combined plot
combined_plot

ggsave("output/PRS_interaction.png", g, width = 15, height = 16, units = "cm")
