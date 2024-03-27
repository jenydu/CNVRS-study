library(dplyr)

pLI_loeuf_dat <- read.csv("data/raw/pLI_loeuf_scores_raw.tsv", sep = '\t',
                          header = TRUE)

pLI_loeuf_dat <- pLI_loeuf_dat %>%
  group_by(gene) %>%
  arrange(p, .by_group = TRUE) %>%  # Arrange by p-value in ascending order
  slice(1)
length(unique(pLI_loeuf_dat$gene))
pLI_loeuf_dat <- select(pLI_loeuf_dat, gene, pLI, oe_lof_upper)
saveRDS(pLI_loeuf_dat, file = "data/pLI_LOEUF_data.rds")

pHI_pTS_dat <- read.csv("data/raw/haplo_triplo_scores_raw.tsv", sep = '\t',
                          header = TRUE)
length(unique(pHI_pTS_dat$X.gene))
saveRDS(pHI_pTS_dat, file = "data/pHaplo_pTriplo_data.rds")
