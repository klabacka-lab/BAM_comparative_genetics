library(rstatix)
library(tidyverse)

entero_bam <- read_csv("entero_align_farmed_stats.csv")
entero_tam <- read_csv("enterobact_tamA_rp_farmed_stats.csv")
proteo_bam <- read_csv("proteo_align_farmed_stats.csv")
proteo_tam <- read_csv("proteobacteria_tamA_rp_farmed_stats.csv")

combined_entero <- bind_rows(
  entero_bam %>%
    select(Position, Proportion) %>%
    mutate(type = "BAM"),
  entero_tam %>%
    select(Position, Proportion) %>%
    mutate(type = "TAM")
)

combined_proteo <- bind_rows(
  proteo_bam %>%
    select(Position, Proportion) %>%
    mutate(type = "BAM"),
  proteo_tam %>%
    select(Position, Proportion) %>%
    mutate(type = "TAM")
)

combined_entero %>%
  group_by(type) %>%
  summarize(mean = mean(Proportion),
            median = median(Proportion),
            sd = sd(Proportion),
            records = n())

combined_proteo %>%
  group_by(type) %>%
  summarize(mean = mean(Proportion),
            median = median(Proportion),
            sd = sd(Proportion),
            records = n())

ggplot(combined_entero, aes(x = Proportion, fill = type)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(x = "Proportion", y = "Count") +
  scale_fill_manual(values = c("BAM" = "blue", "TAM" = "red")) +
  theme_minimal()


ggplot(combined_proteo, aes(x = Proportion, fill = type)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(x = "Proportion", y = "Count") +
  scale_fill_manual(values = c("BAM" = "blue", "TAM" = "red")) +
  theme_minimal()