library(ggplot2)

fel <- read.csv("HyPhy_FEL_output_dataframe.csv")

# Calculate omega
fel$omega <- with(fel, ifelse(alpha == 0, NA, beta / alpha))

# Flag significant sites
fel$Significant <- fel$p_value < 0.05

fel$E.coli_position <- seq_len(nrow(fel))

ggplot(fel, aes(x = E.coli_position, y = omega, color = Significant)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "black",
                                "TRUE" = "red"),
                     guide = "none") +
  labs(x = "Codon position",
       y = expression(omega~"(dN/dS)")) +
  theme_classic()
