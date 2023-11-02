setwd('../BAM_comparative_genetics')

# Read in the data
dat <- read.csv("enterobacteralesOA_farmed_stats.csv")
domain_comparison <- lm(Proportion~Domain, dat)

# Filter out data with no domain
filtered_dat <- subset(dat, Domain != "")
filtered_dat <- subset(filtered, Domain != "SIGNAL_SEQ")

# Plot the data
p <- ggplot(filtered_dat, aes(x=Domain, y=Proportion)) + geom_boxplot()
plot1 <- p + geom_jitter(shape=16, position=position_jitter(0.2))
plot2 <- ggplot(filtered_dat, aes(x=Domain, y=Proportion)) + geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)

# Stack plots
library(gridExtra)
grid.arrange(plot1, plot2, nrow=2)


