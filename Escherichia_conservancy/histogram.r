#install.packages("ggplot2",repos = "http://cran.us.r-project.org")
library("ggplot2")

stats <- read.csv("Gene_Proportions.csv")
meanProp <- stats$Mean.Prop
medProp <-stats$Med.Prop

prop_of_interest = 0.99

hist_plot <-ggplot(stats,aes(x=Mean.Prop)) +
  geom_histogram(aes(fill = Mean.Prop > .995),
                 binwidth = 0.002) +
  geom_vline(aes(xintercept=prop_of_interest))

hist_plot

