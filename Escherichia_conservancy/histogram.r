#install.packages("ggplot2",repos = "http://cran.us.r-project.org")
library("ggplot2")

stats <- read.csv("Gene_Proportions.csv")
meanProp <- stats$Mean.Prop
medProp <-stats$Med.Prop

prop_of_interest = 0.9991

tamA_proportion = 0.9889

hist_plot <-ggplot(stats,aes(x=Mean.Prop)) +
  geom_histogram(aes(fill = Mean.Prop > .9815),
                 binwidth = 0.0015, show.legend = F) +
		 theme_classic() + 
		 scale_color_manual(values = c("gray", "black")) +
		 scale_fill_manual(values = c("gray", "black")) +
      	         labs(x = "Sequence Conservation", y = "Number of Proteins") +
		 theme(axis.text = element_text(size = 18),  
                 axis.title = element_text(size = 20)) +
  geom_vline(aes(xintercept=tamA_proportion), color = "blue") +
  geom_vline(aes(xintercept=prop_of_interest), color = "red")
  
  

hist_plot

n <- nrow(subset(stats, Mean.Prop > 0.9815))
prop_greater <- n / nrow(stats)
