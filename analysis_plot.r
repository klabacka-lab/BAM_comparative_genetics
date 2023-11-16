# Shebang line

library(fs)
library(ggplot2)
library(gridExtra)
library(readr)

make_plots <- function(filename, conserved_value){

    # Create output filenames for both plots
    base_name <- tools::file_path_sans_ext(filename)
    hist_filename <- paste(base_name, sep="", "_histogram.png")
    box_filename <- paste(base_name, sep="", "_box.png")


    # Filter data for histogram plot
    data <- read_csv(filename, col_names = TRUE)
    data$Proportion <- as.numeric(data$Proportion) * 100

    # Filter data for box and whisker plot
    domain_compariosn <- lm(Proportion~Domain, data)
    filtered_dat <- subset(data, Domain!="")

    # Create histogram
    hist_plot <- ggplot(data, aes(x=Position, y=Proportion)) +
        geom_col(aes(fill = Proportion > conserved_value)) +
        theme(legend.position = "none") +
        scale_fill_manual(values = c("gray32", "tomato3"))

    # Create box and whisker plot
    box_plot <- ggplot(filtered_dat, aes(x=Domain, y=Proportion)) + geom_boxplot()
    plot1 <- box_plot + geom_jitter(shape=16, position=position_jitter(0.2))
    plot2 <- ggplot(filtered_dat, aes(x=Domain, y=Proportion)) +
        geom_violin()
    combined_plots <- grid.arrange(plot1, plot2, nrow=2)

    # Save histogram plot
    ggsave(filename = hist_filename,
           plot = hist_plot,
           width = 20, height = 4, dpi = 250, units = "in", device = "png",
           limitsize = FALSE)

    # Save box and whisker plot
    ggsave(filename = box_filename,
           plot = combined_plots,
           width = 10, height = 8, units = "in", device = "png")
}

main <- function()
    setwd("../BAM_comparative_genetics/alignments/statistics/farmed_stats")
    file_paths <- fs::dir_ls(path = "../BAM_comparative_genetics/alignments/statistics/farmed_stats", regexp = "*stats.csv$")
    conserved_value <- 98
    for (filename in file_paths){
	filename
        make_plots(filename, conserved_value)
    }

