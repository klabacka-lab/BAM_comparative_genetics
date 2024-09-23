# 1) Loops through csv files in /data.
# 2) Saves generated charts in /picture_data

# This should allow us to quickly make changes to all charts rather than having to 
# remake them one at a time if we want to tweak something

library(ggplot2)
library(readr)
library(fs)


# Conserved value and greater are highlighted on the graph
conserved_value <- 98


file_paths <-fs::dir_ls("data")
make_plt <- function(file_name,conserved_value){

  # Generating file name for the saved .png
  pic_name <- paste('picture',substr(file_name, 1, nchar(file_name)-3),sep = '_')
  pic_name <- paste(pic_name,'png',sep = '')


  # R doesn't recognise 'N/A'. This just replaces them with NA
  data <- read_csv(file_name,col_names = TRUE)
  data[data == 'N/A'] <- NA
  data <-na.omit(data)

  # Couldn't get chart to work with proportion values less than 1
  data$Proportion <- as.numeric(data$Proportion) * 100
  
  # I did my best here, but I have no clue how to make this thing look great
  g<-ggplot(data,aes(x=Position,y=Proportion))+
    geom_col(aes(fill = Proportion > conserved_value))+
    theme(legend.position = "none") +
    scale_fill_manual(values= c('gray32','tomato3'))
  
  # Set the dimentions and resolution of the plots here
  ggsave(filename = pic_name,
         plot = g,
         width = 20, height = 4, dpi = 250, units = "in",device = "png",
         limitsize = FALSE)  
}

for (file_name in file_paths){
  make_plt(file_name,conserved_value)
}




