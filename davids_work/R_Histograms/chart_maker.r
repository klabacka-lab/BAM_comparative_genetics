# 1) Loops through csv files in /data.
# 2) Saves generated charts in /picture_data

# This should allow us to quickly make changes to all charts rather than having to 
# remake them one at a time if we want to tweak something

library(ggplot2)
library(readr)
library(fs)



file_paths <-fs::dir_ls("data")
file_paths

conserved_value <- 98

make_plt <- function(file_name,conserved_value){
  pic_name <- paste('picture',substr(file_name, 1, nchar(file_name)-3),sep = '_')
  pic_name <- paste(pic_name,'png',sep = '')
  
  data <- read_csv(file_name,col_names = TRUE)
  data[data == 'N/A'] <- NA
  data <-na.omit(data)
  data$Proportion <- as.numeric(data$Proportion) * 100
  
  g<-ggplot(data,aes(x=Position,y=Proportion))+
    geom_col(aes(fill = Proportion > conserved_value))+
    theme(legend.position = "none") +
    scale_fill_manual(values= c('gray32','tomato3'))
  
  ggsave(filename = pic_name,
         plot = g,
         width = 20, height = 4, dpi = 150, units = "in",device = "png",
         limitsize = FALSE)  
}

for (file_name in file_paths){
  make_plt(file_name,conserved_value)
}




