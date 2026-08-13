# Read the file
lines <- readLines("HyPhy_FEL_output_dataframe.txt")

# Keep only the data rows
lines <- lines[grepl("^\\|\\s*[0-9]", lines)]

# Remove leading/trailing |
lines <- gsub("^\\||\\|$", "", lines)

# Split each row on |
dat <- do.call(rbind, strsplit(lines, "\\|"))

# Trim whitespace
dat <- apply(dat, 2, trimws)

# Convert to data frame
fel <- data.frame(
  Codon = as.integer(dat[,1]),
  Partition = as.integer(dat[,2]),
  alpha = as.numeric(dat[,3]),
  beta = as.numeric(dat[,4]),
  LRT = as.numeric(dat[,5]),
  Selection = sub("\\..*", "", dat[,6]),                # Neg or Pos
  p_value = as.numeric(sub(".*=\\s*", "", dat[,6])),
  stringsAsFactors = FALSE
)

# View the result
head(fel)

# Save to CSV
write.csv(fel, "HyPhy_FEL_output_dataframe.csv", row.names = FALSE)

