#!/usr/bin/env Rscript

# Define the file path
file_path <- "/mnt/storage/dept/medonc/beroukhim/mia/J-GaTSV/svm_train/example_test_scaled.rds"

# Check if the file exists before proceeding
if (file.exists(file_path)) {
  cat("Reading RDS file...\n")
  data <- readRDS(file_path)
  
  # Determine the number of rows
  num_rows <- nrow(data)
  
  # Add the 'sv_class' column with random 0s and 1s
  # sample() randomly selects from c(0, 1) with replacement
  data$sv_class <- sample(c(0, 1), size = num_rows, replace = TRUE)
  
  # Save the modified data back to the same location
  cat("Adding 'sv_class' column and saving file...\n")
  saveRDS(data, file_path)
  
  cat("Success! Column added and file updated.\n")
} else {
  stop(paste("Error: File not found at", file_path))
}