# 1. Read the file
data <- readRDS("/mnt/storage/dept/medonc/beroukhim/mia/J-GaTSV/outputs/example_test_scaled.rds")

# 2. Convert to CSV
write.csv(data, "/mnt/storage/dept/medonc/beroukhim/mia/J-GaTSV/outputs/example_test_scaled.csv", row.names = FALSE)