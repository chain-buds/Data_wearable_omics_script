omics_sample_info <- read.csv2(file = "Data_wearable_omics/omics_data/sample_info.csv",sep = ",")
omics_variable_info <- read.csv2(file = "Data_wearable_omics/omics_data/variable_info.csv",sep = ",")
omics_expression_data <- read_csv("Data_wearable_omics/omics_data/expression_data.csv")

# 确保 sample_id 是字符型以便匹配
sample_ids <- as.character(colnames(omics_expression_data))

# 创建映射向量
id_to_time <- setNames(omics_sample_info$accurate_time, omics_sample_info$sample_id)

# 替换列名
colnames(omics_expression_data) <- id_to_time[sample_ids]

omics_all <- cbind(omics_variable_info[, 1, drop = FALSE], omics_expression_data)

write.csv(omics_all, file = "omics_all.csv", row.names = FALSE)



# 创建输出文件夹
output_dir <- "work_dir"

# 提取时间点（列名除 variable_id）
time_points <- colnames(omics_all)[-1]

# 初始化 mapping 数据框
mapping <- data.frame(variable_id = character(), filename = character(), stringsAsFactors = FALSE)

# 循环写入每个变量的 CSV
for (i in seq_len(nrow(omics_all))) {
  variable_id <- as.character(omics_all$variable_id[i])
  values <- as.numeric(omics_all[i, -1])
  
  df_out <- data.frame(
    time = time_points,
    value = values
  )
  colnames(df_out)[2] <- variable_id
  
  # 清洗 variable_id 生成合法文件名
  safe_filename <- gsub("[^A-Za-z0-9_]", "_", variable_id)
  full_filename <- paste0(safe_filename, ".csv")
  file_path <- file.path(output_dir, full_filename)
  
  # 写入 CSV 文件
  write.csv(df_out, file = file_path, row.names = FALSE)
  
  # 添加到 mapping 表
  mapping <- rbind(mapping, data.frame(
    variable_id = variable_id,
    filename = full_filename,
    stringsAsFactors = FALSE
  ))
}

# 保存 mapping.csv 到输出文件夹
write.csv(mapping, file = file.path(output_dir, "mapping.csv"), row.names = FALSE)
