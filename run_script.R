# mol_1_name,mol_2_name,max_cor,global_cor,all_cor,all_cor_p,shift_time_numeric,shift_time
# cortisol_1,cytokine_1,,,,,,

library(progress)
library(dplyr)

folder_path <- "omics_run_test/"
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
# 提取不带路径的文件名
file_names <- basename(file_list)  # 例如 "sample1.csv"

# 去掉扩展名（".csv"）
clean_names <- tools::file_path_sans_ext(file_names)  # 例如 "sample1"

# 给 data_list 命名
names(data_list) <- clean_names
#将list中每个数据框的time列全部转化为POSIXct格式
data_list <- lapply(data_list, function(df) {
  df$time <- as.POSIXct(df$time)
  return(df)
})

# 计算总组合数（C(n,2)）
n_molecules <- length(data_list)
total_pairs <- choose(n_molecules, 2)

# 初始化进度条
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) | ETA: :eta | Pair [:mol1,:mol2]",
  total = total_pairs,
  clear = FALSE,
  width = 80
)

# 外部：先建空 list
result_list <- list()
pair_count <- 0  # 当前进度计数器


start_time <- Sys.time()

for (i in 1:16) {
  for (j in (i+1):17) {
    
    # 更新进度条前，设置当前分子名（会显示）
    pb$tick(tokens = list(mol1 = i, mol2 = j))
    pair_count <- pair_count + 1
    
    #提取分子名
    mol_1_name <- names(data_list[i])
    mol_2_name <- names(data_list[j])
    #给calculate_lagged_correlation用到的参数赋值
    mol_1 <- data_list[[mol_1_name]][[-1]]
    time_1 <- data_list[[mol_1_name]]$time
    mol_2 <- data_list[[mol_2_name]][[-1]]
    time_2 <- data_list[[mol_2_name]]$time
    #开始计算
    result <-
      calculate_lagged_correlation(
        x = mol_1,
        y = mol_2,
        time1 = time_1,
        time2 = time_2,
        time_tol = 0.2,
        step = 2 / 60,
        min_matched_sample = 10,
        threads = 16,
        cor_method = "spearman"
      )
    #从result中提取有用的值并进行处理
    max_cor <- extract_max_cor(object = result)
    global_cor <- extract_global_cor(object = result)
    all_cor <- extract_all_cor(result)
    all_cor_p <- extract_all_cor_p(result)
    shift_time_numeric <- extract_shift_time(result, numeric = TRUE)
    shift_time <- extract_shift_time(result, numeric = FALSE)
    
    max_cor <- as.data.frame(max_cor)
    global_cor <- as.data.frame(global_cor)
    all_cor <- as.data.frame(all_cor)
    all_cor_p <- as.data.frame(all_cor_p)
    shift_time <- as.data.frame(shift_time)
    shift_time_numeric <- as.data.frame(shift_time_numeric)
    #创建一个空的临时list用于存放向量结果
    result_df_list <- list()
    
    result_df_list[["max_cor"]] <- max_cor
    result_df_list[["global_cor"]] <- global_cor
    result_df_list[["all_cor"]] <- all_cor
    result_df_list[["all_cor_p"]] <- all_cor_p
    result_df_list[["shift_time"]] <- shift_time
    result_df_list[["shift_time_numeric"]] <- shift_time_numeric
    
    # 构造单行结果
    single_result <- list(
      mol_1_name = mol_1_name,
      mol_2_name = mol_2_name
    )
    
    # 添加 max_cor
    max_cor <- result_df_list[["max_cor"]]
    max_cor_named <- setNames(as.numeric(max_cor[, 1]), paste0(rownames(max_cor), "_max_cor"))
    single_result <- c(single_result, max_cor_named)
    
    # 添加 global_cor
    global_cor <- result_df_list[["global_cor"]]
    global_cor_named <- setNames(as.numeric(global_cor[, 1]), paste0(rownames(global_cor), "_global_cor"))
    single_result <- c(single_result, global_cor_named)
    
    # 添加 all_cor
    all_cor <- result_df_list[["all_cor"]]
    all_cor_named <- setNames(as.numeric(all_cor[, 1]), paste0(rownames(all_cor), "_all_cor"))
    single_result <- c(single_result, all_cor_named)
    
    # 添加 all_cor_p
    all_cor_p <- result_df_list[["all_cor_p"]]
    all_cor_p_named <- setNames(as.numeric(all_cor_p[, 1]), paste0(rownames(all_cor_p), "_all_cor_p"))
    single_result <- c(single_result, all_cor_p_named)
    
    # shift_time（逗号连接）
    shift_time_vec <- result_df_list[["shift_time"]][, 1]
    single_result[["shift_time"]] <- paste(shift_time_vec, collapse = ",")
    
    # shift_time_numeric（逗号连接）
    shift_time_numeric_vec <- result_df_list[["shift_time_numeric"]][, 1]
    single_result[["shift_time_numeric"]] <- paste(shift_time_numeric_vec, collapse = ",")
    
    # 转换为单行数据框，保留原始列名
    single_result_df <- as.data.frame(single_result, stringsAsFactors = FALSE, check.names = FALSE)
    
    result_list[[pair_count]] <- single_result_df

  }
}

result_df <- bind_rows(result_list)

write.csv(result_df, file = "result.csv", sep = ",",row.names = FALSE)

end_time <- Sys.time()
elapsed <- end_time - start_time
cat("运行结束，耗时：", round(elapsed, 2), "秒\n")



# 
# mol_1_name <- names(data_list[1])
# mol_2_name <- names(data_list[2])
# 
# mol_1 <- data_list[[mol_1_name]][[-1]]
# time_1 <- data_list[[mol_1_name]]$time
# mol_2 <- data_list[[mol_2_name]][[-1]]
# time_2 <- data_list[[mol_2_name]]$time
# 
# result <-
#   calculate_lagged_correlation(
#     x = mol_1,
#     y = mol_2,
#     time1 = time_1,
#     time2 = time_2,
#     time_tol = 0.2,
#     step = 2/60,
#     min_matched_sample = 10,
#     threads = 16,
#     cor_method = "spearman"
#   )
# 
# max_cor <- extract_max_cor(object = result)
# global_cor <- extract_global_cor(object = result)
# all_cor <- extract_all_cor(result)
# all_cor_p <- extract_all_cor_p(result)
# shift_time_numeric <- extract_shift_time(result, numeric = TRUE)
# shift_time <- extract_shift_time(result, numeric = FALSE)
# 
# max_cor <- as.data.frame(max_cor)
# global_cor <- as.data.frame(global_cor)
# all_cor <- as.data.frame(all_cor)
# all_cor_p <- as.data.frame(all_cor_p)
# shift_time <- as.data.frame(shift_time)
# shift_time_numeric <- as.data.frame(shift_time_numeric)
# #创建一个空的临时list用于存放向量结果
# result_df_list <- list()
# 
# result_df_list[["max_cor"]] <- max_cor
# result_df_list[["global_cor"]] <- global_cor
# result_df_list[["all_cor"]] <- all_cor
# result_df_list[["all_cor_p"]] <- all_cor_p
# result_df_list[["shift_time"]] <- shift_time
# result_df_list[["shift_time_numeric"]] <- shift_time_numeric
# 
# 
# 
# result_df <- result_df %>%
#   add_row()
