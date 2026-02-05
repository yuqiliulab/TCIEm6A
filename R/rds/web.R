#!/usr/bin/env Rscript

# ========================================================
# 1. 环境设置与依赖加载
# ========================================================
suppressMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicScores)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(jsonlite)
  library(rtracklayer) # 用于读取 RDS/GTF 等
  library(pbapply)
  library(tidyr)      # <--- 必须添加这一行，提供 unnest 函数
  library(purrr)
  # 显式覆盖冲突函数
  mutate <- dplyr::mutate
  select <- dplyr::select
  filter <- dplyr::filter
})

# --- 调试信息 (保留服务器习惯) ---
system("whoami > /www/wwwroot/rnamd.org/YvqiWeb/data/output/debug_user.txt")
system("which conda > /www/wwwroot/rnamd.org/YvqiWeb/data/output/debug_conda.txt")

# ========================================================
# 2. 参数解析与路径定义
# ========================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript script.R <input_rds_path> <jobID> <output_dir> <mode>")
}

seq_file_path <- normalizePath(args[1])
jobID         <- args[2]
output_dir    <- normalizePath(args[3])
input_mode    <- args[4] # --- 修改点 2：获取模式参数 ---

# --- 修改点 3：验证模式参数是否合法 ---
valid_modes <- c("mode1", "mode2", "mode3")
if (!input_mode %in% valid_modes) {
  stop(paste("Error: Invalid mode selected.", input_mode, "is not allowed. Choose from:", paste(valid_modes, collapse=", ")))
}

message(paste(">> 检测到运行模式:", input_mode))
# --- 服务器路径定义 ---
BASE_DIR   <- "/www/wwwroot/rnamd.org/YvqiWeb"
DATA_DIR   <- file.path(BASE_DIR, "data")
RDS_DIR    <- file.path(DATA_DIR, "rds")
SCRIPT_DIR <- file.path(BASE_DIR, "server_script")    # 假设 server_function.R 在此

# 创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
trans_out_dir <- file.path(output_dir, "transdecoder")
if (!dir.exists(trans_out_dir)) dir.create(trans_out_dir, recursive = TRUE)

# ========================================================
# 3. 加载核心函数与参考数据
# ========================================================

# 加载辅助函数 (假设 server_function.R 存在于服务器)
server_func_path <- file.path(SCRIPT_DIR, "server_function.R")
if(file.exists(server_func_path)){
  source(server_func_path)
} else {
  stop(paste("找不到 server_function.R:", server_func_path))
}

message(">> 正在加载参考数据 (RDS)...")

# 加载 m6A 参考集 (hg38_CL_Tech)
m6a_gr_path <- file.path(RDS_DIR, "hg38_CL_Tech.rds") 
# 或者如果是你在Target代码中用到的 hg38_CL_Tech1.rds，请修改文件名
hg38_CL_Tech <- readRDS(m6a_gr_path)
m6a_gr <- hg38_CL_Tech[, c(1, 8)] # 根据你的逻辑提取列
names(mcols(m6a_gr)) <- c('m6A_source','m6A_ID')

# 加载特征注释所需的参考集
RBP_site_gr_38 <- readRDS(file.path(RDS_DIR, "RBP_site_gr_38.rds"))
miRNA_ALL_gr   <- readRDS(file.path(RDS_DIR, "miRNA_ALL_gr.rds"))
hg38_splicing_sites_100bp_gr <- readRDS(file.path(RDS_DIR, "hg38_splicing_sites_100bp_gr.rds"))
ires_granges   <- readRDS(file.path(RDS_DIR, "ires_granges.rds"))

# ========================================================
# 修正后的 Conservation Scores 加载逻辑
# ========================================================
phastcons <- NULL
if(!exists("phastcons") || is.null(phastcons)) {
  message(">> 正在尝试加载 phastCons100way.UCSC.hg38...")
  
  # 方案 A: 尝试从 R 包加载 (如果你安装了该 Data Package)
  tryCatch({
    library(phastCons100way.UCSC.hg38)
    phastcons <- phastCons100way.UCSC.hg38
    message("   - 已通过 R 包加载成功")
  }, error = function(e) {
    
    # 方案 B: 如果没有包，尝试从 RDS 文件加载 (请确保路径正确)
    # 你可以手动指定一个路径，或者暂时设为 NULL
    phastcons_rds <- file.path(RDS_DIR, "phastcons_hg38.rds") 
    
    if(file.exists(phastcons_rds)) {
      phastcons <- readRDS(phastcons_rds)
      message("   - 已通过本地 RDS 加载成功")
    } else {
      message("Warning: 未找到 phastcons 数据源。保守性评分将设为默认值 0。")
      phastcons <- NULL
    }
  })
}
# ========================================================
# 4. 定义核心业务函数 (直接嵌入以防路径问题)
# ========================================================

get_m6a_sites_by_mode <- function(mode = c("mode1", "mode2", "mode3"), result_gr, seq_df, m6a_gr = NULL, target_1mismatch_kmers = NULL) {
  mode <- match.arg(mode)
  m6a_sites <- GRanges()
  
  # Mode 1: Known m6A
  if (mode == "mode1") {
    if (is.null(m6a_gr)) stop("mode1 需要 m6a_gr")
    m6a_sites <- get_relative_m6a_pos(m6a_gr, result_gr) # 需确保此函数在 server_function.R 中
  }
  
  # Mode 2: DRACH Motif
  if (mode == "mode2") {
    drach_gr <- find_and_map_drach(result_gr, seq_df) # 需确保此函数在 server_function.R 中
    if (length(drach_gr) > 0) {
      m6a_sites <- GRanges(
        seqnames = seqnames(drach_gr),
        ranges   = IRanges(start = drach_gr$genomic_motif_start + 2, width = 1),
        strand   = strand(drach_gr),
        query_id = drach_gr$id,
        relative_pos = drach_gr$motif_q_start + 2
      )
    }
  }
  
  # Mode 3: Near DRACH
  if (mode == "mode3") {
    near_gr <- find_and_map_near_drach(result_gr, seq_df, target_1mismatch_kmers)
    if (length(near_gr) > 0) {
      near_gr <- annotate_drach_mutations(near_gr)
      near_gr <- add_query_context(near_gr)
      m6a_sites <- GRanges(
        seqnames = seqnames(near_gr),
        ranges   = IRanges(start = near_gr$hg38_A_pos, width = 1),
        strand   = strand(near_gr),
        query_id = near_gr$id,
        relative_pos = near_gr$query_A_pos,
        mut_pos_genome = near_gr$hg38_mut_pos,
        mut_relative_pos_circRNA = near_gr$query_mut_pos,
        mut_info = near_gr$mutation_code
      )
    }
  }
  
  # 去重逻辑
  if (length(m6a_sites) > 0) {
    m6a_sites$check <- paste0(m6a_sites$query_id, ":", m6a_sites$relative_pos)
    if (mode %in% c("mode2", "mode3") && !is.null(m6a_gr)) {
      m6a_known <- get_relative_m6a_pos(m6a_gr, result_gr)
      if(length(m6a_known) > 0){
        known_checks <- paste0(m6a_known$query_id, ":", m6a_known$relative_pos)
        m6a_sites <- m6a_sites[!(m6a_sites$check %in% known_checks)]
      }
    }
  }
  return(m6a_sites)
}

annotate_candidate_features <- function(candidate_sites, RBP_ref, miRNA_ref, splicing_ref, phastcons_ref) {
  if (length(candidate_sites) == 0) return(candidate_sites)
  
  # 1. 重叠计数
  message(">> 计算 Feature Overlaps...")
  # --- A. 基础计数 (miRNA 和 Splicing 较小，可以先算) ---
  mcols(candidate_sites)$miRNA <- countOverlaps(candidate_sites, miRNA_ref)
  mcols(candidate_sites)$SplicingSite_Num <- countOverlaps(candidate_sites, splicing_ref)
  
  # --- B. 核心优化：逐个计算 RBP ---
  rbp_list <- c('HNRNPC','YTHDF3','YTHDC1','YTHDF1','YTHDF2','METTL3','METTL14','METTL16','WTAP','ALKBH5','FTO')
  
  # 初始化 RBP_Num 为 0
  mcols(candidate_sites)$RBP_Num <- 0
  
  for (rbp in rbp_list) {
    message(paste("   - 正在处理 RBP:", rbp))
    # 仅提取当前这一个 RBP 的子集
    rbp_sub <- RBP_ref[mcols(RBP_ref)$RBP_name == rbp]
    
    if(length(rbp_sub) > 0) {
      # 计算是否重叠 (0/1)
      is_overlap <- as.integer(countOverlaps(candidate_sites, rbp_sub) > 0)
      mcols(candidate_sites)[[rbp]] <- is_overlap
      # 累加到总数
      mcols(candidate_sites)$RBP_Num <- mcols(candidate_sites)$RBP_Num + is_overlap
    } else {
      mcols(candidate_sites)[[rbp]] <- 0
    }
    
    # 关键：每算完一个 RBP，立刻从内存抹除子集并回收
    rm(rbp_sub)
    gc() 
  }
  
  # 2. 保守性与序列
  if(!is.null(phastcons_ref)){
    mcols(candidate_sites)$phastcons <- gscores(phastcons_ref, candidate_sites)$default
    mcols(candidate_sites)$phastcons[is.na(mcols(candidate_sites)$phastcons)] <- 0
  } else {
    mcols(candidate_sites)$phastcons <- 0
  }
  
  expanded <- resize(candidate_sites, width = 41, fix = "center")
  mcols(candidate_sites)$seq_count <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, expanded))
  
  # 3. 清洗与归一化
  char_cols <- c("ORF", "IRES")
  num_cols <- c("ORF_Count", "ORF_Length", "m6a_in_ORF", "Dist_to_ORF", 
                "IRES_Count", "IRES_Length", "m6a_in_IRES", "Dist_to_IRES")
  
  for(col in char_cols) {
    if (col %in% colnames(mcols(candidate_sites))) {
      val <- mcols(candidate_sites)[[col]]
      mcols(candidate_sites)[[col]] <- ifelse(!is.na(val) & val != "" & val != "None", 1, 0)
    } else { mcols(candidate_sites)[[col]] <- 0 }
  }
  
  for(col in num_cols) {
    if (col %in% colnames(mcols(candidate_sites))) {
      mcols(candidate_sites)[[col]][is.na(mcols(candidate_sites)[[col]])] <- 0
    } else { mcols(candidate_sites)[[col]] <- 0 }
  }
  
  # Normalize
  count_cols <- c("ORF_Count","IRES_Count","RBP_Num","SplicingSite_Num")
  for (col in count_cols) mcols(candidate_sites)[[col]] <- pmin(mcols(candidate_sites)[[col]] / 4, 1)
  
  len_cols <- c("IRES_Length", "ORF_Length")
  for (col in len_cols) mcols(candidate_sites)[[col]] <- pmin(mcols(candidate_sites)[[col]] / 1000, 1)
  
  return(candidate_sites)
}

# ========================================================
# 5. 主流程：数据处理
# ========================================================
message(">> 正在处理输入数据...")

# --- 定义读取函数 ---
load_and_normalize_input <- function(file_path) {
  
  # 获取文件扩展名 (转为小写)
  ext <- tools::file_ext(file_path) %>% tolower()
  
  df_out <- NULL
  
  # -------------------------------------------------------
  # 场景 A: 用户上传了 FASTA (.fa, .fasta)
  # -------------------------------------------------------
  if (ext %in% c("fa", "fasta")) {
    message(">> 检测到输入为 FASTA 格式")
    # 读取 FASTA
    fa <- readDNAStringSet(file_path)
    
    # 解析 Header 获取坐标 (假设 Header 格式为: chr1:100-200(+))
    # 如果 Header 仅仅是 ID，没有坐标，后续步骤会报错，这里需要做正则提取
    ids <- names(fa)
    seqs <- as.character(fa)
    
    # 尝试提取坐标信息
    meta <- str_match(ids, "(chr[0-9XYM]+)[:|]([0-9]+)[-|_]([0-9]+).*?([+-])")
    
    # 构造 DataFrame
    df_out <- data.frame(
      id = ids,
      chr = meta[,2],
      genomic_start = as.integer(meta[,3]),
      genomic_end   = as.integer(meta[,4]),
      strand = meta[,5],
      query_sequence = seqs,
      stringsAsFactors = FALSE
    )
    
    # 如果正则提取失败（比如 Header 只有 >seq1），给默认值防止程序直接崩溃
    # 注意：如果缺失坐标，后续 Conservation/RBP 分析将无法进行
    if (any(is.na(df_out$chr))) {
      warning("部分 FASTA Header 不包含标准坐标 (chr:start-end(strand))，后续基因组映射可能失败。")
    }
    
    # -------------------------------------------------------
    # 场景 B: 用户上传了 RDS (保留旧功能)
    # -------------------------------------------------------
  } else if (ext == "rds") {
    message(">> 检测到输入为 RDS 格式")
    df_out <- readRDS(file_path)
    
    # -------------------------------------------------------
    # 场景 C: 用户在文本框输入 (.txt, .csv, .tsv)
    # -------------------------------------------------------
  } else {
    message(">> 检测到输入为 文本/表格 格式")
    # 尝试读取，支持逗号或制表符分隔
    # fill=TRUE 允许行长度不一致，header=TRUE 读取表头
    raw_df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    
    # 如果分隔符不对，尝试逗号
    if (ncol(raw_df) <= 1) {
      raw_df <- read.table(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
    }
    
    # --- 列名标准化 (Mapping) ---
    # 这一步是为了兼容用户输入的各种奇奇怪怪的列名
    colnames(raw_df) <- tolower(colnames(raw_df))
    
    df_out <- data.frame(stringsAsFactors = FALSE)
    
    # 1. 处理 ID / Seqnames
    if("id" %in% colnames(raw_df)) { df_out$id <- raw_df$id }
    else if("seqnames" %in% colnames(raw_df)) { df_out$id <- raw_df$seqnames } # 有时用户把 seqnames 当 ID
    else { df_out$id <- paste0("input_", 1:nrow(raw_df)) }
    
    # 2. 处理 Chromosome
    if("chr" %in% colnames(raw_df)) { df_out$chr <- raw_df$chr }
    else if("seqnames" %in% colnames(raw_df)) { df_out$chr <- raw_df$seqnames }
    
    # 3. 处理 Start/End
    if("start" %in% colnames(raw_df)) { df_out$genomic_start <- raw_df$start }
    else if("ranges" %in% colnames(raw_df)) { 
      # 尝试解析 ranges 列 (例如 "100-200")
      ranges <- str_split(raw_df$ranges, "[-|]", simplify = TRUE)
      df_out$genomic_start <- as.integer(ranges[,1])
      df_out$genomic_end <- as.integer(ranges[,2])
    }
    
    if("end" %in% colnames(raw_df)) { df_out$genomic_end <- raw_df$end }
    
    # 4. 处理 Strand
    if("strand" %in% colnames(raw_df)) { df_out$strand <- raw_df$strand }
    else { df_out$strand <- "+" } # 默认正链
    
    # 5. 处理 Sequence (最关键)
    if("query_sequence" %in% colnames(raw_df)) { df_out$query_sequence <- raw_df$query_sequence }
    else if("sequence" %in% colnames(raw_df)) { df_out$query_sequence <- raw_df$sequence }
    else if("seq" %in% colnames(raw_df)) { df_out$query_sequence <- raw_df$seq }
  }
  
  return(df_out)
}

# --- 调用读取函数 ---
# 无论前端传文件名叫什么，这里都统一读取
data <- load_and_normalize_input(seq_file_path)

# --- 数据完整性检查 ---
if (is.null(data) || nrow(data) == 0) stop("输入文件为空或格式无法识别。")
if (!"query_sequence" %in% colnames(data)) stop("输入数据缺失序列列 (sequence/query_sequence)。")

# 确保 ID 存在
if(!"id" %in% colnames(data)) data$id <- paste0('input_', seq_along(data[,1]))

# 如果是纯文本表格进来，可能缺失 result_table 的其他列
# 但只要有坐标和序列，process_granges_data 应该能处理（或者直接构造 GRanges）

# 构造 seq_df (供 Mode 2/3 使用)
seq_df <- data %>% dplyr::select(id, query_sequence)

# 构造 result_gr
# 如果输入提供了坐标，则使用坐标；如果没提供（只有序列），则可能无法进行 RBP 分析
if (all(c("chr", "genomic_start", "genomic_end") %in% colnames(data))) {
  result_gr <- makeGRangesFromDataFrame(
    data, 
    keep.extra.columns = TRUE, 
    ignore.strand = FALSE,
    seqnames.field = "chr", 
    start.field = "genomic_start", 
    end.field = "genomic_end", 
    strand.field = "strand"
  )
} else {
  message("Warning: 输入数据缺少基因组坐标信息，仅依靠序列进行分析。部分功能(RBP/Conservation)将被跳过。")
  # 构造一个假的 GRanges 防止报错，或者你需要在这里决定是否终止
  # 这里假设必须要坐标
  stop("错误: 输入数据必须包含基因组坐标 (chr, start, end) 用于特征注释。")
}
# 准备 seq_df (Mode 2/3 需要)
# 注意：find_and_map_drach 可能需要 seq_df 有特定列 (sequence, id 等)
# seq_df <- data # 假设输入 data 就包含 sequence 列
# ---- 关键修复：确保 seq_df 是 data.frame ----
if (inherits(data, "GRanges")) {
  seq_df <- as.data.frame(mcols(data)) %>%
    dplyr::select(id, query_sequence)
} else {
  seq_df <- data %>%
    as.data.frame() %>%
    dplyr::select(id, query_sequence)
}
if (!"q_start" %in% colnames(mcols(result_gr))) {
  mcols(result_gr)$q_start <- 1
}
if (!"q_end" %in% colnames(mcols(result_gr))) {
  # q_end 等于序列本身的长度
  mcols(result_gr)$q_end <- nchar(data$query_sequence)
}

# 确保 seq_df 包含 id 和 query_sequence
if (inherits(data, "GRanges")) {
  seq_df <- as.data.frame(mcols(data)) %>% dplyr::select(id, query_sequence)
} else {
  seq_df <- data %>% as.data.frame() %>% dplyr::select(id, query_sequence)
}
# --- [Step A] 获取 m6A 位点 ---
# 这里默认使用 Mode 1，你可以根据需要改为 "mode2" 或从 args 获取
CURRENT_MODE <- input_mode
message(paste(">> 获取 m6A 位点, Mode:", CURRENT_MODE))

m6a_sites <- get_m6a_sites_by_mode(
  mode = CURRENT_MODE,
  result_gr = result_gr,
  seq_df = seq_df,
  m6a_gr = m6a_gr,
  target_1mismatch_kmers = target_1mismatch_kmers # Mode 3 需要，如果是 Mode 1/2 可为 NULL
)



# ========================================================
# 6. ORF 分析 (TransDecoder on Server)
# ========================================================
message(">> 开始 ORF 分析 (TransDecoder)...")

# 1. 从基因组提取序列并写入临时 FASTA
circ_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, result_gr)
names(circ_sequences) <- paste0(seqnames(result_gr), ":", start(result_gr), "-", end(result_gr), "(", strand(result_gr), ")")

temp_fasta_path <- file.path(trans_out_dir, paste0(jobID, "_circ_seqs.fa"))
writeXStringSet(circ_sequences, filepath = temp_fasta_path, format = "fasta", width = 80)

# 2. 运行 TransDecoder (Conda 环境)
transdecoder_output_dir <- file.path(trans_out_dir, paste0(jobID, "_td_output"))
cmd <- paste0(
  "conda run -n transdecoder_env TransDecoder.LongOrfs ",
  "-t ", temp_fasta_path, " ",
  "--output_dir ", transdecoder_output_dir
)
message("Running: ", cmd)
system(cmd)

# 3. 解析结果 (.cds)
orf_fasta_path <- file.path(transdecoder_output_dir, "longest_orfs.cds")

if(file.exists(orf_fasta_path)) {
  orf_fasta <- readDNAStringSet(orf_fasta_path)
  
  # 解析 ORF 信息
  orf_info <- tibble(name = names(orf_fasta)) %>%
    mutate(
      id = str_extract(name, "chr[0-9XYM]+:[0-9]+-[0-9]+\\([+-]\\)"),
      last_field = str_extract(name, "[^ ]+$"),
      start_rel = as.integer(str_extract(last_field, "(?<=:)\\d+(?=-)")),
      end_rel   = as.integer(str_extract(last_field, "(?<=-)\\d+(?=\\()")),
      orf_strand = str_extract(last_field, "[+-](?=\\)$)")
    )
  
  # 聚合结果 (仅保留 sense 链预测)
  orf_aggregated <- orf_info %>%
    filter(orf_strand == "+") %>%
    group_by(id) %>%
    summarise(
      ORF_Rel = paste0(start_rel, "-", end_rel, collapse = ";"),
      ORF_Count = n(),
      .groups = "drop"
    )
  
  # 映射回 result_gr
  mcols(result_gr)$ORF_Rel <- NA_character_
  mcols(result_gr)$ORF_Count <- NA_integer_
  
  # 匹配 ID (FASTA header vs GRanges ID)
  # result_gr 的 ID 必须与 names(circ_sequences) 一致
  # 我们上面构造了 ID: chr:start-end(strand)
  result_gr_ids <- paste0(seqnames(result_gr), ":", start(result_gr), "-", end(result_gr), "(", strand(result_gr), ")")
  
  match_idx <- match(result_gr_ids, orf_aggregated$id)
  
  mcols(result_gr)$ORF_Rel[!is.na(match_idx)] <- orf_aggregated$ORF_Rel[na.omit(match_idx)]
  mcols(result_gr)$ORF_Count[!is.na(match_idx)] <- orf_aggregated$ORF_Count[na.omit(match_idx)]
  
  # 计算 Genomic ORF 坐标 (复用你的逻辑)
  message(">> 计算 ORF 绝对坐标...")
  orf_absolute <- vapply(seq_along(result_gr), function(i) {
    orf_entry <- mcols(result_gr)$ORF_Rel[i]
    if (is.na(orf_entry)) return(NA_character_)
    
    circ_start <- start(result_gr[i])
    circ_end   <- end(result_gr[i])
    strand_val <- as.character(strand(result_gr[i]))
    
    parts <- str_split(orf_entry, ";")[[1]]
    abs_ranges <- sapply(parts, function(rg) {
      pos <- as.numeric(str_split(rg, "-")[[1]])
      if (strand_val == "+") {
        a_start <- circ_start + pos[1] - 1
        a_end   <- circ_start + pos[2] - 1
      } else {
        a_start <- circ_end - pos[2] + 1
        a_end   <- circ_end - pos[1] + 1
      }
      paste0(a_start, "-", a_end)
    })
    paste0(abs_ranges, collapse = ";")
  }, character(1))
  
  mcols(result_gr)$ORF_Genomic <- orf_absolute
  
  # 计算长度
  mcols(result_gr)$ORF_Length <- sapply(mcols(result_gr)$ORF_Genomic, function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    parts <- str_split(x, ";")[[1]]
    lens <- sapply(parts, function(rg) {
      coords <- as.numeric(str_split(rg, "-", simplify = TRUE))
      abs(coords[2] - coords[1]) + 1
    })
    mean(lens)
  })
  
} else {
  message("Warning: TransDecoder output not found.")
  mcols(result_gr)$ORF_Genomic <- NA
  mcols(result_gr)$ORF_Length <- 0
  mcols(result_gr)$ORF_Count <- 0
}

# --- 将 ORF 信息映射到 m6A 位点 ---
message(">> 映射 ORF 到 m6A 位点...")
hits <- findOverlaps(m6a_sites, result_gr)
m6a_sites$ORF <- NA; m6a_sites$ORF_Count <- 0; m6a_sites$ORF_Length <- 0
m6a_sites$ORF[hits@from] <- mcols(result_gr)$ORF_Genomic[hits@to]
m6a_sites$ORF_Count[hits@from] <- mcols(result_gr)$ORF_Count[hits@to]
m6a_sites$ORF_Length[hits@from] <- mcols(result_gr)$ORF_Length[hits@to]

# 计算 Dist_to_ORF 和 m6a_in_ORF
m6a_sites$m6a_in_ORF <- 0
m6a_sites$Dist_to_ORF <- 0

# (保留你的循环逻辑)
for (i in seq_along(m6a_sites)) {
  ORF_raw <- m6a_sites$ORF[i]
  if(is.na(ORF_raw) || ORF_raw == "0" || ORF_raw == "") next
  
  m6a_pos <- start(m6a_sites[i])
  ranges <- strsplit(ORF_raw, ";")[[1]]
  
  for(r in ranges){
    coords <- as.integer(strsplit(r, "-")[[1]])
    if(length(coords) != 2) next
    
    if(m6a_pos >= coords[1] && m6a_pos <= coords[2]){
      m6a_sites$m6a_in_ORF[i] <- 1
      
      d_s <- m6a_pos - coords[1]
      d_e <- coords[2] - m6a_pos
      m6a_sites$Dist_to_ORF[i] <- d_s / (d_s + d_e)
      break 
    }
  }
}

# ========================================================
# 7. IRES 分析
# ========================================================
message(">> 开始 IRES 分析...")
# 找出 circRNA 内的 IRES
overlaps <- findOverlaps(ires_granges, result_gr, type = "within")
IRES_list <- vector("list", length(result_gr))

q_hits <- queryHits(overlaps)
s_hits <- subjectHits(overlaps)

for (i in seq_along(q_hits)) {
  IRES_list[[s_hits[i]]] <- c(IRES_list[[s_hits[i]]], paste0(start(ires_granges[q_hits[i]]), "-", end(ires_granges[q_hits[i]])))
}

mcols(result_gr)$IRES <- sapply(IRES_list, function(x) if (length(x) == 0) NA_character_ else paste(x, collapse = ";"))
mcols(result_gr)$IRES_Count <- sapply(IRES_list, length)
mcols(result_gr)$IRES_Length <- sapply(IRES_list, function(x) {
  if (length(x) == 0) return(0)
  mean(sapply(strsplit(x, "-"), function(y) as.integer(y[2]) - as.integer(y[1]) + 1))
})

# 映射到 m6A 位点
hits <- findOverlaps(m6a_sites, result_gr)
m6a_sites$IRES <- NA; m6a_sites$IRES_Count <- 0; m6a_sites$IRES_Length <- 0
m6a_sites$IRES[hits@from] <- mcols(result_gr)$IRES[hits@to]
m6a_sites$IRES_Count[hits@from] <- mcols(result_gr)$IRES_Count[hits@to]
m6a_sites$IRES_Length[hits@from] <- mcols(result_gr)$IRES_Length[hits@to]

# 计算 Dist 和 m6a_in_IRES (逻辑同 ORF)
m6a_sites$m6a_in_IRES <- 0
m6a_sites$Dist_to_IRES <- 0

for (i in seq_along(m6a_sites)) {
  IRES_raw <- m6a_sites$IRES[i]
  if(is.na(IRES_raw) || IRES_raw == "0" || IRES_raw == "") next
  
  m6a_pos <- start(m6a_sites[i])
  ranges <- strsplit(IRES_raw, ";")[[1]]
  
  for(r in ranges){
    coords <- as.integer(strsplit(r, "-")[[1]])
    if(length(coords) != 2) next
    if(m6a_pos >= coords[1] && m6a_pos <= coords[2]){
      m6a_sites$m6a_in_IRES[i] <- 1
      d_s <- m6a_pos - coords[1]
      d_e <- coords[2] - m6a_pos
      m6a_sites$Dist_to_IRES[i] <- d_s / (d_s + d_e)
      break
    }
  }
}

# ========================================================
# 8. 综合注释与导出
# ========================================================

message(">> 执行最终特征注释与导出...")
rm(ires_granges, hg38_CL_Tech)
gc()

# 2. 核心优化：只筛选出当前任务涉及的染色体，极大减小参考集规模
message(">> 正在精简 RBP 参考集以节省内存...")
current_chrs <- unique(as.character(seqnames(m6a_sites)))
RBP_site_gr_38 <- RBP_site_gr_38[seqnames(RBP_site_gr_38) %in% current_chrs]
# 强制收缩内存占用
RBP_site_gr_38 <- keepSeqlevels(RBP_site_gr_38, current_chrs, pruning.mode="coarse")
gc()
# 运行特征注释
final_gr <- annotate_candidate_features(
  m6a_sites, 
  RBP_ref = RBP_site_gr_38, 
  miRNA_ref = miRNA_ALL_gr, 
  splicing_ref = hg38_splicing_sites_100bp_gr,
  phastcons_ref = phastcons
)

# 准备导出 DataFrame
cols_to_extract <- c(
  "IRES","IRES_Count","IRES_Length","m6a_in_IRES","Dist_to_IRES",
  "ORF","ORF_Count","ORF_Length","m6a_in_ORF","Dist_to_ORF",
  "HNRNPC","YTHDC1","YTHDF1","YTHDF2","METTL3","METTL14","METTL16",
  "WTAP","ALKBH5","FTO","YTHDF3","miRNA","phastcons",
  "RBP_Num","SplicingSite_Num","seq_count"
)

# # 确保列存在
# missing_cols <- setdiff(cols_to_extract, colnames(mcols(final_gr)))
# if(length(missing_cols) > 0) {
#   for(c in missing_cols) mcols(final_gr)[[c]] <- 0
# }
# 
# df_main <- data.frame(
#   seqnames = as.character(seqnames(final_gr)),
#   ranges   = paste(start(final_gr), end(final_gr), sep = "-"),
#   strand   = as.character(strand(final_gr))
# )
# 
# df_meta <- as.data.frame(mcols(final_gr)[, cols_to_extract])
df_export <- as.data.frame(final_gr)


# 写入文件
csv_file  <- file.path(output_dir, paste0(jobID, "_features.csv"))
json_file <- file.path(output_dir, paste0(jobID, "_features.json"))

write.csv(df_export, csv_file, row.names = FALSE)
write_json(df_export, json_file)

cat("✓ 处理完成，结果已保存:\n", csv_file, "\n", json_file, "\n")