####
if (!requireNamespace("stringdist", quietly = TRUE)) {
  install.packages("stringdist")
}

get_spliced_coordinates <- function(query_seq, subject_dna, region_start, region_end, strand) {
  
  s_char_full <- as.character(subject_dna)
  q_char_full <- query_seq
  q_len <- nchar(q_char_full)
  
  blocks <- data.frame()
  curr_q <- 1
  
  while(curr_q <= q_len) {
    seed_len <- min(20, q_len - curr_q + 1)
    seed <- substr(q_char_full, curr_q, curr_q + seed_len - 1)
    
    matches <- matchPattern(seed, subject_dna)
    
    if (length(matches) == 0) {
      warning(paste("Warning: Break at query pos", curr_q))
      break
    }
    
    s_start_rel <- start(matches)[1]
    
    max_ext_len <- min(q_len - curr_q + 1, length(subject_dna) - s_start_rel + 1)
    sub_q <- substr(q_char_full, curr_q, curr_q + max_ext_len - 1)
    sub_s <- substr(s_char_full, s_start_rel, s_start_rel + max_ext_len - 1)
    
    vec_q <- charToRaw(sub_q)
    vec_s <- charToRaw(sub_s)
    mismatch_idx <- which(vec_q != vec_s)
    
    if(length(mismatch_idx) > 0) match_len <- mismatch_idx[1] - 1
    else match_len <- max_ext_len
    
    s_end_rel <- s_start_rel + match_len - 1
    
    if (strand == "-") {
      g_start <- region_end - s_end_rel + 1
      g_end   <- region_end - s_start_rel + 1
    } else {
      g_start <- region_start + s_start_rel - 1
      g_end   <- region_start + s_end_rel - 1
    }
    
    blocks <- rbind(blocks, data.frame(
      q_start = curr_q,
      q_end = curr_q + match_len - 1,
      genomic_start = g_start,
      genomic_end = g_end,
      len = match_len
    ))
    
    curr_q <- curr_q + match_len
  }
  return(blocks)
}

####
process_granges_data <- function(gr_obj) {
  
  message("正在从 BSgenome 提取基因组序列...")
  
  gr_plus <- gr_obj
  strand(gr_plus) <- "+" 
  genome_seqs_plus <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_plus)
  
  message("开始处理 ", length(gr_obj), " 条数据...")
  
  
  results_list <- pblapply(seq_along(gr_obj), function(i) {
    
    curr_id <- mcols(gr_obj)$id[i]
    curr_q_seq <- mcols(gr_obj)$query_sequence[i]
    curr_strand <- as.character(strand(gr_obj)[i])
    curr_region_start <- start(gr_obj)[i]
    curr_region_end <- end(gr_obj)[i]
    curr_chr <- as.character(seqnames(gr_obj)[i])
    
    raw_genome_seq <- genome_seqs_plus[[i]]
    
    if (curr_strand == "-") {
      subject_seq <- reverseComplement(raw_genome_seq)
    } else {
      subject_seq <- raw_genome_seq
    }
    
    blocks <- get_spliced_coordinates(
      query_seq = curr_q_seq, 
      subject_dna = subject_seq, 
      region_start = curr_region_start, 
      region_end = curr_region_end, 
      strand = curr_strand
    )
    
    if(nrow(blocks) > 0) {
      blocks$id <- curr_id
      blocks$chr <- curr_chr
      blocks$strand <- curr_strand
      blocks$hg38_coord <- paste0(blocks$chr, ":", blocks$genomic_start, "-", blocks$genomic_end)
    }
    
    return(blocks)
  })
  
  final_df <- do.call(rbind, results_list)
  final_df <- final_df %>% dplyr::select(id, chr, strand, hg38_coord, len, q_start, q_end, everything())
  
  return(final_df)
}

####
get_relative_m6a_pos <- function(m6a_gr, query_gr) {
  
  # 1. 寻找重叠
  # queryHits 是 m6a_gr 的索引，subjectHits 是 query_gr 的索引
  hits <- findOverlaps(m6a_gr, query_gr)
  
  if (length(hits) == 0) {
    message("未发现重叠位点。")
    return(GRanges())
  }
  
  # 2. 提取匹配的对象
  m6a_matched <- m6a_gr[queryHits(hits)]
  query_matched <- query_gr[subjectHits(hits)]
  
  # 3. 获取计算所需坐标
  m6a_coord <- start(m6a_matched)
  q_start  <- start(query_matched)
  q_end    <- end(query_matched)
  q_strand <- as.character(strand(query_matched))
  
  # 4. 计算相对位置
  # 正链：从起始端往后算 (m6A - Start + 1)
  # 负链：从终止端往前算 (End - m6A + 1)
  # "*" 无明确链向时，默认按正链处理
  rel_pos <- ifelse(q_strand == "-", 
                    q_end - m6a_coord + 1, 
                    m6a_coord - q_start + 1)
  
  # 5. 整合结果
  # 保留原始 m6A 的元数据
  res <- m6a_matched
  
  # 添加来自 query_gr 的 ID 信息（假设 query_gr 中有 id 列）
  if ("id" %in% colnames(mcols(query_gr))) {
    mcols(res)$query_id <- mcols(query_matched)$id
  }
  
  # 添加计算出的相对位置
  mcols(res)$relative_pos <- rel_pos
  
  return(res)
}

####
# find_and_map_drach <- function(result_gr, seq_df) {
#   
#   # 1. 转换 result_gr 为 Data Frame 以便处理
#   exons_df <- as.data.frame(result_gr) %>%
#     dplyr::select(seqnames, start, end, strand, id, q_start, q_end) %>%
#     dplyr::rename(exon_g_start = start, exon_g_end = end) 
#   
#   # 2. 在序列中查找 DRACH motif (已修复报错)
#   message("1. 正在序列中查找 DRACH motif...")
#   
#   motif_pattern <- "[AGT][AG]AC[ACT]" 
#   
#   drach_in_seq <- seq_df %>%
#     mutate(
#       # 原始查找结果 (List of Matrices)
#       matches_raw = str_locate_all(query_sequence, motif_pattern),
#       
#       # --- 关键修复 ---
#       # 强制将矩阵列表转换为数据框列表，确保有 "start" 和 "end" 列
#       matches_df = map(matches_raw, as.data.frame)
#       # ----------------
#     ) %>%
#     # 展开列表
#     unnest(matches_df) %>%
#     # 现在可以安全地选择列了
#     dplyr::select(id, motif_q_start = start, motif_q_end = end) %>%
#     # 生成唯一 motif ID
#     mutate(motif_id = paste0(id, "_motif_", row_number()))
#   
#   message(paste0("   共找到 ", nrow(drach_in_seq), " 个 motif。"))
find_and_map_drach <- function(result_gr, seq_df) {
  
  # 1. 准备外显子信息
  exons_df <- as.data.frame(result_gr) %>%
    dplyr::select(seqnames, start, end, strand, id, q_start, q_end) %>%
    dplyr::rename(exon_g_start = start, exon_g_end = end) 
  
  message("1. 正在序列中查找 DRACH motif...")
  motif_pattern <- "[AGT][AG]AC[ACT]" 
  
  # 2. 核心修复逻辑
  drach_in_seq <- seq_df %>%
    as_tibble() %>%
    # --- 关键步骤：只保留 id 和序列，剔除原始的 start/end 以免冲突 ---
    dplyr::select(id, query_sequence) %>% 
    mutate(
      matches_raw = str_locate_all(query_sequence, motif_pattern),
      matches_df = map(matches_raw, as.data.frame)
    ) %>%
    # 显式调用 tidyr 并展开
    tidyr::unnest(matches_df) %>% 
    # 此时的 start/end 是来自 str_locate_all 的相对位置
    dplyr::select(id, motif_q_start = start, motif_q_end = end) %>%
    dplyr::mutate(motif_id = paste0(id, "_motif_", row_number()))
  
  message(paste0("   共找到 ", nrow(drach_in_seq), " 个 motif。"))
  # 3. 将 motif 映射回基因组坐标
  message("2. 正在将 motif 映射回基因组坐标...")
  
  mapped_df <- drach_in_seq %>%
    # 将 motif 数据与外显子数据根据 ID 合并
    inner_join(exons_df, by = "id", relationship = "many-to-many") %>%
    # 筛选：找到该 motif 具体落在哪一个外显子区间内
    filter(motif_q_start >= q_start & motif_q_end <= q_end) %>%
    # 计算基因组坐标
    mutate(
      genomic_motif_start = if_else(
        strand == "-",
        # 负链逻辑
        exon_g_end - (motif_q_end - q_start),
        # 正链逻辑
        exon_g_start + (motif_q_start - q_start)
      ),
      genomic_motif_end = if_else(
        strand == "-",
        exon_g_end - (motif_q_start - q_start),
        exon_g_start + (motif_q_end - q_start)
      )
    )
  
  # 4. 转换为 GRanges 输出
  # 注意：GRanges 要求 start <= end。
  # 负链计算结果中 start 可能会比 end 大（如果手动算反了），
  # 但上面的公式算出来 genomic_motif_start 是对应 sequence 5' 端的物理位置。
  # 在构建 GRanges 时，必须保证 start 是较小值，end 是较大值。
  
  final_df <- mapped_df %>%
    mutate(
      final_start = pmin(genomic_motif_start, genomic_motif_end),
      final_end   = pmax(genomic_motif_start, genomic_motif_end)
    )
  
  drach_gr <- makeGRangesFromDataFrame(
    final_df,
    seqnames.field = "seqnames",
    start.field = "final_start",
    end.field = "final_end",
    strand.field = "strand",
    keep.extra.columns = TRUE 
  )
  
  return(drach_gr)
}

####
valid_bases <- list(D=c("A","G","T"), R=c("A","G"), A=c("A"), C=c("C"), H=c("A","C","T"))
real_drach_seqs <- expand.grid(valid_bases) %>% apply(1, paste0, collapse = "")
bases <- c("A", "C", "G", "T")
all_5mers <- expand.grid(p1=bases, p2=bases, p3=bases, p4=bases, p5=bases) %>% apply(1, paste0, collapse = "")
dna_all <- DNAStringSet(all_5mers)         
dna_real <- DNAStringSet(real_drach_seqs)  
combined_set <- c(dna_all, dna_real)
full_dist_matrix <- as.matrix(stringDist(combined_set, method = "hamming"))
n_all <- length(dna_all)
n_real <- length(dna_real)
sub_mat <- full_dist_matrix[1:n_all, (n_all + 1):(n_all + n_real)]
min_dists <- apply(sub_mat, 1, min)
target_1mismatch_kmers <- all_5mers[min_dists == 1]

####
find_and_map_near_drach <- function(result_gr, seq_df, target_kmers) {
  if (!is.data.frame(seq_df)) {
    stop("seq_df 必须是 data.frame / tibble，当前类型是: ",
         paste(class(seq_df), collapse = ", "))
  }
  
  if (!all(c("id", "query_sequence") %in% colnames(seq_df))) {
    stop("seq_df 必须包含列: id, query_sequence")
  }
  
  if (nrow(seq_df) == 0) {
    stop("seq_df 为空，无法扫描 motif")
  }
  # 1. 准备基因组外显子信息
  exons_df <- as.data.frame(result_gr) %>%
    dplyr::select(seqnames, start, end, strand, id, q_start, q_end) %>%
    dplyr::rename(exon_g_start = start, exon_g_end = end)
  
  message("1. 正在序列中扫描 1-mismatch motif (滑动窗口法)...")
  
  # 为了提高速度，我们使用 base R 的 substring 进行向量化滑动窗口
  # 这种方法比循环快得多
  
  scan_results_list <- lapply(1:nrow(seq_df), function(i) {
    curr_id <- seq_df$id[i]
    curr_seq <- seq_df$query_sequence[i]
    seq_len <- nchar(curr_seq)
    
    if(seq_len < 5) return(NULL)
    
    # 生成所有 5bp 切片
    starts <- 1:(seq_len - 4)
    ends <- 5:seq_len
    windows <- substring(curr_seq, starts, ends)
    
    # 检查哪些切片在我们的白名单里
    # %in% 极快，因为利用了哈希表
    hit_idx <- which(windows %in% target_kmers)
    
    if(length(hit_idx) == 0) return(NULL)
    
    # 返回命中结果
    data.frame(
      id = curr_id,
      motif_seq = windows[hit_idx],
      motif_q_start = starts[hit_idx],
      motif_q_end = ends[hit_idx]
    )
  })
  
  # 合并所有扫描结果
  near_drach_in_seq <- do.call(rbind, scan_results_list)
  
  message(paste0("   共在序列中找到 ", nrow(near_drach_in_seq), " 个潜在位点。"))
  
  # 2. 映射回基因组 (复用之前的逻辑)
  message("2. 正在映射回基因组坐标...")
  
  mapped_df <- near_drach_in_seq %>%
    # 生成唯一ID
    mutate(motif_id = paste0(id, "_mutmotif_", row_number())) %>%
    # 联立外显子表
    inner_join(exons_df, by = "id", relationship = "many-to-many") %>%
    # 过滤：motif 必须完整落在一个外显子内 (不能跨越剪接位点)
    filter(motif_q_start >= q_start & motif_q_end <= q_end) %>%
    mutate(
      # 坐标换算 (正/负链逻辑)
      genomic_motif_start = if_else(
        strand == "-",
        exon_g_end - (motif_q_end - q_start),
        exon_g_start + (motif_q_start - q_start)
      ),
      genomic_motif_end = if_else(
        strand == "-",
        exon_g_end - (motif_q_start - q_start),
        exon_g_start + (motif_q_end - q_start)
      )
    )
  
  # 3. 构建 GRanges 输出
  # 同样需要处理 start/end 大小问题
  final_df <- mapped_df %>%
    mutate(
      final_start = pmin(genomic_motif_start, genomic_motif_end),
      final_end   = pmax(genomic_motif_start, genomic_motif_end)
    )
  
  near_drach_gr <- makeGRangesFromDataFrame(
    final_df,
    seqnames.field = "seqnames",
    start.field = "final_start",
    end.field = "final_end",
    strand.field = "strand",
    keep.extra.columns = TRUE 
  )
  
  return(near_drach_gr)
}

####
get_mutation_info_table <- function(observed_seqs, real_targets) {
  
  # 去重，只对唯一的序列进行计算
  unique_seqs <- unique(observed_seqs)
  dna_obs <- DNAStringSet(unique_seqs)
  
  # 计算距离矩阵 (Observed x Real)
  # 既然数据量不大，直接用 stringDist 没问题
  # 注意：这里我们合并计算以规避之前的报错，或者因为只有 18 个 target，我们可以循环
  
  info_list <- lapply(unique_seqs, function(seq_str) {
    # 将当前序列与所有 18 个标准 DRACH 比较
    # stringDist 支持 character vector 对 character vector
    dists <- stringdist::stringdist(seq_str, real_targets, method = "hamming")
    
    # 找到距离为 1 的那个目标 (可能匹配多个，取第一个即可，或者取最有可能的)
    match_idx <- which(dists == 1)[1]
    
    if (is.na(match_idx)) return(NULL) # 理论上不应该发生，因为输入就是筛选过的
    
    target_seq <- real_targets[match_idx]
    
    # 逐位比较，找到不同的位置
    v_obs <- charToRaw(seq_str)
    v_tgt <- charToRaw(target_seq)
    diff_pos <- which(v_obs != v_tgt)[1] # 突变的相对位置 (1-5)
    
    return(data.frame(
      motif_seq = seq_str,
      target_drach = target_seq,
      rel_mut_index = diff_pos # 1-based index (1,2,3,4,5)
    ))
  })
  
  do.call(rbind, info_list)
}

####
annotate_drach_mutations <- function(near_drach_gr) {
  
  message("1. 正在分析 motif 的突变模式...")
  
  # A. 提取 GRanges 中的序列列
  # 假设你的 near_drach_gr 的 metadata 里有一列叫 motif_seq
  if(is.null(near_drach_gr$motif_seq)) stop("GRanges 中必须包含 'motif_seq' 列")
  
  current_seqs <- as.character(near_drach_gr$motif_seq)
  
  # B. 生成对照表 (Seq -> Mutation Index)
  lookup_table <- get_mutation_info_table(current_seqs, real_drach_seqs)
  
  # C. 将对照表合并回 GRanges 的 Metadata
  # 为了保持 GRanges 顺序不变，我们使用 match
  match_indices <- match(near_drach_gr$motif_seq, lookup_table$motif_seq)
  
  near_drach_gr$target_drach  <- lookup_table$target_drach[match_indices]
  near_drach_gr$rel_mut_index <- lookup_table$rel_mut_index[match_indices]
  
  message("2. 正在计算 hg38 绝对坐标...")
  
  # D. 计算基因组坐标
  # -----------------------------------------------------------
  # 逻辑：
  # 正链 (+): Start 是 1, End 是 5。坐标 = Start + (Index - 1)
  # 负链 (-): End 是 1, Start 是 5。坐标 = End - (Index - 1)
  # -----------------------------------------------------------
  
  strands <- as.character(strand(near_drach_gr))
  starts  <- start(near_drach_gr)
  ends    <- end(near_drach_gr)
  mut_idxs <- near_drach_gr$rel_mut_index
  
  # 1. 计算突变位点 (Mutation Site) 的坐标
  near_drach_gr$hg38_mut_pos <- ifelse(
    strands == "-",
    ends - (mut_idxs - 1),   # 负链：从 End 往回数
    starts + (mut_idxs - 1)  # 正链：从 Start 往后数
  )
  
  # 2. 计算 A 位点 (Index = 3) 的坐标
  # 无论突变在哪里，DRACH 里的 A 始终在第 3 位
  near_drach_gr$hg38_A_pos <- ifelse(
    strands == "-",
    ends - (3 - 1),          # 负链：End - 2
    starts + (3 - 1)         # 正链：Start + 2
  )
  
  return(near_drach_gr)
}

####
add_query_context <- function(gr_obj) {
  
  message("正在计算 Query 序列上的绝对坐标和碱基变化...")
  
  # 1. 提取必要信息
  # motif_q_start: 之前计算的，motif 在拼接序列中的起始位置 (1-based)
  # rel_mut_index: 之前计算的，突变在 motif 内部的相对位置 (1-5)
  
  q_starts <- gr_obj$motif_q_start
  mut_idxs <- gr_obj$rel_mut_index
  
  # 2. 计算 Query 序列上的绝对坐标
  # -------------------------------------------------------------
  # 突变位点坐标 = Motif起始 + 相对偏移 - 1
  gr_obj$query_mut_pos <- q_starts + mut_idxs - 1
  
  # DRACH 'A' 位点坐标 = Motif起始 + 3 - 1
  # (不管突变在哪里，核心 A 永远在第 3 位)
  gr_obj$query_A_pos <- q_starts + 2 
  
  # 3. 提取碱基变化信息 (Ref -> Alt)
  # -------------------------------------------------------------
  # Ref: 当前序列 (motif_seq) 在突变位置的碱基
  # Alt: 目标序列 (target_drach) 在突变位置的碱基
  
  # 使用 substr 向量化操作，速度很快
  gr_obj$ref_base <- substr(gr_obj$motif_seq, mut_idxs, mut_idxs)
  gr_obj$alt_base <- substr(gr_obj$target_drach, mut_idxs, mut_idxs)
  
  # 生成一个易读的 mutation 字符串，例如 "C>T"
  gr_obj$mutation_code <- paste0(gr_obj$ref_base, ">", gr_obj$alt_base)
  
  return(gr_obj)
}
