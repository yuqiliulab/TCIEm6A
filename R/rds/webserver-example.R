##====================================================##
## web-server example
##====================================================##
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(pbapply)
source('/home/share/bowen/circm6a/server_function.R')

#### 读取example序列
data <- readRDS('/home/share/bowen/circm6a/server_example_input.rds')
data$id <- paste0('input sequence ',seq_along(data))

result_table <- process_granges_data(data)

result_gr <- makeGRangesFromDataFrame(
  result_table,
  keep.extra.columns = TRUE,      
  ignore.strand = FALSE,          
  seqnames.field = "chr",         
  start.field = "genomic_start",  
  end.field = "genomic_end",      
  strand.field = "strand"       
)


###########################################
#### mode 1
###########################################
m6a_gr <- readRDS("/home/share/bowen/circm6a/hg38_CL_Tech.rds")
#names(mcols(m6a_gr))
m6a_gr <- m6a_gr[,c(1,8)]
names(mcols(m6a_gr)) <- c('m6A_source','m6A_ID')

mode_1_result <- get_relative_m6a_pos(m6a_gr, result_gr)
mode_1_result$check <- paste0(mode_1_result$query_id,':',mode_1_result$relative_pos)


###########################################
#### mode 2
###########################################
seq_df <- mcols(data) %>% 
  as.data.frame() %>% 
  dplyr::select(id, query_sequence) %>% 
  distinct()

drach_gr <- find_and_map_drach(result_gr, seq_df)
drach_gr$A_base <- drach_gr$motif_q_start + 2

new_m6a <- GRanges(seqnames = seqnames(drach_gr),
                   IRanges(start = drach_gr$genomic_motif_start + 2,
                           width = 1),
                   strand = strand(drach_gr),
                   relative_pos = drach_gr$A_base,
                   query_id = drach_gr$id)

new_m6a$check <- paste0(new_m6a$query_id,':',new_m6a$relative_pos)
indx <- setdiff(new_m6a$check,mode_1_result$check)
match <- match(new_m6a$check,indx)
indx2 <- which(match != 'NA')

mode_2_result <- new_m6a[indx2]


###########################################
#### mode 3
###########################################
near_drach_gr <- find_and_map_near_drach(result_gr, seq_df, target_1mismatch_kmers)

valid_bases <- list(D=c("A","G","T"), R=c("A","G"), A=c("A"), C=c("C"), H=c("A","C","T"))
real_drach_seqs <- expand.grid(valid_bases) %>% apply(1, paste0, collapse = "")
dna_real <- DNAStringSet(real_drach_seqs)

final_gr <- annotate_drach_mutations(near_drach_gr)

final_gr_complete <- add_query_context(final_gr)

result_df_final <- as.data.frame(final_gr_complete) %>%
  dplyr::select(
    seqnames,
    id,                    # 序列ID
    strand,                # 基因组方向
    
    # --- 序列信息 ---
    motif_seq,             # 当前序列 (如 ATACA)
    target_drach,          # 目标序列 (如 AAACA)
    mutation_code,         # 突变类型 (如 T>A)
    
    # --- 基因组坐标 (h38) ---
    hg38_mut_pos,          # 突变位点的基因组坐标
    hg38_A_pos,            # 变成A的位置的基因组坐标
    
    # --- Query (转录本) 坐标 ---
    query_mut_pos,         # 突变位点在拼接序列中的第几个碱基
    query_A_pos            # A 在拼接序列中的第几个碱基
  )

m6a_new_2 <- GRanges(seqnames = result_df_final$seqnames,
                     IRanges(start = result_df_final$hg38_A_pos,
                             width = 1),
                     strand = result_df_final$strand,
                     relative_pos = result_df_final$query_A_pos,
                     query_id = result_df_final$id,
                     mut_pos_genome = result_df_final$hg38_mut_pos,
                     mut_relative_pos_circRNA = result_df_final$query_mut_pos,
                     mut_info = result_df_final$mutation_code)

m6a_new_2$check <- paste0(m6a_new_2$query_id,':',m6a_new_2$relative_pos)
indx <- setdiff(m6a_new_2$check,mode_1_result$check)
match <- match(m6a_new_2$check,indx)
indx2 <- which(match != 'NA')

mode_3_result <- m6a_new_2[indx2]

