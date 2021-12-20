library(tidyverse)
#Set path to file that has unique, fixed SNP positions for samples
# HiC_scaffold_1   302546 SMi  
# HiC_scaffold_1   303409 SMi  
# HiC_scaffold_1   303422 SMi
FILEPATH <- "snp.unique.ann.fixed.positions.txt"
#Karyotype file (to get lengths of chromosomes)
KARY <- "kary_V3_only.txt"
#quantile cut-off
QCO <- 0.80
snp_pos <- read_tsv(FILEPATH, col_names = c("chr","pos","loc"))

#Get length of 8 scaffolds from circos karyotype file
len <- read_delim(KARY, delim = " ", col_names = F) %>% select(3,6)
colnames(len) <- c("chr","len")

#Sum SNP positions over blocks of 1Mb across the chromosomes
blocks_1Mb <- snp_pos %>% 
  inner_join(len, by = "chr") %>% 
  rowwise() %>% 
  mutate(start = as.integer(pos/1e6)*1e6, end = (as.integer(pos/1e6)+1)*1e6) %>% 
  mutate(end=if_else(end > len, len, end)) %>% 
  group_by(chr, loc, start, end) %>% 
  summarise(snpc=n()) %>% 
  ungroup()

#create some summary stats
sum_stats <- blocks_1Mb %>% group_by(chr,loc) %>% summarise(med=median(snpc), mean=mean(snpc), sd=sd(snpc), min=min(snpc), max=max(snpc), quants=quantile(snpc, probs=QCO))
write_tsv(sum_stats, "summary_stats_unique.ann.fixed_counts.tsv", col_names = T)

#Select all 1Mb blocks with density higher than the QCO-quantile
select_high_dens <- 
  blocks_1Mb %>% 
  inner_join(sum_stats, by = c("chr","loc")) %>% 
  select(-starts_with("m"), -sd) %>% 
  mutate(high_dense=if_else(snpc > quants, T, F)) %>% 
  filter(high_dense)

#quantify how much of the whole genome the high dense clusters represent and how much that is in percent of all snps
high_dens_distr <- 
  blocks_1Mb %>% 
  inner_join(sum_stats, by = c("chr","loc")) %>% 
  select(-starts_with("m"), -sd) %>% 
  mutate(high_dense=if_else(snpc > quants, T, F)) %>% 
  group_by(loc, high_dense) %>% summarise(total_snp=sum(snpc))
    
block_cnt <- blocks_1Mb %>% 
    inner_join(sum_stats, by = c("chr","loc")) %>% 
    select(-starts_with("m"), -sd) %>% 
    mutate(high_dense=if_else(snpc > quants, T, F)) %>% 
    group_by(loc, high_dense) %>% count(high_dense)
  
SNP_distr_table <- full_join(high_dens_distr,block_cnt, by = c("loc","high_dense")) %>%
  pivot_wider(names_from = high_dense, values_from = c(total_snp, n)) %>% 
  mutate(total_1Mb=n_FALSE+n_TRUE, allSNP=total_snp_FALSE+total_snp_TRUE) %>% 
  mutate(percHD=n_TRUE/total_1Mb*100, percLD=n_FALSE/total_1Mb*100, SNPpercHD=total_snp_TRUE/allSNP*100, SNPpercLD=total_snp_FALSE/allSNP*100)

#Write bed files
select_high_dens %>% filter(loc == "SMi") %>% select(-loc,-quants,-high_dense) %>% write_tsv(paste0("SMi_high_dens_",QCO,"_1Mb.bed"), col_names = F)
select_high_dens %>% filter(loc == "SS4") %>% select(-loc,-quants,-high_dense) %>% write_tsv(paste0("SS4_high_dens_",QCO,"_1Mb.bed"), col_names = F)
select_high_dens %>% filter(loc == "SMs") %>% select(-loc,-quants,-high_dense) %>% write_tsv(paste0("SMs_high_dens_",QCO,"_1Mb.bed"), col_names = F)
select_high_dens %>% filter(loc == "SZ1") %>% select(-loc,-quants,-high_dense) %>% write_tsv(paste0("SZ1_high_dens_",QCO,"_1Mb.bed"), col_names = F)

#############create SNP file for circos###############
gene_cov <- read_tsv("Shae.V3_final.coverage.txt", col_names = c("chr","start","end", "cov")) 

g_cnt_1Mb <- gene_cov %>% 
  inner_join(len, by = "chr") %>% 
  rowwise() %>% 
  mutate(left_end = end %% 1e6, start = as.integer(start/1e6)*1e6, end = (as.integer(end/1e6)+1)*1e6) %>% 
  mutate(end=if_else(left_end == 0, ((end/1e6)-1)*1e6, end)) %>% 
  mutate(end=if_else(end > len, len, end)) %>%
  select(-left_end, -len) %>% 
  group_by(chr, start, end) %>% 
  summarise(gc=sum(cov)) %>% 
  ungroup()

gc_snpc_tbl_im <- full_join(g_cnt_1Mb, blocks_1Mb, by = c("chr","start","end")) %>% filter(!is.na(loc)) 
dist_loc <- blocks_1Mb %>% distinct(loc,chr)

add_missed <- full_join(g_cnt_1Mb, blocks_1Mb, by = c("chr","start","end")) %>% filter(is.na(loc)) %>% select(-loc) %>% left_join(dist_loc, by = "chr") %>% 
  select(1:4,6,5) %>% replace_na(list(snpc = 0))

gc_snpc_tbl <- rbind(gc_snpc_tbl_im,add_missed) %>% arrange(loc, chr, start)

gc_snpc_tbl %>% write_tsv("gc_snpc_tbl_1Mb.tsv", col_names = F)
