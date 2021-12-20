library(tidyverse)
#install.packages("ExPosition")
library(ExPosition)
#devtools::install_github("jbengler/tidyheatmaps")
library(tidyheatmaps)
tpm_table <- read_tsv("transcription_levels_10_06_21.tsv", col_names = c("rnalib", "trans_id", "tpm")) %>% spread(rnalib,tpm)

#take median for samples with multiple libraries and remove additional mixed adult and 2012 male/female libraries
subs <- tpm_table %>% select(-starts_with("Ad"), -oldFemale, -oldMale) %>% 
  rowwise() %>% 
mutate(medEgg=median(c_across(cols = starts_with("Egg"))),
       medF=median(c_across(cols = starts_with("f"))),
       medM=median(c_across(cols = starts_with("m")))
       ) %>% select(trans_id, medEgg, oldEgg, B1, C1, S1, medM, medF) %>% 
  filter(!if_all(where(is.numeric), ~ . == 0))

#add gene id
gid_subs <- subs %>% ungroup() %>% mutate(gid = str_split(trans_id, "\\.", 2)) %>% 
  unnest(cols = c(gid)) %>% 
  filter(str_detect(gid,"^MS3"))

#median TPM
gid_subs %>% select(-trans_id, -gid) %>%  summarise_all(median)
# # A tibble: 1 x 7
# medEgg oldEgg    B1    C1    S1  medM  medF
# <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#   1   15.0   9.95  5.95  5.74  26.3  9.81  7.79

#top TPM
gid_subs %>% select(-trans_id, -gid) %>%  summarise_all(max)
# # A tibble: 1 x 7
# medEgg oldEgg     B1      C1    S1   medM   medF
# <dbl>  <dbl>  <dbl>   <dbl> <dbl>  <dbl>  <dbl>
#   1 27401.  7451. 10103. 114119. 6054. 25077. 18915.

#average number of isos per gene
gid_subs %>% select(medF, gid) %>% filter(medF > 0.5) %>% count(gid) %>% summarise(mean=mean(n))

#summary counts per stage
gid_subs %>% filter(medEgg > 0.5) %>% count()
#11106
gid_subs %>% filter(medEgg > 0.5) %>% distinct(gid) %>% count()
#8153

#oldEgg
gid_subs %>% filter(oldEgg > 0.5) %>% count()
#9584
gid_subs %>% filter(oldEgg > 0.5) %>% distinct(gid) %>% count()
# 7446

#B1
gid_subs %>% filter(B1 > 0.5) %>% count()
# 7990
gid_subs %>% filter(B1 > 0.5) %>% distinct(gid) %>% count()
# 6506

#C1
gid_subs %>% filter(C1 > 0.5) %>% count()
#9280
gid_subs %>% filter(C1 > 0.5) %>% distinct(gid) %>% count()
#7202

#S1
gid_subs %>% filter(S1 > 0.5) %>% count()
#10657
gid_subs %>% filter(S1 > 0.5) %>% distinct(gid) %>% count()
#7696

#medM
gid_subs %>% filter(medM > 0.5) %>% count()
#11935
gid_subs %>% filter(medM > 0.5) %>% distinct(gid) %>% count()
#8182

#medF
gid_subs %>% filter(medF > 0.5) %>% count()
#11714
gid_subs %>% filter(medF > 0.5) %>% distinct(gid) %>% count()
#8112

#cluster transcription profiles
tpm_m <- rowNorms(as.matrix(column_to_rownames(subs, "trans_id")), type = "z")
dist_tpm <- dist(tpm_m, method="euclidean")
clust_tpm <- hclust(dist_tpm, method="ward.D")
cluster7 <- cutree(clust_tpm, k=7)  
cluster8 <- cutree(clust_tpm, k=8)  
cluster9 <- cutree(clust_tpm, k=9)  
cluster10 <- cutree(clust_tpm, k=10)  
cluster11 <- cutree(clust_tpm, k=11)  
cluster12 <- cutree(clust_tpm, k=12)

clust_subs <- as_tibble(cbind(subs, cluster7, cluster8, cluster9, cluster10, cluster11, cluster12)) %>% arrange(match(cluster7, c(3,2,6,5,1,7,4)), desc(cluster8), desc(cluster9), desc(cluster10), desc(cluster11), desc(cluster12)) %>% select(1:9) %>% rename(cluster = cluster7)

#create heatmap from clusters
clust_subs %>%  mutate(cluster=as.character(cluster)) %>%
  gather("sample", "tpm", 2:8) %>% 
  tidy_heatmap(rows=trans_id,
               cellwidth = 20,
               columns = sample,
               values = tpm,
               scale = "row",
               gaps_row = cluster,
               gaps_col = NULL,
               color_legend_n = 7,
               colors = c("gray30","#ffffff","goldenrod2"),
               annotation_row = cluster,
               show_rownames = F,
               height = 4,
               width = 3,
               filename = "direct_out_heatmap.pdf")

#Create file with cluster ids
clust_subs %>% select(trans_id, cluster) %>% write_tsv("clusters_tpm.tsv", col_names = F)

clust_subs <- 
clust_subs %>% 
mutate(gid = str_split(trans_id, "\\.", 2)) %>% 
unnest(cols = c(gid)) %>% 
filter(str_detect(gid,"^MS3"))

#Fix cluster names (reorder)
clust_lkp <- tibble(cluster = c(
  3,
  2,
  6,
  5,
  1,
  7,
  4
),

cleanclust = c(
  "1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7"
)
)

fix_clust <- clust_subs %>% 
left_join(clust_lkp, by = "cluster") %>%
  select(-cluster) %>% 
  rename(clust = cleanclust)

t_fix_clust <- fix_clust %>% select(1:8) %>% gather("stage", "tpm", 2:8) 

#top 1% transcription
fix_clust %>% summarise(across(2:8, ~quantile(.x, probs=0.99))) %>%
  gather("stage", "q90") %>% 
  right_join(t_fix_clust, by = "stage") %>% 
  filter(tpm > q90) %>% 
  #select(-q90) %>% 
  pivot_wider(names_from = stage, values_from = tpm) %>%
  summarise(across(3:9, ~median(.x, na.rm = T)))

#################################################

#isoforms with microexons
mexon_isos <- read_delim("micro_exon_isos.tsv", col_names = c("m_exon","trans_id"), delim = " ") %>% 
  mutate(m_exon=as.integer(str_trim(m_exon, side = "left"))) %>% 
  select(trans_id, m_exon) 

#number of exons per isoform
exon_isos <- read_delim("exon_isos.tsv", col_names = c("exon","trans_id"), delim = " ") %>% 
  mutate(exon=as.integer(str_trim(exon, side = "left"))) %>% 
  select(trans_id, exon) 

#extract info for two genes with many isoforms and microexons
fix_clust %>% left_join(mexon_isos, by = "trans_id") %>%
  mutate(m_exon = replace_na(m_exon, 0)) %>% 
  left_join(exon_isos, by = "trans_id") %>%
  filter(gid == "MS3_00004678" | gid =="MS3_00008061") %>% 
  write_tsv("MS3_00004678_MS3_00008061_exon_cluster_TPM_info.tsv")

isoforms_diff_clust <- clust_subs %>% distinct(cluster,gid) %>% count(gid) %>% filter(n > 1)
#2648

#number of cluster memberships vs. number of microexons (correlation)
clust_me_corr <- clust_subs %>% distinct(cluster,gid) %>% count(gid) %>% left_join(me_genes, by = "gid") %>% mutate(meCnt=replace_na(meCnt, 0)) %>% rename(numClust = n)
library(ggpubr)
clust_me_corr %>% ggplot(aes(x = numClust, y = meCnt)) + geom_point() + geom_smooth(method=lm) +  stat_cor(method = "pearson")

#awk -F"\t" '$3 == "exon"' Shae.V3_final.gff3 | sed 's/;Parent=/\t/' | cut -f4,5,10  | sed 's/\./\t/' | cut -f1-3 | sort | uniq | cut -f3 | sort | uniq -c | sed 's/ \+//' | awk 'BEGIN{print "gid\texCnt"}{print $2"\t"$1}' > exonCounts.tsv

#number of cluster memberships vs. number of exons (correlation)
exo_genes <- read_tsv("exonCounts.tsv")
clust_exo_corr <- clust_subs %>% distinct(cluster,gid) %>% count(gid) %>% left_join(exo_genes, by = "gid") %>% mutate(exCnt=replace_na(exCnt, 0)) %>% rename(numClust = n)
clust_exo_corr %>% ggplot(aes(x = numClust, y = exCnt)) + geom_point() + geom_smooth(method=lm) +  stat_cor(method = "pearson")