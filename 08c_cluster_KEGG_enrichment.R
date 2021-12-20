library(tidyverse)
library(qvalue)

#Read in pathway enrichment info for clusters and set qvalue cutoff
pathway <- 
  read_tsv("pathway_enrichment_all.tsv", col_names = F) %>% 
  select(1:4,7:8) %>% 
  setNames(c("pathId", "pVal", "set", "all", "pw", "sam")) %>% 
  group_by(sam) %>% 
  mutate(qVal = qvalue(pVal)$qvalues) %>% 
  filter(qVal < 0.05) %>% 
  separate(pw, c("lv1","lv2","lv3"), sep = ";")

#fix cluster names
clust_lkp <- tibble(sam = c(
  "pw3",
  "pw2",
  "pw6",
  "pw5",
  "pw1",
  "pw7",
  "pw4"
),

clust = c(
  "1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7"
)
)

#enriched pathway counts
pathway %>% 
  group_by(sam, lv2) %>% 
  summarise(allpw=sum(set)) %>% 
  arrange(desc(allpw)) %>% 
  group_by(sam) %>% 
  top_n(2) %>%
  left_join(clust_lkp, by = "sam") %>% 
  ungroup() %>% 
  select(-sam) %>% 
  arrange(lv2, desc(allpw)) %>% 
  count(lv2)

pathway %>% write_tsv("pathway_enrichment_qval.tsv", col_names = T)

#Q-value cut-off for BRITE enrichment
brite_tpm <- read_tsv("all_BRITE_enrichment.tsv", col_names = F) %>% 
  setNames(c("briteId", "descr", "set", "all", "pVal", "sam")) %>% 
  group_by(sam) %>% 
  mutate(qVal = qvalue(pVal)$qvalues) %>% 
  filter(qVal < 0.05)

brite_tpm %>% 
group_by(sam, briteId, descr) %>% 
  summarise(allpw=sum(set)) %>% 
  arrange(desc(allpw)) %>% 
  group_by(sam) %>% 
  top_n(2) %>%
  mutate(sam=tolower(sam)) %>% 
  left_join(clust_lkp, by = "sam") %>% 
  select(-sam) %>% 
  arrange(descr, desc(allpw))
  
brite_tpm %>% group_by(sam, descr) %>% count() %>% arrange(sam, desc(n)) %>% 
  write_tsv("brite_counts.tsv", col_names = T)

brite_tpm %>% write_tsv("brite_enrichment_tpm_qval.tsv", col_names = T)