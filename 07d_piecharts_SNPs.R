library(tidyverse)
#These are the summarised counts of how many (cnt) of the 500 2kb sections hit which species in the database
SZ1 <- read_tsv("SZ1_high_dens_0.8_1Mb.bed.fa.windows.bed.cov90.result.txt", col_names = c("chr","range","spec","cnt")) %>% separate(range, c("start","end")) %>% mutate(sample="SZ1")
SS4 <- read_tsv("SS4_high_dens_0.8_1Mb.bed.fa.windows.bed.cov90.result.txt", col_names = c("chr","range","spec","cnt")) %>% separate(range, c("start","end")) %>% mutate(sample="SS4")
SMs <- read_tsv("SMs_high_dens_0.8_1Mb.bed.fa.windows.bed.cov90.result.txt", col_names = c("chr","range","spec","cnt")) %>% separate(range, c("start","end")) %>% mutate(sample="SMs")
SMi <- read_tsv("SMi_high_dens_0.8_1Mb.bed.fa.windows.bed.cov90.result.txt", col_names = c("chr","range","spec","cnt")) %>% separate(range, c("start","end")) %>% mutate(sample="SMi")

chr_map <- tibble(
chr = c(
"HiC_scaffold_1",
"HiC_scaffold_4",
"HiC_scaffold_3",
"HiC_scaffold_5",
"HiC_scaffold_7",
"HiC_scaffold_6",
"HiC_scaffold_8",
"HiC_scaffold_2"),

chrn = c(
"Chr_1",
"Chr_2",
"Chr_3",
"Chr_4",
"Chr_5",
"Chr_6",
"Chr_7",
"Chr_ZW"))

color_mapping <- tibble(spec = c(
  "HiC",   
  "SBOVIS",
  "SCUD",  
  "SMRZ",  
  "SMTD"  
),

col = c(
  "#F0F0F0",
  "#A6611A",
  "#DFC27D",
  "#80CDC1",
  "#018571"
)
)

all_sams <- rbind(SZ1, SS4, SMs, SMi) %>% group_by(chr, start, end, sample) %>% mutate(id = cur_group_id()) %>% ungroup() %>% left_join(color_mapping, by = "spec") %>% left_join(chr_map, by = "chr") %>% mutate(chr=chrn) %>% select(-chrn)
totals <- rbind(SZ1, SS4, SMs, SMi) %>% group_by(chr, start, end, sample) %>% mutate(id = cur_group_id()) %>% group_by(id) %>% summarise(total=sum(cnt)) %>% ungroup()
all_sams <- all_sams %>% left_join(totals, by = "id")

gt20 <- all_sams %>% filter(spec == "HiC") %>% filter(cnt/total < 0.8) %>% select(id) %>% mutate(gt20 = T)

all_sams <- all_sams %>% 
  right_join(gt20, by = "id")

#Create table:
pchart_table <- all_sams %>% select(-total, -gt20, -col) %>% 
  spread(spec, cnt) %>% 
  mutate_at(vars(6:10), ~replace_na(.,0)) %>% 
  mutate(chr = str_extract(chr, "[^Chr_]")) %>% 
  select(-id)

pchart_perc <- all_sams %>% select(-gt20, -col, -id) %>% 
  spread(spec, cnt) %>% 
  mutate_at(vars(6:10), ~replace_na(.,0)) %>% 
  mutate(chr = str_extract(chr, "[^Chr_]")) %>% 
  mutate_at(vars(6:10), ~. / total *100) %>% 
  select(chr, start, end, sample, 6:10)

pchart_table %>% left_join(pchart_perc, by = c("chr", "start", "end", "sample")) %>% write_tsv("piechart_table.tsv")


lapply(split(all_sams, all_sams$id), function(d) {
  c <- paste(unique(d$chr), unique(d$start), unique(d$end), unique(d$sample), sep = "_")
  col <- as.character(d$col)
  #names(col) <- as.character(all_sams$col)
  d %>% 
    ggplot(aes(x = "", y=cnt, fill = factor(spec))) + 
    geom_bar(width = 1, stat = "identity") +
    theme(axis.line = element_blank(), 
          plot.title = element_text(hjust=0.5)) + 
    labs(fill="spec", 
         x=NULL, 
         y=NULL) +
        coord_polar(theta = "y", start=0) +
    theme_void() +
    geom_vline(xintercept = 1.5, size = 0.3, color="#CCCCCC") +
    ##CCCCCC" "#969696" 
    scale_fill_manual(values=col)
  ggsave(paste0("piecharts_final/plot_",c,".svg"))
})
