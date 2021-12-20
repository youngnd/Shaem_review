#library(BiocManager)
#install("ballgown")
library(tidyverse)
library(ggsignif)
ipro <- read_tsv("suppl_anno/Shae.V3_IPRO.tsv.cleancol", col_names = c(
  "qseqid",
  "length",
  "db",
  "db_id",
  "db_descr",
  "qstart",
  "qend",
  "evalue",
  "ipro_id",
  "ipro_descr",
  "go"))

tpm_shv3 <- read_tsv("transcription_levels_10_06_21.tsv", col_names = c("rnalib", "trans_id", "tpm")) %>% spread(rnalib,tpm)

library(ballgown)
mf_st <- ballgown(dataDir="male_female/", samplePattern = "*", meas='all')

#create groups for differential expression analysis
pData(mf_st) = data.frame(id=sampleNames(mf_st), group=rep(c(1,0), each=6))

pData(mf_st)

male_female_de = stattest(mf_st, feature='transcript', meas='FPKM', covariate='group', getFC = T)
library(tidyverse)

# If you are doing two-group comparisons, you can use the `getFC` option in stattest to get the estimated fold change,
# which indicates direction. (fold change < 1 indicates downregulation, fold change > 1 indicates upregulation,
# and the reference group is determined by whatever category from your `group` variable the `lm` function in R
# would consider the reference category. Generally this is the category that comes first alphabetically.
# If you used 0/1 labels, the reference category will be category 0, so FC > 1 means upregulated in the 1 group.

#female = 1
#male = 0

diff_expr_f <- male_female_de %>% filter(qval < 0.05) %>% filter(pval < 0.05) %>% filter(fc >= 2)
diff_expr_m <- male_female_de %>% filter(qval < 0.05) %>% filter(pval < 0.05) %>% filter(fc <= 0.5)

mf_up_transcripts <- as_tibble(rbind((diff_expr_f %>% mutate(sex="female")), (diff_expr_m %>% mutate(sex="male"))) %>% select(id, sex, fc, qval))
tid_tname <- as_tibble(cbind(ballgown::expr(mf_st)$trans[,6],ballgown::expr(mf_st)$trans[,1]))
colnames(tid_tname) <- c("tid","id")

# make volcano plot
as_tibble(left_join(male_female_de, tid_tname, by = "id") %>% select(tid, fc, qval)) %>% 
filter(!is.na(qval)) %>% 
  mutate(color = if_else(fc >= 2 & qval < 0.05, "red", "black")) %>% 
  mutate(color = if_else(fc <= 0.5 & qval < 0.05, "blue", color)) %>%
  ggplot(aes(x=log2(fc), y=-log10(qval), color=color)) + geom_point(shape = 20) +
  scale_color_manual(values=c("grey", "blue", "red")) +
  theme_minimal() +
  geom_vline(xintercept=c(log2(0.5), log2(2)), col="black", size = 0.25) +
  geom_hline(yintercept=-log10(0.05), col="black", size = 0.25) +
  xlab("log2(FC)") +
  ylab("-log10(q)") +
  theme(legend.position = "none")
  ggsave("volcano_mf.svg")

male_up_out <- left_join(mf_up_transcripts, tid_tname, by = "id") %>% filter(sex == "male") %>% mutate(fc=1/fc) %>% select(tid,fc,qval)
female_up_out <- left_join(mf_up_transcripts, tid_tname, by = "id") %>% filter(sex == "female") %>% select(tid,fc,qval)
write_tsv(male_up_out, file = "male_up_trans.tsv", col_names = F) 
write_tsv(female_up_out, file = "female_up_trans.tsv", col_names = F)

#suppl tables
read_tsv("male_up_trans.tsv", col_names = c("trans_id", "FC", "qval")) %>%
  inner_join(tpm_shv3, by = "trans_id") %>%
  select(trans_id, FC, qval, starts_with("m")) %>%
  write_tsv("suppl_anno/male_DE_suppl.tsv", col_names = T)

read_tsv("female_up_trans.tsv", col_names = c("trans_id", "FC", "qval")) %>%
  inner_join(tpm_shv3, by = "trans_id") %>%
  select(trans_id, FC, qval, starts_with("f")) %>%
  write_tsv("suppl_anno/female_DE_suppl.tsv", col_names = T)

#run pathway enrichment
#08a_transcription.sh

library(qvalue)
male_PW <- read_tsv("pwMale.res", col_names = T)
male_PW$Qval <- qvalue(male_PW$Pval)$qvalues
sigPW_male <- male_PW %>% filter(Qval < 0.05)

female_PW <- read_tsv("pwFemale.res", col_names = T)
female_PW$Qval <- qvalue(female_PW$Pval)$qvalues
sigPW_female <- female_PW %>% filter(Qval < 0.05)

male_BRITE <- read_tsv("m_BRITE_enrichment.tsv", col_names = c("brite", "descr", "deCnt", "allCnt", "pval"))
male_BRITE$qval <- qvalue(male_BRITE$pval)$qvalues
sigBRITE_male <- male_BRITE %>% filter(qval < 0.05)

female_BRITE <- read_tsv("f_BRITE_enrichment.tsv", col_names = c("brite", "descr", "deCnt", "allCnt", "pval"))
female_BRITE$qval <- qvalue(female_BRITE$pval)$qvalues
sigBRITE_female <- female_BRITE %>% filter(qval < 0.05)

write_tsv(sigPW_female, file = "sigPW_female.tsv", col_names = T)
write_tsv(sigPW_male, file = "sigPW_male.tsv", col_names = T)
write_tsv(sigBRITE_female, file = "sigBRITE_female.tsv", col_names = F)
write_tsv(sigBRITE_male, file = "sigBRITE_male.tsv", col_names = F)

#scaffold enrichment
#See 08a_transcription.sh
all_scaff <- read_tsv("all_IDs_trans_scaffold.tsv", col_names = c("scaff","trans_id"))
m_scaff <- read_tsv("DE_male_scaffold_counts.tsv", col_names = c("male","scaff"))
f_scaff <- read_tsv("DE_female_scaffold_counts.tsv", col_names = c("female","scaff"))

DE_M = 1963
DE_F = 1512

pvalScaff <- full_join(m_scaff,f_scaff) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(mnot=DE_M-male) %>%
  mutate(fnot=DE_F-female) %>% 
  rowwise() %>% 
mutate(pvalm=fisher.test(matrix(c(male,female,mnot,fnot),nrow=2,ncol=2),alternative="greater")$p.value) %>%
  mutate(pvalf=fisher.test(matrix(c(male,female,mnot,fnot),nrow=2,ncol=2),alternative="less")$p.value) %>%
  ungroup() %>% 
mutate(qvalm=qvalue(pvalm)$qvalues) %>% 
  mutate(qvalf=qvalue(pvalf)$qvalues)

filter(pvalScaff,if_any(starts_with("qval"), ~ . < 0.05)) %>% select(starts_with("q"),starts_with("scaff"))

#are any of the genes on scaffold 194 transcribed in  males?
scaff194 <- all_scaff %>% filter(scaff == "HiC_scaffold_194")

tpm_shv3 %>% select(trans_id, starts_with("m"),oldMale) %>% 
  inner_join(scaff194, by = "trans_id") %>% 
  summarise_at(vars(2:7), median) %>% 
  rowwise() %>%  max()

#scaff194 annotation
# all 
inner_join(scaff194, ipro, by = c("trans_id" = "qseqid")) %>%  filter(db == "Pfam") %>% select(trans_id, db_descr, ipro_descr) %>% mutate(gid = str_split(trans_id, "\\.", 2)) %>% 
  unnest(cols = c(gid)) %>% 
  filter(str_detect(gid,"^MS3")) %>%  distinct(gid)
#9

#figure enrichment scaffolds
#Barplot Fig 3c
bar_mf_de_scaff <-   tibble(
  chr = c(
  "Chromosome 1",
  "Chromosome 2",
  "Chromosome 3",
  "Chromosome 4",
  "Chromsome 5",
  "Chromsome 6",
  "Chromosome 7",
  "Chromosome ZW",
  "Scaffold 194"),
  
  fracfem = c(
  24.86772487,
  12.83068783,
  14.88095238,
  11.24338624,
  4.894179894,
  5.886243386,
  3.637566138,
  18.45238095,
  0.595238095
  ),
  
  fracmal = c(
  20.83545593,
  10.39225675,
  9.628120224,
  9.169638309,
  3.565970453,
  6.979113602,
  4.075394804,
  34.13143148,
  0
  )
) %>% 
  gather(sex, perc, fracmal, fracfem) %>% 
  mutate(sex = fct_reorder(sex, sex, .desc = TRUE)) %>% 
    mutate(chr=as_factor(chr))

bar_mf_de_scaff %>% 
  ggplot(aes(x = chr, y = perc, fill=sex)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  geom_signif(
    y_position = c(25.5, 15.5, 34.7, 1.2), xmin = c(0.75, 2.75, 7.75, 8.75), xmax = c(1.25, 3.25, 8.25, 9.25),
    annotation = "*", tip_length = 0.01,
    size = 1, textsize = 6
  ) +
  scale_fill_manual(values = c("blue", "red")) + 
  theme_minimal(base_size = 14) + 
  ylab("Percentage of DE transcripts") + 
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust=0.5, vjust = 0.55, size=14))
ggsave("barplot_de_scaffs.svg", width = 8, height = 5)