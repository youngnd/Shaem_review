### R script 07b_select_high_dens_SNPs.R

cut -f1-4 gc_snpc_tbl_1Mb.tsv | sort | uniq | sort -k1,1 -k2,2n > 1Mb_gene_cnt.tsv
awk -F"[ \t]" 'BEGIN{while(getline < "kary_V3_only.txt"){col[$3]=$7}}{print $0"\tfill_color="col[$1]",color="col[$1]}' 1Mb_gene_cnt.tsv  > colored_1Mb_gene_cnt.tsv

cut -f5 gc_snpc_tbl_1Mb.tsv | sort | uniq | while read line; do awk -v var=$line -F"\t" '{if($5 == var){print $1"\t"$3"\t"$5"\t"$6}}' gc_snpc_tbl_1Mb.tsv > "$line".density.tsv; done

#run circos
#/home/andreas/DNAZOO_Shae/circos/SNP_NDY/circos_final.conf

### R script 07d_piecharts_SNPs.R