#filtering and trimming
fastp -i short_reads_R1.fq -I short_reads_R2.fq -o fastp.short_reads_R1.fq -O fastp.short_reads_R2.fq -m --merged_out fastp_merged_short_reads.fq --unpaired1 fastp_unpaired_SE_short_reads.fq --unpaired2 fastp_unpaired_SE_short_reads.fq 2> fastpShaemShort.log 1> /dev/null

centrifuge -x nt -1 fastp.short_reads_R1.fq -2 fastp.short_reads_R2.fq -U fastp_merged_short_reads.fq,fastp_unpaired_SE_short_reads.fq --report-file fastp.Centrifuge.Shae.report.txt -S fastp.Shae.classification.txt -p 20 > centr.log 2>&1 &
centrifuge-kreport -x nt fastp.Shae.classification.txt > fastp.Shae.kraken > kraken.log 2>&1 &
awk -F"\t" '{if($5=="6157"){print $5; getline; while(!($4 == "P")){print $5; getline;}}}END{print "0"}' fastp.Shae.kraken > platyhelm_unclass_IDs.ls
awk -F"\t" 'BEGIN{while(getline < "platyhelm_unclass_IDs.ls"){plaun[$0]=$0}; getline;}{if($3 in plaun){print $1}}' fastp.Shae.classification.txt | sort | uniq > fastp_centrifuge_filtered_read_ids.ls
awk -F"\t" '{if($5 == "6157" || $5 == "0"){print $1"\t"$2"\t"$5"\t"$6}}' fastp.Shae.kraken > fastp.centrifuge_Shae_kraken_summary.txt
awk -F"\t" '{if($4 == "P" || $4 == "U"){print $0}}' fastp.Shae.kraken | sort -k1,1nr | grep -v "0.00" > fastp.phyla_Shae_centrifuge.txt

#remove ".1" at the end to get both read pairs later
sed 's/\.1$//' fastp_centrifuge_filtered_read_ids.ls > fixed.ids

#Now remove all reads that aren't either unclassified or Platyhelminthes
fq2tab fastp.short_reads_R1.fq | awk 'BEGIN{while(getline < "fastp_centrifuge_filtered_read_ids.ls"){good["@"$0]=$0}}{if($1 in good){print $0}}' | awk -F"\t" '{print $1; print $2; print $3; print $4}' > centr_fastp.short_reads_R1_Shae.fq
fq2tab fastp.short_reads_R2.fq | awk 'BEGIN{while(getline < "fixed.ids"){good["@"$0".1"]="y"; good["@"$0".2"]="y"; good["@"$0]="y"}}{if($1 in good){print $0}}' | awk -F"\t" '{print $1; print $2; print $3; print $4}' > fixed_centr_fastp.short_reads_R2_Shae.fq
mv fixed_centr_fastp.short_reads_R2_Shae.fq centr_fastp.short_reads_R2_Shae.fq 
cat fastp_merged_short_reads.fq fastp_unpaired_SE_short_reads.fq | fq2tab | awk 'BEGIN{while(getline < "fixed.ids"){good["@"$0".1"]="y"; good["@"$0".2"]="y"; good["@"$0]="y"}}{if($1 in good){print $0}}' | awk -F"\t" '{print $1; print $2; print $3; print $4}' > centr_fastp_SE_Shae.fq

#fix the SRR readname error
fq2tab centr_fastp.short_reads_R2_Shae.fq | awk -F"\t" '{if($1 ~ /^@SRR/){sub(/length=[0-9]+$/," 2:N:0",$1); sub(/\.2 /,":",$1)}; for(i=1; i <=4; i++){print $i}}' > SRR_centr_fastp.short_reads_R2_Shae.fq
fq2tab centr_fastp.short_reads_R1_Shae.fq | awk -F"\t" '{if($1 ~ /^@SRR/){sub(/length=[0-9]+$/," 1:N:0",$1); sub(/\.1 /,":",$1)}; for(i=1; i <=4; i++){print $i}}' > SRR_centr_fastp.short_reads_R1_Shae.fq
#check that all reads are properly paired
cat SRR_centr_fastp.short_reads_R1_Shae.fq SRR_centr_fastp.short_reads_R2_Shae.fq | fq2tab | cut -f1 | awk '{print $1}' | sort | uniq -c | awk '{print $1}' | sort | uniq -c

sed 's/  2:N:0$/\/2/' SRR_centr_fastp.short_reads_R2_Shae.fq > 2_SRR_centr_fastp.short_reads_R2_Shae.fq
sed 's/  1:N:0$/\/1/' SRR_centr_fastp.short_reads_R1_Shae.fq > 1_SRR_centr_fastp.short_reads_R1_Shae.fq

grep "@SRR" -A4 SRR_centr_fastp.short_reads_R1_Shae.fq | head -2000 > test_R1.fq
grep "@SRR" -A4 SRR_centr_fastp.short_reads_R2_Shae.fq | head -2000 > test_R2.fq
interleave-reads.py --gzip -o test_il.gz test_R1.fq test_R2.fq

#khmer 3.0.0a3
interleave-reads.py --gzip -o centr_fastp.short_reads_IL_Shae.fq.gz 1_SRR_centr_fastp.short_reads_R1_Shae.fq 2_SRR_centr_fastp.short_reads_R2_Shae.fq 1> il.log 2> il.err
normalize-by-median.py -k 26 -p -C 20 -M 600G -o fastp_centr_K26_khmer_Shae.fq --unpaired-reads centr_fastp_SE_Shae.fq centr_fastp.short_reads_IL_Shae.fq.gz 1>fastp_centr_K26_khmer.out 2>fastp_centr_K26_khmer.err
extract-paired-reads.py -p fastp_centr_K26_khmer_Shae_paired.fq -s fastp_centr_K26_khmer_Shae_single.fq fastp_centr_K26_khmer_Shae.fq

# wrote to: fastp_centr_K26_khmer_Shae_paired.fq and fastp_centr_K26_khmer_Shae_single.fq
seqtk seq fastp_centr_K26_khmer_Shae_paired.fq -1 > fastp_centr_K26_khmer_Shae_R1.fq
seqtk seq fastp_centr_K26_khmer_Shae_paired.fq -2 > fastp_centr_K26_khmer_Shae_R2.fq