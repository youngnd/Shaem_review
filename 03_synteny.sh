#circos_final.sh
#AJS 21.06.2021

#Create list of longest isoforms for new gene set
fa2tab ../../../updated/final_07_06_21_prot.fasta | awk -F"\t" '{split($1,seq,"."); print $0"\t"seq[1] }'  | awk -F"\t" '{if(length($2) > length(longest[$3])){longest[$3]=$2; longid[$3]=$1}}END{for(i in longid){print longid[i]; print longest[i]}}' > longest_iso_shae.fasta

#Run orthofinder to create SCOs
conda activate orthof
orthofinder -t 48 -S diamond -n LI_f -f long_iso_final 2>longiso_f.err 1>longiso_f.out &


ls | while read line; do sed 's/, /,/g' $line | cut -f2- | grep -v longest_iso | awk -F"\t" '{split($1,a,","); split($2,b,","); if(length(a)==1 && length(b)==1){print $0}}' > SCOs_"$line"; done

#create coordinate files for linked SCOs
cat schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3 Shae.V3_final.gff3 | awk -F"\t" '$3 == "mRNA"' | cut -f1,4,5,9 | sed 's/;/\t/g' | cut -f1-4 | sed 's/transcript://' | sed 's/ID=//' | awk -F"\t" 'BEGIN{while(getline < "SCOs_longest_iso_shae__v__longest_iso_sman.tsv"){sha[$1]=$2}}{coord[$4]=$1"\t"$2"\t"$3}END{for(i in sha){print coord[i]"\t"coord[sha[i]]"\t"i"\t"sha[i]}}' | sort | cut -f1-6 > links_SCOs_sman_shae.tsv
cat schistosoma_bovis.PRJNA451066.WBPS15.annotations.gff3 Shae.V3_final.gff3| awk -F"\t" '$3 == "mRNA"' | cut -f1,4,5,9 | sed 's/;/\t/g' | cut -f1-4 | sed 's/transcript://' | sed 's/ID=//' | awk -F"\t" 'BEGIN{while(getline < "SCOs_longest_iso_shae__v__longest_iso_sbovis.tsv"){sha[$1]=$2}}{coord[$4]=$1"\t"$2"\t"$3}END{for(i in sha){print coord[i]"\t"coord[sha[i]]"\t"i"\t"sha[i]}}' | sort | cut -f1-6 > links_SCOs_sbov_shae.tsv
cat schistosoma_japonicum.PRJNA520774.WBPS15.annotations.gff3 Shae.V3_final.gff3| awk -F"\t" '$3 == "mRNA"' | cut -f1,4,5,9 | sed 's/;/\t/g' | cut -f1-4 | sed 's/transcript://' | sed 's/ID=//' | awk -F"\t" 'BEGIN{while(getline < "SCOs_longest_iso_shae__v__longest_iso_sjap.tsv"){sha[$1]=$2}}{coord[$4]=$1"\t"$2"\t"$3}END{for(i in sha){print coord[i]"\t"coord[sha[i]]"\t"i"\t"sha[i]}}' | sort | cut -f1-6 > links_SCOs_sjap_shae.tsv
sed -i 's/transcript_//' SCOs_longest_iso_shae__v__longest_iso_shae9314.tsv 
cat schistosoma_haematobium.PRJNA78265.WBPS15.annotations.gff3 Shae.V3_final.gff3 | awk -F"\t" '$3 == "mRNA"' | cut -f1,4,5,9 | sed 's/;/\t/g' | cut -f1-4 | sed 's/transcript://' | sed 's/ID=//' | awk -F"\t" 'BEGIN{while(getline < "SCOs_longest_iso_shae__v__longest_iso_shae9314.tsv"){sha[$1]=$2}}{coord[$4]=$1"\t"$2"\t"$3}END{for(i in sha){print coord[i]"\t"coord[sha[i]]"\t"i"\t"sha[i]}}' | sort | cut -f1-6 > links_SCOs_V2_shae.tsv

source /home/andreas/miniconda2/etc/profile.d/conda.sh

#define circos parameters
bun_mem=5
gap=1e6
bun_size=1e4

#bundle links
dirCircos=/home/andreas/bin/circos-tools-0.23
conda activate circos
ls links_* | while read line; do $dirCircos/tools/bundlelinks/bin/bundlelinks -min_bundle_membership $bun_mem -max_gap $gap -min_bundle_size $bun_size -links $line > bundled_"$line"; done

#define karyotype colors depending on linked chromosomes
#ls bundled_links* | while read line; do id=$(echo $line | awk -F"_" '{print $3"_"$4}'); awk '{print $1"\t"$4}' $line | awk '{count[$2]++;cline[$0]++;line[$0]=$2}END{for(i in line){print cline[i]" "i" "count[line[i]]}}' | sed 's/ /\t/g' | awk -F"\t" '{print $2"\t"$3"\t"$1/$4}' | sort -k2,2 -k3,3nr | awk -F"\t" 'BEGIN{while(getline < "col_scheme_circos_shv3.tsv"){dark[$1]=$2; light[$1]=$3}}{if(!($2 in found)){found[$2]=$2; print $1"\t"dark[$1]"skip\t"$2"\t"light[$1]"\t0"}else{print $1"\t"dark[$1]"skip\t"$2"\t"dark[$1]"\t5"}}' | awk -F"\t" '$5 == 0' | sed 's/skip/\n/' | sed 's/^\t\+//' | cut -f1,2  > kary_cols_"$line"; done
ls bundled_links_SCOs_* | while read line; do awk 'OFS="\t"{print $1,$4,$6-$5}' $line | awk -F"\t" '{occ[$1"\t"$2]=occ[$1"\t"$2]+$3}END{for(i in occ){print i"\t"occ[i]}}' | sort -t$'\t' -k2,2 -k3,3nr | awk -F"\t" 'BEGIN{while(getline < "col_scheme_circos_shv3.tsv"){dark[$1]=$2; light[$1]=$3}}{if(!($2 in found)){found[$2]=$2; print $1"\t"dark[$1]"skip\t"$2"\t"light[$1]"\t0"}else{print $1"\t"dark[$1]"skip\t"$2"\t"dark[$1]"\t5"}}' | awk -F"\t" '$5 == 0' | sed 's/skip/\n/' | sed 's/^\t\+//' | cut -f1,2  > kary_cols_"$line"; done

#Create karyotype files
cat 163_wgene_hypo_genome_ShV3.softmasked.fasta schistosoma_bovis.PRJNA451066.WBPS15.genomic_softmasked.fa | fa2len | sed 's/^>//' | awk '{print $1" "$NF}' | awk 'BEGIN{while(getline < "bundled_links_SCOs_sbov_shae.tsv"){id[$1]=$1; id[$4]=$4}}{if($1 in id){print $0}}' | awk 'BEGIN{while(getline < "kary_cols_bundled_links_SCOs_sbov_shae.tsv"){col[$1]=$2}}{print "chr - "$1" "$1" 0 "$2" "col[$1]}' > kary_sbov_shae.txt
awk '{print $1"\t"$4}' bundled_links_SCOs_sbov_shae.tsv | sed 's/\t/\n/g' | sort | uniq | while read line; do grep -w $line kary_sbov_shae.txt; done > tmp; mv tmp kary_sbov_shae.txt

cat 163_wgene_hypo_genome_ShV3.softmasked.fasta schistosoma_japonicum.PRJNA520774.WBPS15.genomic_softmasked.fa | fa2len | sed 's/^>//' | awk '{print $1" "$NF}' | awk 'BEGIN{while(getline < "bundled_links_SCOs_sjap_shae.tsv"){id[$1]=$1; id[$4]=$4}}{if($1 in id){print $0}}' | awk 'BEGIN{while(getline < "kary_cols_bundled_links_SCOs_sjap_shae.tsv"){col[$1]=$2}}{print "chr - "$1" "$1" 0 "$2" "col[$1]}' > kary_sjap_shae.txt
awk '{print $1"\t"$4}' bundled_links_SCOs_sjap_shae.tsv | sed 's/\t/\n/g' | sort | uniq | while read line; do grep -w $line kary_sjap_shae.txt; done > tmp; mv tmp kary_sjap_shae.txt

cat 163_wgene_hypo_genome_ShV3.softmasked.fasta schistosoma_mansoni.PRJEA36577.WBPS15.genomic_softmasked.fa | fa2len | sed 's/^>//' | awk '{print $1" "$NF}' | awk 'BEGIN{while(getline < "bundled_links_SCOs_sman_shae.tsv"){id[$1]=$1; id[$4]=$4}}{if($1 in id){print $0}}' | awk 'BEGIN{while(getline < "kary_cols_bundled_links_SCOs_sman_shae.tsv"){col[$1]=$2}}{print "chr - "$1" "$1" 0 "$2" "col[$1]}' > kary_sman_shae.txt
awk '{print $1"\t"$4}' bundled_links_SCOs_sman_shae.tsv | sed 's/\t/\n/g' | sort | uniq | while read line; do grep -w $line kary_sman_shae.txt; done > tmp; mv tmp kary_sman_shae.txt

cat 163_wgene_hypo_genome_ShV3.softmasked.fasta schistosoma_haematobium.PRJNA78265.WBPS15.genomic_softmasked.fa | fa2len | sed 's/^>//' | awk '{print $1" "$NF}' | awk 'BEGIN{while(getline < "bundled_links_SCOs_V2_shae.tsv"){id[$1]=$1; id[$4]=$4}}{if($1 in id){print $0}}' | awk 'BEGIN{while(getline < "kary_cols_bundled_links_SCOs_V2_shae.tsv"){col[$1]=$2}}{print "chr - "$1" "$1" 0 "$2" "col[$1]}' > kary_V2_shae.txt
awk '{print $1"\t"$4}' bundled_links_SCOs_V2_shae.tsv | sed 's/\t/\n/g' | sort | uniq | while read line; do grep -w $line kary_V2_shae.txt; done > tmp; mv tmp kary_V2_shae.txt


#order based on the way you color, just order by HiC chr order, then spacing as usual
ls bundled_links_SCOs_* | while read line; do id=$(echo $line | awk -F"_" '{print $4}'); awk 'OFS="\t"{print $1,$3,$4,$6-$5}' $line | awk -F"\t" '{occ[$1"\t"$3]=occ[$1"\t"$3]+$4; if($2 > long[$1"\t"$3]){long[$1"\t"$3]=$2}}END{for(i in occ){print i"\t"occ[i]"\t"long[i]}}' | sort -t$'\t' -k2,2 -k3,3nr | awk -F"\t" 'OFS="\t"{if(!($2 in found)){found[$2]=$2; print $1,$2,$4}}' | sort -t$'\t' -k1,1 -k3,3nr | awk -F"\t" 'BEGIN{while(getline < "ShV3_fixed_order.txt"){ref[$0]=++cnt}}{scaff[$1]=scaff[$1]","$2}END{printf "chromosomes_order = "; for(i in scaff){print i""scaff[i]"\t"ref[i]}}' | sort -t$'\t' -k2,2n | cut -f1 | sed ':a;N;$!ba;s/\n/,/g' > chr.order_"$id".txt; done

#create spacing file to leave bigger distance between chromosome and associated linked scaffolds
ls chr.order_*.txt | while read line; do id=$(echo $line | awk -F"[_.]" '{print $3}'); sed 's/chromosomes_order = //' $line | sed 's/,HiC/\nHiC/g' | awk -F"," '{if(NR == 1){first=$1; last=$NF}else{print $1"\t"last; last=$NF}; }END{print first"\t"last}' > spacing_"$id".tsv; done

#coloring by size not number of connections
#size-based
ls bundled_links_SCOs_* | while read line; do id=$(echo $line | awk -F"_" '{print $4}'); awk 'OFS="\t"{print $1,$4,$6-$5}' $line | awk -F"\t" '{occ[$1"\t"$2]=occ[$1"\t"$2]+$3}END{for(i in occ){print i"\t"occ[i]}}' | sort -t$'\t' -k2,2 -k3,3nr | awk -F"\t" 'BEGIN{while(getline < "col_scheme_circos_shv3.tsv"){dark[$1]=$2; light[$1]=$3}}{if(!($2 in found)){found[$2]=$2; print $1"\t"$2"\t"light[$1]"\t0"}else{print $1"\t"$2"\t"dark[$1]"\t5"}}' | while read chr1 chr2 col z; do echo "<rule>"; echo "condition = between("$chr1","$chr2")"; echo "color = "$col""; echo "z = "$z""; echo "thickness  = 5"; echo "</rule>"; done > col_rules_"$id".txt; done

#Create final circos configuration files with coloring
awk 'BEGIN{while(getline < "col_rules_sbov.txt"){c++; a[c]=$0}}{if($0 == "</rules>"){for(i=1; i <= length(a); i++){print a[i]}; print $0}else{print $0}}' colored_links_template.conf | cat chr.order_sbov.txt - | awk 'BEGIN{while(getline < "spacing_sbov.tsv"){link[$1]=$2}}{if($0 == "<spacing>"){print $0; getline; print $0; for(i in link){print "<pairwise "i" "link[i]">"; print "    spacing = 6r"; print "  </pairwise>"}}else{print $0}}' | sed 's/^KARYOTYPE$/karyotype = kary_sbov_shae.txt/; s/^BUNDLES$/file =  bundled_links_SCOs_sbov_shae.tsv/' > colored_links_sbov.conf
awk 'BEGIN{while(getline < "col_rules_sman.txt"){c++; a[c]=$0}}{if($0 == "</rules>"){for(i=1; i <= length(a); i++){print a[i]}; print $0}else{print $0}}' colored_links_template.conf | cat chr.order_sman.txt - | awk 'BEGIN{while(getline < "spacing_sman.tsv"){link[$1]=$2}}{if($0 == "<spacing>"){print $0; getline; print $0; for(i in link){print "<pairwise "i" "link[i]">"; print "    spacing = 6r"; print "  </pairwise>"}}else{print $0}}' | sed 's/^KARYOTYPE$/karyotype = kary_sman_shae.txt/; s/^BUNDLES$/file =  bundled_links_SCOs_sman_shae.tsv/' > colored_links_sman.conf
awk 'BEGIN{while(getline < "col_rules_sjap.txt"){c++; a[c]=$0}}{if($0 == "</rules>"){for(i=1; i <= length(a); i++){print a[i]}; print $0}else{print $0}}' colored_links_template.conf | cat chr.order_sjap.txt - | awk 'BEGIN{while(getline < "spacing_sjap.tsv"){link[$1]=$2}}{if($0 == "<spacing>"){print $0; getline; print $0; for(i in link){print "<pairwise "i" "link[i]">"; print "    spacing = 6r"; print "  </pairwise>"}}else{print $0}}' | sed 's/^KARYOTYPE$/karyotype = kary_sjap_shae.txt/; s/^BUNDLES$/file =  bundled_links_SCOs_sjap_shae.tsv/' > colored_links_sjap.conf
awk 'BEGIN{while(getline < "col_rules_V2.txt"){c++; a[c]=$0}}{if($0 == "</rules>"){for(i=1; i <= length(a); i++){print a[i]}; print $0}else{print $0}}' colored_links_template.conf | cat chr.order_V2.txt - | awk 'BEGIN{while(getline < "spacing_V2.tsv"){link[$1]=$2}}{if($0 == "<spacing>"){print $0; getline; print $0; for(i in link){print "<pairwise "i" "link[i]">"; print "    spacing = 6r"; print "  </pairwise>"}}else{print $0}}' | sed 's/^KARYOTYPE$/karyotype = kary_V2_shae.txt/; s/^BUNDLES$/file =  bundled_links_SCOs_V2_shae.tsv/' > colored_links_V2.conf

#Run circos
ls colored_links_* | grep -v "template" | while read line; do id=$(echo $line | awk -F"_" '{print $3}'); circos -configfile $line 2> "$id".err 1> "$id".log; mv circos.png circos_"$id".png; mv circos.svg circos_"$id".svg; done
