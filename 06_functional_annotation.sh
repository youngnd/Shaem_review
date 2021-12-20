#Create final protein file
gffread -g ../163_wgene_hypo_genome_ShV3.softmasked.fasta nr_9431_14700.gff3 -y nr_9431_14700.pep

#BUSCO
busco -m protein -c 24 -i nr_9431_14700.pep -o BUSCO_nr_9431_14700 -l metazoa_odb10 --update-data 2>BUSCO_nr_9431_14700.err 1>BUSCO_nr_9431_14700.log

#EggNog
conda create -n eggnog python=2.7
wget https://github.com/eggnogdb/eggnog-mapper/archive/1.0.3.tar.gz
tar -xzvf 1.0.3.tar.gz
./download_eggnog_data.py -y euk 1>dl_en.log 2>dl_en.err &
python ../bin/eggnog-mapper-1.0.3/emapper.py -i nr_9431_14700.pep --output EN_nr_9431_14700 -d euk --usemem --cpu 24 -m diamond 2>EN_nr_9431_14700.err 1>EN_nr_9431_14700.out &


#Run InterProScan
module load interproscan/5.44-79.0
module load python3/3.7.7
interproscan.sh -i nr_9431_14700.pep -f TSV,XML --goterms -dp -iprlookup -exclappl MobiDBLite,Coils,Hamap -cpu 72 1> nr_9431_14700_ipro.log 2> nr_9431_14700_ipro.err &
#done

#annotation of the gff
#Prepare database first
#NDY: modified the file:  ~/miniconda3/envs/funannotate/lib/python2.7/site-packages/funannotate/annotate.py so that we keep proteins with 30% protein id and 30% length. initial cutoff was 60% for both and that was only 400 proteins
funannotate setup -i all -d /home/andreas/ShV3/FUN/db --update
funannotate setup -b metazoa -d /home/andreas/ShV3/FUN/db
funannotate annotate -d /home/andreas/ShV3/FUN/db --gff nr_9431_14700.gff3 --fasta 163_wgene_hypo_genome_ShV3.softmasked.fasta -s "Schistosoma haematobium" -o anno_res --busco_db metazoa --iprscan nr_9431_14700.pep.xml --eggnog EN_nr_9431_14700.emapper.annotations --cpus 48 --sbt ../NCBI/template.sbt 2>anno.err 1>anno.log &

#Fix some product names
#All manually fixed product names ok, now run tbl2asn to create sqn for submission
~/bin/linux64.tbl2asn -p . -t ../../../NCBI/template.sbt -M n -Z discrepancy_manual_fix -a r10k -l paired-ends -n "Schistosoma haematobium"
# All fixed
# Final submission on 28.05.2021
# SUB9749424
# /home/andreas/ShV3/FUN/anno_res/annotate_results/Schistosoma_haematobium_v3_28_05_21_submit.sqn

#04.06.2021
#NCBI errors:

# [1] You are using the wrong Locus Tag Prefix 'Sha'. 
# The registered LTP for this BioProject/BioSample pair 
# is MS3. If you are not tracking the annotation from 
# the previous version 2 (which is using 7 digits IDs 
# after the underscore), make sure to have 8 digits 
# identifiers after the underscore . For example, like 
# this:

#   locus_tag MS3_00000001
#   locus_tag MS3_00000002
#   ...


# [2] 2 coding regions have mismatching mRNA. Make sure the 
# product names are identical between the mRNA and the 
# corresponding CDS features.

# CDS ribosome biogenesis protein tsr1  [lcl|HiC_scaffold_2:61617472-61617699, 61618604-61618731, 61619607-61620192, 61621590-61621766, 61622880-61623068, 61623279-61623608] Sha_04077
# CDS CH  [lcl|HiC_scaffold_2:78889001-78889061, 78894160-78894221, 78895738-78895833, 78900446-78900529, 78902107-78902199, 78905782-78905897, 78907841-78907861, 78909850-78910046, 78910760-78910857, 78914133-78914249, 78917511-78917915, 78918874-78919245, 78921641-78921907, 78923760-78923936, 78934980-78935222, 78942006-78942296, 78945348-78945506, 78947617-78947796, 78949846-78949977, 78961118-78961366, 78968261-78968413, 78972838-78973230, 78974597-78975112, 78978202-78982707, 78983756-78983773, 78985818-78986441, 78987708-78987860, 78991789-78991953, 78995087-78995398, 78999272-78999496, 79000292-79000435, 79001633-79002127, 79003542-79003945, 79005079-79005223, 79007244-79007417, 79009049-79009078, 79014006-79014209, 79019816-79019985, 79022092-79022463, 79024175-79024450, 79026041-79026077]  Sha_10501

#Change Sha_ IDs to MS3
sed -i 's/Sha_/MS3_000/g' man_Schistosoma_haematobium.tbl
grep locus_tag bkp_man_Schistosoma_haematobium.tbl | cut -f5 | awk '{print $0"\t"$0}' | sed 's/Sha_/MS3_000/1' > ../gene_Sha_MS3_conversion.txt

#Fix two product names
#renaming protein sequences
fasta_formatter < Schistosoma_haematobium.proteins.fa | awk '{print $1}'  | awk -F"\t" 'BEGIN{while(getline < "gene_Sha_MS3_conversion.txt"){ms3[$1]=$2}}{for(i in ms3){gsub(ms3[i],i)}; print $0}' > final_07_06_21_prot.fasta &
awk -F"\t" '$3 == "gene"' Schistosoma_haematobium.gff3 | cut -f1 | sort | uniq -c > gene_density_scaffolds.tsv
sed 's/;$//' Schistosoma_haematobium.gff3 | sed 's/Sha_/MS3_000/g' | grep -v "#" | gt gff3 -sort -tidy -retainids - | sed 's/\tfunannotate\t/\tShae.V3\t/' > Shae.V3_final.gff3

#Run eggnog mapper (for suppl table)
conda activate py27
python ~/bin/eggnog-mapper-1.0.3/emapper.py -i final_07_06_21_prot.fasta --output EggNog_final/ -d euk --usemem --cpu 48 -m diamond 1>en.log 2>en.err &
cd EggNog_final
mv .emapper.annotations eggnog.emapper.annotations
mv .emapper.seed_orthologs eggnog.emapper.seed_orthologs
awk -F"\t" 'BEGIN{getline; getline; getline}{if(NF == 1){}else{print $0}}' eggnog.emapper.annotations | cut -f1-7,13 > 14700_eggnog.tsv

#kegg annotation
diamond blastp -q final_07_06_21_prot.fasta --db kegg_euk.dmnd -p 48 --out kegg_blast.blast6 --outfmt 6 --evalue 1e-8 1>kegg.log 2>kegg.err &

#run keggAnno
keggAnno -i kegg_blast.blast6 -d keggAnno -r 0 > keggAnno.log 2>&1 &
grep -v "^##" kegg_blast.blast6.kegg.pathway.details | cut -f1,7 | sed 's/[[:digit:]]\{5\} //g' | awk -F"\t" '{split($2,levels,/;/);if(levels[1]!="Human Diseases"){for(i in levels){print $1"\tlevel"i"\t"levels[i];}}}' | sort | uniq | awk -F"\t" '{hit[$1]=1; if(out[$1";"$2]!=""){out[$1";"$2]=out[$1";"$2]";"$3}else{out[$1";"$2]=$3}}END{for(i in hit){print i"\t"out[i";level1"]"\t"out[i";level2"]"\t"out[i";level3"]}}' | sed 's/ - yeast//g' | sed 's/ - fly//g' | sed 's/ - plant//g' > Sh_kegg_pathway_anno.txt
grep -v "^##" kegg_blast.blast6.kegg.pathway.details  | cut -f1,4,5 | sed 's/ \[EC/\t/g' | cut -f1-3 | sed "s/'\+$//" | sort | uniq > KO_sh.tsv
awk -F"\t" 'BEGIN{while(getline < "KO_sh.tsv"){ko[$1]=$0}}{print ko[$1]"\t"$0}' Sh_kegg_pathway_anno.txt | cut -f1-3,5- > KEGG_anno_Sh.tsv
cut -f1,3,4 kegg_blast.blast6.kegg.brite.details | grep -v "##" | sed 's/ \[EC/\t/g' | cut -f1-3 | sed "s/'\+$//" | sort | uniq > allKO.tsv
awk -F"\t" 'BEGIN{while(getline < "KEGG_anno_Sh.tsv"){a[$1""$2]=$0}}{k=$1""$2; if(k in a){print a[$1""$2]}else{print $0}}' allKO.tsv > KEGG_anno_Sh_KO.tsv
sed 's/ - [^-;\t]\+//g' KEGG_anno_Sh_KO.tsv.cleancol > cleanedPWs_KEGG_anno_Sh_KO.tsv.cleancol
cut -f1,6 cleanedPWs_KEGG_anno_Sh_KO.tsv.cleancol | awk -F"\t" '{split($2,a,";"); for(i in a){print $1"\t"a[i]}}' | sed 's/\./\t/1' | cut -f1,3| sort | uniq > all_gene_pathways.tsv
#3203 genes with PW annotation
awk -F"\t" '{printf $1"\t"$2"\t"$3"\t"; split($4,fo,";"); for(i in fo){fou[fo[i]]=fo[i]}; for(k in fou){printf fou[k]";"}; printf "\t"; split($5,fi,";"); for(l in fi){fiu[fi[l]]=fi[l]}; for(m in fiu){printf fiu[m]";"}; printf "\t"; split($6,si,";"); for(n in si){siu[si[n]]=si[n]}; for(o in siu){printf siu[o]";"}; print ""; delete fo; delete fi; delete si; delete fou; delete fiu; delete siu;}' cleanedPWs_KEGG_anno_Sh_KO.tsv.cleancol | sed 's/;\t/\t/g; s/;$//' > UniqCleanedPWs_KEGG_anno_Sh_KO.tsv.cleancol
grep "KO|09180 Brite Hierarchies"  kegg_blast.blast6.kegg.brite.details | cut -f1,7 | sed 's/|/\t/g' | cut -f1,5 | sed 's/^ \+//' | sed 's/ /\t/1' | sed 's/ \[/\t/' | cut -f1,3 > KEGG_BRITE.tsv
awk -F"\t" '{if($6 == "Enzymes"){print $1"\t"$7}}' kegg_blast.blast6.kegg.brite.details | awk -F"\t" '$2 != ""' | sed 's/|/\t/g' | sed 's/^ \+//' | cut -f1,3 | sed 's/ /\t/g' | cut -f1,3 > Enzymes_KEGG.tsv
awk -F"\t" '{a[$1]=a[$1]", "$2}END{for(i in a){print i"\t"a[i]}}' KEGG_BRITE.tsv | sed 's/\t, /\t/1' > tmp; mv tmp KEGG_BRITE.tsv
awk -F"\t" '{a[$1]=a[$1]", "$2}END{for(i in a){print i"\t"a[i]}}' Enzymes_KEGG.tsv | sed 's/\t, /\t/1' > tmp; mv tmp Enzymes_KEGG.tsv
awk -F"\t" 'BEGIN{while(getline < "KEGG_BRITE.tsv"){a[$1]=$2}; while(getline < "Enzymes_KEGG.tsv"){b[$1]=$2}}{print $0"\t"a[$1]"\t"b[$1]}' UniqCleanedPWs_KEGG_anno_Sh_KO.tsv.cleancol > final_KEGG_BRITE_ENZ.tsv
awk -F"\t" '{printf $0; split($7, a, ", "); printf "\t"a[1]; split($8, b, ", "); print "\t"b[1]}' final_KEGG_BRITE_ENZ.tsv | cut -f1-6,9,10 > tmp; mv tmp final_KEGG_BRITE_ENZ.tsv

#total of annotated proteins
cat EggNog_final/14700_eggnog.tsv KEGG_final/keggAnno/final_KEGG_BRITE_ENZ.tsv Shae.V3_IPRO.tsv.cleancol | cut -f1 | grep "^MS3" | sort | uniq | wc -l

#version 2.2.6
#Make sure fasta files have .fasta ending
orthofinder -t 24 -S diamond -n 3specV3 -f sman_sjap_sbov_shv3/ 2>3specV3.err 1>3specV3.out &
cut -f2- Orthogroups.tsv | sed '/^$/d' | grep -v "final_07_06_21_prot" | sed 's/, /|/g' | awk -F"\t" '{split($3,Sha,"|"); for(i in Sha){printf Sha[i]"\t"; for(k in Sha){if(i == k){}else{printf Sha[k]"|"}}print "\t"$1"\t"$2"\t"$4}}' | sed 's/|\t/\t/' | sed 's/|/,/g' | awk -F"\t" '{print $0; orth[$1]=$1}END{while(getline < "ShV3.ids"){if($0 in orth){}else{print $0}}}' | sort -n | sed '/^$/d' > orthoShae.V3.tsv
grep -e "DC041_" -e "EWB00_" -e "Smp_" Orthogroups.tsv | grep "MS3" > 8370_groups_MS3.tsv

awk -F"\t" '{if($3 != ""){print $0}}' orthoShae.V3.tsv.cleancol | cut -f1 | sed 's/\./\t/' |cut -f1 | sort | uniq | wc -l
#8462
awk -F"\t" '{if($4 != ""){print $0}}' orthoShae.V3.tsv.cleancol | cut -f1 | sed 's/\./\t/' |cut -f1 | sort | uniq | wc -l
#7953
awk -F"\t" '{if($5 != ""){print $0}}' orthoShae.V3.tsv.cleancol | cut -f1 | sed 's/\./\t/' |cut -f1 | sort | uniq | wc -l
#9246