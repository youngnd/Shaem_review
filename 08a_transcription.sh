#Short-read RNA-Seq Mapping
ls *.fq | awk -F"[_.]" '{print $1}' | sort | uniq | while read line; do echo "hisat2 --dta -p 24 -x ShV3 -1 "$line"_R1.fq -2 "$line"_R2.fq -U "$line"_SE.fq --summary-file hisat."$line".summary 2> /dev/null | samtools sort -@ 24 -o "$line".aln.sorted.bam"; done > runAlignments.sh

cat runAlignments.sh
# source /home/andreas/miniconda2/etc/profile.d/conda.sh
# conda activate hisat2
# hisat2 --dta -p 24 -x ShV3 -1 Ad1_R1.fq -2 Ad1_R2.fq -U Ad1_SE.fq --summary-file hisat.Ad1.summary 2> /dev/null | samtools sort -@ 24 -o Ad1.aln.sorted.bam
# hisat2 --dta -p 24 -x ShV3 -1 Ad2_R1.fq -2 Ad2_R2.fq -U Ad2_SE.fq --summary-file hisat.Ad2.summary 2> /dev/null | samtools sort -@ 24 -o Ad2.aln.sorted.bam
# hisat2 --dta -p 24 -x ShV3 -1 Ad3_R1.fq -2 Ad3_R2.fq -U Ad3_SE.fq --summary-file hisat.Ad3.summary 2> /dev/null | samtools sort -@ 24 -o Ad3.aln.sorted.bam
# hisat2 --dta -p 24 -x ShV3 -1 B1_R1.fq -2 B1_R2.fq -U B1_SE.fq --summary-file hisat.B1.summary 2> /dev/null | samtools sort -@ 24 -o B1.aln.sorted.bam
# hisat2 --dta -p 24 -x ShV3 -1 C1_R1.fq -2 C1_R2.fq -U C1_SE.fq --summary-file hisat.C1.summary 2> /dev/null | samtools sort -@ 24 -o C1.aln.sorted.bam
#...

#Create file with StringTie transcription analysis command for each library
ls *.bam | while read line; do spec=$(echo $line | awk -F"." '{print $1}'); mkdir $spec; echo "stringtie $line -eB -G Shae.V3_final.gff3 -p 48 -o "$spec"/"$spec"_stringtie.gtf > /dev/null 2>&1"; done >> runStringTie.sh

cat runStringTie.sh
# source /home/andreas/miniconda3/etc/profile.d/conda.sh
# conda activate stringtie

# stringtie Ad1.aln.sorted.bam -eB -G Shae.V3_final.gff3 -p 48 -o Ad1/Ad1_stringtie.gtf > /dev/null 2>&1
# stringtie Ad2.aln.sorted.bam -eB -G Shae.V3_final.gff3 -p 48 -o Ad2/Ad2_stringtie.gtf > /dev/null 2>&1
#...

#stringTie TPM transcript levels
find . -iname "*.gtf" | while read line; do sam=$(echo $line | awk -F"/" '{print $2}'); sed 's/;$//' $line | awk -v sample=$sam -F"\t" '{if($3 == "transcript"){split($9,trans,";"); print sample"\t"trans[2]"\t"trans[5]}}' | sed 's/ transcript_id "//; s/ TPM "//; s/"//g'; done > transcription_levels.tsv

#for cluster transcript analysis
#find all genes that contain micro exons (<= 54bp)
awk -F"\t" '$3 == "exon" && $5-$4 <= 54' Shae.V3_final.gff3 | sed 's/=/\t/g' | cut -f11 | sed 's/\./\t/' | cut -f1 | sort | uniq -c | sed 's/^ \+//' | sed 's/ /\t/g' | awk -F"\t" '{print $2"\t"$1}' | sort -k2,2nr |  awk -F"\t" 'BEGIN{print "gid\tmeCnt"}{print $0}' > micro_exon_genes.tsv

#isoform microecxon and exon counts
awk -F"\t" '$3 == "exon" && $5-$4 <= 54' Shae.V3_final.gff3 | sed 's/=/\t/g' | cut -f11 | sort | uniq -c > micro_exon_isos.tsv 
awk -F"\t" '$3 == "exon"' Shae.V3_final.gff3 | sed 's/=/\t/g' | cut -f11 | sort | uniq -c > exon_isos.tsv 

##R script 08b_clustering_TPM_ShV3.R

#Calculate pathway enrichment for individual clusters
cut -f2 clusters_tpm.tsv | sort | uniq | while read line; do awk -F"\t" -v clust=$line '{if($2 == clust){print $1"\t2\t0.05"}}' clusters_tpm.tsv > "$line"_clust.tsv; python /home/pakorhon/Codebase/diffExpKeggPwEnrichment.py -d KeggPw"$line" -i "$line"_clust.tsv -r 0.05 -b ../kegg_blast.blast6 -s + | sort -k 2,2 -gr | grep -v "Human Dis" > pw"$line".res; done 1> tpm_kegg.out 2> tpm_kegg.err &
ls pw*.res | while read line; do grep "^path:ko" $line | awk -F"\t" -v sam=$line '{print $0"\t"sam}' >> pathway_enrichment_all.tsv; done
sed -i 's/.res$//' pathway_enrichment_all.tsv

#automate BRITE processing
ls | grep KeggPw | while read line; do python /home/pakorhon/Codebase/diffExpKeggBrEnrichment.py -i "$line"/kegg_blast.blast6.kegg.brite.freqs -r "$line"/kegg_blast.blast6.kegg.brite.hierarchy.with.gene.ids -s "$line"/-.kegg.brite.hierarchy.with.gene.ids -n 10 > "$line"/br.res
awk -F"\t" 'NF == 5' "$line"/br.res | sed 's/^\s\+//' | sed 's/|/\t/' > "$line"/"$line"BRITE_enrichment.tsv; done

##R script 08c_cluster_KEGG_enrichment.R

##R script 08d_ballgown_MF_DE.R

#male-female pathway enrichment (see 08d_ballgown_MF_DE.R)
python /home/pakorhon/Codebase/diffExpKeggPwEnrichment.py -d KeggPwMale -i male_up_trans.tsv -r 0.05 -b kegg_blast.blast6 -s + | sort -k 2,2 -gr | grep -v "Human Dis" > pwMale.res
python /home/pakorhon/Codebase/diffExpKeggPwEnrichment.py -d KeggPwFemale -i female_up_trans.tsv -r 0.05 -b kegg_blast.blast6 -s + | sort -k 2,2 -gr | grep -v "Human Dis" > pwFemale.res

python /home/pakorhon/Codebase/diffExpKeggBrEnrichment.py -i KeggPwMale/kegg_blast.blast6.kegg.brite.freqs -r  KeggPwMale/kegg_blast.blast6.kegg.brite.hierarchy.with.gene.ids -s KeggPwMale/-.kegg.brite.hierarchy.with.gene.ids -n 10 > KeggPwMale/br.res
awk -F"\t" 'NF == 5' KeggPwMale/br.res | sed 's/^\s\+//' | sed 's/|/\t/' > KeggPwMale/BRITE_enrichment.tsv

python /home/pakorhon/Codebase/diffExpKeggBrEnrichment.py -i KeggPwFemale/kegg_blast.blast6.kegg.brite.freqs -r  KeggPwFemale/kegg_blast.blast6.kegg.brite.hierarchy.with.gene.ids -s KeggPwFemale/-.kegg.brite.hierarchy.with.gene.ids -n 10 > KeggPwFemale/br.res
awk -F"\t" 'NF == 5' KeggPwFemale/br.res | sed 's/^\s\+//' | sed 's/|/\t/' > KeggPwFemale/BRITE_enrichment.tsv

#scaffold enrichment analysis (see 08d_ballgown_MF_DE.R)
cut -f1 male_up_trans.tsv | while read line; do grep -w $line Shae.V3_final.gff3; done | awk -F"\t" '$3 == "mRNA"' | cut -f1 | sort | uniq -c | sed 's/^\s\+//' | sed 's/ /\t/' | sort -nr > DE_male_scaffold_counts.tsv
cut -f1 female_up_trans.tsv | while read line; do grep -w $line Shae.V3_final.gff3; done | awk -F"\t" '$3 == "mRNA"' | cut -f1 | sort | uniq -c | sed 's/^\s\+//' | sed 's/ /\t/' | sort -nr > DE_female_scaffold_counts.tsv
awk -F"\t" '$3 == "mRNA"' Shae.V3_final.gff3 | cut -f1 | sort | uniq -c | sed 's/^\s\+//' | sed 's/ /\t/' | sort -nr > all_trans_scaffold.tsv
awk -F"\t" '$3 == "mRNA"' Shae.V3_final.gff3 | cut -f1,9 | sed 's/;/\t/' | cut -f1,2 | sed 's/ID=//' > all_IDs_trans_scaffold.tsv

#For figure for DE genes on scaffolds
head -9 DE_female_scaffold_counts.tsv | awk -F"\t" 'BEGIN{print "female\tscaffold\tlength"}{print $0}' > fem_scaff_DE.tsv
fa2len 163_wgene_hypo_genome_ShV3.softmasked.fasta | sed 's/^>//' > scaffold_lengths.tsv