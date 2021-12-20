############Run LiftOver to transfer ShV2 genes onto new assembly#################
conda activate ucsc
faToTwoBit fixed_scaffolds_666_NNN.fa.masked.fasta ShV2.2bit
twoBitInfo ShV2.2bit stdout | sort -k2,2nr > target.chrom.sizes
faToTwoBit hypo_genome_ShV3.fa query.2bit
twoBitInfo query.2bit stdout | sort -k2,2nr > query.chrom.sizes
 
blat ShV2.2bit /dev/null /dev/null -tileSize=11 -makeOoc=target.ooc -repMatch=300
 
rm -rf run.blat/
rm -rf run.chain/
rm -rf doCleanup.csh doLoad.csh fb.target.chain.QueryLink.txt chainQuery.bb chainQueryLink.bb targetToQuery.over.chain.gz
 
/usr/bin/time -v doSameSpeciesLiftOver.pl -verbose 2 -buildDir=`pwd` -ooc=`pwd`/target.ooc -fileServer=localhost -localTmp="/dev/shm" -bigClusterHub=localhost -dbHost=localhost -workhorse=localhost -target2Bit=`pwd`/ShV2.2bit -targetSizes=`pwd`/target.chrom.sizes -query2Bit=`pwd`/query.2bit -querySize=`pwd`/query.chrom.sizes target query > blat.out 2> blat.err
 
zcat target.fna.gz > target.fna
gffread shv2.introns_noComments.gff3 -g fixed_scaffolds_666_NNN.fa.masked.fasta -x target.cds.fasta
ldHgGene -nobin -out=target.gp -requireCDS ignored RefSeq target.gff
liftOver -genePred target.gp targetToQuery.over.chain.gz query.gp unmapped


awk -F"\t" '{split($1,a,"[;=]"); if(a[3] == "Parent"){pos=4}else{pos=2}; cnt[a[pos]]++; $1=a[pos]"."cnt[a[pos]]; for(i=1; i<NF; i++){printf $i"\t"}; print $NF}' query.gp > fixed_query.gp
genePredToGtf file fixed_query.gp query.gtf
sed -i 's|\(gene_id \"mRNA[^.]\+\).[^;]\+;|\1";|' query.gtf
gffread query.gtf -L -g hypo_genome_ShV3.fa > query.gff
sed -i 's|fixed_query.gp|liftOver|' query.gff


###############Read mapping and transcript inference using StringTie#############
#map long reads against genome using minimap2
LONGR=merged_long_reads.fq
minimap2 -t 24 -ax splice -uf -k14 hypo_genome_ShV3.fa $LONGR | samtools sort -@ 24 -o long_reads.mm2_ShV3.sorted.bam
python ~/miniconda2/envs/flair/bin/bin/bam2Bed12.py -i long_reads.mm2_ShV3.sorted.bam > long_reads.mm2_ShV3.sorted.bed12

#Map Short reads (needed for both StringTie and polishing long reads using FLAIR)
hisat2-build hypo_genome_ShV3.fa ShV3
hisat2 --dta -p 24 -x ShV3 -1 fastp_centr_K26_khmer_Shae_R1.fq -2 fastp_centr_K26_khmer_Shae_R2.fq -U fastp_centr_K26_khmer_Shae_single.fq --summary-file hisat.summary 2> /dev/null | samtools sort -o short_reads.hisat2_ShV3.sorted.bam

#recreate stringent isoforms
samtools index -b long_reads.mm2_ShV3.sorted.bam 
samtools index -b short_reads.hisat2_ShV3.sorted.bam


#Create high-confidence long-read transcripts using FLAIR and both mapped short- and long-reads
python ~/miniconda2/envs/flair/bin/bin/junctions_from_sam.py -s short_reads.hisat2_ShV3.sorted.bam -n short_reads 2>junc.err 1>junc.out

awk '{if($5 > 2){print}}' short_reads_junctions.bed > short_reads_junctions.filtered.bed
sed -i 's/^chr//' short_reads_junctions.filtered.bed

python ~/miniconda2/envs/flair/bin/flair.py correct -t 24 -q long_reads.mm2_ShV3.sorted.bed12 -g hypo_genome_ShV3.fa -j short_reads_junctions.filtered.bed 2> flair_correct.err 1> flair_correct.out

python ~/miniconda2/envs/flair/bin/flair.py collapse -t 24 -n best_only --max_ends 1 --filter nosubset -o stringent --stringent -q flair_all_corrected.bed -g hypo_genome_ShV3.fa -r merged_long_reads.fq 2> flair_collapse_stringent.err 1> flair_collapse_stringent.out

python ~/miniconda2/envs/flair/bin/bin/psl_to_gtf.py stringent.isoforms.bed --force > stringent.isoforms.gtf

#Make gff3 file from gtf
gffread -E stringent.isoforms.gtf -o stringent.isoforms.gff3
sed -i /gffread/d stringent.isoforms.gff3
#Use gt interfeat to add introns
gt interfeat -o introns_stringent.isoforms.gff3 stringent.isoforms.gff3 
grep "geneID" introns_stringent.isoforms.gff3 | cut -f9 | sed 's/;/\t/g' | sed 's/ID=/\t/g' | cut -f2,4 > geneTrans.lkp
fa2tab stringent.isoforms.fa | sed 's/;/:/' | awk -F"\t" 'BEGIN{while(getline < "geneTrans.lkp"){lkp[">"$2]=">"$1}}{print lkp[$1]; print $2}' > stringent.isoforms_transfixed.fa

#merge long-read evidence with LiftOver GFF
mv query.gff LiftOverShV3.gff
stringtie --merge -G LiftOverShV3.gff -l MLRSHV3 -o merged_LR.gtf stringent.isoforms.gtf > LR_merge.log 2>&1

#Run TransDecoder
conda activate transdec
/data/miniconda3/envs/transdec/bin/gtf_genome_to_cdna_fasta.pl merged_LR.gtf hypo_genome_ShV3.fa > merged_LR.fasta
/data/miniconda3/envs/transdec/bin/gtf_to_alignment_gff3.pl merged_LR.gtf > merged_LR.gff3
TransDecoder.LongOrfs -t merged_LR.fasta
TransDecoder.Predict -t merged_LR.fasta
/data/miniconda3/envs/transdec/bin/cdna_alignment_orf_to_genome_orf.pl merged_LR.fasta.transdecoder.gff3 merged_LR.gff3 merged_LR.fasta > merged_LR.fasta.transdecoder.genome.gff3

gffread merged_LR.fasta.transdecoder.genome.gff3 -g hypo_genome_ShV3.fa -y proteins_merged_LR.fasta


###############SHORT READS###########
stringtie ../short_reads.hisat2_ShV3.sorted.bam -G ../LiftOverShV3.gff -l SRSH -o SR_predict_stringtie.gtf -p 24 -m 90 -A abundances.tab -C ref_coverage.gtf --conservative > SR_predict.log 2>&1
stringtie --merge -G ../LiftOverShV3.gff -l MSRSH -o merged_SR.gtf SR_predict_stringtie.gtf -p 24 > SR_merge.log 2>&1

#Run TransDecoder
conda activate transdec
/data/miniconda3/envs/transdec/bin/gtf_genome_to_cdna_fasta.pl merged_SR.gtf ../hypo_genome_ShV3.fa > merged_SR.fasta
/data/miniconda3/envs/transdec/bin/gtf_to_alignment_gff3.pl merged_SR.gtf > merged_SR.gff3
TransDecoder.LongOrfs -t merged_SR.fasta
TransDecoder.Predict -t merged_SR.fasta

/data/miniconda3/envs/transdec/bin/cdna_alignment_orf_to_genome_orf.pl merged_SR.fasta.transdecoder.gff3 merged_SR.gff3 merged_SR.fasta > merged_SR.fasta.transdecoder.genome.gff3
gffread merged_SR.fasta.transdecoder.genome.gff3 -g ../hypo_genome_ShV3.fa -y proteins_merged_SR.fasta


#MERGE ShortRead and LongRead
stringtie --merge -G merged_SR.fasta.transdecoder.genome.gff3 -l SLmerge -o merged_both_SR_LR.gtf ../merged_LR.fasta.transdecoder.genome.gff3 > SL_merge.log 2>&1

conda activate transdec
/data/miniconda3/envs/transdec/bin/gtf_genome_to_cdna_fasta.pl merged_both_SR_LR.gtf ../fixed_scaffolds_666_NNN.fa > merged_both_SR_LR.fasta
/data/miniconda3/envs/transdec/bin/gtf_to_alignment_gff3.pl merged_both_SR_LR.gtf > merged_both_SR_LR.gff3
TransDecoder.LongOrfs -t merged_both_SR_LR.fasta
TransDecoder.Predict -t merged_both_SR_LR.fasta
/data/miniconda3/envs/transdec/bin/cdna_alignment_orf_to_genome_orf.pl merged_both_SR_LR.fasta.transdecoder.gff3 merged_both_SR_LR.gff3 merged_both_SR_LR.fasta > merged_both_SR_LR.fasta.transdecoder.genome.gff3
gffread merged_both_SR_LR.fasta.transdecoder.genome.gff3 -g ../fixed_scaffolds_666_NNN.fa -y proteins_merged_both_SR_LR.fasta

busco -c 20 -m protein -i proteins_merged_both_SR_LR.fasta -o metazoa_both_SR_LR_StringTie_busco -l metazoa_odb10 --update-data

#Run StringTie for final gene set again with long read data to infer transcription levels and keep only most hihgly transcribed isoforms
stringtie ../../long_reads.mm2_ShV3.sorted.bam -eB -G ../merged_both_SR_LR.fasta.transdecoder.genome.gff3 -p 24 -o LR_transcripts_cov_tpm.gtf > /dev/null 2>&1
awk -F"\t" '$3 == "transcript"' LR_transcripts_cov_tpm.gtf | cut -f9 | sed 's/[ ;]/\t/g' | sed 's/"//g' | cut -f2,5,8,11,14 > transcript_gene_levels_LR.tsv
#In R
# lr <- read_tsv(file = "Downloads/transcript_gene_levels_LR.tsv", col_names = T)
# lr %>% filter(TPM > 1) %>% ggplot(aes(x= TPM)) + geom_histogram() + xlim(0,50) + theme_minimal() + xlab("Long-read TPM") + ylab("Number of transcripts with TPM > 1 (n = 8238/9932 genes)") + ggsave("Downloads/long_read_tpm.pdf")

#map long reads back to transcripts and infer length coverage
gffread merged_both_SR_LR.fasta.transdecoder.genome.gff3 -g ../hypo_genome_ShV3.fa -w SR_LR_transcripts.fa
minimap2 -t 24 -ax map-ont -N 5 --secondary=no -uf -k14 SR_LR_transcripts.fa ../merged_long_reads.fq | samtools sort -@ 24 -o long_reads_transcripts_ShV3.sorted.bam
fa2len SR_LR_transcripts.fa | sed 's/ /\t/' | cut -f1,3 > transcripts.txt
genomeCoverageBed -ibam long_reads_transcripts_ShV3.sorted.bam -g transcripts.txt -d > long_reads_transcripts_ShV3.coverage.txt
awk -F"\t" '{len[$1]++; if($3 > 0){cov[$1]++}}END{for(i in len){print i"\t"cov[i]/len[i]}}' long_reads_transcripts_ShV3.coverage.txt > perc_bases_covered_lr.tsv
#in R
# base_cov_lr %>% mutate(perc100=perc*100) %>% ggplot(aes(x = perc100)) + geom_histogram() + theme_minimal() + xlab("Percent of transcript covered by long reads") + ylab("Number of transcripts") + ggsave("Downloads/long_read_coverage.pdf")

#Run StringTie for final gene set again with short read data to infer transcription levels and keep only most hihgly transcribed isoforms
stringtie ../short_reads.hisat2_ShV3.sorted.bam -eB -G merged_both_SR_LR.fasta.transdecoder.genome.gff3 -p 24 -o transcripts_cov_tpm.gtf > /dev/null 2>&1 &
awk -F"\t" '$3 == "transcript"' transcripts_cov_tpm.gtf | cut -f9 | sed 's/[ ;]/\t/g' | sed 's/"//g' | cut -f2,5,8,11,14 > transcript_gene_levels_SR.tsv
#in R
# sr <- read_tsv(file = "Downloads/transcript_gene_levels_SR.tsv", col_names = T)
# sr %>% filter(TPM > 1) %>% ggplot(aes(x= TPM)) + geom_histogram() + xlim(0,100) + theme_minimal() + xlab("Short-read TPM") + ylab("Number of transcripts with TPM > 1 (n = 9780/9932 genes)") + ggsave("Downloads/short_read_tpm.pdf")

#Map protein sequences from S. mansoni and S. japonicum
fa2tab Shaem_PR_Sjap_Sman_100.fasta | grep -v ">mRNA" | awk -F"\t" '{print $1; print $2}' > Sman_Sjap.fasta
#Run FastExonerate on Sman_Sjap.fasta to create SmanSjapFastExon.gff
#Run StringTie merge from final set and Protein GFF
stringtie --merge -G SRLR_final.gff -l PROTSL -o PROT_SR_LR.gtf SmanSjapFastExon.gff > PROT_SL_merge.log 2>&1
/data/miniconda3/envs/transdec/bin/gtf_genome_to_cdna_fasta.pl PROT_SR_LR.gtf hypo_genome_ShV3.fa > PROT_SR_LR.fasta
/data/miniconda3/envs/transdec/bin/gtf_to_alignment_gff3.pl PROT_SR_LR.gtf > PROT_SR_LR.gff3
TransDecoder.LongOrfs -t PROT_SR_LR.fasta
TransDecoder.Predict -t PROT_SR_LR.fasta
/data/miniconda3/envs/transdec/bin/cdna_alignment_orf_to_genome_orf.pl PROT_SR_LR.fasta.transdecoder.gff3 PROT_SR_LR.gff3 PROT_SR_LR.fasta > PROT_SR_LR.fasta.transdecoder.genome.gff3
gffread PROT_SR_LR.fasta.transdecoder.genome.gff3 -g hypo_genome_ShV3.fa -y proteins_PROT_SR_LR.fasta

#################ISOFORM SELECTION################
#number of unique sequences
fa2tab SR_LR_transcripts.fa | cut -f2 | sort | uniq | wc -l
#33772

#Many isoforms predicted by StringTie are exact duplicates remove those first then remap short and long reads
fa2tab SR_LR_transcripts.fa | sed 's/ /\t/g' | cut -f1,3 | sed 's/^>//' | awk -F"\t" 'BEGIN{while(getline < "gene_trans.tsv"){gene[$1]=$2}}{print $1"\t"gene[$1]"\t"$2}' | awk -F"\t" '{iso[$2"\t"$3]=iso[$2"\t"$3]","$1}END{for(i in iso){print i"\t"iso[i]}}' | cut -f1,3 | sed 's/\t,/\t/' > duplicate_isoforms.tsv
sed 's/,/\t/1' duplicate_isoforms.tsv | cut -f1,2 > select_dupl_isoforms.ls

fa2tab SR_LR_transcripts.fa | sed 's/ /\t/g' | cut -f1,3 | sed 's/^>//' | awk -F"\t" 'BEGIN{while(getline < "select_dupl_isoforms.ls"){sel[$0]=$0}}{if($1 in sel){print ">"$1; print $2}}' > cleaned_isoforms.fa
fa2tab SR_LR_transcripts.fa | sed 's/ /\t/g' | cut -f1,3 | sed 's/^>//' | awk -F"\t" 'BEGIN{while(getline < "select_dupl_isoforms.ls"){sel[$2]=$1}}{if($1 in sel){print sel[$1]"\t"$1"\t"$2"\t"length($2)}}' > cleaned_isoforms.tsv

#Check if there are also duplicate genes (in the same location, we will not throw out genes predicted in different locations on the genome that produce the same transcript, they may be real)
 awk -F"\t" '$3 == "gene"' merged_both_SR_LR.fasta.transdecoder.genome.gff3 | cut -f1,4,5 | sort | uniq -c | grep "^ *2" | sed 's/ \+/\t/g' | cut -f3- | awk -F"\t" '{print $1"\ttransdecoder\tgene\t"$2"\t"$3}' > duplicate_genes.tsv

#remove SLmerge.4409
#remove SLmerge.197
awk -F"\t" '{split($9,att,"[=;]"); if($3 == "gene" && (att[2] == "SLmerge.4409" || att[2] == "SLmerge.197")){getline; while($3 != "gene"){getline}; print $0}else{print $0}}' merged_both_SR_LR.fasta.transdecoder.genome.gff3 > dupl_genes_removed.gff


#Remove isoforms fully contained in another isoform
#by definition a longer isoform cannot be fully contained in a shorter isoform so if you order them you can check for containment
cat cleaned_isoforms.tsv | while read gene trans seq len;  do echo $trans; grep -w "^$gene" cleaned_isoforms.tsv | grep "$seq" | wc -l; done | awk '{printf $0; getline; print "\t"$0}' > contained_isos.tsv
awk -F"\t" '$2 == 1' contained_isos.tsv | cut -f1 > isos_to_keep.ls
grep -v "#" full_table.tsv | cut -f3 | sed '/^$/d' | sort | uniq > ../../BUSCO_isos.ls
cat isos_to_keep.ls BUSCO_isos.ls | sort | uniq > isos_plus_BUSCO.ls
sed 's/[=;]/\t/g' gene_lines.gff | cut -f10 | sort | uniq > genes_to_keep.ls
awk -F"\t" 'BEGIN{while(getline < "isos_to_keep.ls"){trans[$1]=$1}; while(getline < "genes_to_keep.ls"){genes[$1]=$1}; print "##gff-version 3"}{split($9,a,"[=;]"); for(i=1; i<length(a); i++){att[a[i]]=a[i+1]}; if(att["ID"] in genes || att["ID"] in trans || att["Parent"] in trans){print $0}}' dupl_genes_removed.gff > isos_genes_kept.gff3

##########Remap short and long reads###############

gffread isos_genes_kept.gff3 -g ../hypo_genome_ShV3.fa -w isos_genes_kept_transcripts.fa
hisat2-build isos_genes_kept_transcripts.fa ISOS
hisat2 -p 24 -x ISOS -1 ../fastp_centr_K26_khmer_Shae_R1.fq -2 ../fastp_centr_K26_khmer_Shae_R2.fq -U ../fastp_centr_K26_khmer_Shae_single.fq --no-spliced-alignment --summary-file ISOS.summary 2> /dev/null | samtools sort -@ 24 -o short_reads.ISOS_no-splice.sorted.bam &
minimap2 -t 24 -ax map-ont -N 5 --secondary=no -uf -k14 isos_genes_kept_transcripts.fa ../merged_long_reads.fq | samtools sort -@ 24 -o long_reads.ISOS.sorted.bam

stringtie long_reads.mm2_ShV3.sorted.bam -eB -G isos_genes_kept.gff3 -p 24 -o ISOS_LR_transcripts_cov_tpm.gtf > /dev/null 2>&1 &
awk -F"\t" '$3 == "transcript"' ISOS_LR_transcripts_cov_tpm.gtf | cut -f9 | sed 's/[ ;]/\t/g' | sed 's/"//g' | cut -f2,5,8,11,14 > ISOS_transcript_gene_levels_LR.tsv

stringtie short_reads.hisat2_ShV3.sorted.bam -eB -G isos_genes_kept.gff3 -p 24 -o ISOS_SR_transcripts_cov_tpm.gtf > /dev/null 2>&1 &
awk -F"\t" '$3 == "transcript"' ISOS_SR_transcripts_cov_tpm.gtf | cut -f9 | sed 's/[ ;]/\t/g' | sed 's/"//g' | cut -f2,5,8,11,14 > ISOS_transcript_gene_levels_SR.tsv

fa2len isos_genes_kept_transcripts.fa | sed 's/ /\t/' | cut -f1,3 > ISOS_transcripts.txt
genomeCoverageBed -ibam long_reads.ISOS.sorted.bam -g ISOS_transcripts.txt -d > ISOS_LR_transcripts_ShV3.coverage.txt
awk -F"\t" '{len[$1]++; if($3 > 0){cov[$1]++}}END{for(i in len){print i"\t"cov[i]/len[i]}}' ISOS_LR_transcripts_ShV3.coverage.txt > ISOS_perc_bases_covered_lr.tsv
genomeCoverageBed -ibam short_reads.ISOS_no-splice.sorted.bam -g ISOS_transcripts.txt -d > ISOS_SR_transcripts_ShV3.coverage.txt
awk -F"\t" '{len[$1]++; if($3 > 0){cov[$1]++}}END{for(i in len){print i"\t"cov[i]/len[i]}}' ISOS_SR_transcripts_ShV3.coverage.txt > ISOS_perc_bases_covered_sr.tsv

fa2len ../proteins_isos_genes_kept.fasta | awk -F"\t" 'BEGIN{print "trans_id\tcds_len"}{print substr($1,2)"\t"$2*3}' > ISOS_CDS_length.tsv
fa2tab ../proteins_isos_genes_kept.fasta | awk -F"\t" 'BEGIN{while(getline < "15404_selected_isoforms.tsv"){trans[">"$1]=$1}}{if($1 in trans){print $1; print $2}}' > 15404_selected_isoforms.fasta

busco -c 24 -m protein -i 15404_selected_isoforms.fasta -o BUSCO_isos_transcribed -l metazoa_odb10 --update-data

grep -v "#" full_table.tsv | cut -f1,3 | awk -F"\t" '$2 != ""' |  grep -w -f /home/centos/store_sh/V3StringTieSh/SR/ISOS/BUSCO_isos_transcribed/run_metazoa_odb10/missing_busco_list.tsv | cut -f2 > ../../keep_missing_BUSCOs.ls

#add isoforms that represent complete BUSCOs in the previous busco but are fragmented now
grep -v "#" full_table.tsv | grep -e "Complete" -e "Duplicated" | cut -f1,3 |  grep -w -f /home/centos/store_sh/V3StringTieSh/SR/ISOS/BUSCO_isos_transcribed/run_metazoa_odb10/fragmented.ls | cut -f2 > ../../keep_fragmented_BUSCOs.ls
cat ISOS/15404_selected_isoforms.tsv keep_* | grep -v "trans_id" > 15555_keep.ls

awk -F"\t" 'BEGIN{while(getline < "15555_keep.ls"){trans[$1]=$1}; while(getline < "genes_to_keep.ls"){genes[$1]=$1}; print "##gff-version 3"}{split($9,a,"[=;]"); for(i=1; i<length(a); i++){att[a[i]]=a[i+1]}; if(att["ID"] in genes || att["ID"] in trans || att["Parent"] in trans){print $0}}' dupl_genes_removed.gff > final_to_keep.gff3

#select one that matches best a CDS of S. mansoni in addition to the ones transcribed
gffread merged_both_SR_LR.fasta.transdecoder.genome.gff3 -g ../hypo_genome_ShV3.fa -x cds_SR_LR.fasta
fa2tab cds_SR_LR.fasta | sed 's/^>//' | awk -F"\t" 'BEGIN{while(getline < "isos_to_keep.ls"){trans[$1]=$1}}{if($1 in trans){print $0}}' | awk -F"\t" 'BEGIN{while(getline < "cleaned_isoforms.tsv"){trans[$2]=$1}}{if($1 in trans){print trans[$1]"\t"$0}}' | awk -F"\t" 'BEGIN{while(getline < "genes_to_keep.ls"){trans[$1]=$1}}{if($1 in trans){print ">"$1"__"$2; print $3}}' > 32119_CDSs_nonred.fa
conda create -n blast -c bioconda blast

makeblastdb -in SMANCDS.fa -dbtype nucl
blastn -query 32119_CDSs_nonred.fa -db SMANCDS.fa -out sman_blastn_32119.blast6 -evalue 1e-30 -num_threads 24 -max_target_seqs 1 -outfmt 6
sed 's/__/\t/1' sman_blastn_32119.blast6 | sort -t$'\t' -k1,1 -k13,13nr | sort -u -k1,1 --merge | cut -f1,2 > best_iso_hits_32119.tsv
cut -f2 best_iso_hits_32119.tsv > ~/store_sh/V3StringTieSh/SR/best_sman_hit.tsv
cat best_sman_hit.tsv ISOS/15404_selected_isoforms.tsv | sort | uniq > transcribed_or_sman.ls

awk -F"\t" 'BEGIN{while(getline < "transcribed_or_sman.ls"){trans[$1]=$1}; while(getline < "genes_to_keep.ls"){genes[$1]=$1}; print "##gff-version 3"}{split($9,a,"[=;]"); for(i=1; i<length(a); i++){att[a[i]]=a[i+1]}; if(att["ID"] in genes || att["ID"] in trans || att["Parent"] in trans){print $0}}' dupl_genes_removed.gff > transcribed_or_sman.gff3

cat transcribed_or_sman.ls keep_* | sort | uniq | grep -v trans > 17685_keep_final.ls

awk -F"\t" 'BEGIN{while(getline < "17685_keep_final.ls"){trans[$1]=$1}; while(getline < "genes_to_keep.ls"){genes[$1]=$1}; print "##gff-version 3"}{split($9,a,"[=;]"); for(i=1; i<length(a); i++){att[a[i]]=a[i+1]}; if(att["ID"] in genes || att["ID"] in trans || att["Parent"] in trans){print $0}}' dupl_genes_removed.gff > 17685_keep_final.gff3
bash clean_gff.sh 17685_keep_final.gff3 > cleaned_17685_final.gff3
gt gff3 -sort -tidy -retainids cleaned_17685_final.gff3 > gt_cleaned_17685_final.gff3 2> gt_report.txt
sed -i 's/\ttransdecoder\t/\tShaeV3\t/' gt_cleaned_17685_final.gff3 
awk -F"\t" '$3 == "mRNA"' gt_cleaned_17685_final.gff3 | sed 's/;/\t/' | cut -f9 | sed 's/^ID=//' | awk '{print $0"\tmRNA"++a}' > mRNA_introns.lkp
gt interfeat -o gt_introns_cleaned_17685_final.gff3 gt_cleaned_17685_final.gff3 
awk -F"\t" '$3 == "intron"' gt_introns_cleaned_17685_final.gff3 | sed 's/\./ShaeV3/1' | sed 's/\tParent=/\t/' | awk -F"\t" 'BEGIN{while(getline < "mRNA_introns.lkp"){lkp[$2]=$1}}{icnt[lkp[$9]]++; $9="ID="lkp[$9]".intron"icnt[lkp[$9]]";Parent="lkp[$9]; for(i=1; i<NF; i++){printf $i"\t"}; print $NF}' >> introns.gff
cat gt_cleaned_17685_final.gff3 introns.gff  | grep -v "#"  | sort -k1 -k4n -k5n > add_introns.gff
gt gff3 -sort -tidy -retainids add_introns.gff > 17685_introns_final.gff3

#There are still isoforms that are exactly the same
#when the protein is an ecxact duplicate, keep the one with longest mRNA of those duplicates only and discard the others
awk -F"\t" '$3 == "mRNA"' 17685_keep_final.gff3 | cut -f9 | sed 's/[;=]/\t/g' | cut -f2,4 > 17685_ids.tsv
fa2len test_mrna.fa | sed 's/>//' | sed 's/ /\t/' | cut -f1,3 > mrna_len.tsv
fa2tab test_proteins.fa | sed 's/>//' | awk -F"\t" 'BEGIN{while(getline < "mrna_len.tsv"){len[$1]=$2}}{print $0"\t"len[$1]}' | awk -F"\t" 'BEGIN{while(getline < "17685_ids.tsv"){gene[$1]=$2}}{print $0"\t"gene[$1]}' | awk -F"\t" '{if($4" "$2 in lkp){if($3 > lkp[$4" "$2]){lkp[$4" "$2]=$3; id[$4" "$2]=$1}}else{lkp[$4" "$2]=$3; id[$4" "$2]=$1;}}END{for(i in id){print id[i]}}' > 16303_keep.ls
awk -F"\t" 'BEGIN{while(getline < "16303_keep.ls"){trans[$1]=$1}; while(getline < "genes_to_keep.ls"){genes[$1]=$1}; print "##gff-version 3"}{split($9,a,"[=;]"); for(i=1; i<length(a); i++){att[a[i]]=a[i+1]}; if(att["ID"] in genes || att["ID"] in trans || att["Parent"] in trans){print $0}}' 17685_keep_final.gff3 > 16303_keep_final.gff3

bash clean_gff.sh 16303_keep_final.gff3 > cleaned_16303_final.gff3
gt gff3 -sort -tidy -retainids cleaned_16303_final.gff3 > gt_cleaned_16303_final.gff3 2> gt_report.txt
sed -i 's/\ttransdecoder\t/\tShaeV3\t/' gt_cleaned_16303_final.gff3 
awk -F"\t" '$3 == "mRNA"' gt_cleaned_16303_final.gff3 | sed 's/;/\t/' | cut -f9 | sed 's/^ID=//' | awk '{print $0"\tmRNA"++a}' > mRNA_introns.lkp
gt interfeat -o gt_introns_cleaned_16303_final.gff3 gt_cleaned_16303_final.gff3 
awk -F"\t" '$3 == "intron"' gt_introns_cleaned_16303_final.gff3 | sed 's/\./ShaeV3/1' | sed 's/\tParent=/\t/' | awk -F"\t" 'BEGIN{while(getline < "mRNA_introns.lkp"){lkp[$2]=$1}}{icnt[lkp[$9]]++; $9="ID="lkp[$9]".intron"icnt[lkp[$9]]";Parent="lkp[$9]; for(i=1; i<NF; i++){printf $i"\t"}; print $NF}' > introns.gff
cat gt_cleaned_16303_final.gff3 introns.gff  | grep -v "#"  | sort -k1 -k4n -k5n > add_introns.gff
gt gff3 -sort -tidy -retainids add_introns.gff > 16303_introns_final.gff3

#Do it the other way round: when the transcript is an exact duplicate pick the one that encocdes the longest protein
fa2len test_proteins.fa | sed 's/>//' | sed 's/ /\t/' | cut -f1,2 > prot_len.tsv
fa2tab test_mrna.fa | sed 's/>//' | sed 's/ /\t/' | cut -f1,3 | awk -F"\t" 'BEGIN{while(getline < "prot_len.tsv"){len[$1]=$2}}{print $0"\t"len[$1]}' | awk -F"\t" 'BEGIN{while(getline < "17685_ids.tsv"){gene[$1]=$2}}{print $0"\t"gene[$1]}' | awk -F"\t" '{if($4" "$2 in lkp){if($3 < lkp[$4" "$2]){lkp[$4" "$2]=$3; print id[$4" "$2]; id[$4" "$2]=$1}}else{lkp[$4" "$2]=$3; id[$4" "$2]=$1;}}' > remove_ident_trans.ls
grep -w -v -f remove_ident_trans.ls 16303_keep.ls  > 16284_keep.ls
awk -F"\t" 'BEGIN{while(getline < "16284_keep.ls"){trans[$1]=$1}; while(getline < "genes_to_keep.ls"){genes[$1]=$1}; print "##gff-version 3"}{split($9,a,"[=;]"); for(i=1; i<length(a); i++){att[a[i]]=a[i+1]}; if(att["ID"] in genes || att["ID"] in trans || att["Parent"] in trans){print $0}}' 17685_keep_final.gff3 > 16284_keep_final.gff3

conda activate gtools
bash clean_gff.sh 16284_keep_final.gff3 > cleaned_16284_final.gff3
gt gff3 -sort -tidy -retainids cleaned_16284_final.gff3 > gt_cleaned_16284_final.gff3 2> gt_report.txt
sed -i 's/\ttransdecoder\t/\tShaeV3\t/' gt_cleaned_16284_final.gff3 
awk -F"\t" '$3 == "mRNA"' gt_cleaned_16284_final.gff3 | sed 's/;/\t/' | cut -f9 | sed 's/^ID=//' | awk '{print $0"\tmRNA"++a}' > mRNA_introns.lkp
gt interfeat -force -o gt_introns_cleaned_16284_final.gff3 gt_cleaned_16284_final.gff3 
awk -F"\t" '$3 == "intron"' gt_introns_cleaned_16284_final.gff3 | sed 's/\./ShaeV3/1' | sed 's/\tParent=/\t/' | awk -F"\t" 'BEGIN{while(getline < "mRNA_introns.lkp"){lkp[$2]=$1}}{icnt[lkp[$9]]++; $9="ID="lkp[$9]".intron"icnt[lkp[$9]]";Parent="lkp[$9]; for(i=1; i<NF; i++){printf $i"\t"}; print $NF}' > introns.gff
cat gt_cleaned_16284_final.gff3 introns.gff  | grep -v "#"  | sort -k1 -k4n -k5n > add_introns.gff
gt gff3 -sort -tidy -retainids add_introns.gff > 16284_introns_final.gff3
conda deactivate

conda activate funannotate
funannotate util gff2tbl -g 16284_introns_final.gff3 -f 163_wgene_hypo_genome_ShV3.softmasked.fasta > 16284_introns_final.tbl
funannotate util tbl2gbk -i 16284_introns_final.tbl -f 163_wgene_hypo_genome_ShV3.softmasked.fasta -s "Schistosoma haematobium" --isolate Egypt --sbt template.sbt -o ShV3 -M n -J -c w -euk -gaps-min 10 -l paired-ends -locus-tag-prefix MS3

#Remove all problematic isoforms and genes and repredict the rest using Augustus and PASA, also use all ones without problems as input for Augustus training
grep "Sha" ShV3/genome.val | sed 's/\[/\n/g; s/\]/\n/g; s/|/\n/g; s/ /\n/g; s/(/\n/g; s/:/\n/g;' | grep Sha | sort | uniq > genes_isos_w_problems.tsv

#Create mRNA sequences from good ones
gffread 16284_introns_final.gff3 -g 163_wgene_hypo_genome_ShV3.softmasked.fasta -w - | awk '{print $1}' | fa2tab | sed 's/^>//' | awk -F"\t" 'BEGIN{while(getline < "genes_isos_w_problems.tsv"){problem[$1]=$1}}{split($1,a,"."); if($1 in problem || a[1] in problem){}else{print ">"$1; print $2}}' > 13318_mRNA_no_probs.fa

#create protein seqs
gffread 16284_introns_final.gff3 -g 163_wgene_hypo_genome_ShV3.softmasked.fasta -y - | awk '{print $1}' | fa2tab | sed 's/^>//' | awk -F"\t" 'BEGIN{while(getline < "genes_isos_w_problems.tsv"){problem[$1]=$1}}{split($1,a,"."); if($1 in problem || a[1] in problem){}else{print ">"$1; print $2}}' > 13318_prots_no_probs.fasta

#Make the gene line coordinates match up with the coordinates of the longest isoform that is present (we may have removed isoforms but haven't adjusted the coordinates in the gene line)
gffread -g 163_wgene_hypo_genome_ShV3.softmasked.fasta 16284_introns_final.gff3 -F --tlf | grep -v -f genes_isos_w_problems.tsv | gffread -g 163_wgene_hypo_genome_ShV3.softmasked.fasta - -F | awk -F"\t" 'BEGIN{while(getline < "16284_introns_final.gff3"){if($3 == "gene"){split($9,a,"="); gene[a[2]]=$0}}}{split($9,g,"="); if($3 == "mRNA" && fnd[g[3]] != 1){fnd[g[3]]=1; print gene[g[3]]; print $0}else{print $0}}' | sed 's/geneID/Parent/g' > 13317_no_probs.gff3

conda activate gtools
bash add_IDs_gff.sh 13317_no_probs.gff3 > cleaned_13317_no_probs.gff3
gt gff3 -sort -tidy -retainids cleaned_13317_no_probs.gff3 > gt_cleaned_13317_no_probs.gff3 2> gt_report.log

awk -F"\t" '$3 == "mRNA"' gt_cleaned_13317_no_probs.gff3 | sed 's/;Parent=/\t/g' | cut -f4,5,10 | awk -F"\t" '{if(min[$3]){if(min[$3] > $1){min[$3]=$1}}else{min[$3]=$1}; if(max[$3]){if(max[$3] < $2){max[$3]=$2}}else{max[$3]=$2};}END{for(i in min){print i"\t"min[i]"\t"max[i]}}'> max_gene_length.tsv
awk -F"\t" 'BEGIN{while(getline < "max_gene_length.tsv"){min[$1]=$2; max[$1]=$3}}{if($3 == "gene"){split($9,g,"="); $4=min[g[2]]; $5=max[g[2]]; printf $1; for(i=2; i<=8; i++){printf "\t"$i}; print "\t"$9}else{print $0}}' gt_cleaned_13317_no_probs.gff3 > fixed_gene_line_13317_no_probs.gff3
#use this for Augustus training


#run liftOver on the original 16284 dataset to create liftover hints for 
#Run liftover
#Create the .gp file so it can be supplied to GenePrediction pipeline using the liftOver input option
conda activate ucsc
ldHgGene -nobin -out=16284.gp -requireCDS ignored ShaeV3 16284_introns_final.gff3
sed 's/Parent=/\t/' 16284.gp | cut -f2- > query.matches.gp

#run training and prediction
#create workDir training.gff (fixed_gene_line_13317_no_probs.gff3)
#--> pipeline finished succesfully (10 hours runtime with liftOver and training)

#Now let's check what the gene models are like
#Do they come up with problems when running funannotate?
grep -v "# " 9822_PKK.gff3 | awk 'BEGIN{getline}{print $0}' > no_comm_9822_PKK.gff3 
#let's add one bad one so funannotate creates the output folder with the stats (if no problems occur no output files are created so we are using a "bad" dummy gene to create the output folder)
cp ../no_comm_9822_PKK.gff3 ../add_1_bad_no_comm_9822_PKK.gff3
grep "SEQ_FEAT.InternalStop" genome.val | sed 's/|/\t/g'  | cut -f5 | sed 's/\./\t/' | cut -f1 | sort | uniq > 20_problematic.ls
grep -f 20_problematic.ls ../15236_isos_retained.gff >> ../add_1_bad_no_comm_9822_PKK.gff3 
funannotate util gff2tbl -g add_1_bad_no_comm_9822_PKK.gff3 -f 163_wgene_hypo_genome_ShV3.softmasked.fasta > PKK_9822.tbl
funannotate util tbl2gbk -i PKK_9822.tbl -f 163_wgene_hypo_genome_ShV3.softmasked.fasta -s "Schistosoma haematobium" --isolate Egypt --sbt template.sbt -o PKK_9822_out

#Run PASA updater to add UTRs and check again

#################Create database and map long-read transcripts##########
docker run --rm -it -d --name=pasa_cont_ajs \
      -v  /home/andreas/PASA_Sh/tmp:/home/andreas/PASA_Sh/tmp \
      -v /home/andreas/PASA_Sh:/home/andreas/PASA_Sh \
       pasapipeline/pasapipeline:latest \
        bash -c 'cd /home/andreas/PASA_Sh \
              && /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
              -c alignAssembly.config -C -R \
              --ALIGNER gmap -g 163_wgene_hypo_genome_ShV3.softmasked.fasta -t stringent.isoforms.fa --CPU 48 2> align.err 1> align.log'
# fa38b36b018f67604e7474a21bc6ce366132dc4c4f83e5425313afaf2b7b4e65

##################Build transcriptome DB################
docker run --rm -it -d --name=pasa_cont_ajs \
      -v  /home/andreas/PASA_Sh/tmp:/home/andreas/PASA_Sh/tmp \
      -v /home/andreas/PASA_Sh:/home/andreas/PASA_Sh \
       pasapipeline/pasapipeline:latest \
        bash -c ' cd /home/andreas/PASA_Sh \
              && /usr/local/src/PASApipeline/scripts/build_comprehensive_transcriptome.dbi \
           -c alignAssembly.config \
           -t stringent.isoforms.fa \
           --min_per_ID 85 \
           --min_per_aligned 30 --CPU 48 2> build_trans.err 1> build_trans.log'
# 008f775aa29419302f1c75194d7e94633accefda1a6ff54529b9b014722f646e

##############Check GFF compatibility###################
docker run --rm -it -v /home/andreas/PASA_Sh/tmp:/home/andreas/PASA_Sh/tmp -v /home/andreas/PASA_Sh:/home/andreas/PASA_Sh pasapipeline/pasapipeline:latest \
bash -c 'cd /home/andreas/PASA_Sh && /usr/local/src/PASApipeline/misc_utilities/pasa_gff3_validator.pl no_comm_9822_PKK.gff3'
#ALL OK

###########Load current annotation#####################
docker run --rm -it -d --name=pasa_cont_ajs \
      -v  /home/andreas/PASA_Sh/tmp:/home/andreas/PASA_Sh/tmp \
      -v /home/andreas/PASA_Sh:/home/andreas/PASA_Sh \
       pasapipeline/pasapipeline:latest \
bash -c 'cd /home/andreas/PASA_Sh && /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g 163_wgene_hypo_genome_ShV3.softmasked.fasta -P no_comm_9822_PKK.gff3 2> load_anno.err 1> load_anno.log'
# 85eee54b8b6ccf174abb071c41b3284944ce092bd149d14d792058921077c22a

#########run updater################################
docker run --rm -it -d --name=pasa_cont_ajs \
      -v  /home/andreas/PASA_Sh/tmp:/home/andreas/PASA_Sh/tmp \
      -v /home/andreas/PASA_Sh:/home/andreas/PASA_Sh \
       pasapipeline/pasapipeline:latest \
       bash -c '
       bash /home/andreas/PASA_Sh/fix_fasta.ph.sh \
        && cd /home/andreas/PASA_Sh \
         && /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c annotCompare.config -A \
        -g 163_wgene_hypo_genome_ShV3.softmasked.fasta \
        -t stringent.isoforms.fa --CPU 48 2> update2.err 1> update2.log'
# 324e0a67322a2c2d017eb28cb8cefb46fb842216ce250227fa1a336d8f52ec5a

#FAILED because one of the sequences was >80% ACGT and therefore recognised as nucleotide then alignment couldn't be done against a protein seq
#Patch this in the docker beforehand
# using this patch
# cat fix_fasta.ph.sh 
# awk '{if($0 ~ /my \$cmd = "\$prog \$file1 \$file2 > \$result_file";/){gsub(/\$prog/,"$prog -p"); print}else{print}}' /usr/local/src/PASApipeline/PerlLib/fasta.ph > /usr/local/src/PASApipeline/PerlLib/tmp
# mv /usr/local/src/PASApipeline/PerlLib/tmp /usr/local/src/PASApipeline/PerlLib/fasta.ph
# head -76 /usr/local/src/PASApipeline/PerlLib/fasta.ph | tail -1 > /home/andreas/PASA_Sh/test.fph.out

# 4c5e336e142ec402bb562d685be0d032ba8b4352489d1843f7abcee11649d668

#Finished fine

conda activate gtools
sed '/^$/d' 9678_14052_PASIPASA.gff3 |grep -v "#" | gt gff3 -sort -tidy -retainids - > gt_cleaned_9678_14052_PASIPASA.gff3 2> gt_clean_PASIPASA.log

funannotate util gff2tbl -g gt_cleaned_9678_14052_PASIPASA.gff3 -f 163_wgene_hypo_genome_ShV3.softmasked.fasta > PASI_PASA.tbl
funannotate util tbl2gbk -i PASI_PASA.tbl -f 163_wgene_hypo_genome_ShV3.softmasked.fasta -s "Schistosoma haematobium" --isolate Egypt --sbt template.sbt -o PASI_PASA_out
#  418 SEQ_FEAT.NotSpliceConsensusAcceptor
   #   16 SEQ_FEAT.NotSpliceConsensusDonor
   #  302 SEQ_FEAT.PartialProblem
   # 1078 SEQ_FEAT.RareSpliceConsensusDonor
   #   76 SEQ_FEAT.ShortExon

#For ones where there is a good isoform don't add additional ones
sed 's/\./\t/g' genes_isos_w_problems.tsv | cut -f1 | sort | uniq > genes_w_probs.ls
#2260
#find if there are isoforms left that are good
#grep -f genes_w_probs.ls fixed_gene_line_13317_no_probs.gff3 | grep gene | wc -l
#662 genes that have at least one good isoform, not adding any new ones
awk -F"\t" 'BEGIN{while(getline < "genes_w_probs.ls"){probs[$1]=$1}}{split($9,g,"[=.]"); if(g[2] in probs){probs[g[2]]=""}}END{for(i in probs){if(probs[i] != ""){print probs[i]}}}' fixed_gene_line_13317_no_probs.gff3 > 1598_genes_w_probs.ls

#Extract all transcript isoforms from the 16284 file 
grep -f 1598_genes_w_probs.ls 16284_introns_final.gff3 | grep gene | wc -l
#1598
grep -f 1598_genes_w_probs.ls 16284_introns_final.gff3 | grep mRNA | wc -l
#2120

grep -f 1598_genes_w_probs.ls 16284_introns_final.gff3 > 2120_problem_isoforms.gff3

#Run trmap against the PASI-PASA one
#remove the "gene" line first, causes segmentation fault running trmap
awk -F"\t" '$3 == "mRNA" || $3 == "exon" || $3 == "CDS"' ../2120_problem_isoforms.gff3 > q2120_trans.gff3
awk -F"\t" '$3 == "mRNA" || $3 == "exon" || $3 == "CDS"' ../gt_cleaned_9678_14052_PASIPASA.gff3 > rPP9678_14052.gff3
trmap -o rPP9678_14052_q2120.out rPP9678_14052.gff3 q2120_trans.gff3

grep ">" rPP9678_14052_q2120.out | awk '{print $1}' | sort | uniq | wc -l
#1705 isoforms have >=1 matching isoform in the PP pipeline

grep ">" rPP9678_14052_q2120.out | awk '{print $1}' | sed 's/\./\t/' | cut -f1 | sort | uniq | wc -l
#1195 genes have >=1 matching gene in the PP pipeline

grep -v ">" rPP9678_14052_q2120.out | cut -f6 | grep "_" > 31_merged.ls

grep -v ">"  rPP9678_14052_q2120.out | cut -f6 | sort | uniq | wc -l
#1539 transcripts represented from PP pipeline

grep -v ">"  rPP9678_14052_q2120.out | cut -f6 | sed 's/\.t1[\.$]/\t/' | cut -f1 | sort | uniq | wc -l
#1210 genes represented from PP pipeline

awk -F"[\t ]" '{if($1 ~ /^>/){q=$1}else{l[q]=l[q]"\t"$6}}END{for(i in l){print i""l[i]} }' rPP9678_14052_q2120.out | cut -f2- | numcol
   # 1275 1
   #  247 2
   #    1 23
   #   97 3
   #   48 4
   #   21 5
   #   11 6
   #    2 7
   #    1 8
   #    2 9

awk -F"[\t ]" '{if($1 ~ /^>/){q=$1}else{print $6"\t"q}}' rPP9678_14052_q2120.out | sed 's/>//' | awk -F"\t" '{split($1,g,"."); split($2,i,"."); print g[1]"\t"$1"\t"i[1]"\t"$2}' > gene_iso_relships.tsv
cut -f1,2 gene_iso_relships.tsv | awk '{if($0 in fnd){}else{fnd[$0]=1; print $0}}' | awk -F"\t" '{if(!g[$1]){g[$1]=10000+len;} len=0; for(i in g){++len}; ++t[$1]; print $0"\tSha_"g[$1]"\tSha_"g[$1]"."t[$1]}' > rename_g.t_Sha.tsv

# merged gene models have different name structure
# #g9396.t1_g9394.t1
awk -F"\t" '{if($2 ~ /_/){orig=$2; gsub(".t1","",$2); print $2"\t"orig"\t"$3"\t"$4}else{print $0}}' rename_g.t_Sha.tsv  > merged_rename_g.t_Sha.tsv

cut -f1 merged_rename_g.t_Sha.tsv | sort | uniq | while read line; do grep -w $line gt_cleaned_9678_14052_PASIPASA.gff3; done | awk -F"\t" '{if($3 == "gene"){split($9,g,"[=;]"); print g[2]}}' > genes_9678_14052.txt

cut -f9 gt_cleaned_9678_14052_PASIPASA.gff3 | grep -v "#" | sed 's/[=;]/\t/g' | awk -F"\t" '{for(i=1; i<=NF; i=i+2){print $i}}' | sort | uniq -c
 #    126 5_prime_partial
 # 264457 ID
 #  23730 Name
 # 254779 Parent

grep -v "#" gt_cleaned_9678_14052_PASIPASA.gff3 | awk -F"\t" '
BEGIN{while(getline < "merged_rename_g.t_Sha.tsv"){g[$1]=$3; t[$2]=$4}}
{
	split($9,attr,"[=;]");
	$2="ShaeV3";
	if($3 == "gene" && attr[2] in g){fnd=1; attr[2]=g[attr[2]]; a=2};
	if($3 == "mRNA" && attr[2] in t){fnd=1; attr[2]=t[attr[2]]; attr[4]=g[attr[4]]; a=4};
	if($3 == "CDS" && attr[4] in t){fnd=1; attr[2]="cds."t[attr[4]]; attr[4]=t[attr[4]]; if(attr[6]=="true"){a=6}else{a=4}};
	if($3 ~ /prime_UTR/ && attr[4] in t){fnd=1; split(attr[2],u,".utr"); attr[2]=t[u[1]]".utr"u[2]; attr[4]=t[attr[4]]; a=4};
	if($3 == "exon" && attr[4] in t){fnd=1; split(attr[2],e,".exon"); attr[2]=t[e[1]]".exon"e[2]; attr[4]=t[attr[4]]; a=4};
	if(fnd == 1){
		fnd=0;
		for(i=1; i<=8; i++){printf $i"\t"};
		printf attr[1]"="attr[2];
		for(k=3; k<a; k=k+2){printf ";"attr[k]"="attr[k+1]};
		print ""}
}' > renamed_PP.gff

grep ">" 13317_final.pep | sed 's/>//' | awk -F"\t" '{split($0,a,"."); print a[1]"\t"$0}' > 13317.ls

#13317_no_probs.gff3 does not have UTRs
#16284 has the wrong exon numbering
awk -F"\t" 'BEGIN{while(getline < "13317.ls"){g[$1]=$1; t[$2]=$2}}{split($9,attr,"[=;]"); for(i in attr){if((attr[i] in g && $3 == "gene") || attr[i] in t){print $0; break}}}' 16284_introns_final.gff3  | grep -v "intron" > 13317_nextAdd.gff3
bash add_IDs_gff.sh 13317_nextAdd.gff3 > cleaned_13317_nextAdd.gff3

sort gt_cleaned_13317_no_probs.gff3 > diff1
gt gff3 -sort -tidy -retainids cleaned_13317_nextAdd.gff3 > gt_addUTRs_13317.gff3
sort gt_addUTRs_13317.gff3 > diff2
grep -v -e "UTR" gt_addUTRs_13317.gff3 | sort > diff3

#########OKAY ALL CLEAN and UTRs added#########

cp gt_addUTRs_13317.gff3 gt_cleaned_13317_no_probs.gff3

cat gt_cleaned_13317_no_probs.gff3 renamed_PP.gff | grep -v "#" > final_ShV3_9542_14856.gff3
gt gff3 -sort -tidy -retainids final_ShV3_9542_14856.gff3 > gt_final_ShV3_9542_14856.gff3 2> gt_report.log

funannotate util gff2tbl -g gt_final_ShV3_9542_14856.gff3 -f 163_wgene_hypo_genome_ShV3.softmasked.fasta > 20_05_21.tbl
funannotate util tbl2gbk -i 20_05_21.tbl -f 163_wgene_hypo_genome_ShV3.softmasked.fasta -s "Schistosoma haematobium" --isolate Egypt --sbt template.sbt -o 20_05_21_out

busco -m protein -c 24 -i 13318_prots_no_probs.fasta -o BUSCO_Shae_clean -l metazoa_odb10 --update-data 2>clean13318.err 1>clean13318.log

#Add in good ones from PP pipeline if they are BUSCOs
fa2tab ../gt_cleaned_9678_14052_PASIPASA.pep | awk -F"\t" 'BEGIN{while(getline < "1539_isos_gt.ls"){a[">"$0]=$0}}{if($1 in a){}else{print $1; print $2}}' > ../minus_1539.pep
busco -m protein -c 24 -i minus_1539.pep -o BUSCO_minus_1539 -l metazoa_odb10 --update-data 2>BUSCO_minus_1539.err 1>BUSCO_minus_1539.log

#find the ones to add but don't add if they overlap
#minus_1539
grep -v "#" full_table.tsv | grep -v "Missing" | cut -f1-3 > 1065_CDF.tsv
mv BUSCO_minus_1539/run_metazoa_odb10/1065_CDF.tsv BUSCO_fixing/
#all that have some sort of "good" BUSCO hit.

#Check how many of those aren't already represented in the final set
grep -v "#" ../BUSCO_ShV3_20_05_21/run_metazoa_odb10/full_table.tsv | grep -v "Missing" | cut -f1-3 > 917_CDF_Sha.tsv
awk -F"\t" 'BEGIN{while(getline < "917_CDF_Sha.tsv"){sha[$1]=$2; anno[$1]="\t"$2"\t"$3}}{if($1 in sha){}else{print $0}}' 1065_CDF.tsv > 27_missing_Sha.tsv
cut -f3 27_missing_Sha.tsv | sed 's/\./\t/' | cut -f1 | sort | uniq > 19_genes_of_missing.tsv
#Add those genes to the final
cut -f3 27_missing_Sha.tsv | awk -F"\t" '{split($1,iso,"."); print iso[1]"\t"$1}' | awk -F"\t" '{if(!g[$1]){g[$1]=20000+len;} len=0; for(i in g){++len}; ++t[$1]; print $0"\tSha_"g[$1]"\tSha_"g[$1]"."t[$1]}' > add27.lkp

grep -v "#" ../mergePP_Sha/gt_cleaned_9678_14052_PASIPASA.gff3 | awk -F"\t" '
BEGIN{while(getline < "add27.lkp"){g[$1]=$3; t[$2]=$4}}
{
  split($9,attr,"[=;]");
  $2="ShaeV3";
  if($3 == "gene" && attr[2] in g){fnd=1; attr[2]=g[attr[2]]; a=2};
  if($3 == "mRNA" && attr[2] in t){fnd=1; attr[2]=t[attr[2]]; attr[4]=g[attr[4]]; a=4};
  if($3 == "CDS" && attr[4] in t){fnd=1; attr[2]="cds."t[attr[4]]; attr[4]=t[attr[4]]; if(attr[6]=="true"){a=6}else{a=4}};
  if($3 ~ /prime_UTR/ && attr[4] in t){fnd=1; split(attr[2],u,".utr"); attr[2]=t[u[1]]".utr"u[2]; attr[4]=t[attr[4]]; a=4};
  if($3 == "exon" && attr[4] in t){fnd=1; split(attr[2],e,".exon"); attr[2]=t[e[1]]".exon"e[2]; attr[4]=t[attr[4]]; a=4};
  if(fnd == 1){
    fnd=0;
    for(i=1; i<=8; i++){printf $i"\t"};
    printf attr[1]"="attr[2];
    for(k=3; k<a; k=k+2){printf ";"attr[k]"="attr[k+1]};
    print ""}
}' > add27_PP.gff

awk -F"\t" 'BEGIN{while(getline < "917_CDF_Sha.tsv"){sha[$1]=$2; anno[$1]="\t"$2"\t"$3}}{if($1 in sha && $2 != "Fragmented" && sha[$1]=="Fragmented"){print $0""anno[$1]}}' 1065_CDF.tsv > 15_fragmented_replace.tsv
#check if fragmented ones are already represented
cut -f1 15_fragmented_replace.tsv | sort | uniq | while read line; do grep -w $line ../BUSCO_ShV3_20_05_21/run_metazoa_odb10/full_table.tsv; done
#the only ones with that BUSCO ID are the fragmented ones
#check for the isoforms from gt that match that Sha ID
cut -f3 15_fragmented_replace.tsv | sed 's/\.t1/\t/' | cut -f1 | sort | uniq | while read line; do grep -w $line add27.lkp; done
cut -f3 15_fragmented_replace.tsv | sed 's/\.t1/\t/' | cut -f1 | sort | uniq | while read line; do grep -w $line ../mergePP_Sha/merged_rename_g.t_Sha.tsv; done
cut -f5 15_fragmented_replace.tsv | sed 's/\./\t/' | cut -f1 | sort | uniq | while read line; do grep -w $line ../BUSCO_ShV3_20_05_21/run_metazoa_odb10/full_table.tsv; done
cut -f5 15_fragmented_replace.tsv | sort | uniq | while read line; do grep -w $line ../BUSCO_ShV3_20_05_21/run_metazoa_odb10/full_table.tsv; done 
#for all 12 Sha genes there is only one isoform that has a Fragmented BUSCO match no others
#replace those entire genes by exactly the ones that do have better (D/C) BUSCO scores (include all isoforms)
#remove 12:
cut -f5 15_fragmented_replace.tsv | sed 's/\./\t/' | cut -f1 | sort | uniq > 12_remove.ls
grep -v -w -f 12_remove.ls gt_final_ShV3_9542_14856.gff3 > 12_removed.gff
cut -f3,5 15_fragmented_replace.tsv | sed 's/\./\t/g' | awk -F"\t" '{print $1"\t"$(NF-1)}' | sort | uniq > name_match_13_missing.lkp

#manually change g9710   Sha_04710 to g9710   Sha_04714 in name_match_13_missing.lkp
#Extract 13 genes
#g2858_g2859 became g2858
#fix that in the lkp table
cut -f1 name_match_13_missing.lkp | while read line; do grep -w $line gt_cleaned_9678_14052_PASIPASA.gff3; done | grep -e "mRNA" | cut -f9 | sed 's/[;=]/\t/g' | cut -f2,4 | awk -F"\t" '{print $2"\t"$1}' > 13_genes_isos_add.tsv
awk -F"\t" 'BEGIN{while(getline < "name_match_13_missing.lkp"){mt[$1]=$2}}{print $0"\t"mt[$1]"\t"mt[$1]"."++g[$1]}' 13_genes_isos_add.tsv > add_13.lkp

grep -v "#" ../mergePP_Sha/gt_cleaned_9678_14052_PASIPASA.gff3 | awk -F"\t" '
BEGIN{while(getline < "add_13.lkp"){g[$1]=$3; t[$2]=$4}}
{
  split($9,attr,"[=;]");
  $2="ShaeV3";
  if($3 == "gene" && attr[2] in g){fnd=1; attr[2]=g[attr[2]]; a=2};
  if($3 == "mRNA" && attr[2] in t){fnd=1; attr[2]=t[attr[2]]; attr[4]=g[attr[4]]; a=4};
  if($3 == "CDS" && attr[4] in t){fnd=1; attr[2]="cds."t[attr[4]]; attr[4]=t[attr[4]]; if(attr[6]=="true"){a=6}else{a=4}};
  if($3 ~ /prime_UTR/ && attr[4] in t){fnd=1; split(attr[2],u,".utr"); attr[2]=t[u[1]]".utr"u[2]; attr[4]=t[attr[4]]; a=4};
  if($3 == "exon" && attr[4] in t){fnd=1; split(attr[2],e,".exon"); attr[2]=t[e[1]]".exon"e[2]; attr[4]=t[attr[4]]; a=4};
  if(fnd == 1){
    fnd=0;
    for(i=1; i<=8; i++){printf $i"\t"};
    printf attr[1]"="attr[2];
    for(k=3; k<a; k=k+2){printf ";"attr[k]"="attr[k+1]};
    print ""}
}' > add13_PP.gff

cat 12_removed.gff add27_PP.gff add13_PP.gff | grep -v "#" > 9562_14885_no_comm.gff3

#run the "make longest gene ID" script. 
#Make the gene line coordinates match up with the coordinates of the longest isoform that is present (we may have removed isoforms but haven't adjusted the coordinates in the gene line)
awk -F"\t" 'NF > 1' 9562_14885_no_comm.gff3 | cut -f1,3-5 | grep -v -e "intron" -e "mRNA" -e "gene" | awk -F"\t" '{scaffs[$1]=$1; if($2 ~ /UTR/ || $2 == "CDS"){for(n=$3; n <= $4; n++){cover[$1"@"n]=1}}; if($2 == "exon"){for(k=$3; k <= $4; k++){exon[$1"@"k]=1}}}END{for(i in exon){split(i,a,"@"); if(!(i in cover)){print a[1]"\t"a[2]}}}' > bases_not_covered.tsv
sort -t$'\t' -k1,1 -k2,2n bases_not_covered.tsv > tmp; mv tmp bases_not_covered.tsv

#let's fix the bases that are not covered
awk -F"\t" 'BEGIN{while(getline < "bases_not_covered.tsv"){lkp[$1"@"$2]=1}}{if($3 == "exon"){for(n=$4; n<=$5; n++){if($1"@"n in lkp){if(n-$4 < $5-n){$4=n+1}else{$5=n-1}}};  printf $1; for(i=2; i<=NF; i++){printf "\t"$i}; print ""}else{print $0}}' 9562_14885_no_comm.gff3 > fixed_uncovered_exons.gff

#check all genes and mRNAs are still there
#should be 9562 genes and 14885 mRNAs - OK
awk -F"\t" '{if($3 == "exon"){split($9,a,";Parent="); split(a[2],g,"."); if(min_iso[a[2]] > $4 || min_iso[a[2]] == 0){min_iso[a[2]]=$4}; if(max_iso[a[2]] < $5){max_iso[a[2]]=$5}; if(min_gene[g[1]] > $4  || min_gene[g[1]] == 0){min_gene[g[1]]=$4}; if(max_gene[g[1]] < $5){max_gene[g[1]]=$5};}}END{for(i in min_iso){print i"\t"min_iso[i]"\t"max_iso[i]}}' fixed_uncovered_exons.gff > min_max_iso_pos.tsv
awk -F"\t" '{if($3 == "exon"){split($9,a,";Parent="); split(a[2],g,"."); if(min_iso[a[2]] > $4 || min_iso[a[2]] == 0){min_iso[a[2]]=$4}; if(max_iso[a[2]] < $5){max_iso[a[2]]=$5}; if(min_gene[g[1]] > $4  || min_gene[g[1]] == 0){min_gene[g[1]]=$4}; if(max_gene[g[1]] < $5){max_gene[g[1]]=$5};}}END{for(i in min_gene){print i"\t"min_gene[i]"\t"max_gene[i]}}' fixed_uncovered_exons.gff > min_max_gene_pos.tsv
awk -F"\t" 'BEGIN{while(getline < "min_max_iso_pos.tsv"){min_iso[$1]=$2; max_iso[$1]=$3};while(getline < "min_max_gene_pos.tsv"){min_gene[$1]=$2; max_gene[$1]=$3};}{split($9,g,"[;=]");if($3 == "mRNA"){if($4 < min_iso[g[2]]){$4=min_iso[g[2]]};if($5 > max_iso[g[2]]){$5=max_iso[g[2]]}};if($3 == "gene"){if($4 < min_gene[g[2]]){$4=min_gene[g[2]]};if($5 > max_gene[g[2]]){$5=max_gene[g[2]]}};printf $1; for(i=2; i<=NF; i++){printf "\t"$i};print ""}' fixed_uncovered_exons.gff > fixed_gene_mRNA.gff
gt gff3 -sort -tidy -retainids fixed_gene_mRNA.gff > gt_9562_14885_no_comm.gff3 2> gt_report.log

#run trmap to find overlapping NOT belonging to the same gene
grep -v gene gt_9562_14885_no_comm.gff3 > tmp1
trmap -o check_overlaps.out tmp1 tmp1
awk -F"[\t ]" '{if($1 ~ /^>/){q=$1}else{l[q]=l[q]"\t"$6}}END{for(i in l){print i""l[i]} }' check_overlaps.out | sed 's/>//' | awk -F"\t" '{split($1,g,"."); for(i=2; i<=NF; i++){split($NF,a,"."); if(a[1] != g[1]){print $0; break}}}' > overlaps.tsv
sed 's/\t/\n/g' overlaps.tsv | sed 's/\./\t/g' | cut -f1 | sort | uniq > overlapping_genes.ls
#381

grep -v ">" check_overlaps.out | cut -f1 | sort | uniq -c
  # 14935 =
  #    43 c
  #    15 e
  #    67 i
  # 10958 j
  #    43 k
  #  1680 m
  #  2941 n
  #   226 o
  #    93 y

awk -F"\t" '$2 == "remove"' overlapping_genes.ls | cut -f1 > 132_remove_overlapping.ls
grep -v -f 132_remove_overlapping.ls mf_gt_9562_14885_no_comm.gff3 > nr_9431_14700.gff3
grep -v -e "#" -e "gene" nr_9431_14700.gff3 > tmp2
trmap -o nr_trmap.out tmp2 tmp2
awk -F"[\t ]" '{if($1 ~ /^>/){q=$1}else{l[q]=l[q]"\t"$6}}END{for(i in l){print i""l[i]} }' nr_trmap.out | sed 's/>//' | awk -F"\t" '{split($1,g,"."); for(i=2; i<=NF; i++){split($NF,a,"."); if(a[1] != g[1]){print $0; break}}}' > remaining_76_overlap.ls

#funannotate tbl2asn check
funannotate util gff2tbl -g nr_9431_14700.gff3 -f 163_wgene_hypo_genome_ShV3.softmasked.fasta > nr_9431_14700.tbl
funannotate util tbl2gbk -i nr_9431_14700.tbl -f 163_wgene_hypo_genome_ShV3.softmasked.fasta -s "Schistosoma haematobium" --isolate Egypt --sbt template.sbt -o nr_9431_14700_out


#gene/protein set stats
#macropus
#/home/andreas/ShV3_anno/final_ShV3/

#genes/mRNA
awk -F"\t" 'NF == 9' Shae.V3_final.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | cut -f1 | sort | uniq -c

#UTRs
awk -F"\t" 'NF == 9' Shae.V3_final.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep prime | cut -f1,4 | sort | uniq | cut -f1 | sort | uniq -c

#gene length
awk -F"\t" 'NF == 9' Shae.V3_final.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep gene | awk -F"\t" '{print $3-$2+1}' > gene_lengths.tsv

#CDS
awk -F"\t" 'NF == 9' Shae.V3_final.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep CDS | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > CDS_lengths.tsv

#exon
awk -F"\t" 'NF == 9' Shae.V3_final.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep exon | awk -F"\t" '{print $3-$2+1}' > exon_lengths.tsv

#mRNA (length of all exons combined per mRNA)
awk -F"\t" 'NF == 9' Shae.V3_final.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep exon | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > mRNA_lengths.tsv

#prot lengths
fa2len final_07_06_21_prot.fasta | cut -f2 > prot_lengths.tsv



#schistosoma_haematobium.PRJNA78265.WBPS15.annotations.gff3

#genes/mRNA
awk -F"\t" 'NF == 9' schistosoma_haematobium.PRJNA78265.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | cut -f1 | sort | uniq -c

#UTRs
awk -F"\t" 'NF == 9' schistosoma_haematobium.PRJNA78265.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=.\+://' | grep prime | sed 's/Parent=transcript://' | cut -f1,4 | sort | uniq | cut -f1 | sort | uniq -c

#gene length
awk -F"\t" 'NF == 9' schistosoma_haematobium.PRJNA78265.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep gene | awk -F"\t" '{print $3-$2+1}' > gene_lengths.tsv

#CDS
awk -F"\t" 'NF == 9' schistosoma_haematobium.PRJNA78265.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep CDS | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > CDS_lengths.tsv

#exon
awk -F"\t" 'NF == 9' schistosoma_haematobium.PRJNA78265.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep exon | awk -F"\t" '{print $3-$2+1}' > exon_lengths.tsv

#mRNA (length of all exons combined per mRNA)
awk -F"\t" 'NF == 9' schistosoma_haematobium.PRJNA78265.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep exon | sed 's/\.[^\.]\+$//' | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > mRNA_lengths.tsv

#prot lengths
fa2len Shae_prots.fasta | cut -f2 > ../final_ShV3/shaev2/prot_lengths.tsv


#schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3

#genes/mRNA
awk -F"\t" 'NF == 9' schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | cut -f1 | sort | uniq -c

#UTRs
awk -F"\t" 'NF == 9' schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=.\+://' | grep prime | sed 's/Parent=transcript://' | cut -f1,4 | sort | uniq | cut -f1 | sort | uniq -c

#gene length
awk -F"\t" 'NF == 9' schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep -w gene | awk -F"\t" '{print $3-$2+1}' > gene_lengths.tsv

#CDS
awk -F"\t" 'NF == 9' schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep CDS | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > CDS_lengths.tsv

#exon
awk -F"\t" 'NF == 9' schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep exon | awk -F"\t" '{print $3-$2+1}' > exon_lengths.tsv

#mRNA (length of all exons combined per mRNA)
awk -F"\t" 'NF == 9' schistosoma_mansoni.PRJEA36577.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep -w exon | sed 's/\.[^\.]\+$//' | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > mRNA_lengths.tsv

#prot lengths
fa2len schistosoma_mansoni.PRJEA36577.WBPS15.protein.fa | cut -f2 > ../final_ShV3/sman/prot_lengths.tsv


#schistosoma_bovis.PRJNA451066.WBPS15.annotations.gff3

#genes/mRNA
awk -F"\t" 'NF == 9' schistosoma_bovis.PRJNA451066.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | cut -f1 | sort | uniq -c

#gene length
awk -F"\t" 'NF == 9' schistosoma_bovis.PRJNA451066.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep -w mRNA | awk -F"\t" '{print $3-$2+1}' > gene_lengths.tsv

#CDS
awk -F"\t" 'NF == 9' schistosoma_bovis.PRJNA451066.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep CDS | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > CDS_lengths.tsv

#exon
awk -F"\t" 'NF == 9' schistosoma_bovis.PRJNA451066.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep exon | awk -F"\t" '{print $3-$2+1}' > exon_lengths.tsv

#mRNA (length of all exons combined per mRNA)
awk -F"\t" 'NF == 9' schistosoma_bovis.PRJNA451066.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep -w exon | sed 's/\.[^\.]\+$//' | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > mRNA_lengths.tsv

#prot lengths
fa2len schistosoma_bovis.PRJNA451066.WBPS15.protein.fa | cut -f2 > ../final_ShV3/sbovis/prot_lengths.tsv


#schistosoma_japonicum.PRJNA520774.WBPS15.annotations.gff3

#genes/mRNA
grep "gene" schistosoma_japonicum.PRJNA520774.WBPS15.annotations.gff3 | grep "biotype=protein_coding" | wc -l
#10089

#UTRs
awk -F"\t" 'NF == 9' schistosoma_japonicum.PRJNA520774.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=.\+://' | grep prime | sed 's/Parent=transcript://' | cut -f1,4 | sort | uniq | cut -f1 | sort | uniq -c

#gene length
grep "gene" schistosoma_japonicum.PRJNA520774.WBPS15.annotations.gff3 | grep "biotype=protein_coding" | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | awk -F"\t" '{print $3-$2+1}' > gene_lengths.tsv

#CDS
awk -F"\t" 'NF == 9' schistosoma_japonicum.PRJNA520774.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep CDS | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > CDS_lengths.tsv

#exon
awk -F"\t" 'NF == 9' schistosoma_japonicum.PRJNA520774.WBPS15.annotations.gff3 | cut -f3,4,5,9 | sed 's/\.[^[:digit:]]/\t/' | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | grep exon | awk -F"\t" '{print $3-$2+1}' > exon_lengths.tsv

#mRNA (length of all exons combined per mRNA)
awk -F"\t" 'NF == 9' schistosoma_japonicum.PRJNA520774.WBPS15.annotations.gff3 | cut -f3,4,5,9  | awk -F"\t" '{if($1 == "mRNA"){split($4,attr,"[:;=]"); mrna[attr[3]]=attr[3]}; if($1 == "exon"){split($4,exattr,"[:;=]"); if(exattr[6] in mrna){print $0}}}' | sed 's/;/\t/' | cut -f1-3,5 | awk -F"\t" '{len[$4]=len[$4]+$3-$2+1}END{for(i in len){print len[i]}}' > mRNA_lengths.tsv

#prot lengths
fa2len schistosoma_japonicum.PRJNA520774.WBPS15.protein.fa | cut -f2 > ../final_ShV3/sjap/prot_lengths.tsv