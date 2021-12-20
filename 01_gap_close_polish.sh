##Data for assembly and polishing
NUMTH=24
DRAFT=SciHae_3.0_HiC.fasta

#R1 and R2 for 500,800 libraries

R1_500=500_paired_1.fq
R2_500=500_paired_2.fq

R1_800=800_paired_1.fq
R2_800=800_paired_2.fq

RTYPE=map-ont
LONGR=shaem.trimmedReads.fasta
 
source /data/miniconda3/etc/profile.d/conda.sh
#conda create -n tgsgc -c bioconda tgsgapcloser
conda activate tgsgc
echo "Running TGSGAPCLOSER $(date)"
tgsgapcloser --tgstype ont --ne --min_match 1000 --scaff $DRAFT --reads $LONGR --thread $NUMTH --output tgsgc_gapclose
conda deactivate
echo "Finished running TGSGAPCLOSER $(date)"


TGS_DRAFT=$(ls tgsgc_gapclose.scaff_seq*)
echo "TGSDRAFT: $TGS_DRAFT"

gunzip $TGS_DRAFT > /dev/null 2>&1

TGS_FASTA=$(ls tgsgc_gapclose.scaff_seq* | grep -v "gz")
echo "TGSFASTA: $TGS_FASTA"

conda activate busco
busco -m genome -i $TGS_FASTA -o TGS_BUSCO -l metazoa_odb10 -c $NUMTH -f --update-data
conda deactivate
echo "Finished running TGS BUSCO $(date)"
 
#Polishing with long and short reads
conda activate minimap2
minimap2 --secondary=no --MD -ax sr -t $NUMTH $TGS_FASTA $R1_500 $R2_500 | samtools view -Sb - > mapped-sr_500.bam
samtools sort -@$NUMTH -o mapped-sr_500.sorted.bam mapped-sr_500.bam
samtools index mapped-sr_500.sorted.bam
rm mapped-sr_500.bam
echo "Finished mapping 500 $(date)"


minimap2 --secondary=no --MD -ax sr -t $NUMTH $TGS_FASTA $R1_800 $R2_800 | samtools view -Sb - > mapped-sr_800.bam
samtools sort -@$NUMTH -o mapped-sr_800.sorted.bam mapped-sr_800.bam
samtools index mapped-sr_800.sorted.bam
rm mapped-sr_800.bam
echo "Finished mapping 800 $(date)"


echo "Merging mapped short-read bam files $(date)"
samtools merge merged_500_800.sorted.bam mapped-sr_500.sorted.bam mapped-sr_800.sorted.bam
samtools index merged_500_800.sorted.bam
echo "Merged mapped short-read bam files $(date)"

##Mapping the long reads to contigs

minimap2 --secondary=no --MD -ax $RTYPE -t $NUMTH $TGS_FASTA $LONGR | samtools view -Sb - > mapped-lg.bam
samtools sort -@$NUMTH -o mapped-lg.sorted.bam mapped-lg.bam
samtools index mapped-lg.sorted.bam
rm mapped-lg.bam
conda deactivate
echo "Finished mapping long reads $(date)"

##Running Hypo:
#conda create -n hypo -c bioconda hypo

##Create a text file containing the names of the short reads files.

echo -e "$R1_500\n$R2_500\n$R1_800\n$R2_800" > il_names.txt

## Coverage estimates from mapping 500 bp and 800 bp lib to gap closed genome

echo "Estimating mean coverage for mapped short-reads $(date)"
conda activate bedtools
bedtools genomecov -ibam merged_500_800.sorted.bam -g $TGS_FASTA -d > coverage_estimates.tsv
COV=$(cut -f3 coverage_estimates.tsv | awk '{sum=sum+$0}END{print sum/NR}')
echo "Mean coverage: $COV"
conda deactivate
echo "Done estimating mean coverage for mapped short-reads $(date)"

##run Hypo (for short reads as well as long reads polishing)
echo "Running hypo $(date)"
conda activate hypo
hypo -d $TGS_FASTA -r @il_names.txt -s 0.38g -c $COV -b merged_500_800.sorted.bam -B mapped-lg.sorted.bam -t $NUMTH -o hypo_genome_ShV3.fa
conda deactivate
echo "Done hypo $(date)"

echo "Running BUSCO HYPO $(date)"
conda activate busco
busco -m genome -i hypo_genome_ShV3.fa -o TGS_HYPO -l metazoa_odb10 -c $NUMTH -f --update-data
conda deactivate
echo "Done BUSCO HYPO $(date)"