###RUNNING REPEATMASKER TO MASK THE GENOME###
export PERL5LIB=/home/andreas/miniconda2/envs/repmod/share/RepeatMasker

RepeatMasker -dir ./customShae -pa 48 -engine ncbi -excln -gccalc -s -no_is -gff -nolow -lib final.library hypo_genome_ShV3.fa

rmOutToGFF3.pl hypo_genome_ShV3.fa.out > hypo_genome_ShV3.fa.out.gff3

#bedtools-2.30.0
bedtools sort -i hypo_genome_ShV3.fa.out.gff3 | bedtools merge -d 1 > hypo_genome_ShV3.fa.masked.gff3
bedtools maskfasta -soft -fi ../hypo_genome_ShV3.fa -bed hypo_genome_ShV3.fa.masked.gff3 -fo hypo_genome_ShV3.softmasked.fasta