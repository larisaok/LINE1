# 14.01.2022
Discussion with Gael:
* decided to use Gencode versions
* use only canoical chromosomes
* how to define gene and filter gtf/gff file according to this strategy?
* organise the control actions - Laurine took the motifs of fasta and made tha logo
* check outliers
* include orientation of the motif in the analysis


# 17.01.2022
Genes - how they are annotated in the human genome
https://genome.ucsc.edu/FAQ/FAQgenes.html
gene model tables in gencode human genome: ncbiRefSeq, refGene, ensGene, knownGene

Ensembl and GENCODE gene models are the same
ensGene are very comfortable to use 
* possible use the otrhologs from the ensemble (were identified using the synteny approach)
* same gene model for all species, similar to the gencode annotation
* wider, than refseq models

## downloading genome

hg38.fa.gz - "Soft-masked" assembly sequence in one file.

rsync options 
--verbose, -v            increase verbosity
--archive, -a            archive mode; equals -rlptgoD (no -H,-A,-X)
--compress, -z           compress file data during the transfer
-P     The -P option is equivalent to --partial --progress.  Its
       purpose is to make it much easier to specify these two
       options for a long transfer that may be interrupted.

```
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz .
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.gc5Base.bw .
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz .
```
## removing alt loci from fasta file
```
lokorokova@professorx:~/INTRANET/cristo-nas-shared/Personal_documents/Larisa/genomes/gencode_genomes> zgrep ">" hg38.fa.gz | head
```

>chr1
>chr10
>chr11
>chr11_KI270721v1_random
>chr12
>chr13
>chr14
>chr14_GL000009v2_random
>chr14_GL000225v1_random
>chr14_KI270722v1_random
>chr14_GL000194v1_random

```
zgrep ">" hg38.fa.gz | wc -l
```
`455`
```
zgrep ">" hg38.fa.gz | grep "alt" | wc -l
```
`261`
```
zgrep ">" hg38.fa.gz | grep "chrUn" | wc -l
```
`127`
```
zgrep ">" hg38.fa.gz | grep "random" | wc -l
```
`42`
```
gunzip hg38.fa.gz
samtools faidx hg38.fa
cut -f 1 hg38.fa.fai | grep -v "alt" | wc -l
```
`194`
```
cut -f 1 hg38.fa.fai | grep -v "alt" > selected_human_chr.list
while read chr
do
#echo $chr
samtools faidx hg38.fa $chr >> hg38_selected_chr.fa
done < selected_human_chr.list
```
## run homer on human genome
```
/home/lokorokova/bin/homer/bin/scanMotifGenomeWide.pl \
    /home/lokorokova/LINE1_project/line1.motif hg38_selected_chr.fa -bed > \
    homer_results_human.bed &
    
# groom homer results
grep -v "total" homer_results_human.bed | grep "^chr" | awk -F "\t" -v OFS="\t" '{print $1,$2-1,$3,$4,$5,$6}' > homer_results_human_corrected.bed
pigz homer_results_human.bed
```

```
head -n 20 homer_results_human.bed
```
`Using Custom Genome
Processing chr1
Reading input files...
24996 total sequences read
1 motifs loaded
Finding instances of 1 motif(s)
|0%                                    50%                                  100%|
    =================================================================================
chr1    11532    11542        7.188263    +
chr1    11538    11548        7.164812    +
chr1    11725    11735        9.443662    -
chr1    11887    11897        5.185282    -
chr1    13907    13917        6.548140    +
chr1    14163    14173        5.933078    -
chr1    14380    14390        4.966773    -
chr1    15466    15476        6.758355    +
chr1    16315    16325        5.276019    +
chr1    16457    16467        5.193864    -`

```
head -n 20 homer_results_human_corrected.bed 
```
`chr1    11531    11542        7.188263    +
chr1    11537    11548        7.164812    +
chr1    11724    11735        9.443662    -
chr1    11886    11897        5.185282    -
chr1    13906    13917        6.548140    +
chr1    14162    14173        5.933078    -
chr1    14379    14390        4.966773    -
chr1    15465    15476        6.758355    +
chr1    16314    16325        5.276019    +
chr1    16456    16467        5.193864    -
chr1    16502    16513        7.815973    +
chr1    19732    19743        5.142425    +
chr1    20114    20125        5.046636    +
chr1    20310    20321        5.569618    +
chr1    20428    20439        5.112428    +
chr1    20495    20506        7.714726    -
chr1    20576    20587        7.649706    +
chr1    22371    22382        5.143932    -
chr1    22548    22559        5.611923    -
chr1    23195    23206        7.641769    +`

```
wc -l homer_results_human_corrected.bed 
```
`28700814 homer_results_human_corrected.bed`

## L1 motif quantification

```
zcat hg38.ensGene.gtf.gz | head
```
`chr1    ensGene    transcript    11869    14409    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328";  gene_name "ENSG00000223972";
chr1    ensGene    exon    11869    12227    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1"; exon_id "ENST00000456328.1"; gene_name "ENSG00000223972";
chr1    ensGene    exon    12613    12721    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "2"; exon_id "ENST00000456328.2"; gene_name "ENSG00000223972";
chr1    ensGene    exon    13221    14409    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "3"; exon_id "ENST00000456328.3"; gene_name "ENSG00000223972";
chr1    ensGene    transcript    12010    13670    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305";  gene_name "ENSG00000223972";
chr1    ensGene    exon    12010    12057    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305"; exon_number "1"; exon_id "ENST00000450305.1"; gene_name "ENSG00000223972";
chr1    ensGene    exon    12179    12227    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305"; exon_number "2"; exon_id "ENST00000450305.2"; gene_name "ENSG00000223972";
chr1    ensGene    exon    12613    12697    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305"; exon_number "3"; exon_id "ENST00000450305.3"; gene_name "ENSG00000223972";
chr1    ensGene    exon    12975    13052    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305"; exon_number "4"; exon_id "ENST00000450305.4"; gene_name "ENSG00000223972";
chr1    ensGene    exon    13221    13374    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305"; exon_number "5"; exon_id "ENST00000450305.5"; gene_name "ENSG00000223972";`

## features in the gtf file 
```
zcat hg38.ensGene.gtf.gz | cut -f 3 | sort -u
```
`3UTR
5UTR
CDS
exon
start_codon
stop_codon
transcript`

```
bedtools sort -i hg38_only_transcripts.gtf > sorted_hg38_only_transcripts.gtf
bedtools sort -i homer_results_human_corrected.bed > sorted_homer_results_homer_corrected.bed
bedtools map -a sorted_hg38_only_transcripts.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -s > count_motif_same_strand.bed
bedtools map -a sorted_hg38_only_transcripts.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -S > count_motif_opposite_strand.bed
```
```head count_motif_same_strand.bed
```
`chr1    ensGene    transcript    11869    14409    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328";  gene_name "ENSG00000223972";    1
chr1    ensGene    transcript    12010    13670    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305";  gene_name "ENSG00000223972";    0
chr1    ensGene    transcript    14404    29570    .    -    .    gene_id "ENSG00000227232"; transcript_id "ENST00000488147";  gene_name "ENSG00000227232";    23
chr1    ensGene    transcript    17369    17436    .    -    .    gene_id "ENSG00000278267"; transcript_id "ENST00000619216";  gene_name "ENSG00000278267";    0
chr1    ensGene    transcript    29554    31097    .    +    .    gene_id "ENSG00000243485"; transcript_id "ENST00000473358";  gene_name "ENSG00000243485";    6
chr1    ensGene    transcript    30267    31109    .    +    .    gene_id "ENSG00000243485"; transcript_id "ENST00000469289";  gene_name "ENSG00000243485";    2
chr1    ensGene    transcript    30366    30503    .    +    .    gene_id "ENSG00000274890"; transcript_id "ENST00000607096";  gene_name "ENSG00000274890";    1
chr1    ensGene    transcript    34554    36081    .    -    .    gene_id "ENSG00000237613"; transcript_id "ENST00000417324";  gene_name "ENSG00000237613";    2
chr1    ensGene    transcript    35245    36073    .    -    .    gene_id "ENSG00000237613"; transcript_id "ENST00000461467";  gene_name "ENSG00000237613";    0
chr1    ensGene    transcript    52473    53312    .    +    .    gene_id "ENSG00000268020"; transcript_id "ENST00000606857";  gene_name "ENSG00000268020";    3`

```head count_motif_opposite_strand.bed
```
`chr1    ensGene    transcript    11869    14409    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328";  gene_name "ENSG00000223972";    3
chr1    ensGene    transcript    12010    13670    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305";  gene_name "ENSG00000223972";    0
chr1    ensGene    transcript    14404    29570    .    -    .    gene_id "ENSG00000227232"; transcript_id "ENST00000488147";  gene_name "ENSG00000227232";    22
chr1    ensGene    transcript    17369    17436    .    -    .    gene_id "ENSG00000278267"; transcript_id "ENST00000619216";  gene_name "ENSG00000278267";    0
chr1    ensGene    transcript    29554    31097    .    +    .    gene_id "ENSG00000243485"; transcript_id "ENST00000473358";  gene_name "ENSG00000243485";    3
chr1    ensGene    transcript    30267    31109    .    +    .    gene_id "ENSG00000243485"; transcript_id "ENST00000469289";  gene_name "ENSG00000243485";    2
chr1    ensGene    transcript    30366    30503    .    +    .    gene_id "ENSG00000274890"; transcript_id "ENST00000607096";  gene_name "ENSG00000274890";    1
chr1    ensGene    transcript    34554    36081    .    -    .    gene_id "ENSG00000237613"; transcript_id "ENST00000417324";  gene_name "ENSG00000237613";    7
chr1    ensGene    transcript    35245    36073    .    -    .    gene_id "ENSG00000237613"; transcript_id "ENST00000461467";  gene_name "ENSG00000237613";    5
chr1    ensGene    transcript    52473    53312    .    +    .    gene_id "ENSG00000268020"; transcript_id "ENST00000606857";  gene_name "ENSG00000268020";    4`

```bedtools map -a sorted_hg38_only_transcripts.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count> count_motif_both_strand.bed
head count_motif_both_strand.bed
```
`chr1    ensGene    transcript    11869    14409    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328";  gene_name "ENSG00000223972";    4
chr1    ensGene    transcript    12010    13670    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305";  gene_name "ENSG00000223972";    0
chr1    ensGene    transcript    14404    29570    .    -    .    gene_id "ENSG00000227232"; transcript_id "ENST00000488147";  gene_name "ENSG00000227232";    45
chr1    ensGene    transcript    17369    17436    .    -    .    gene_id "ENSG00000278267"; transcript_id "ENST00000619216";  gene_name "ENSG00000278267";    0
chr1    ensGene    transcript    29554    31097    .    +    .    gene_id "ENSG00000243485"; transcript_id "ENST00000473358";  gene_name "ENSG00000243485";    9
chr1    ensGene    transcript    30267    31109    .    +    .    gene_id "ENSG00000243485"; transcript_id "ENST00000469289";  gene_name "ENSG00000243485";    4
chr1    ensGene    transcript    30366    30503    .    +    .    gene_id "ENSG00000274890"; transcript_id "ENST00000607096";  gene_name "ENSG00000274890";    2
chr1    ensGene    transcript    34554    36081    .    -    .    gene_id "ENSG00000237613"; transcript_id "ENST00000417324";  gene_name "ENSG00000237613";    9
chr1    ensGene    transcript    35245    36073    .    -    .    gene_id "ENSG00000237613"; transcript_id "ENST00000461467";  gene_name "ENSG00000237613";    5
chr1    ensGene    transcript    52473    53312    .    +    .    gene_id "ENSG00000268020"; transcript_id "ENST00000606857";  gene_name "ENSG00000268020";    7`

```
paste count_motif_same_strand.bed count_motif_opposite_strand.bed | cut -f 1,2,3,4,5,6,7,8,9,10,20 > count_motif_human.bed
head count_motif_human.bed
```
`chr1    ensGene    transcript    11869    14409    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328";  gene_name "ENSG00000223972";    1    3
chr1    ensGene    transcript    12010    13670    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000450305";  gene_name "ENSG00000223972";    0    0
chr1    ensGene    transcript    14404    29570    .    -    .    gene_id "ENSG00000227232"; transcript_id "ENST00000488147";  gene_name "ENSG00000227232";    23    22
chr1    ensGene    transcript    17369    17436    .    -    .    gene_id "ENSG00000278267"; transcript_id "ENST00000619216";  gene_name "ENSG00000278267";    0    0
chr1    ensGene    transcript    29554    31097    .    +    .    gene_id "ENSG00000243485"; transcript_id "ENST00000473358";  gene_name "ENSG00000243485";    6    3
chr1    ensGene    transcript    30267    31109    .    +    .    gene_id "ENSG00000243485"; transcript_id "ENST00000469289";  gene_name "ENSG00000243485";    2    2
chr1    ensGene    transcript    30366    30503    .    +    .    gene_id "ENSG00000274890"; transcript_id "ENST00000607096";  gene_name "ENSG00000274890";    1    1
chr1    ensGene    transcript    34554    36081    .    -    .    gene_id "ENSG00000237613"; transcript_id "ENST00000417324";  gene_name "ENSG00000237613";    2    7
chr1    ensGene    transcript    35245    36073    .    -    .    gene_id "ENSG00000237613"; transcript_id "ENST00000461467";  gene_name "ENSG00000237613";    0    5
chr1    ensGene    transcript    52473    53312    .    +    .    gene_id "ENSG00000268020"; transcript_id "ENST00000606857";  gene_name "ENSG00000268020";    3    4`


# 18/01/22
Starting from the beginning using ensemble reference files

```
wget http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz
wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 
```
```
grep ">" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa | head -n 30
```
>1 dna_sm:chromosome chromosome:GRCh38:1:1:248956422:1 REF
>10 dna_sm:chromosome chromosome:GRCh38:10:1:133797422:1 REF
>11 dna_sm:chromosome chromosome:GRCh38:11:1:135086622:1 REF
>12 dna_sm:chromosome chromosome:GRCh38:12:1:133275309:1 REF
>13 dna_sm:chromosome chromosome:GRCh38:13:1:114364328:1 REF
>14 dna_sm:chromosome chromosome:GRCh38:14:1:107043718:1 REF
>15 dna_sm:chromosome chromosome:GRCh38:15:1:101991189:1 REF
>16 dna_sm:chromosome chromosome:GRCh38:16:1:90338345:1 REF
>17 dna_sm:chromosome chromosome:GRCh38:17:1:83257441:1 REF
>18 dna_sm:chromosome chromosome:GRCh38:18:1:80373285:1 REF
>19 dna_sm:chromosome chromosome:GRCh38:19:1:58617616:1 REF
>2 dna_sm:chromosome chromosome:GRCh38:2:1:242193529:1 REF
>20 dna_sm:chromosome chromosome:GRCh38:20:1:64444167:1 REF
>21 dna_sm:chromosome chromosome:GRCh38:21:1:46709983:1 REF
>22 dna_sm:chromosome chromosome:GRCh38:22:1:50818468:1 REF
>3 dna_sm:chromosome chromosome:GRCh38:3:1:198295559:1 REF
>4 dna_sm:chromosome chromosome:GRCh38:4:1:190214555:1 REF
>5 dna_sm:chromosome chromosome:GRCh38:5:1:181538259:1 REF
>6 dna_sm:chromosome chromosome:GRCh38:6:1:170805979:1 REF
>7 dna_sm:chromosome chromosome:GRCh38:7:1:159345973:1 REF
>8 dna_sm:chromosome chromosome:GRCh38:8:1:145138636:1 REF
>9 dna_sm:chromosome chromosome:GRCh38:9:1:138394717:1 REF
>MT dna_sm:chromosome chromosome:GRCh38:MT:1:16569:1 REF
>X dna_sm:chromosome chromosome:GRCh38:X:1:156040895:1 REF
>Y dna_sm:chromosome chromosome:GRCh38:Y:2781480:56887902:1 REF
>KI270728.1 dna_sm:scaffold scaffold:GRCh38:KI270728.1:1:1872759:1 REF
>KI270727.1 dna_sm:scaffold scaffold:GRCh38:KI270727.1:1:448248:1 REF
>KI270442.1 dna_sm:scaffold scaffold:GRCh38:KI270442.1:1:392061:1 REF
>KI270729.1 dna_sm:scaffold scaffold:GRCh38:KI270729.1:1:280839:1 REF
>GL000225.1 dna_sm:scaffold scaffold:GRCh38:GL000225.1:1:211173:1 REF
>KI270743.1 dna_sm:scaffold scaffold:GRCh38:KI270743.1:1:210658:1 REF

Primary assembly contains all toplevel sequence regions excluding haplotypes and patches. 

## run the homer
```
# rename chromosomes in fasta file 
sed -i 's/\s.*$//' Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

/home/lokorokova/bin/homer/bin/scanMotifGenomeWide.pl \
    /home/lokorokova/LINE1_project/line1.motif Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed -p 20 > \
    homer_results_human.bed &

# groom homer results
awk -F "\t" -v OFS="\t" '$5>=4 {print $1,$2-1,$3,$4,$5,$6}' homer_results_human.bed > homer_results_human_corrected.bed
bedtools sort -i homer_results_human_corrected.bed > sorted_homer_results_homer_corrected.bed
pigz homer_results_human.bed
pigz homer_results_human_corrected.bed
```
## exploring ensemble gtf
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4919035/ The Ensembl gene annotation system
https://www.youtube.com/watch?v=7S1-nqYRtrY - homology in Ensemble

```
zcat Homo_sapiens.GRCh38.105.gtf.gz | grep -v "^#" | cut -f 3 | sort -u
```
`CDS
exon
five_prime_utr
gene
Selenocysteine
start_codon
stop_codon
three_prime_utr
transcript`

## counting L1 motifs per feature

```
## filter the gtf - remove start_codon, stop_codon, Selenocysteine
zcat Homo_sapiens.GRCh38.105.gtf.gz | awk '$3 != "start_codon" && $3 != "stop_codon" && $3 != "Selenocysteine"' > filtered_Homo_sapiens.GRCh38.105.gtf

bedtools sort -i filtered_Homo_sapiens.GRCh38.105.gtf > sorted_filtered_Homo_sapiens.GRCh38.105.gtf

bedtools map -a sorted_filtered_Homo_sapiens.GRCh38.105.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -s > count_motif_same_strand.bed
bedtools map -a sorted_filtered_Homo_sapiens.GRCh38.105.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -S > count_motif_opposite_strand.bed


paste count_motif_same_strand.bed count_motif_opposite_strand.bed | cut -f 1,2,3,4,5,6,7,8,9,10,20 > count_motif_human.bed
head count_motif_human.bed

rm count_motif_same_strand.bed
rm count_motif_opposite_strand.bed
```

## control

``` 
cat homer_results_human.bed | awk -F "\t" -v OFS="\t" '$5>=4.96256 {print}' | wc -l
```
28681846
```
cat homer_results_human.bed | awk -F "\t" -v OFS="\t" '$5>=0 {print}' | wc -l
```
28681846

Nucleotide sequence extraction

`bedtools getfasta -s -fi ./Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
-bed ./homer_results_human_corrected.bed -tab > extracted_motifs.tsv` 

-s    Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.

`awk -F"\t" -v OFS="\t" '{split($1,a,"("); split(a[1],b,"-"); split(a[2],c,")"); split(b[1],d,":"); print(d[1],d[2],b[2],c[1],toupper($2))}' extracted_motifs.tsv > extracted_motifs.bed`

Generation of the weblogo

`cut -f 5 extracted_motifs.bed > motifL1.seq.fa`

`weblogo -c classic -s large -i -3 < motifL1.seq.fa > sequence.WebLogo.motifL1.eps`

## new gtf, prepared in R for bedtools map

```{r}
gtf <- readGFF("../genomes/Homo_sapiens.GRCh38.105.gtf")
gtf$rownumber <- rownames(gtf)
gtf_for_bedtools_map <- gtf[, c(1:8,27)]
fwrite(gtf_for_bedtools_map, "../genomes/Homo_sapiens.GRCh38.for_bedtoolmap.gtf", sep = "\t", col.names = FALSE)
```

```{bash}
bedtools sort -i Homo_sapiens.GRCh38.for_bedtoolmap.gtf > sorted_Homo_sapiens.GRCh38.for_bedtoolmap.gtf

bedtools map -a sorted_Homo_sapiens.GRCh38.for_bedtoolmap.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -s > count_motif_same_strand.bed

bedtools map -a sorted_Homo_sapiens.GRCh38.for_bedtoolmap.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -S > count_motif_opposite_strand.bed

paste count_motif_same_strand.bed count_motif_opposite_strand.bed | cut -f 1,2,3,4,5,6,7,8,9,10,20 > count_motif_human.bed
head count_motif_human.bed

rm count_motif_same_strand.bed
rm count_motif_opposite_strand.bed
```

## Add AT-content to gtf file
bedtools nuc -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed Homo_sapiens.GRCh38.for_bedtoolmap.gtf > Homo_sapiens.GRCh38_AT_content.gtf

## Extract introns from gtf - something went wrong
zcat Homo_sapiens.GRCh38.105.gtf.gz | awk '$3 == "transcript"' > Homo_sapiens.GRCh38.105_transcripts.gtf

zcat Homo_sapiens.GRCh38.105.gtf.gz | awk '$3 == "exon"' > Homo_sapiens.GRCh38.105_exons.gtf

bedtools subtract -a Homo_sapiens.GRCh38.105_transcripts.gtf -b Homo_sapiens.GRCh38.105_exons.gtf | head
`1    ensembl_havana    transcript    1211833    1211941    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000379236"; transcript_version "4"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS11"; tag "basic"; transcript_support_level "1 (assigned to previous version 3)";
1    ensembl_havana    transcript    1212139    1212637    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000379236"; transcript_version "4"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS11"; tag "basic"; transcript_support_level "1 (assigned to previous version 3)";
1    ensembl_havana    transcript    1212705    1212991    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000379236"; transcript_version "4"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS11"; tag "basic"; transcript_support_level "1 (assigned to previous version 3)";
1    ensembl_havana    transcript    1213786    1213982    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000379236"; transcript_version "4"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS11"; tag "basic"; transcript_support_level "1 (assigned to previous version 3)";
1    havana    transcript    1211833    1211941    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000497869"; transcript_version "5"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-203"; transcript_source "havana"; transcript_biotype "retained_intron"; transcript_support_level "2";
1    havana    transcript    1212139    1212637    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000497869"; transcript_version "5"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-203"; transcript_source "havana"; transcript_biotype "retained_intron"; transcript_support_level "2";
1    havana    transcript    1212705    1212991    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000497869"; transcript_version "5"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-203"; transcript_source "havana"; transcript_biotype "retained_intron"; transcript_support_level "2";
1    havana    transcript    1213786    1213982    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000497869"; transcript_version "5"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-203"; transcript_source "havana"; transcript_biotype "retained_intron"; transcript_support_level "2";
1    havana    transcript    1212139    1212637    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000453580"; transcript_version "1"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; transcript_support_level "3";
1    havana    transcript    1212705    1212991    .    -    .    gene_id "ENSG00000186827"; gene_version "11"; transcript_id "ENST00000453580"; transcript_version "1"; gene_name "TNFRSF4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF4-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; transcript_support_level "3";
`
## something went wrong file sorted_Homo_sapiens.GRCh38.intron_for_bedtoolmap.gtf looks like a lot of introns were lost
bedtools subtract -a Homo_sapiens.GRCh38.105_transcripts.gtf -b Homo_sapiens.GRCh38.105_exons.gtf > Homo_sapiens.GRCh38.105_introns.gtf

## preparing gtf for bedtools count - something went wrong file sorted_Homo_sapiens.GRCh38.intron_for_bedtoolmap.gtf looks like a lot of introns were lost

```{r}
intron_gtf <- readGFF("../genomes/Homo_sapiens.GRCh38.105_introns.gtf")
intron_gtf$rownumber <- rownames(intron_gtf)
intron_gtf_for_bedtools_map <- intron_gtf[, c(1:8, 22)]
fwrite(intron_gtf_for_bedtools_map, "../genomes/Homo_sapiens.GRCh38.intron_for_bedtoolmap.gtf", sep = "\t", col.names = FALSE)
```


```
bedtools sort -i Homo_sapiens.GRCh38.intron_for_bedtoolmap.gtf > sorted_Homo_sapiens.GRCh38.intron_for_bedtoolmap.gtf

bedtools map -a sorted_Homo_sapiens.GRCh38.intron_for_bedtoolmap.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -s > count_motif_introns_same_strand.bed

bedtools map -a sorted_Homo_sapiens.GRCh38.intron_for_bedtoolmap.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -S > count_motif_introns_opposite_strand.bed

paste count_motif_introns_same_strand.bed count_motif_introns_opposite_strand.bed | cut -f 1,2,3,4,5,6,7,8,9,10,19,20 > count_motif_introns_human.bed

rm count_motif_introns_same_strand.bed
rm count_motif_introns_opposite_strand.bed

## Add AT-content to gtf file
bedtools nuc -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed sorted_Homo_sapiens.GRCh38.intron_for_bedtoolmap.gtf > Homo_sapiens.GRCh38_AT_content_introns.gtf
```
## preparing gtf for bedtools count - only longest protein coding transcript

```{r}
fwrite(only_longest_transcripts[, c(1:8,27) ], 'only_longest_transcript_transcripts.gtf')
longest_transcripts <- only_longest_transcripts$transcript_id

fwrite(protein_coding_genes[feature == 'exon' & transcript_id %in% longest_transcripts, c(1:9)], 'only_longest_transcript_exons.gtf')
```

```
bedtools sort -i only_longest_transcript_exons.gtf > sorted_only_longest_transcript_exons.gtf
bedtools sort -i only_longest_transcript_transcripts.gtf > sorted_only_longest_transcript_transcripts.gtf

bedtools subtract -a sorted_only_longest_transcript_transcripts.gtf -b sorted_only_longest_transcript_exons.gtf > only_longest_transcripts_introns.gtf

bedtools sort -i only_longest_transcripts_introns.gtf > sorted_only_longest_transcripts_introns.gtf

bedtools map -a sorted_only_longest_transcripts_introns.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -s > count_motif_introns_same_strand.bed

bedtools map -a sorted_only_longest_transcripts_introns.gtf -b sorted_homer_results_homer_corrected.bed -c 6 -o count -S > count_motif_introns_opposite_strand.bed

paste count_motif_introns_same_strand.bed count_motif_introns_opposite_strand.bed | cut -f 1,2,3,4,5,6,7,8,9,10,19,20 > count_motif_introns_human.bed

rm count_motif_introns_same_strand.bed
rm count_motif_introns_opposite_strand.bed

## Add AT-content to gtf file
bedtools nuc -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed sorted_only_longest_transcripts_introns.gtf > Homo_sapiens.GRCh38_AT_content_introns.gtf
```

G-banding, G banding or Giemsa banding is a technique used in cytogenetics to produce a visible karyotype by staining condensed chromosomes. It is useful for identifying genetic diseases through the photographic representation of the entire chromosome complement.[1] The metaphase chromosomes are treated with trypsin (to partially digest the chromosome) and stained with Giemsa stain. Heterochromatic regions, which tend to be rich with adenine and thymine (AT-rich) DNA and relatively gene-poor, stain more darkly in G-banding. In contrast, less condensed chromatin (Euchromatin)—which tends to be rich with guanine and cytosine (GC-rich) and more transcriptionally active—incorporates less Giemsa stain, and these regions appear as light bands in G-banding.

gpos100
gpos25
gpos50
gpos75

"acen" is centromeric;
"stalk" refers to the short arm of acrocentric chromosomes chr13,14,15,21,22;
"gvar" bands tend to be heterochomatin, either pericentric or telomeric.

variable-length regions (gvar), or tightly-constricted regions (stalk).


```
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/hg38.fa.gz .

```

## Preparing files made in R for IGV

`cat U_zones_for_bedtools_count.gtf | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,".",$6, ".", $7}' > U_zones_for_IGV.gtf`

## Metaplots in deeptools

### Creation of bigWig file for L1 motif

- creation of the file

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$5)}' homer_results_human_corrected.bed > positions.motifL1.bedgraph`

`sort -k1,1 -k2,2n positions.motifL1.bedgraph > sorted.positions.motifL1.bedgraph`

overlapping motifs causes problems!
`bedtools merge -i sorted.positions.motifL1.bedgraph > not_overlapped_sorted.positions.motifL1.bed`
`bedtools map -a not_overlapped_sorted.positions.motifL1.bed -b sorted.positions.motifL1.bedgraph -c 4 -o mean > not_overlapped_sorted.positions.motifL1.bedgraph`

`cut -f1,2 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai > hg38.chrom.sizes`

`bedGraphToBigWig not_overlapped_sorted.positions.motifL1.bedgraph hg38.chrom.sizes bigWig.positions.motifL1.bw`


### Creation of bigWig AT%

#### make 5b wide bins and create a BED file with windows of wished width
`bedtools makewindows -g hg38.chrom.sizes -w 5 > hg38_5bp.bed`

#### compute base frequencies across all bins with BedTools nuc
`bedtools nuc -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    -bed hg38_5bp.bed  \
    > hg38_AT_content.txt`
    
$head hg38_AT_content.txt
```
    #1_usercol    2_usercol    3_usercol    4_pct_at    5_pct_gc    6_num_A    7_num_C    8_num_G    9_num_T    10_num_N    11_num_oth    12_seq_len
    1    0    5    0.000000    0.000000    0    0    0    0    5    0    5
    1    5    10    0.000000    0.000000    0    0    0    0    5    0    5
    1    10    15    0.000000    0.000000    0    0    0    0    5    0    5
    1    15    20    0.000000    0.000000    0    0    0    0    5    0    5
    1    20    25    0.000000    0.000000    0    0    0    0    5    0    5
    1    25    30    0.000000    0.000000    0    0    0    0    5    0    5
    1    30    35    0.000000    0.000000    0    0    0    0    5    0    5
    1    35    40    0.000000    0.000000    0    0    0    0    5    0    5
    1    40    45    0.000000    0.000000    0    0    0    0    5    0    5
```
    
`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$4)}' hg38_AT_content.txt > hg38_AT_content.bedgraph`

`sort -k1,1 -k2,2n hg38_AT_content.bedgraph > sorted.hg38_AT_content.bedgraph`

`bedGraphToBigWig sorted.hg38_AT_content.bedgraph hg38.chrom.sizes bigWig.hg38_AT_content.bw`




## Conservation of region

`wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw`

`bigWigToBedGraph hg38.phyloP100way.bw hg38.phyloP100way.bedgraph` 

`sort -k1,1 -k2,2n all_exons_PSI_data_derived.bed > sorted_all_exons_PSI_data_derived.bed`

`awk -F "\t" -v OFS="\t" {'print "chr"$1,$2,$3,$4,$5,$6'} sorted_all_exons_PSI_data_derived.bed`

`awk -F "\t" -v OFS="\t" {'print "chr"$1,$2,$3,$4,$5,$6'} sorted_all_exons_PSI_data_derived.bed > test.bed`

`bedtools map -a test.bed -b hg38.phyloP100way.bedgraph -c 4 -o mean > all_exons_PSI_data_derived_conservation_level.bed`

## Different exons

`computeMatrix scale-regions -S bigWig.hg38_AT_content.bw bigWig.positions.motifL1.bw \
                              -R all_exons_PSI_data_derived.bed \
                              --binSize 40 \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 100 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o matrix_all_exons.mat.gz`

`plotProfile -m matrix_all_exons.mat.gz \
              -out exons.png \
              --numPlotsPerRow 2 \
              --plotTitle "All exons"`


## Test metaplots in deeptools for chr 22
`samtools faidx chr22.fa`
`cut -f1,2 chr22.fa.fai > hg38_chr22_chrom.sizes`

### Creation of bigWig file for L1 motif

- creation of the file
extract L1 motifs on chr22

`awk '$1=="22" {print $1"\t"$2"\t"$3"\t"$4"\t"".""\t"$5}' homer_results_human_corrected.bed > chr22_L1motif.bed`

`bedtools genomecov -i chr22_L1motif.bed -g hg38_chr22_chrom.sizes -bga -strand + > chr22_L1motif_plus_strand.bedgraph`
`bedtools genomecov -i chr22_L1motif.bed -g hg38_chr22_chrom.sizes -bga -strand - > chr22_L1motif_minus_strand.bedgraph`
`bedGraphToBigWig chr22_L1motif_plus_strand.bedgraph hg38_chr22_chrom.sizes chr22_L1motif_plus_strand.bw`
`bedGraphToBigWig chr22_L1motif_minus_strand.bedgraph hg38_chr22_chrom.sizes chr22_L1motif_minus_strand.bw`

### Creation of bigWig AT%

#### make 5b wide bins and create a BED file with windows of wished width
`bedtools makewindows -g hg38_chr22_chrom.sizes -w 5 > hg38_22chr_5bp.bed`

#### compute base frequencies across all bins with BedTools nuc
`bedtools nuc -fi chr22.fa \
    -bed hg38_22chr_5bp.bed  \
    > hg38_chr22_AT_content.txt`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$4)}' hg38_chr22_AT_content.txt > hg38_chr22_AT_content.bedgraph`
`bedGraphToBigWig hg38_chr22_AT_content.bedgraph hg38_chr22_chrom.sizes chr22_AT_content.bw`

### Extracting exons 
`awk '$1=="22" {print ($1"\t"$4-1"\t"$5"\t"$7"\t"".""\t"$6)}' only_longest_transcript_exons.gtf > exons_chr22_from_only_longest_transcript.bed`

`awk '$6=="+"' exons_chr22_from_only_longest_transcript.bed > exons_chr22_plus_strand_from_only_longest_transcript.bed`

`awk '$6=="-"' exons_chr22_from_only_longest_transcript.bed > exons_chr22_minus_strand_from_only_longest_transcript.bed`

## Creating plots


`computeMatrix scale-regions -S chr22_L1motif_plus_strand.bw chr22_L1motif_minus_strand.bw chr22_AT_content.bw\
                              -R exons_chr22_plus_strand_from_only_longest_transcript.bed exons_chr22_minus_strand_from_only_longest_transcript.bed \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 500 \
                              --regionBodyLength 500 \
                              --afterRegionStartLength 500 \
                              -o matrix_exons_chr22.mat.gz`

`plotProfile -m matrix_exons_chr22.mat.gz \
              -out exons.png \
              --numPlotsPerRow 3 \
              --plotTitle "Exons_22chr"`
              
`computeMatrix scale-regions -S chr22_AT_content.bw\
                              -R exons_chr22_plus_strand_from_only_longest_transcript.bed exons_chr22_minus_strand_from_only_longest_transcript.bed \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 500 \
                              --regionBodyLength 500 \
                              --afterRegionStartLength 500 \
                              -o matrix_exons_chr22_AT.mat.gz`
                              
`computeMatrix scale-regions -S chr22_L1motif_plus_strand.bw chr22_L1motif_minus_strand.bw \
                              -R exons_chr22_plus_strand_from_only_longest_transcript.bed exons_chr22_minus_strand_from_only_longest_transcript.bed \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 500 \
                              --regionBodyLength 500 \
                              --afterRegionStartLength 500 \
                              -o matrix_exons_chr22_L1.mat.gz`

`plotProfile -m matrix_exons_chr22_AT.mat.gz \
              -out exons_AT.png \
              --numPlotsPerRow 3 \
              --kmeans 4 \
              --plotTitle "Exons_22chr"`
              
`plotProfile -m matrix_exons_chr22_L1.mat.gz \
              -out exons_L1.png \
              --numPlotsPerRow 2 \
              --kmeans 4 \
              --plotTitle "Exons_22chr"`



## Adding RE coverage into table

We process RepeatMasker annotation data (as described in:
http://www.repeatmasker.org/webrepeatmaskerhelp.html) converting
them into the first six UCSC BED columns, as follows:

- Query sequence            <-->   chromosome (1st column)
- Query start -  1          <-->   start (2nd column)
- Query end                 <-->   stop (3rd column)
- Match repeat name         <-->   id (4th column)
- SW score                  <-->   score (5th column)
- strand                    <-->   strand (6th column)

`
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz .
gunzip hg38.fa.out.gz
convert2bed -i rmsk < hg38.fa.out > hg38.repeats.bed
##rename chromosome
awk '{gsub(/chr1/,"1")\
;gsub(/chr10/,"10");gsub(/chr11/,"11");gsub(/chr12/,"12")\
;gsub(/chr13/,"13");gsub(/chr14/,"14");gsub(/chr15/,"15")\
;gsub(/chr16/,"16");gsub(/chr17/,"17")\
;gsub(/chr18/,"18");gsub(/chr19/,"19")\
;gsub(/chr2/,"2")\
;gsub(/chr20/,"20");gsub(/chr21/,"21")\
;gsub(/chr22/,"22");gsub(/chr3/,"3");gsub(/chr4/,"4");gsub(/chr5/,"5")\
;gsub(/chr6/,"6");gsub(/chr7/,"7")\
;gsub(/chr8/,"8");gsub(/chr9/,"9");gsub(/chrX/,"x")\
;gsub(/chrY/,"Y");gsub(/chrM/,"MT")\
;gsub(/chr14_GL000009v2_random/,"GL000009.2")\
;gsub(/chr14_GL000194v1_random/,"GL000194.1")\
;gsub(/chrUn_GL000195v1/,"GL000195.1")\
;gsub(/chr17_GL000205v2_random/,"GL000205.2")\
;gsub(/chrUn_GL000213v1/,"GL000213.1")\
;gsub(/chrUn_GL000216v2/,"GL000216.2")\
;gsub(/chrUn_GL000218v1/,"GL000218.1")\
;gsub(/chrUn_GL000219v1/,"GL000219.1")\
;gsub(/chrUn_GL000220v1/,"GL000220.1")\
;gsub(/chr14_GL000225v1_random/,"GL000225.1")\
;gsub(/chrUn_KI270442v1/,"KI270442.1")\
;gsub(/chr1_KI270711v1_random/,"KI270711.1")\
;gsub(/chr1_KI270713v1_random/,"KI270713.1")\
;gsub(/chr11_KI270721v1_random/,"KI270721.1")\
;gsub(/chr14_KI270726v1_random/,"KI270726.1")\
;gsub(/chr15_KI270727v1_random/,"KI270727.1")\
;gsub(/chr16_KI270728v1_random/,"KI270728.1")\
;gsub(/chr22_KI270731v1_random/,"KI270731.1")\
;gsub(/chr22_KI270733v1_random/,"KI270733.1")\
;gsub(/chr22_KI270734v1_random/,"KI270734.1")\
;gsub(/chrUn_KI270744v1/,"KI270744.1")\
;gsub(/chrUn_KI270750v1/,"KI270750.1")}1' hg38.repeats.bed > hg38.repeats_chr.bed
cut -f 1,2,3,4,5,6 hg38.repeats_chr.bed > hg38.repeats.6col.bed
bedtools coverage -a sorted_Homo_sapiens.GRCh38.for_bedtoolmap.gtf -b hg38.repeats.6col.bed > Homo_sapiens.GRCh38_total_REcoverage.gtf
##extract L1 repeats
awk -v OFS="\t" '$11=="LINE/L1"' hg38.repeats_chr.bed | cut -f 1,2,3,4,5,6 > hg38_repeats_L1.bed
bedtools coverage -a sorted_Homo_sapiens.GRCh38.for_bedtoolmap.gtf -b hg38_repeats_L1.bed > Homo_sapiens.GRCh38_L1_REcoverage.gtf
##extract SINE repeats
awk -v OFS="\t" '$11=="SINE*"' hg38.repeats_chr.bed | cut -f 1,2,3,4,5,6 > hg38_repeats_SINEs.bed
bedtools coverage -a sorted_Homo_sapiens.GRCh38.for_bedtoolmap.gtf -b hg38_repeats_SINEs.bed > Homo_sapiens.GRCh38_SINE_REcoverage.gtf
##extract Simple repeats
awk -v OFS="\t" '$11=="Simple_repeat"' hg38.repeats_chr.bed | cut -f 1,2,3,4,5,6 > hg38_repeats_simpleR.bed
bedtools coverage -a sorted_Homo_sapiens.GRCh38.for_bedtoolmap.gtf -b hg38_repeats_simpleR.bed > Homo_sapiens.GRCh38_simpleR_REcoverage.gtf
`

cut -f 11 hg38.repeats.bed | sort -u
DNA
DNA?
DNA/hAT
DNA/hAT?
DNA/hAT-Ac
DNA/hAT-Blackjack
DNA/hAT-Charlie
DNA/hAT-Tag1
DNA/hAT-Tip100
DNA/hAT-Tip100?
DNA/Merlin
DNA/MULE-MuDR
DNA/PIF-Harbinger
DNA/PiggyBac
DNA/PiggyBac?
DNA/TcMar
DNA/TcMar?
DNA/TcMar-Mariner
DNA/TcMar-Pogo
DNA/TcMar-Tc2
DNA/TcMar-Tigger
LINE/CR1
LINE/Dong-R4
LINE/L1
LINE/L2
LINE/Penelope
LINE/RTE-BovB
LINE/RTE-X
Low_complexity
LTR
LTR?
LTR/ERV1
LTR/ERV1?
LTR/ERVK
LTR/ERVL
LTR/ERVL?
LTR/ERVL-MaLR
LTR/Gypsy
LTR/Gypsy?
RC?/Helitron?
RC/Helitron
Retroposon/SVA
RNA
rRNA
Satellite
Satellite/acro
Satellite/centr
Satellite/telo
scRNA
Simple_repeat
SINE?
SINE/5S-Deu-L2
SINE/Alu
SINE/MIR
SINE?/tRNA
SINE/tRNA
SINE/tRNA-Deu
SINE/tRNA-RTE
snRNA
srpRNA
tRNA
Unknown


## Making metaplots for introns

`awk '{print $1"\t"$2"\t"$3"\t"$4"\t"".""\t"$5}' sorted_homer_results_homer_corrected.bed > L1motif.bed`

`bedtools genomecov -i L1motif.bed -g ../hg38.chrom.sizes -bga -strand + > L1motif_plus_strand.bedgraph`
`bedSort L1motif_plus_strand.bedgraph L1motif_plus_strand.bedgraph`
`bedtools genomecov -i L1motif.bed -g ../hg38.chrom.sizes -bga -strand - > L1motif_minus_strand.bedgraph`
`bedSort L1motif_minus_strand.bedgraph L1motif_minus_strand.bedgraph`
`bedGraphToBigWig L1motif_plus_strand.bedgraph ../hg38.chrom.sizes L1motif_plus_strand.bw`
`bedGraphToBigWig L1motif_minus_strand.bedgraph ../hg38.chrom.sizes L1motif_minus_strand.bw`

`for i in *.bed
do
BED=${i%%.bed}
##compute matrix for AT content
computeMatrix scale-regions -S ../genomes/for_metaplots/bigWig.hg38_AT_content.bw \
                              -R $BED".bed" \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 100 \
                              --regionBodyLength 1000 \
                              --afterRegionStartLength 100 \
                              -o "matrix_"$BED"_AT.mat.gz"
##compute matrix for L1motif
computeMatrix scale-regions -S ../genomes/for_metaplots/L1motif_plus_strand.bw \
                                                                ../genomes/for_metaplots/L1motif_minus_strand.bw \
                              -R  $BED".bed" \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 100 \
                              --regionBodyLength 1000 \
                              --afterRegionStartLength 100 \
                              -o "matrix_"$BED"_L1.mat.gz"
##make plots
plotProfile -m "matrix_"$BED"_AT.mat.gz" \
              -out $BED"_AT.png" \
              --startLabel "Intron start" \
              --endLabel "Intron end" \
              --kmeans 20 
plotProfile -m "matrix_"$BED"_L1.mat.gz" \
              -out $BED"_L1.png" \
              --startLabel "Intron start" \
              --endLabel "Intron end" \
              --kmeans 20
done
`

`for i in *.bed
do
BED=${i%%.bed}
##make plots
plotProfile -m "matrix_"$BED"_AT.mat.gz" \
              -out $BED"_AT.png" \
              --startLabel "RS" \
              --endLabel "RE" \
              --legendLocation none \
              --kmeans 10 \
              --outFileNameData $BED"_L1.data"
plotProfile -m "matrix_"$BED"_L1.mat.gz" \
              -out $BED"_L1.png" \
              --startLabel "RS" \
              --numPlotsPerRow 1 \
              --legendLocation none \
              --endLabel "RE" \
              --kmeans 10 \
              --outFileNameData $BED"_L1.data"
done
`
                              
## All intron analysis

`for i in all_tails_for_bedtools_count.gtf all_promoters_for_bedtools_count.gtf
do
echo "Processing..."$i
GTF=${i%%.gtf}
sortBed -i $GTF".gtf" > $GTF"_sorted.gtf" 
echo "Sorting done..."$i
rm $GTF".gtf"
mv $GTF"_sorted.gtf" $GTF".gtf"
bedtools map -a $GTF".gtf" -b ~/INTRANET/cristo-nas-shared/Personal_documents/Larisa/genomes/sorted_homer_results_homer_corrected.bed  -c 6 -o count -s > $GTF"_count_motif_same_strand.bed"
bedtools map -a $GTF".gtf" -b ~/INTRANET/cristo-nas-shared/Personal_documents/Larisa/genomes/sorted_homer_results_homer_corrected.bed -c 6 -o count -S > $GTF"_count_motif_opposite_strand.bed"
paste $GTF"_count_motif_same_strand.bed" $GTF"_count_motif_opposite_strand.bed" | cut -f 1,2,3,4,5,6,7,8,9,10,19,20 > $GTF"_count_motif_human.bed"
echo "Motif counting done..."$i
rm $GTF"_count_motif_same_strand.bed"
rm $GTF"_count_motif_opposite_strand.bed"
##add AT content
bedtools nuc -fi ~/INTRANET/cristo-nas-shared/Personal_documents/Larisa/genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed $GTF".gtf" > $GTF"_ATcontent.bed"
echo "AT content done..."$i
done
`

## Compute distance to L1 in introns


`
awk '{gsub(/chr1/,"1")\
> ;gsub(/chr10/,"10");gsub(/chr11/,"11");gsub(/chr12/,"12")\
> ;gsub(/chr13/,"13");gsub(/chr14/,"14");gsub(/chr15/,"15")\
> ;gsub(/chr16/,"16");gsub(/chr17/,"17")\
> ;gsub(/chr18/,"18");gsub(/chr19/,"19")\
> ;gsub(/chr2/,"2")\
> ;gsub(/chr20/,"20");gsub(/chr21/,"21")\
> ;gsub(/chr22/,"22");gsub(/chr3/,"3");gsub(/chr4/,"4");gsub(/chr5/,"5")\
> ;gsub(/chr6/,"6");gsub(/chr7/,"7")\
> ;gsub(/chr8/,"8");gsub(/chr9/,"9");gsub(/chrX/,"x")\
> ;gsub(/chrY/,"Y");gsub(/chrM/,"MT")\
> ;gsub(/chr14_GL000009v2_random/,"GL000009.2")\
> ;gsub(/chr14_GL000194v1_random/,"GL000194.1")\
> ;gsub(/chrUn_GL000195v1/,"GL000195.1")\
> ;gsub(/chr17_GL000205v2_random/,"GL000205.2")\
> ;gsub(/chrUn_GL000213v1/,"GL000213.1")\
> ;gsub(/chrUn_GL000216v2/,"GL000216.2")\
> ;gsub(/chrUn_GL000218v1/,"GL000218.1")\
> ;gsub(/chrUn_GL000219v1/,"GL000219.1")\
> ;gsub(/chrUn_GL000220v1/,"GL000220.1")\
> ;gsub(/chr14_GL000225v1_random/,"GL000225.1")\
> ;gsub(/chrUn_KI270442v1/,"KI270442.1")\
> ;gsub(/chr1_KI270711v1_random/,"KI270711.1")\
> ;gsub(/chr1_KI270713v1_random/,"KI270713.1")\
> ;gsub(/chr11_KI270721v1_random/,"KI270721.1")\
> ;gsub(/chr14_KI270726v1_random/,"KI270726.1")\
> ;gsub(/chr15_KI270727v1_random/,"KI270727.1")\
> ;gsub(/chr16_KI270728v1_random/,"KI270728.1")\
> ;gsub(/chr22_KI270731v1_random/,"KI270731.1")\
> ;gsub(/chr22_KI270733v1_random/,"KI270733.1")\
> ;gsub(/chr22_KI270734v1_random/,"KI270734.1")\
> ;gsub(/chrUn_KI270744v1/,"KI270744.1")\
> ;gsub(/chrUn_KI270750v1/,"KI270750.1")}1' rm_out_L1.bed > rm_out_L1_chr.bed
bedtools intersect -a all_introns_for_bedtools_count.gtf -b rm_out_L1_chr.bed -wa -wb > introns_intersect_L1.bed
`



### Creation of bigWig for L1, SINEs
```
awk '{gsub(/chr1/,"1")\
;gsub(/chr10/,"10");gsub(/chr11/,"11");gsub(/chr12/,"12")\
;gsub(/chr13/,"13");gsub(/chr14/,"14");gsub(/chr15/,"15")\
;gsub(/chr16/,"16");gsub(/chr17/,"17")\
;gsub(/chr18/,"18");gsub(/chr19/,"19")\
;gsub(/chr2/,"2")\
;gsub(/chr20/,"20");gsub(/chr21/,"21")\
;gsub(/chr22/,"22");gsub(/chr3/,"3");gsub(/chr4/,"4");gsub(/chr5/,"5")\
;gsub(/chr6/,"6");gsub(/chr7/,"7")\
;gsub(/chr8/,"8");gsub(/chr9/,"9");gsub(/chrX/,"X")\
;gsub(/chrY/,"Y");gsub(/chrM/,"MT")\
;gsub(/chr14_GL000009v2_random/,"GL000009.2")\
;gsub(/chr14_GL000194v1_random/,"GL000194.1")\
;gsub(/chrUn_GL000195v1/,"GL000195.1")\
;gsub(/chr17_GL000205v2_random/,"GL000205.2")\
;gsub(/chrUn_GL000213v1/,"GL000213.1")\
;gsub(/chrUn_GL000216v2/,"GL000216.2")\
;gsub(/chrUn_GL000218v1/,"GL000218.1")\
;gsub(/chrUn_GL000219v1/,"GL000219.1")\
;gsub(/chrUn_GL000220v1/,"GL000220.1")\
;gsub(/chr14_GL000225v1_random/,"GL000225.1")\
;gsub(/chrUn_KI270442v1/,"KI270442.1")\
;gsub(/chr1_KI270711v1_random/,"KI270711.1")\
;gsub(/chr1_KI270713v1_random/,"KI270713.1")\
;gsub(/chr11_KI270721v1_random/,"KI270721.1")\
;gsub(/chr14_KI270726v1_random/,"KI270726.1")\
;gsub(/chr15_KI270727v1_random/,"KI270727.1")\
;gsub(/chr16_KI270728v1_random/,"KI270728.1")\
;gsub(/chr22_KI270731v1_random/,"KI270731.1")\
;gsub(/chr22_KI270733v1_random/,"KI270733.1")\
;gsub(/chr22_KI270734v1_random/,"KI270734.1")\
;gsub(/chrUn_KI270744v1/,"KI270744.1")\
;gsub(/chrUn_KI270750v1/,"KI270750.1")}1' hg38.repeats.bed > hg38_repeats_chr.bed
```
`awk -v OFS="\t" '$11 ~/LINE\/L1/ && $1 !~/.alt/ && $1 !~/.random/ && $1 !~/chrUn./' hg38_repeats_chr.bed | cut -f 1,2,3,4,5,6 > hg38_repeats_L1.bed`

`bedtools genomecov -i hg38_repeats_L1.bed -g hg38.chrom.sizes -bga -strand + > hg38_L1_plus.bedgraph`
`bedtools genomecov -i hg38_repeats_L1.bed -g hg38.chrom.sizes -bga -strand - > hg38_L1_minus.bedgraph`
`bedSort hg38_L1_plus.bedgraph hg38_L1_plus.bedgraph`
`bedSort hg38_L1_minus.bedgraph hg38_L1_minus.bedgraph`
`bedGraphToBigWig hg38_L1_plus.bedgraph hg38.chrom.sizes hg38_L1_plus.bw`
`bedGraphToBigWig hg38_L1_minus.bedgraph hg38.chrom.sizes hg38_L1_minus.bw`

`awk -v OFS="\t" '$11 ~/SINE./ && $1 !~/.alt/ && $1 !~/.random/ && $1 !~/chrUn./' hg38_repeats_chr.bed | cut -f 1,2,3,4,5,6 > hg38_repeats_SINEs.bed`

`bedtools genomecov -i hg38_repeats_SINEs.bed -g hg38.chrom.sizes -bga -strand + > hg38_SINEs_plus.bedgraph`
`bedtools genomecov -i hg38_repeats_SINEs.bed -g hg38.chrom.sizes -bga -strand - > hg38_SINEs_minus.bedgraph`
`bedSort hg38_SINEs_plus.bedgraph hg38_SINEs_plus.bedgraph`
`bedSort hg38_SINEs_minus.bedgraph hg38_SINEs_minus.bedgraph`
`bedGraphToBigWig hg38_SINEs_plus.bedgraph hg38.chrom.sizes hg38_SINEs_plus.bw`
`bedGraphToBigWig hg38_SINEs_minus.bedgraph hg38.chrom.sizes hg38_SINEs_minus.bw`


## Metaplots of all exons of only longest transcript
`
for i in only_longest_transcripts_intermedite_exons.bed only_longest_transcripts_first_exons.bed only_longest_transcripts_last_exons.bed only_longest_transcripts_monoexonic_exons.bed
do
echo "Processing .. " $i
BED=${i%%.bed}
##compute matrix for AT content
echo "compute matrix for AT content.. "$i
computeMatrix scale-regions -S ../genomes/for_metaplots/bigWig.hg38_AT_content.bw \
                              -R $BED".bed" \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 2000 \
                              --regionBodyLength 300 \
                              --afterRegionStartLength 2000 \
                              -o "matrix_"$BED"_AT.mat.gz"
echo "Plotting heatmap for AT content..."
plotHeatmap -m "matrix_"$BED"_AT.mat.gz" \
                         -out $BED"_ATcontent.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                         --endLabel "RE"
##compute matrix for L1motif
echo "compute matrix for L1 motif.. "$i
computeMatrix scale-regions -S ../genomes/for_metaplots/L1motif_plus_strand.bw \
                                                                ../genomes/for_metaplots/L1motif_minus_strand.bw \
                              -R  $BED".bed" \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 2000 \
                              --regionBodyLength 300 \
                              --afterRegionStartLength 2000 \
                              -o "matrix_"$BED"_L1motif.mat.gz"
echo "plotting heatmap for L1 motif .. "$i
plotHeatmap -m "matrix_"$BED"_L1motif.mat.gz" \
                         -out $BED"_L1motif.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                         --endLabel "RE"
##compute matrix for L1
echo "compute matrix for L1.. "$i
computeMatrix scale-regions -S ../genomes/hg38_L1_plus.bw \
                                                                ../genomes/hg38_L1_minus.bw \
                                            -R  $i \
                                            --binSize 5 \
                                            --missingDataAsZero \
                                            --beforeRegionStartLength 2000 \
                                            --regionBodyLength 300 \
                                            --afterRegionStartLength 2000 \
                                            -o "matrix_"$BED"_L1.mat.gz"
echo "plotting heatmap for L1 .. "$i
plotHeatmap -m "matrix_"$BED"_L1.mat.gz" \
                         -out $BED"_L1.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                          --endLabel "RE"
##compute matrix for SINEs
echo "compute matrix for SINEs.. "$i
computeMatrix scale-regions -S ../genomes/hg38_SINEs_plus.bw \
                                                                ../genomes/hg38_SINEs_minus.bw \
                                                          -R  $i \
                                                          --binSize 5 \
                                                         --missingDataAsZero \
                                                         --beforeRegionStartLength 2000 \
                                                         --regionBodyLength 300 \
                                                         --afterRegionStartLength 2000 \
                                                         -o "matrix_"$BED"_SINEs.mat.gz"
echo "plotting heatmap for SINEs .. "$i
plotHeatmap -m "matrix_"$BED"_SINEs.mat.gz" \
                         -out $BED"_SINEs.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                        --endLabel "RE"                        
done
`

## Metaplots for understand density of exons and genes around L1
`
for i in *.bed
do
echo "Processing .. " $i
BED=${i%%.bed}
echo "Calculating genomecov .. "
bedtools genomecov -i $i -g ../genomes/hg38.chrom.sizes -bga -strand + > $BED"_plus.bedgraph"
bedtools genomecov -i $i -g ../genomes/hg38.chrom.sizes -bga -strand - > $BED"_minus.bedgraph"
echo "Sorting .. "
bedSort $BED"_plus.bedgraph" $BED"_plus.bedgraph"
bedSort $BED"_minus.bedgraph" $BED"_minus.bedgraph"
echo "Making bigwigs .. "
bedGraphToBigWig $BED"_plus.bedgraph" ../genomes/hg38.chrom.sizes $BED"_plus.bw"
bedGraphToBigWig $BED"_minus.bedgraph" ../genomes/hg38.chrom.sizes $BED"_minus.bw"
done
`

`
##compute matrix
echo "compute matrix .. "
computeMatrix scale-regions -S only_longest_transcripts_first_exons_minus.bw \
                                                                only_longest_transcripts_first_exons_plus.bw \
                                                                only_longest_transcripts_intermedite_exons_minus.bw \
                                                                only_longest_transcripts_intermedite_exons_plus.bw \
                                                                only_longest_transcripts_last_exons_minus.bw \
                                                                only_longest_transcripts_last_exons_plus.bw \
                                                                only_longest_transcripts_monoexonic_exons_minus.bw \
                                                                only_longest_transcripts_monoexonic_exons_plus.bw \
                                                                only_longest_transcripts_transcripts_minus.bw \
                                                                only_longest_transcripts_transcripts_plus.bw \
                                                          -R  ../genomes/hg38_L1_plus.bedgraph \
                                                                    ../genomes/hg38_L1_minus.bedgraph \
                                                          --binSize 5 \
                                                         --missingDataAsZero \
                                                         --beforeRegionStartLength 20000 \
                                                         --regionBodyLength 1000 \
                                                         --afterRegionStartLength 20000 \
                                                         -o matrix_exons_genes_around_L1.mat.gz
echo "plotting heatmap .. "
plotHeatmap -m matrix_exons_genes_around_L1.mat.gz \
                         -out exons_genes_around_L1.png \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                        --endLabel "RE"                        
`

`
##compute matrix
echo "compute matrix .. "
computeMatrix scale-regions -S ../genomes/for_metaplots/bigWig.hg38_AT_content.bw \
                                                          -R  ../genomes/hg38_L1_plus.bedgraph \
                                                                    ../genomes/hg38_L1_minus.bedgraph \
                                                          --binSize 5 \
                                                         --missingDataAsZero \
                                                         --beforeRegionStartLength 20000 \
                                                         --regionBodyLength 1000 \
                                                         --afterRegionStartLength 20000 \
                                                         -o matrix_AT_around_L1.mat.gz
echo "plotting heatmap .. "
plotHeatmap -m matrix_AT_around_L1.mat.gz \
                         -out AT_around_L1.png \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                        --endLabel "RE"                        
`
`
##compute matrix
echo "compute matrix .. "
computeMatrix scale-regions -S ../genomes/for_metaplots/L1motif_minus_strand.bw \
                                                                ../genomes/for_metaplots/L1motif_plus_strand.bw \
                                                          -R  ../genomes/hg38_L1_plus.bedgraph \
                                                                    ../genomes/hg38_L1_minus.bedgraph \
                                                          --binSize 5 \
                                                         --missingDataAsZero \
                                                         --beforeRegionStartLength 20000 \
                                                         --regionBodyLength 1000 \
                                                         --afterRegionStartLength 20000 \
                                                         -o matrix_L1motif_around_L1.mat.gz
echo "plotting heatmap .. "
plotHeatmap -m matrix_L1motif_around_L1.mat.gz \
                         -out L1motif_around_L1.png \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                        --endLabel "RE"                        
`
`
awk -v OFS="\t" '$5 == "-"' only_longest_transcripts_intermedite_exons.bed > only_longest_transcripts_intermediate_exons_minus.bed
awk -v OFS="\t" '$5 == "+"' only_longest_transcripts_intermedite_exons.bed > only_longest_transcripts_intermediate_exons_plus.bed
BED="only_longest_transcripts_intermediate_exons_by_strand"
##compute matrix for AT content
echo "compute matrix for AT content.. "
computeMatrix scale-regions -S ../genomes/for_metaplots/bigWig.hg38_AT_content.bw \
                              -R only_longest_transcripts_intermediate_exons_plus.bed \
                                    only_longest_transcripts_intermediate_exons_minus.bed  \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 2000 \
                              --regionBodyLength 300 \
                              --afterRegionStartLength 2000 \
                              -o "matrix_"$BED"_AT.mat.gz"
echo "Plotting heatmap for AT content..."
plotHeatmap -m "matrix_"$BED"_AT.mat.gz" \
                         -out $BED"_ATcontent.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                         --endLabel "RE"
##compute matrix for L1motif
echo "compute matrix for L1 motif.. "
computeMatrix scale-regions -S ../genomes/for_metaplots/L1motif_plus_strand.bw \
                                                                ../genomes/for_metaplots/L1motif_minus_strand.bw \
                             -R only_longest_transcripts_intermediate_exons_plus.bed \
                                  only_longest_transcripts_intermediate_exons_minus.bed  \
                              --binSize 5 \
                              --missingDataAsZero \
                              --beforeRegionStartLength 2000 \
                              --regionBodyLength 300 \
                              --afterRegionStartLength 2000 \
                              -o "matrix_"$BED"_L1motif.mat.gz"
echo "plotting heatmap for L1 motif .. "
plotHeatmap -m "matrix_"$BED"_L1motif.mat.gz" \
                         -out $BED"_L1motif.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                         --endLabel "RE"
##compute matrix for L1
echo "compute matrix for L1.. "$i
computeMatrix scale-regions -S ../genomes/hg38_L1_plus.bw \
                                                                ../genomes/hg38_L1_minus.bw \
                                           -R only_longest_transcripts_intermediate_exons_plus.bed \
                                                 only_longest_transcripts_intermediate_exons_minus.bed  \
                                            --binSize 5 \
                                            --missingDataAsZero \
                                            --beforeRegionStartLength 2000 \
                                            --regionBodyLength 300 \
                                            --afterRegionStartLength 2000 \
                                            -o "matrix_"$BED"_L1.mat.gz"
echo "plotting heatmap for L1 .. "
plotHeatmap -m "matrix_"$BED"_L1.mat.gz" \
                         -out $BED"_L1.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                          --endLabel "RE"
##compute matrix for SINEs
echo "compute matrix for SINEs.. "
computeMatrix scale-regions -S ../genomes/hg38_SINEs_plus.bw \
                                                                ../genomes/hg38_SINEs_minus.bw \
                                                          -R only_longest_transcripts_intermediate_exons_plus.bed \
                                                             only_longest_transcripts_intermediate_exons_minus.bed  \
                                                          --binSize 5 \
                                                         --missingDataAsZero \
                                                         --beforeRegionStartLength 2000 \
                                                         --regionBodyLength 300 \
                                                         --afterRegionStartLength 2000 \
                                                         -o "matrix_"$BED"_SINEs.mat.gz"
echo "plotting heatmap for SINEs .. "
plotHeatmap -m "matrix_"$BED"_SINEs.mat.gz" \
                         -out $BED"_SINEs.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                        --endLabel "RE"                        
`

## Genome wide L1 motif density in 1 kb window

#### make 1kb wide bins and create a BED file with windows of wished width
`bedtools makewindows -g hg38.chrom.sizes -w 1000 > hg38_1kb.bed`

#### compute base frequencies across all bins with BedTools nuc
`bedtools nuc -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    -bed hg38_1kb.bed  \
    > hg38_AT_content.bed`
    
    
`bedtools map -a hg38_1kb.bed -b sorted_homer_results_homer_corrected.bed -c 6 -o count -s > count_motif_1kb_same_strand.bed
bedtools map -a hg38_1kb.bed -b sorted_homer_results_homer_corrected.bed -c 6 -o count -S > count_motif_1kb_opposite_strand.bed
paste count_motif_1kb_same_strand.bed count_motif_1kb_opposite_strand.bed | cut -f 1,2,3,4,5,6,7,8,9,10,19,20 > count_motif_1kb_human.bed`
    
## Density of Alu and L1HS elements around exons

`echo "Calculating bw files for Alu..." 
awk -v OFS="\t" '$11 ~/SINE\/Alu/ && $1 !~/.alt/ && $1 !~/.random/ && $1 !~/chrUn./' hg38_repeats_chr.bed | cut -f 1,2,3,4,5,6 > hg38_repeats_Alu.bed
bedtools genomecov -i hg38_repeats_Alu.bed -g hg38.chrom.sizes -bga -strand + > hg38_Alu_plus.bedgraph
bedtools genomecov -i hg38_repeats_Alu.bed -g hg38.chrom.sizes -bga -strand - > hg38_Alu_minus.bedgraph
bedSort hg38_Alu_plus.bedgraph hg38_Alu_plus.bedgraph
bedSort hg38_Alu_minus.bedgraph hg38_Alu_minus.bedgraph
bedGraphToBigWig hg38_Alu_plus.bedgraph hg38.chrom.sizes hg38_Alu_plus.bw
bedGraphToBigWig hg38_Alu_minus.bedgraph hg38.chrom.sizes hg38_Alu_minus.bw`


`awk -v OFS="\t" '$11 ~/LINE\/L1/ && $1 !~/.alt/ && $1 !~/.random/ && $1 !~/chrUn./ && $4 =="L1HS"' hg38_repeats_chr.bed | cut -f 1,2,3,4,5,6 > hg38_repeats_L1HS.bed
bedtools genomecov -i hg38_repeats_L1HS.bed -g hg38.chrom.sizes -bga -strand + > hg38_L1HS_plus.bedgraph
bedtools genomecov -i hg38_repeats_L1HS.bed -g hg38.chrom.sizes -bga -strand - > hg38_L1HS_minus.bedgraph
bedSort hg38_L1HS_plus.bedgraph hg38_L1HS_plus.bedgraph
bedSort hg38_L1HS_minus.bedgraph hg38_L1HS_minus.bedgraph
bedGraphToBigWig hg38_L1HS_plus.bedgraph hg38.chrom.sizes hg38_L1HS_plus.bw
bedGraphToBigWig hg38_L1HS_minus.bedgraph hg38.chrom.sizes hg38_L1HS_minus.bw`


`##compute matrix for Alu
BED="only_longest_transcripts_intermediate_exons_by_strand"
echo "compute matrix for Alu.. "
computeMatrix scale-regions -S ../genomes/hg38_Alu_plus.bw \
                                                                ../genomes/hg38_Alu_minus.bw \
                                                          -R only_longest_transcripts_intermediate_exons_plus.bed \
                                                             only_longest_transcripts_intermediate_exons_minus.bed  \
                                                          --binSize 5 \
                                                         --missingDataAsZero \
                                                         --beforeRegionStartLength 2000 \
                                                         --regionBodyLength 300 \
                                                         --afterRegionStartLength 2000 \
                                                         -o "matrix_"$BED"_Alu.mat.gz"
echo "plotting heatmap for Alu .. "
plotHeatmap -m "matrix_"$BED"_Alu.mat.gz" \
                         -out $BED"_Alu.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                        --endLabel "RE"
##compute matrix for L1
echo "compute matrix for L1HS.. "$i
computeMatrix scale-regions -S ../genomes/hg38_L1HS_plus.bw \
                                                                ../genomes/hg38_L1HS_minus.bw \
                                           -R only_longest_transcripts_intermediate_exons_plus.bed \
                                                 only_longest_transcripts_intermediate_exons_minus.bed  \
                                            --binSize 5 \
                                            --missingDataAsZero \
                                            --beforeRegionStartLength 10000 \
                                            --regionBodyLength 300 \
                                            --afterRegionStartLength 10000 \
                                            -o "matrix_"$BED"_L1HS.mat.gz"
echo "plotting heatmap for L1HS .. "
plotHeatmap -m "matrix_"$BED"_L1HS.mat.gz" \
                         -out $BED"_L1HS.png" \
                         --startLabel "RS" \
                         --legendLocation none \
                         --labelRotation 45 \
                          --endLabel "RE"`
                        
                        

bedtools window -a only_longest_transcripts_intermedite_exons.bed -b ../genomes/hg38_repeats_SINEs.bed -bed -w 2000 > only_longest_transcripts_intermedite_exons_window2000_SINEs.bed

bedtools window -a only_longest_transcripts_intermedite_exons.bed -b ../genomes/hg38_repeats_L1.bed -bed -w 2000 > only_longest_transcripts_intermedite_exons_window2000_LINEs.bed

bedtools window -a only_longest_transcripts_intermedite_exons.bed -b ../genomes/homer_results_human_corrected.bed.gz -bed -w 2000 > only_longest_transcripts_intermedite_exons_window2000_L1motifs.bed
                        
                        
bedSort only_longest_transcripts_intermedite_exons.bed only_longest_transcripts_intermedite_exons.bed                          
bedtools closest -D a -io -a only_longest_transcripts_intermedite_exons.bed -b ../genomes/hg38_repeats_SINEs.bed > only_longest_transcripts_intermedite_exons_dist_SINEs.bed
bedtools closest -D a -io -a only_longest_transcripts_intermedite_exons.bed -b ../genomes/hg38_repeats_Alu.bed > only_longest_transcripts_intermedite_exons_dist_Alu.bed
bedtools closest -D a -io -a only_longest_transcripts_intermedite_exons.bed -b ../genomes/hg38_repeats_L1.bed > only_longest_transcripts_intermedite_exons_dist_L1.bed
bedtools closest -D a -io -a only_longest_transcripts_intermedite_exons.bed -b ../genomes/hg38_repeats_L1HS.bed > only_longest_transcripts_intermedite_exons_dist_L1HS.bed
bedtools closest -D a -io -a only_longest_transcripts_intermedite_exons.bed -b ../genomes/homer_results_human_corrected.bed.gz > only_longest_transcripts_intermedite_exons_dist_L1motif.bed

## 1000 genomes L1 polymorphic
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END<POLARITY; If there is only 5' OR 3' support for this call, will be NULL NULL for START and END">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT

zcat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | grep -v "##" | cut -f 1,2,3,4,5,6,7,8 | awk -v OFS="\t" '$5 ~/INS:ME/' > polymorphic_ME_1000G.tsv

` head polymorphic_ME_1000G.tsv 
1    645710    ALU_umary_ALU_2    A    <INS:ME:ALU>    .    .    AC=35;AF=0.00698882;AFR_AF=0;AMR_AF=0.0072;AN=5008;CS=ALU_umary;EAS_AF=0.0069;EUR_AF=0.0189;MEINFO=AluYa4_5,1,223,-;NS=2504;SAS_AF=0.0041;SITEPOST=0.9998;SVLEN=222;SVTYPE=ALU;TSD=null
1    812283    L1_umary_LINE1_1    G    <INS:ME:LINE1>    .    .    AC=58;AF=0.0115815;AFR_AF=0.0098;AMR_AF=0.0187;AN=5008;CS=L1_umary;EAS_AF=0.0109;EUR_AF=0.0179;MEINFO=LINE1,2926,3363,+;NS=2504;SAS_AF=0.0031;SITEPOST=1;SVLEN=437;SVTYPE=LINE1;TSD=null
1    813866    L1_umary_LINE1_2    C    <INS:ME:LINE1>    .    .    AC=9;AF=0.00179712;AFR_AF=0.003;AMR_AF=0;AN=5008;CS=L1_umary;EAS_AF=0;EUR_AF=0.001;MEINFO=LINE1,4049,4300,+;NS=2504;SAS_AF=0.0041;SITEPOST=1;SVLEN=251;SVTYPE=LINE1;TSD=null
1    1517860    SVA_umary_SVA_1    A    <INS:ME:SVA>    .    .    AC=3;AF=0.00059904;AFR_AF=0.0015;AMR_AF=0;AN=5008;CS=SVA_umary;EAS_AF=0.001;EUR_AF=0;MEINFO=SVA,44,394,+;NS=2504;SAS_AF=0;SITEPOST=0.9899;SVLEN=350;SVTYPE=SVA;TSD=null
1    3011887    ALU_umary_ALU_3    A    <INS:ME:ALU>    .    .    AC=63;AF=0.0125799;AFR_AF=0.0439;AMR_AF=0.0072;AN=5008;CS=ALU_umary;EAS_AF=0;EUR_AF=0;MEINFO=AluYa5,1,281,-;NS=2504;SAS_AF=0;SITEPOST=1;SVLEN=280;SVTYPE=ALU;TSD=AGAAAGTGGAGTA
1    3703963    ALU_umary_ALU_4    A    <INS:ME:ALU>    .    .    AC=1;AF=0.00019968;AFR_AF=0.0008;AMR_AF=0;AN=5008;CS=ALU_umary;EAS_AF=0;EUR_AF=0;MEINFO=AluUndef,1,281,-;NS=2504;SAS_AF=0;SITEPOST=1;SVLEN=280;SVTYPE=ALU;TSD=null
1    3995268    ALU_umary_ALU_5    N    <INS:ME:ALU>    .    .    AC=3934;AF=0.785543;AFR_AF=0.7806;AMR_AF=0.804;AN=5008;CS=ALU_umary;EAS_AF=0.7609;EUR_AF=0.7992;MEINFO=AluUndef,89,280,-;NS=2504;SAS_AF=0.7904;SITEPOST=0.9891;SVLEN=191;SVTYPE=ALU;TSD=null
1    4288465    ALU_umary_ALU_6    T    <INS:ME:ALU>    .    .    AC=496;AF=0.0990415;AFR_AF=0.1861;AMR_AF=0.0677;AN=5008;CS=ALU_umary;EAS_AF=0.0079;EUR_AF=0.0954;MEINFO=AluYa5,1,281,-;NS=2504;SAS_AF=0.1012;SITEPOST=1;SVLEN=280;SVTYPE=ALU;TSD=TTTTCTCACCCTTCT
1    4532042    ALU_umary_ALU_7    G    <INS:ME:ALU>    .    .    AC=62;AF=0.0123802;AFR_AF=0.0454;AMR_AF=0.0029;AN=5008;CS=ALU_umary;EAS_AF=0;EUR_AF=0;MEINFO=AluUndef,1,281,-;NS=2504;SAS_AF=0;SITEPOST=1;SVLEN=280;SVTYPE=ALU;TSD=null
1    5047126    ALU_umary_ALU_8    T    <INS:ME:ALU>    .    .    AC=15;AF=0.00299521;AFR_AF=0.0106;AMR_AF=0.0014;AN=5008;CS=ALU_umary;EAS_AF=0;EUR_AF=0;MEINFO=AluYa4_5,4,281,-;NS=2504;SAS_AF=0;SITEPOST=1;SVLEN=277;SVTYPE=ALU;TSD=null`


