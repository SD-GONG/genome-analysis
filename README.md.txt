### nextdenovo 
ls ../all_reads/20220419-UNL273-P6-PAK06069.pass.fastq.gz > ./input.fofn
<run.cfg>
[General]
job_type = local
job_prefix = nextDenovo
task = all # 'all', 'correct', 'assemble'
rewrite = yes # yes/no
deltmp = yes
parallel_jobs = 22
input_type = raw
read_type = clr
input_fofn = ./input.fofn
workdir = ./01_rundir

[correct_option]
read_cutoff = 1k
genome_size = 1.1g
pa_correction = 5
sort_options = -m 30g -t 35
minimap2_options_raw = -t 20
correction_options = -p 35

[assemble_option]
minimap2_options_cns = -t 20
nextgraph_options = -a 1

see https://nextdenovo.readthedocs.io/en/latest/OPTION.html for a detailed introduction about all the parameters

<run nextdenovo>
../NextDenovo/nextDenovo run.cfg

### nextpolish
ls ../clean_1.fq.gz  ../clean_2.fq.gz > nextpolish.sgs.fofn
ls ../CLR_subreads.fasta > nextpolish.lgs.fofn

<run.cfg>
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = ./nd.asm.fasta
genome_size = auto
workdir = ./01_rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./nextpolish.sgs.fofn
sgs_options = -max_depth 100 -bwa

[lgs_option]
lgs_fofn = ./nextpolish.lgs.fofn
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-pb

see https://nextpolish.readthedocs.io/en/latest/OPTION.html for a detailed introduction about all the parameters

<run nextPolish>
~/gene_assembly/own_job/correct/nextpolish/NextPolish/nextPolish run.cfg

###juicer+3d-dna
<juicer>
bwa index genome.nextpolish.fasta
python ~/gene_assembly/3d-dna/juicer-1.6/misc/generate_site_positions.py DpnII genome genome.nextpolish.fasta
awk 'BEGIN{OFS="\t"}{print $1, $NF}' genome_DpnII.txt > genome.chrom.sizes
~/gene_assembly/3d-dna/juicer-1.6/scripts/juicer.sh -g genome -s DpnII -z reference/genome.nextpolish.fasta -y reference/genome_DpnII.txt -p reference/genome.chrom.sizes -D ~/gene_assembly/3d-dna/juicer-1.6 -t 120
<3d-dna>
~/software/3d-dna/3d-dna-master/run-asm-pipeline.sh -r 2 ./reference/genome.nextpolish.fasta ../aligned/merged_nodups.txt > 3d.log


###purge_haplotigs
mv genome.nextpolish.rawchrom.fasta > GG.contig.fasta
minimap2 -ax map-pb ./GG.contig.fa ./CLR_subreads.fasta | samtools view -hF 256 - |samtools sort -@ 8 -m 1G -o aligned.bam -T tmp.ali
purge_haplotigs  readhist  -b aligned.bam  -g ./GG.contig.fa  -t 20
purge_haplotigs contigcov -i aligned.bam.gencov -o coverage_stats.csv -l 20 -m 100 -h 190
purge_haplotigs purge -g ./GG.contig.fa -c coverage_stats.csv -b aligned.bam -t 60 -a 80



###juicer+3d-dna+juicebox
mv curated.fasta purge_haplotigs.fasta
<juicer>
mkdir reference
cd reference
bwa index purge_haplotigs.fasta
python ~/gene_assembly/3d-dna/juicer-1.6/misc/generate_site_positions.py DpnII genome purge_haplotigs.fasta
awk 'BEGIN{OFS="\t"}{print $1, $NF}' genome_DpnII.txt > genome.chrom.sizes
cd ..
~/gene_assembly/3d-dna/juicer-1.6/scripts/juicer.sh -g genome -s DpnII -z reference/purge_haplotigs.fasta -y reference/genome_DpnII.txt -p reference/genome.chrom.sizes -D ~/gene_assembly/3d-dna/juicer-1.6 -t 120
<3d-dna>
~/software/3d-dna/3d-dna-master/run-asm-pipeline.sh -r 0 ./reference/purge_haplotigs.fasta ./aligned/merged_nodups.txt > 3d.log
<juicebox>
sh ~/software/3d-dna/3d-dna-master/run-asm-pipeline-post-review.sh -r purge_haplotigs-review.assembly reference/purge_haplotigs.fasta aligned/merged_nodups.txt

### TGS-gapcloser
tgsgapcloser --scaff ./YG.chromosome.fa --reads ~/gene_assembly/QC/NanoFilt/NanoFilt_50k.fa --output yege --racon ~/miniconda3/envs/tgsgapcloser/bin/racon --r_round 2  >pipe.log 2>pipe.err

### LR-gapcloser
sh ~/software/LR_Gapcloser/LR_Gapcloser-master2/src/LR_Gapcloser.sh -i TGS-result.fa -l ../clean.NanoFilt_10k.fa -t 60 -s n -o ./minimap_result

### pilon

bwa index LR-result.fa
bwa mem -t 60 LR-result.fa ../clean_1.fq.gz ../clean_2.fq.gz | samtools view -Sb > align1.bam
samtools sort -@ 60 -o align1.sorted.bam align1.bam
samtools index -@ 60 align1.sorted.bam
picard MarkDuplicates REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=align1.sorted.bam OUTPUT=align1.sorted.markup.bam METRICS_FILE=align1.sorted.markup.bam.metrics
samtools view -@ 30 -q 30 -b align1.sorted.markup.bam > align1.sorted.markup.filter.bam
samtools index -@ 30 align1.sorted.markup.filter.bam
java -Xmx512G -jar ~/gene_assembly/correct/pilon/pilon-1.23.jar --genome LR-result.fa --frags align1.sorted.markup.filter.bam --fix bases,amb --output pilon_polished1 >pilon.log

bwa index  pilon_polished1.fasta
bwa mem -t 60 pilon_polished1.fasta ../clean_1.fq.gz ../clean_2.fq.gz| samtools view -Sb > align2.bam
samtools sort -@ 60 -o align2.sorted.bam align2.bam
samtools index -@ 60 align2.sorted.bam
picard MarkDuplicates REMOVE_DUPLICATES= false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=align2.sorted.bam OUTPUT=align2.sorted.markup.bam METRICS_FILE=align2.sorted.markup.bam.metrics
samtools view -@ 30 -q 30 -b align2.sorted.markup.bam > align2.sorted.markup.filter.bam
samtools index -@ 30 align2.sorted.markup.filter.bam
java -Xmx512G -jar ~/gene_assembly/correct/pilon/pilon-1.23.jar --genome pilon_polished1.fasta --frags align2.sorted.markup.filter.bam --fix bases,amb --output pilon_polished2  >pilon.log


### plot heatmap

bowtie2-build genome.fa genome
python ~/software/Hic-pro2/HiC-Pro_3.1.0/bin/utils/digest_genome.py genome.fa -r dpnii -o genome_dpnii.bed
samtools faidx genome.fa
awk '{print $1"\t"$2}' genome.fa.fai > genome.size
HIC-pro -c config-hicpro.txt -o analysisi -i data
hicConvertFormat -m sample1_100000_iced.matrix -o result --inputFormat hicpro --outputFormat h5 --bedFileHicpro sample1_100000_abs.bed
hicPlotMatrix -m result.h5 --rotationX 270 --log1p -out YlOrBr --colorMap YlOrBr --dpi 800



### SNV calling

bwa index -p yege ../YG.genome.fasta
bwa mem -t 60 yege ~/pacbio_data/clean_1.fq.gz ~/pacbio_data/clean_2.fq.gz | samtools view -Sb > align1.bam
samtools sort -@ 60 -o align1.sorted.bam align1.bam
samtools index -@ 60 align1.sorted.bam
picard MarkDuplicates REMOVE_DUPLICATES=true INPUT=align1.sorted.bamOUTPUT=align1.sorted.markup.bam METRICS_FILE=align1.sorted.markup.bam.metrics 
bcftools mpileup  align1.sorted.markup.bam -f ../YG.genome.fasta | bcftools call -mv -o raw.vcf 
vcftools --vcf raw.vcf --max-missing 0.9  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --minQ 30 --out  result.vcf


###busco

busco -i ../genome.fasta -o busco -l ~/data/eudicots_odb10 -m genome -c 20 --offline


### orthofinder

orthofinder -f orthofinder/ -t 60 -S blast

### species tree + divergence time

$vi bash.sh
for i in *.fa
do  
  muscle -in $i -out $i.1
  Gblocks $i.1 -b4=5 -b5=h -t=p -e=.2
  seqkit sort $i.1.2 > $i.1.2.3
  seqkit seq $i.1.2.3 -w 0 > $i.1.2.3.4
done
$sh bash.sh


mkdir new
mv *.4 new/
cd new/
paste -d " " *.4 > all.fa
sed 's/ //g' all.fa > all0.fa

raxmlHPC-AVX2 -T 60 -f a -x 123 -p 123 -N 100 -m PROTGAMMAIJTTF -o O.sativa -s all0.fa -n all.tre

$vi bash4.sh
for $i in *.fa
  muscle -in $i -out $i.1
  grep '>' $i|cut -c2- > $i.1.2
  ./faSomeRecords all_cds.fa $i.1.2 $i.1.2.3
  pal2nal.pl -nogap -nomismatch $i.1 $i.1.2.3 -output fasta > $i.1.2.3.4
  seqkit sort $i.1.2.3.4 > $i.1.2.3.4.5
  seqkit seq $i.1.2.3.4.5 -w 0 > $i.1.2.3.4.5.6
$sh bash4.sh

paste -d " " *.6 > all_cds.fa
sed 's/ //g' all_cds.fa > all0_cds.fa

java -cp EasyCodeML.jar SeqFormatConvert.seqFactory.SeqConverter -i ~/mcmctree/cds -o ~/mcmctree/cds -oF phylip 

mcmctree mcmctree.ct1
mkdir new/
mv mcmctree.ct1 all_final.phy out.BV RAxML_bestTree.all.tre new  
cd new
mv out.BV in.BV
mcmctree mcmctree.ct2 


###cafe5

cafe5 -i gene_families_filter.txt -t tree.txt -p -k 3 -o k3

###JCVI

python -m jcvi.compara.catalog ortholog P.lobata P.thomsonii --no_strip_names
python -m jcvi.compara.catalog ortholog P.lobata P.montana --no_strip_names
python -m jcvi.compara.catalog ortholog P.thomsonii P.montana --no_strip_names

