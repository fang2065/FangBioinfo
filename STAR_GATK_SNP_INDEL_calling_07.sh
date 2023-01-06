#!/bin/bash
#SBATCH --ntasks-per-node=20
#SBATCH -p siberia 
#SBATCH --time=0-23:00:00
#SBATCH --error=zf.%J.error
#SBATCH --output=zf.%J.out

SAMPLES="07"

for SAMPLE in $SAMPLES; 
do  
    /home/USSR/zf250/software/STAR-2.7.10a/bin/Linux_x86_64_static/STAR \
         --genomeDir /home/USSR/zf250/index/star_index \
         --runThreadN 10 \
         --sjdbOverhang 149 \
         --outFileNamePrefix /home/USSR/zf250/bam/${SAMPLE}_ \
         --outSAMtype BAM Unsorted \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --twopassMode Basic \
         --readFilesCommand gunzip -c \
         --readFilesIn /home/USSR/zf250/input/baseline/${SAMPLE}_00_R1.fastq.gz /home/USSR/zf250/input/baseline/${SAMPLE}_00_R2.fastq.gz


    /home/USSR/zf250/software/miniconda/pkgs/sambamba-0.6.6-2/bin/sambamba sort \
         -o /home/USSR/zf250/bam/${SAMPLE}_sorted.bam \
         /home/USSR/zf250/bam/${SAMPLE}_Aligned.out.bam

     rm /home/USSR/zf250/bam/${SAMPLE}_Aligned.out.bam
     rm -r /home/USSR/zf250/bam/${SAMPLE}__STARpass1
     rm -r /home/USSR/zf250/bam/${SAMPLE}__STARtmp
     rm -r /home/USSR/zf250/bam/${SAMPLE}__STARgenome

    /home/USSR/zf250/software/miniconda/pkgs/sambamba-0.6.6-2/bin/sambamba markdup \
         --overflow-list-size 600000 \
         --tmpdir='./' \
         -r /home/USSR/zf250/bam/${SAMPLE}_sorted.bam \
         /home/USSR/zf250/bam/${SAMPLE}_rmd.bam

     rm /home/USSR/zf250/bam/${SAMPLE}_sorted.bam
     rm /home/USSR/zf250/bam/${SAMPLE}_sorted.bam.bai


    /home/USSR/zf250/software/gatk-4.2.6.1/gatk SplitNCigarReads -R /home/USSR/zf250/reference/GRCh38.p13.genome.fa \
         -I /home/USSR/zf250/bam/${SAMPLE}_rmd.bam \
         -O /home/USSR/zf250/bam/${SAMPLE}_rmd_split.bam
     
     rm /home/USSR/zf250/bam/${SAMPLE}_rmd.bam
     rm /home/USSR/zf250/bam/${SAMPLE}_rmd.bam.bai


     /home/USSR/zf250/software/gatk-4.2.6.1/gatk AddOrReplaceReadGroups -I /home/USSR/zf250/bam/${SAMPLE}_rmd_split.bam \
         -O  /home/USSR/zf250/bam/${SAMPLE}_rmd_split_add.bam \
         -LB ${SAMPLE} \
         -PL ILLUMINA \
         -PU ${SAMPLE} \
         -SM ${SAMPLE}

     rm /home/USSR/zf250/bam/${SAMPLE}_rmd_split.bam
     rm /home/USSR/zf250/bam/${SAMPLE}_rmd_split.bam.bai


    /home/USSR/zf250/software/gatk-4.2.6.1/gatk BaseRecalibrator \
         -I /home/USSR/zf250/bam/${SAMPLE}_rmd_split_add.bam \
         -R /home/USSR/zf250/reference/GRCh38.p13.genome.fa \
         --known-sites /home/USSR/zf250/reference/GATK/Homo_sapiens_assembly38.dbsnp138.vcf \
         --known-sites /home/USSR/zf250/reference/GATK/Homo_sapiens_assembly38.known_indels.vcf.gz \
         --known-sites /home/USSR/zf250/reference/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
         -O /home/USSR/zf250/bam/${SAMPLE}_recal.table


    /home/USSR/zf250/software/gatk-4.2.6.1/gatk ApplyBQSR \
         --bqsr-recal-file /home/USSR/zf250/bam/${SAMPLE}_recal.table \
         -R /home/USSR/zf250/reference/GRCh38.p13.genome.fa \
         -I /home/USSR/zf250/bam/${SAMPLE}_rmd_split_add.bam \
         -O /home/USSR/zf250/bam/${SAMPLE}_recal.bam
     
     rm /home/USSR/zf250/bam/${SAMPLE}_rmd_split_add.bam
     rm /home/USSR/zf250/bam/${SAMPLE}_rmd_split_add.bam.bai

    /home/USSR/zf250/software/gatk-4.2.6.1/gatk HaplotypeCaller \
         --native-pair-hmm-threads 10 \
         -R /home/USSR/zf250/reference/GRCh38.p13.genome.fa \
         --dbsnp /home/USSR/zf250/reference/GATK/Homo_sapiens_assembly38.dbsnp138.vcf \
         -I /home/USSR/zf250/bam/${SAMPLE}_recal.bam \
         -O /home/USSR/zf250/vcf/${SAMPLE}.vcf \
         --standard-min-confidence-threshold-for-calling 30 \
         --dont-use-soft-clipped-bases

     rm /home/USSR/zf250/bam/${SAMPLE}_recal.bam
     rm /home/USSR/zf250/bam/${SAMPLE}_recal.bam.bai

     /home/USSR/zf250/software/gatk-4.2.6.1/gatk VariantFiltration \
          -R /home/USSR/zf250/reference/GRCh38.p13.genome.fa \
          -V /home/USSR/zf250/vcf/${SAMPLE}.vcf \
          -O /home/USSR/zf250/vcf/${SAMPLE}_filtered.vcf \
          --window 35 \
          --cluster 3 \
          --filter-name "FS" \
          --filter "FS > 30.0" \
          --filter-name "QD" \
          --filter "QD < 2.0"
          
     rm /home/USSR/zf250/vcf/${SAMPLE}.vcf

     /home/USSR/zf250/software/gatk-4.2.6.1/gatk SelectVariants \
          --select-type-to-include SNP \
          -R /home/USSR/zf250/reference/GRCh38.p13.genome.fa \
          -V /home/USSR/zf250/vcf/${SAMPLE}_filtered.vcf \
          -O /home/USSR/zf250/vcf/${SAMPLE}_fil_snp.vcf

     
     /home/USSR/zf250/software/gatk-4.2.6.1/gatk SelectVariants \
          --select-type-to-include INDEL \
          -R /home/USSR/zf250/reference/GRCh38.p13.genome.fa \
          -V /home/USSR/zf250/vcf/${SAMPLE}_filtered.vcf \
          -O /home/USSR/zf250/vcf/${SAMPLE}_fil_indel.vcf
     
     rm /home/USSR/zf250/vcf/${SAMPLE}_filtered.vcf

     /home/USSR/zf250/software/annovar/convert2annovar.pl \
         --format vcf4 \
         /home/USSR/zf250/vcf/${SAMPLE}_fil_snp.vcf > /home/USSR/zf250/vcf/${SAMPLE}_fil_snp.avinput

     /home/USSR/zf250/software/annovar/table_annovar.pl \
         /home/USSR/zf250/vcf/${SAMPLE}_fil_snp.avinput \
         /home/USSR/zf250/software/annovar/humandb/ \
         --buildver hg38 \
         -out /home/USSR/zf250/vcf/${SAMPLE}_fil_snp.anno \
         -protocol refGene \
         -operation g \
         -remove \
         -nastring '.'

     /home/USSR/zf250/software/annovar/convert2annovar.pl \
         --format vcf4 \
         /home/USSR/zf250/vcf/${SAMPLE}_fil_indel.vcf > /home/USSR/zf250/vcf/${SAMPLE}_fil_indel.avinput

     /home/USSR/zf250/software/annovar/table_annovar.pl \
         /home/USSR/zf250/vcf/${SAMPLE}_fil_indel.avinput \
         /home/USSR/zf250/software/annovar/humandb/ \
         --buildver hg38 \
         -out /home/USSR/zf250/vcf/${SAMPLE}_fil_indel.anno \
         -protocol refGene \
         -operation g \
         -remove \
         -nastring '.'
done

