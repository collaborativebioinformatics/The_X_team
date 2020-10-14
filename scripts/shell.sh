minimap2 -d human_hs37d5.mmi human_hs37d5.fasta

minimap2 -ax asm5 human_hs37d5.fasta GCA_004796485.2_HG002_CCS_canu_maternal_1.1_genomic.fa > GCA_004796485.2_HG002_CCS_canu_maternal_1.1_genomic.sam

minimap2 -ax asm5 human_hs37d5.fasta GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic.fa > GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic.sam

minimap2 -ax asm5 human_hs37d5.fasta HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap1.fasta > HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap1.sam

minimap2 -ax asm5 human_hs37d5.fasta HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap2.fasta > HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap2.sam

samtools view -bS GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic.sam | samtools sort - > GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic.bam

bedtools bamtobed -i GCA_004796485.2_HG002_CCS_canu_maternal_1.1_genomic.bam > GCA_004796485.2_HG002_CCS_canu_maternal_1.1_genomic_bed.bam

bedtools bamtobed -i GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic.bam > GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic_bed.bam

bedtools bamtobed -i HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap1.bam > HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap1_bed.bam

bedtools bamtobed -i HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap2.bam > HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap2_bed.bam