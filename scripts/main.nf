/*
 * -------------------------------------------------
 *   Nextflow main script
 *   Author: Chunxiao Liao
 * -------------------------------------------------
 * Pipeline management script.
 */

oneRef = file(params.oneRef)
vcfChannel = Channel.fromPath(params.vcfPath)
contigChannel = Channel.fromPath(params.contigPath)


process map_assemblies {

}

//to run Rocio's script to generate a query
process generate_queries_from_vcf{
  tag "${replicateId}"
  publishDir "{params.outputDir}/generateQueries", mode: 'copy'

  input:
    set replicateId, file(vcf) from vcfChannel
    file ref from oneRef

  output:
    set replicateId, file("cuteSV_v108_hg37_hg002_g4015_8reads_variants_v2.fasta") into mapQueriesReceiver

  script:
    """
    gunzip $vcf
    python params.program/GetVariantRef.py --vcf cuteSV_v108_hg37_hg002_g4015_8reads.vcf --ref $ref > cuteSV_v108_hg37_hg002_g4015_8reads_variants_v2.fasta
    python params.program/parseHapBAMs.py --hap1 The_X_team:/minimap2_output/SV2/SV2_paternal.sorted.bam --hap2 The_X_team:/minimap2_output/SV2/SV2_maternal.sorted.bam --varfasta cuteSV_v108_hg37_hg002_g4015_8reads_variants_v2.fasta --fp fp-noBed_default_noPASS_refPASS.vcf > table-Len-Allen-all-v2.tsv
    """ 
}

// priya's code
process map_queries{
  tag "${replicateId}"
  publishDir "{params.outputDir}/mapQueries", mode: 'copy'

  input:
    set replicateId, file(sv4.fasta) from mapQueriesReceiver
    file ref from oneRef

  output:
    set replicateId, file("${replicateId}.bam") into mapQueryReceiver

  script:
    """
    minimap2 -d ref.mmi $ref 
    minimap2 -a ref.mmi $sv4.fasta | samtools view -bS - | samtools sort - > ${replicateId}.bam
    """
}

// just merge results
process merge_map_query {
    input:
        file "*" from mapQueryReceiver.collect()

    output:
        file "mergeMapQuery" into mergeMapQueryReceiver

    script:
    """
    mkdir mergeMapQuery;
    scp *.bam mergeMapQuery;
    """
}

process compare_boundaries{

}

//to evaluate the scores of each bam (e.g. a haplotype 1 bam and a haplotype 2 bam)
process score{
  tag "${replicateId}"
  publishDir "{params.outputDir}/score", mode: 'copy'

  input:
    set replicateId, file(bam1), file(bam2) from mergeMapQueryReceiver

  output:
    set replicateId, file("MaternalPaternal.score") into scoreReceiver

  script:
    """
    python ScoreAlignment.py $bam1 $bam2 maternal.score paternal.score

    ./CreateJointSupportTab.sh maternal.score paternal.score MaternalPaternal.score
    """
  
}

// so then you would run Rocio's script to annotate Truvari
process annotation_support{

  tag "${replicateId}"
  publishDir "{params.outputDir}/annotation", mode: 'copy'

  input:
    set replicateId, file(score) from scoreReceiver

  output:
    set replicateId, file("MaternalPaternal-Extraannot-withHAPcoords.tab") into scoreReceiver

  script:
  """"
  python $params.program/AnnotateJointSupportTab.py -tab MaternalPaternal.tab -truvariFP fp-noBed_default_noPASS_refPASS.vcf -truvariTP tp-noBed_default_noPASS_refPASS.vcf -hap1 SV4_paternal_primary.sorted.bam -hap2 SV4_maternal_primary.sorted.bam > MaternalPaternal-Extraannot-withHAPcoords.tab
  """"
}

process map_contig_ref{
    tag "${replicateId}"
    publishDir "{params.outputDir}/mapContigRef", mode: 'copy'

    input:
    set replicateId, file(contig) from contigChannel
    file ref from oneRef

    output:
    set replicateId, file("${replicateId}_filtered.bed") into mapContigRefReceiver

    script:
    """
    minimap2 -ax asm5 $ref $contig | samtools view -bS - | samtools sort - > ${replicateId}.bam
    bedtools bamtobed -i ${replicateId}.bam > ${replicateId}.bed
    python $params.program/bedParser.py ${replicateId}.bed ${replicateId}_filtered.bed
    """
}
