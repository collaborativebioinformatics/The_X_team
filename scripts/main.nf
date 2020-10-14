oneRef = file(params.oneRef)
if(params.r && params.longMode){
    Channel
    .fromPath( params.longPath )
    .map { file -> tuple(file.baseName, file) }
    .ifEmpty { error "Oops! Cannot find any file matching: ${params.ref}"}
    .into { read_ch_pbhifi; read_ch_NGMLR; read_ch_NGMLR_CLR; read_ch_preprocessSniffles; read_ch_preprocessFai; read_CLR_ch_minimap; read_CCS_ch_minimap}
}

/*Given a SV vcf as input
5:22
we want to run Rocio's script to generate a query
5:22
then map the query to both haplotypes to make a bam
5:22
then run my script to evaluate the scores of each bam (e.g. a haplotype 1 bam and a haplotype 2 bam)
5:23
Rocio's script is GetVariantRef.py
5:23
min is ScoreAlignment.py*/

process map_assemblies {

}

process generate_queries_from_vcf{
  (GetVariantRef.py)
}

process map_queries{
  (minimap2)
}

process compare_boundaries{

}

process score{
  ScoreAlignment.py
}

process minimap2_pacbio_CLR{
    tag "${replicateId}"
    publishDir "{params.outputDir}/minimap2PacbioCLR", mode: 'copy'

    input:
    set replicateId, file(reads) from read_CLR_ch_minimap.mix(pb_optional_CLR_ch_minimap)
    file ref from oneRef

    output:
    set replicateId, file("${replicateId}.bam") into minimap2PacbioCLRReceiver

    script:
    """
    minimap2 -ax map-pb $ref $reads | samtools view -bS - | samtools sort - > ${replicateId}.bam
    """
}
