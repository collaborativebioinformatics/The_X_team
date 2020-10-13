oneRef = file(params.oneRef)
if(params.r && params.longMode){
    Channel
    .fromPath( params.longPath )
    .map { file -> tuple(file.baseName, file) }
    .ifEmpty { error "Oops! Cannot find any file matching: ${params.ref}"}
    .into { read_ch_pbhifi; read_ch_NGMLR; read_ch_NGMLR_CLR; read_ch_preprocessSniffles; read_ch_preprocessFai; read_CLR_ch_minimap; read_CCS_ch_minimap}
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
