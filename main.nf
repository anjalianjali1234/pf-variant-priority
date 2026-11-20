#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * pf-variant-priority
 * Targeted Nanopore amplicon reads -> annotated protein variants -> ESM scoring ->
 * calibration -> ranked watchlist.
 *
 * Alignment and variant calling follow the approach of de Cesare et al. (Nat Commun 2024):
 * minimap2 map-ont, then Clair3 with an ONT model. Consequence annotation uses
 * bcftools csq, which needs a GFF with CDS features matching the reference.
 */

params.reads     = null
params.reference = null
params.gff       = null
params.panel     = "${projectDir}/data/pfmdr15_panel.tsv"
params.proteins  = "${projectDir}/data/panel_proteins.faa"
params.validated = "${projectDir}/data/who_validated_k13.tsv"
params.config    = "${projectDir}/conf/priority.yml"
params.outdir    = "results"
params.clair3_model = "r1041_e82_400bps_sup_v430"
params.esm_checkpoint = "facebook/esm2_t33_650M_UR50D"
params.min_depth = 20

def requireParam(value, name) {
    if (!value) { exit 1, "Missing required parameter --${name}. See README.md" }
    return file(value, checkIfExists: true)
}

process ALIGN {
    tag   "${sample_id}"
    cpus  4
    label 'align'

    input:
    tuple val(sample_id), path(reads)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads} \\
      | samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
}

process CALL_VARIANTS {
    tag   "${sample_id}"
    cpus  4
    label 'call'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path fai

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    run_clair3.sh \\
      --bam_fn=${bam} \\
      --ref_fn=${reference} \\
      --threads=${task.cpus} \\
      --platform=ont \\
      --model_path=\${CONDA_PREFIX}/bin/models/${params.clair3_model} \\
      --include_all_ctgs \\
      --no_phasing_for_fa \\
      --output=clair3_out
    cp clair3_out/merge_output.vcf.gz ${sample_id}.vcf.gz
    """
}

process ANNOTATE {
    tag       "${sample_id}"
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)
    path reference
    path gff

    output:
    path "${sample_id}.csq.vcf.gz"

    script:
    """
    bcftools norm -f ${reference} -m -any ${vcf} -Ou \\
      | bcftools csq -f ${reference} -g ${gff} --force -Oz -o ${sample_id}.csq.vcf.gz
    bcftools index ${sample_id}.csq.vcf.gz
    """
}

process COLLECT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path vcfs
    path panel

    output:
    path "variant_table.tsv"

    script:
    """
    collect_variants.py --vcf '*.csq.vcf.gz' --panel ${panel} \\
        --min-depth ${params.min_depth} --out variant_table.tsv
    """
}

process SCORE {
    label      'gpu'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path variants
    path proteins

    output:
    path "scored.tsv"

    script:
    """
    score_variants.py --variants ${variants} --fasta ${proteins} \\
        --checkpoint ${params.esm_checkpoint} --out scored.tsv
    """
}

process BENCHMARK {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scored
    path validated

    output:
    path "benchmark/*"

    script:
    """
    benchmark.py --scored ${scored} --validated ${validated} --out benchmark/
    """
}

process PRIORITISE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scored
    path cfg
    path validated

    output:
    path "watchlist*.tsv"

    script:
    """
    prioritise.py --scored ${scored} --config ${cfg} \\
        --validated ${validated} --out watchlist.tsv
    """
}

workflow {
    reference = requireParam(params.reference, 'reference')
    gff       = requireParam(params.gff, 'gff')
    requireParam(params.reads, 'reads')

    fai       = file("${reference}.fai", checkIfExists: true)
    panel     = file(params.panel,     checkIfExists: true)
    proteins  = file(params.proteins,  checkIfExists: true)
    validated = file(params.validated, checkIfExists: true)
    cfg       = file(params.config,    checkIfExists: true)

    read_ch = Channel
        .fromPath(params.reads, checkIfExists: true)
        .map { f -> tuple(f.simpleName, f) }

    bams     = ALIGN(read_ch, reference)
    vcfs     = CALL_VARIANTS(bams, reference, fai)
    csq      = ANNOTATE(vcfs, reference, gff)
    variants = COLLECT(csq.collect(), panel)
    scored   = SCORE(variants, proteins)

    BENCHMARK(scored, validated)
    PRIORITISE(scored, cfg, validated)
}

workflow.onComplete {
    log.info """
    ------------------------------------------------------------
    pf-variant-priority finished: ${workflow.success ? 'OK' : 'FAILED'}
    duration : ${workflow.duration}
    results  : ${params.outdir}
    Read ${params.outdir}/benchmark/benchmark.json BEFORE the watchlist.
    ------------------------------------------------------------
    """.stripIndent()
}
