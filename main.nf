#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }          from './modules/snippy.nf'
include { snippy }         from './modules/snippy.nf'
include { count_variants } from './modules/snippy.nf'
include { qualimap_bamqc } from './modules/snippy.nf'
include { samtools_mpileup } from './modules/snippy.nf'
include { plot_coverage } from './modules/snippy.nf'
include { convert_gb_to_fa } from './modules/snippy.nf'

workflow {
  ch_ref = Channel.fromPath( "${params.ref}", type: 'file')

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  main:
    fastp(ch_fastq)
    snippy(fastp.out.reads.combine(ch_ref))
    qualimap_bamqc(snippy.out.alignment)
    count_variants(snippy.out.variants_csv.combine(ch_ref))

    if (file(params.ref).name =~ /.+\.[Gg]b.*/){
      ch_ref_fa = convert_gb_to_fa(ch_ref)
    }else {
      ch_ref_fa = ch_ref
    }

    samtools_mpileup(snippy.out.alignment.combine(ch_ref_fa))
    plot_coverage(samtools_mpileup.out.depths.combine(ch_ref_fa))
}
