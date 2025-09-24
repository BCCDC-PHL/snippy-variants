#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }                      from './modules/snippy.nf'
include { snippy }                     from './modules/snippy.nf'
include { count_variants }             from './modules/snippy.nf'
include { qualimap_bamqc }             from './modules/snippy.nf'
include { pipeline_provenance }        from './modules/provenance.nf'
include { collect_provenance }         from './modules/provenance.nf'
include { hash_files }                 from './modules/hash_files.nf'

workflow {
  ch_workflow_metadata = Channel.value([
	  workflow.sessionId,
	  workflow.runName,
	  workflow.manifest.name,
	  workflow.manifest.version,
	  workflow.start,
  ])

  ch_ref = Channel.fromPath( "${params.ref}", type: 'file')

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  main:
    hash_files(ch_fastq.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq-input")))

    fastp(ch_fastq)
    snippy(fastp.out.reads.combine(ch_ref))
    qualimap_bamqc(snippy.out.alignment)
    count_variants(snippy.out.variants_csv.combine(ch_ref))


    // Collect Provenance
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // At each step, we add another provenance file to the list using the << operator...
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_sample_ids = ch_fastq.map{ it -> it[0] }
    ch_provenance = ch_sample_ids
    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->     [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it ->     [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it ->          [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(snippy.out.provenance).map{ it ->     [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(qualimap_bamqc.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)

}
