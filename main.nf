#!/usr/bin/env nextflow 

nextflow.enable.dsl=2


params.input_folder = "${workflow.launchDir}"

params.ouput_folder = "${params.input_folder}/trimmed"



workflow {

  if (params.se) {

    fq = Channel.fromFilePairs("${params.input_folder}/*.fq.gz", size: 1, flat:true)
    .concat(Channel.fromFilePairs("${params.input_folder}/*.fastq.gz", size: 1, flat:true))

    fq | fastp_se

    fastp_se.out.fq_trimmed.collectFile(newLine: true, storeDir: "${params.input_folder}") { ID, fq ->
      [ "r_file.txt", "$ID\t${params.ouput_folder}/$fq"]
      }
    
    fastp_se.out.fastp_json.collect() | multi_QC

    } else {

    fq = Channel.fromFilePairs("${params.input_folder}/*_{1,2}.fq.gz", flat: true)
    .concat(Channel.fromFilePairs("${params.input_folder}/*_{R1,R2}.fq.gz", flat: true))
    .concat(Channel.fromFilePairs("${params.input_folder}/*_{1P,2P}.fq.gz", flat: true))
    .concat(Channel.fromFilePairs("${params.input_folder}/*_{1,2}.fastq.gz", flat: true))
    .concat(Channel.fromFilePairs("${params.input_folder}/*_{R1,R2}.fastq.gz", flat: true))
    .concat(Channel.fromFilePairs("${params.input_folder}/*_{1P,2P}.fastq.gz", flat: true))

    fq | fastp_pe
    
    fastp_pe.out.fq_trimmed.collectFile(newLine: true, storeDir: "${params.input_folder}") { ID, fq1, fq2 ->
      [ "r_file.txt", "$ID\t${params.ouput_folder}/$fq1\t${params.ouput_folder}/$fq2"]
      }
    
    fastp_pe.out.fastp_json.collect() | multi_QC

    }

//  (params.se ? fastp_se : fastp_pe) | 
}




process fastp_pe {

    publishDir "${params.ouput_folder}", mode: 'copy', pattern: "*.fq.gz"

    publishDir "${params.ouput_folder}/multi_QC", mode: 'copy', pattern: "*_fastp.html"

    input:
      tuple val(sampleID), path(fq1), path(fq2) 

    output:
      tuple val(sampleID), val("${sampleID}_R1_trimmed.fq.gz"), val("${sampleID}_R2_trimmed.fq.gz"), emit: fq_trimmed
      path "*_trimmed.fq.gz"
      path "*_fastp.json", emit: fastp_json
      path "*_fastp.html"

    """
    fastp --in1 $fq1 --in2 $fq2 \\
          --out1 ${sampleID}_R1_trimmed.fq.gz --out2 ${sampleID}_R2_trimmed.fq.gz \\
          --length_required 20 \\
          --detect_adapter_for_pe \\
          -j ${sampleID}_fastp.json -h ${sampleID}_fastp.html
    """
}





process fastp_se {

    publishDir "${params.ouput_folder}", mode: 'copy', pattern: "*.fq.gz"

    publishDir "${params.ouput_folder}/multi_QC", mode: 'copy', pattern: "*_fastp.html"

    input:
      tuple val(sampleID), path(fq) 

    output:
      tuple val(sampleID), val("${sampleID}_trimmed.fq.gz"), emit: fq_trimmed
      path "*_trimmed.fq.gz"
      path "*_fastp.json", emit: fastp_json
      path "*_fastp.html"

    """
    fastp --in1 $fq \\
          --out1 ${sampleID}_trimmed.fq.gz \\
          --length_required 20 \\
          -j ${sampleID}_fastp.json -h ${sampleID}_fastp.html
    """
}



process multi_QC {

    publishDir "${params.ouput_folder}/multi_QC", mode: 'move'

    input:
    path(json) 

    output:
    path "*.html"
    path "multiqc_data/*"
    
    """
    multiqc . --data-format tsv
    """
  }
