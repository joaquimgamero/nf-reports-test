params.outdir = "results"
params.delay = 10

Channel
    .value("${projectDir}/resources/MultiQC Report.html")
    .into { multiqc_html_ch }

Channel
    .from(1..4)
    .into { step_ch }

Channel
    .of("${projectDir}/resources/bin_depths_summary.tsv", "${projectDir}/resources/bin_summary.tsv", "${projectDir}/resources/CAPES_S7.log",
    "${projectDir}/resources/execution_trace_2021-07-29_07-12-59.txt", "${projectDir}/resources/kraken2_report.txt",
    "${projectDir}/resources/pipeline_dag_2021-07-29_07-12-59.svg", "${projectDir}/resources/SPAdesHybrid-CAPES_S7-binDepths.heatmap.png",
    "${projectDir}/resources/taxonomy.krona.html", "${projectDir}/resources/transposed_report.tex", "${projectDir}/resources/transposed_report.tsv",
    "${projectDir}/resources/transposed_report.txt", "${projectDir}/resources/genome.fasta", "${projectDir}/resources/all_sites.fas",
    "${projectDir}/resources/baits.bed", "${projectDir}/resources/genome.dict", "${projectDir}/resources/genome.fasta.fai",
    "${projectDir}/resources/genome.gff3", "${projectDir}/resources/genome.gtf", "${projectDir}/resources/genome.sizes",
    "${projectDir}/resources/proteome.fasta", "${projectDir}/resources/test.baserecalibrator.table", "${projectDir}/resources/test.bed",
    "${projectDir}/resources/test.bedgraph", "${projectDir}/resources/test.bigwig", "${projectDir}/resources/test.paired_end.bam",
    "${projectDir}/resources/test.single_end.bam.readlist.txt", "${projectDir}/resources/test.vcf", "${projectDir}/resources/test.vcf.gz.tbi",
    "${projectDir}/resources/test2.targets.tsv.gz", "${projectDir}/resources/transcriptome.paf", "${projectDir}/resources/report.pdf",
    "${projectDir}/resources/nfcore_chipseq110_samplesheet_test_full_6cols.csv")
    .set { resources_ch }

process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path multiqc_html from multiqc_html_ch
    val step from step_ch

    output:
    path("step_*/*.html") into reports

    script:
    """
    echo "Copying MultiQC reports"
    mkdir step_$step
    cp $multiqc_html step_$step/
    """
}

workflow REPORTS {
    take:
    val step from step_ch

    main:
    multiqc_html_ch.collect().set{ resources_ch }
    resources_ch.view()

    MULTIQC(step)

    emit:
    val reports from reports.collect()
}

workflow {
    REPORTS()
}
