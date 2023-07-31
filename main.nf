params.outdir = "results"
params.delay = 10

Channel
    .value("${projectDir}/resources/MultiQC Report.html")
    .into { multiqc_html_ch }

Channel
    .from(1..4)
    .into { step_ch }

process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path multiqc_html from multiqc_html_ch.first()
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

workflow {
    path bin_depths_summary_tsv from Channel.of("${projectDir}/resources/bin_depths_summary.tsv")
    path bin_summary_tsv from Channel.of("${projectDir}/resources/bin_summary.tsv")
    path CAPES_S7_log from Channel.of("${projectDir}/resources/CAPES_S7.log")
    path execution_trace_txt from Channel.of("${projectDir}/resources/execution_trace_2021-07-29_07-12-59.txt")
    path kraken2_report_txt from Channel.of("${projectDir}/resources/kraken2_report.txt")
    path pipeline_dag_svg from Channel.of("${projectDir}/resources/pipeline_dag_2021-07-29_07-12-59.svg")
    path binDepths_png from Channel.of("${projectDir}/resources/SPAdesHybrid-CAPES_S7-binDepths.heatmap.png")
    path krona_html from Channel.of("${projectDir}/resources/taxonomy.krona.html")
    path transposed_tex from Channel.of("${projectDir}/resources/transposed_report.tex")
    path transposed_tsv from Channel.of("${projectDir}/resources/transposed_report.tsv")
    path transposed_txt from Channel.of("${projectDir}/resources/transposed_report.txt")
    path genome_fasta from Channel.of("${projectDir}/resources/genome.fasta")
    path all_sites_fas from Channel.of("${projectDir}/resources/all_sites.fas")
    path baits_bed from Channel.of("${projectDir}/resources/baits.bed")
    path genome_dict from Channel.of("${projectDir}/resources/genome.dict")
    path genome_fasta_fai from Channel.of("${projectDir}/resources/genome.fasta.fai")
    path genome_gff3 from Channel.of("${projectDir}/resources/genome.gff3")
    path genome_gtf from Channel.of("${projectDir}/resources/genome.gtf")
    path genome_sizes from Channel.of("${projectDir}/resources/genome.sizes")
    path proteome_fasta from Channel.of("${projectDir}/resources/proteome.fasta")
    path test_baserecalibrator_table from Channel.of("${projectDir}/resources/test.baserecalibrator.table")
    path test_bed from Channel.of("${projectDir}/resources/test.bed")
    path test_bedgraph from Channel.of("${projectDir}/resources/test.bedgraph")
    path test_bigwig from Channel.of("${projectDir}/resources/test.bigwig")
    path test_paired_end_bam from Channel.of("${projectDir}/resources/test.paired_end.bam")
    path test_single_end_bam_readlist_txt from Channel.of("${projectDir}/resources/test.single_end.bam.readlist.txt")
    path test_vcf from Channel.of("${projectDir}/resources/test.vcf")
    path test_vcf_gz_tbi from Channel.of("${projectDir}/resources/test.vcf.gz.tbi")
    path test2_targets_tsv_gz from Channel.of("${projectDir}/resources/test2.targets.tsv.gz")
    path transcriptome_paf from Channel.of("${projectDir}/resources/transcriptome.paf")
    path report_pdf from Channel.of("${projectDir}/resources/report.pdf")
    path samplesheet_csv from Channel.of("${projectDir}/resources/nfcore_chipseq110_samplesheet_test_full_6cols.csv")

    tuple(
        path bin_depths_summary_tsv,
        path bin_summary_tsv,
        path CAPES_S7_log,
        path execution_trace_txt,
        path kraken2_report_txt,
        path pipeline_dag_svg,
        path binDepths_png,
        path krona_html,
        path transposed_tex,
        path transposed_tsv,
        path transposed_txt,
        path genome_fasta,
        path all_sites_fas,
        path baits_bed,
        path genome_dict,
        path genome_fasta_fai,
        path genome_gff3,
        path genome_gtf,
        path genome_sizes,
        path proteome_fasta,
        path test_baserecalibrator_table,
        path test_bed,
        path test_bedgraph,
        path test_bigwig,
        path test_paired_end_bam,
        path test_single_end_bam_readlist_txt,
        path test_vcf,
        path test_vcf_gz_tbi,
        path test2_targets_tsv_gz,
        path transcriptome_paf,
        path report_pdf,
        path samplesheet_csv
    ).set { all_resources_ch }

    REPORTS(reports.collect(), all_resources_ch)
}

workflow REPORTS {
    take:
    val multiqc_reports
    tuple path bin_depths_summary_tsv,
        path bin_summary_tsv,
        path CAPES_S7_log,
        path execution_trace_txt,
        path kraken2_report_txt,
        path pipeline_dag_svg,
        path binDepths_png,
        path krona_html,
        path transposed_tex,
        path transposed_tsv,
        path transposed_txt,
        path genome_fasta,
        path all_sites_fas,
        path baits_bed,
        path genome_dict,
        path genome_fasta_fai,
        path genome_gff3,
        path genome_gtf,
        path genome_sizes,
        path proteome_fasta,
        path test_baserecalibrator_table,
        path test_bed,
        path test_bedgraph,
        path test_bigwig,
        path test_paired_end_bam,
        path test_single_end_bam_readlist_txt,
        path test_vcf,
        path test_vcf_gz_tbi,
        path test2_targets_tsv_gz,
        path transcriptome_paf,
        path report_pdf,
        path samplesheet_csv

    main:
    """
    echo "Sleeping ${params.delay} seconds"
    sleep ${params.delay}
    echo "Copying all resource files to results directory for testing!"
    """
}
