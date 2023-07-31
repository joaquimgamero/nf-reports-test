params.outdir = "results"
params.delay = 300

// Channels
Channel.from(1..4).map{ [it, 'step_number'] }.set { steps_ch }
Channel.value("${projectDir}/resources/MultiQC Report.html").set { multiqc_html_ch }
Channel.of("${projectDir}/resources/bin_depths_summary.tsv").set { bin_depths_summary_tsv_ch }
Channel.of("${projectDir}/resources/bin_summary.tsv").set { bin_summary_tsv_ch }
Channel.of("${projectDir}/resources/CAPES_S7.log").set { CAPES_S7_log_ch }
Channel.of("${projectDir}/resources/execution_trace_2021-07-29_07-12-59.txt").set { execution_trace_txt_ch }
Channel.of("${projectDir}/resources/kraken2_report.txt").set { kraken2_report_txt_ch }
Channel.of("${projectDir}/resources/pipeline_dag_2021-07-29_07-12-59.svg").set { pipeline_dag_svg_ch }
Channel.of("${projectDir}/resources/SPAdesHybrid-CAPES_S7-binDepths.heatmap.png").set { binDepths_png_ch }
Channel.of("${projectDir}/resources/taxonomy.krona.html").set { krona_html_ch }
Channel.of("${projectDir}/resources/transposed_report.tex").set { transposed_tex_ch }
Channel.of("${projectDir}/resources/transposed_report.tsv").set { transposed_tsv_ch }
Channel.of("${projectDir}/resources/transposed_report.txt").set { transposed_txt_ch }
Channel.of("${projectDir}/resources/genome.fasta").set { genome_fasta_ch }
Channel.of("${projectDir}/resources/all_sites.fas").set { all_sites_fas_ch }
Channel.of("${projectDir}/resources/baits.bed").set { baits_bed_ch }
Channel.of("${projectDir}/resources/genome.dict").set { genome_dict_ch }
Channel.of("${projectDir}/resources/genome.fasta.fai").set { genome_fasta_fai_ch }
Channel.of("${projectDir}/resources/genome.gff3").set { genome_gff3_ch }
Channel.of("${projectDir}/resources/genome.gtf").set { genome_gtf_ch }
Channel.of("${projectDir}/resources/genome.sizes").set { genome_sizes_ch }
Channel.of("${projectDir}/resources/proteome.fasta").set { proteome_fasta_ch }
Channel.of("${projectDir}/resources/test.baserecalibrator.table").set { test_baserecalibrator_table_ch }
Channel.of("${projectDir}/resources/test.bed").set { test_bed_ch }
Channel.of("${projectDir}/resources/test.bedgraph").set { test_bedgraph_ch }
Channel.of("${projectDir}/resources/test.bigwig").set { test_bigwig_ch }
Channel.of("${projectDir}/resources/test.paired_end.bam").set { test_paired_end_bam_ch }
Channel.of("${projectDir}/resources/test.single_end.bam.readlist.txt").set { test_single_end_bam_readlist_txt_ch }
Channel.of("${projectDir}/resources/test.vcf").set { test_vcf_ch }
Channel.of("${projectDir}/resources/test.vcf.gz.tbi").set { test_vcf_gz_tbi_ch }
Channel.of("${projectDir}/resources/test2.targets.tsv.gz").set { test2_targets_tsv_gz_ch }
Channel.of("${projectDir}/resources/transcriptome.paf").set { transcriptome_paf_ch }
Channel.of("${projectDir}/resources/report.pdf").set { report_pdf_ch }
Channel.of("${projectDir}/resources/nfcore_chipseq110_samplesheet_test_full_6cols.csv").set { samplesheet_csv_ch }

process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val step
    path multiqc_html

    output:
    path("step_*/*.html"), emit: reports

    script:
    """
    echo "Copying MultiQC reports"
    mkdir step_$step
    cp $multiqc_html step_$step/
    """
}

process REPORTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path bin_depths_summary_tsv
    path bin_summary_tsv
    path CAPES_S7_log
    path execution_trace_txt
    path kraken2_report_txt
    path pipeline_dag_svg
    path binDepths_png
    path krona_html
    path transposed_tex
    path transposed_tsv
    path transposed_txt
    path genome_fasta
    path all_sites_fas
    path baits_bed
    path genome_dict
    path genome_fasta_fai
    path genome_gff3
    path genome_gtf
    path genome_sizes
    path proteome_fasta
    path test_baserecalibrator_table
    path test_bed
    path test_bedgraph
    path test_bigwig
    path test_paired_end_bam
    path test_single_end_bam_readlist_txt
    path test_vcf
    path test_vcf_gz_tbi
    path test2_targets_tsv_gz
    path transcriptome_paf
    path report_pdf
    path samplesheet_csv

    output:
    path("*"), emit: publish_ch

    script:
    """
    echo "Sleeping ${params.delay} seconds"
    sleep ${params.delay}
    echo "Copying all resource files to results directory for testing!"
    """
}

workflow {
    steps_ch
        .combine(multiqc_html_ch)
        .set { multiqc_inputs }

    MULTIQC(multiqc_inputs)

    bin_depths_summary_tsv_ch
        .combine(bin_summary_tsv_ch)
        .combine(CAPES_S7_log_ch)
        .combine(execution_trace_txt_ch)
        .combine(kraken2_report_txt_ch)
        .combine(pipeline_dag_svg_ch)
        .combine(binDepths_png_ch)
        .combine(krona_html_ch)
        .combine(transposed_tex_ch)
        .combine(transposed_tsv_ch)
        .combine(transposed_txt_ch)
        .combine(genome_fasta_ch)
        .combine(all_sites_fas_ch)
        .combine(baits_bed_ch)
        .combine(genome_dict_ch)
        .combine(genome_fasta_fai_ch)
        .combine(genome_gff3_ch)
        .combine(genome_gtf_ch)
        .combine(genome_sizes_ch)
        .combine(proteome_fasta_ch)
        .combine(test_baserecalibrator_table_ch)
        .combine(test_bed_ch)
        .combine(test_bedgraph_ch)
        .combine(test_bigwig_ch)
        .combine(test_paired_end_bam_ch)
        .combine(test_single_end_bam_readlist_txt_ch)
        .combine(test_vcf_ch)
        .combine(test_vcf_gz_tbi_ch)
        .combine(test2_targets_tsv_gz_ch)
        .combine(transcriptome_paf_ch)
        .combine(report_pdf_ch)
        .combine(samplesheet_csv_ch)
        .set { reports_inputs }

    REPORTS(reports_inputs)
}
