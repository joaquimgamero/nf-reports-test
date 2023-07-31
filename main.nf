params.outdir = "results"
params.delay = 10

Channel
    .value("${projectDir}/resources/MultiQC Report.html")
    .into { multiqc_html_ch }

Channel
    .from(1..4)
    .into { step_ch }

process GenerateReports {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path multiqc_html from multiqc_html_ch.first()
    val step from step_ch

    output:
    path("step_*/*.html") into reports

    script:
    """
    echo "Generating report for step $step"
    mkdir step_$step
    cp $multiqc_html step_$step/
    """
}

process Delay {
    input:
    val step from step_ch

    script:
    """
    echo "Delaying for step $step"
    sleep ${params.delay}
    """
}

workflow {
    GenerateReports()
    Delay()
}
