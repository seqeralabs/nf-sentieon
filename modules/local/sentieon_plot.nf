process SENTIEON_PLOT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path('*.pdf'), emit: pdf
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sentieon plot \\
        $args \\
        $input_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
