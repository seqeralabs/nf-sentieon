process SENTIEON_PLOT {
    tag "$meta.id"
    label 'process_low'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path('*.pdf'), emit: pdf
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sentieon_exe = params.sentieon_install_dir ? "${params.sentieon_install_dir}/sentieon" : 'sentieon'
    """
    source sentieon_init.sh SENTIEON_LICENSE_BASE64

    $sentieon_exe \\
        plot \\
        $args \\
        $input_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$($sentieon_exe driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
