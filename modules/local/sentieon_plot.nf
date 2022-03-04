process SENTIEON_PLOT {
    tag "$meta.id"
    label 'process_low'
    label 'sentieon'

    secret 'sentieon_license_text'

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
    set -eu
    export SENTIEON_LICENSE=\$(mktemp)
    echo -e "\$sentieon_license_text" > \$SENTIEON_LICENSE

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
