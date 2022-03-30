process SENTIEON_BWAMEM {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(reads)
    path fasta
    path fai
    path index

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.bai'), emit: bai
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ?: "\'@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\'"
    def sentieon_exe = params.sentieon_install_dir ? "${params.sentieon_install_dir}/sentieon" : 'sentieon'
    """
    source sentieon_init.sh SENTIEON_LICENSE_BASE64

    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    $sentieon_exe \\
        bwa \\
        mem \\
        -R $read_group \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        $args \\
        | $sentieon_exe \\
            util \\
            sort \\
            -r $fasta \\
            -o $prefix \\
            -t $task.cpus \\
            $args2 \\
            --sam2bam \\
            -i -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$($sentieon_exe driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$($sentieon_exe bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
