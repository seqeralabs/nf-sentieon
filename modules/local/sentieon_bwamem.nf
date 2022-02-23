process SENTIEON_BWAMEM {
    tag "$meta.id"
    label 'process_high'

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
    def read_group = "-R \'@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\'"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    sentieon bwa mem \\
        $read_group \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        $args \\
        | sentieon util sort \\
            -r $fasta \\
            -o $prefix \\
            -t $task.cpus \\
            $args2 \\
            --sam2bam \\
            -i -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
