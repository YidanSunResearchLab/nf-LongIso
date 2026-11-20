/*
 * Normalised design table (sample,condition[,covariates]) for the R steps. The design matrix is
 * checked here (bin/check_design.R), before alignment, so that a confounded design fails early.
 *
 * Built by a process (not collectFile) so that the file is cached: a collectFile output is
 * re-created on every launch, which would make -resume re-run every step that reads it.
 */

process MAKE_DESIGN {
    tag "design"

    input:
    val lines           // header + one line per sample: sample,condition[,covariates]; values are validated
    val reference       // reference condition

    output:
    path "design.csv"  , emit: csv
    path "versions.yml", emit: versions

    script:
    def quoted = lines.collect { r -> "'${r}'" }.join(' ')
    """
    printf '%s\\n' ${quoted} > design.csv
    check_design.R --design=design.csv --design-formula='${params.design}' --reference=${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
    END_VERSIONS
    """

    stub:
    def quoted = lines.collect { r -> "'${r}'" }.join(' ')
    """
    printf '%s\\n' ${quoted} > design.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
